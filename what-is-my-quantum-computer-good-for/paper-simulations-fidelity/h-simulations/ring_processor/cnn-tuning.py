import numpy as np
import json

import tensorflow as tf
import kerastuner as kt
from keras.losses import MeanSquaredError
from keras.callbacks import EarlyStopping
from tensorflow.keras import layers, models

import ml

import argparse

dtypes = ['train','test','validate']

exp_num = 0
exp_path = f'./experiment_{exp_num}/'
proc_path = exp_path + '/processed_inputs_and_outputs/'
model_path = exp_path + '/models/'
pred_path = exp_path + '/predictions/'
history_path = exp_path + '/training_histories/'

print(f'You are working with experiment {exp_path}.')

dtypes = ['train','test','validate']

parser = argparse.ArgumentParser()
parser.add_argument('--TRIALS', type=int, default = 100)
parser.add_argument('--EPOCHS', type=int, default = 100)

args = parser.parse_args()
N_TRIALS = args.TRIALS
EPOCHS = args.EPOCHS
print('We are using the following args: {}'.format(args))

# a is the actual values while b is the predictions
def bce(a, b):
    return tf.keras.losses.BinaryCrossentropy(from_logits = False)(a,b).numpy()

def mse(a, b):
    return np.mean((a - b) ** 2)

def pick_layers(hp, model, layer):
    if layer != 0:
        shape = model.shape
    #print(shape[0],shape[1])
  
    if layer == 0:
        k_depth = hp.Int('k_depth_{}'.format(layer), min_value = 3, max_value = 10)
        k_width = hp.Int('k_width_{}'.format(layer), min_value = 1, max_value = 4)
        p_width = hp.Int('p_width_{}'.format(layer), min_value = 1, max_value = 2)
        p_depth = hp.Int('p_depth_{}'.format(layer), min_value = 1, max_value = 30)
    else:
        k_depth = hp.Int('k_depth_{}'.format(layer), min_value = 1, max_value = min(3,shape[2]))
        k_width = hp.Int('k_width_{}'.format(layer), min_value = 1, max_value = shape[1])
        p_width = hp.Int('p_width_{}'.format(layer), min_value = 1, max_value = 2)
        p_depth = hp.Int('p_depth_{}'.format(layer), min_value = 1, max_value = 6)
    
    return k_width, k_depth, p_width, p_depth
    

def tune_model(hp):
    input_layer = layers.Input(shape=input_shape)
    
    # pooling_windows = [(1,5), (5, 1), (5,5)]
    
    NUM_CONV_LAYERS = hp.Int('num_conv_layers', min_value = 1, max_value = 2)
    NUM_EXTRA_DENSE_LAYERS = hp.Int('num_dense_layers', min_value = 1, max_value = 5) 
    dropout = hp.Boolean('dropout')
    
    for i in range(NUM_CONV_LAYERS):
        filters = hp.Choice('filters_{}'.format(i), values = [5, 10, 15], ordered = True)
        if i == 0:
            k_width, k_depth, p_width, p_depth = pick_layers(hp, None, i)
            cnn_model = layers.Conv2D(filters = filters, kernel_size = (k_width, k_depth), padding = 'same', activation = 'relu')(input_layer)
        else:
            k_width, k_depth, p_width, p_depth = pick_layers(hp, cnn_model, i)
            cnn_model = layers.Conv2D(filters = filters, kernel_size = (k_width, k_depth), padding = 'same', activation = 'relu')(cnn_model)
        cnn_model = layers.AveragePooling2D((p_width, p_depth))(cnn_model)
        
    cnn_model = layers.Flatten()(cnn_model)
    
    initial_units = hp.Choice('units_0', values = [10,50,100], ordered = True)
    # initial_units = hp.Choice('units_0', min_value = 256, max_value = 512, step = 16)
    merged_model = layers.Dense(units = initial_units, activation = 'relu')(cnn_model)
    if dropout is True:    
        merged_model = layers.Dropout(rate = .1)(merged_model)

    for i in range(1,NUM_EXTRA_DENSE_LAYERS):
        merged_model = layers.Dense(units = 10, activation = 'relu')(merged_model)
        # merged_model = layers.Dense(units = int(d/2), activateion = 'relu')(merged_model)
        if dropout is True and (i+1 != NUM_EXTRA_DENSE_LAYERS):    
            merged_model = layers.Dropout(rate = .1)(merged_model)

    merged_model = layers.Dense(units = 1, activation = 'relu')(merged_model)
    
    model = models.Model(inputs = [input_layer, ], outputs=merged_model)
    
    model.compile(optimizer='Adam',
              loss=MeanSquaredError(),
              metrics=[tf.keras.metrics.MeanSquaredError(name="MSE")])
    
    return model

if __name__=="__main__":

    ### Load up the circuits and the (scaled) infidelities ###
    circuits = np.load(proc_path+'/processed_high_fidelity_circuits.npz')
    infidelities = np.load(proc_path+'/processed_infidelities.npz')  

    ### Load up the meta information ###
    with open(exp_path+'/meta.json', 'r') as f:
        meta = json.load(f)

    num_qubits = meta['num_qubits']
    max_error_weight = meta['max_error_weight']
    num_hops = meta['num_hops']
    num_channels = meta['num_channels']
    if meta['geometry'] == 'ring':
        adj_matrix = ml.newtools.ring_adj_matrix(num_qubits)
        laplace = ml.newtools.laplace_from_qubit_graph(adj_matrix)
    else:
        print(f'You have not implemented the {geometry} geometry.')

    len_encoding = num_qubits * num_channels

    ### Re-shape the circuits
    x_circs = {dt: circuits[dt][:, :, :len_encoding] for dt in dtypes}
    x_c_t = {dt: np.transpose(x_circs[dt], axes = (0, 2, 1)) for dt in dtypes}
    x_c_t_r = {dt: x_c_t[dt].reshape(x_c_t[dt].shape[0], num_qubits, num_channels, -1) for dt in dtypes}
    
    cnn_circs = {dt: np.transpose(x_c_t_r[dt], axes = (0, 1, 3, 2)) for dt in dtypes}
    cnn_fids = infidelities 

    np.savez_compressed(proc_path + '/processed_cnn_high_fidelity_circuits.npz', **cnn_circs)

    input_shape = cnn_circs['train'][0].shape

    ### Initialize the tuner ###
            
    tuner = kt.BayesianOptimization(hypermodel = tune_model,
                                objective='val_loss',
                                max_trials=N_TRIALS, num_initial_points=1, directory=exp_path,
                                project_name='tuning_trials', overwrite=True)
        
    callback = EarlyStopping(monitor='val_loss', patience=10, min_delta = 50)

    ### Run the tuner ###
    tuner.search(x = cnn_circs['train'], y = cnn_fids['train'], epochs = EPOCHS, batch_size = 32, 
             validation_data = (cnn_circs['validate'], cnn_fids['validate']), callbacks = [callback], verbose = 1)
    
    print('The search process is over.')

    ### Extract the best hyperparameters and build the best model ###      
    best_hyperparameters = tuner.get_best_hyperparameters()[0]
    model = tuner.hypermodel.build(best_hyperparameters)

    ### Train the best model ###
    
    callback = EarlyStopping(monitor='val_loss', patience=10, min_delta = 50)
            
    history = model.fit(cnn_circs['train'], cnn_fids['train'], epochs=EPOCHS, 
                        validation_data=(cnn_circs['validate'], cnn_fids['validate']), callbacks = [callback])

    ### Save the trained model, its predictions, and its training history    
    model.save(model_path + '/tuned_cnn')

    predictions = {dtype: model.predict(cnn_circs[dtype]) for dtype in dtypes}
    np.savez_compressed(pred_path + '/processed_cnn_predictions.npz', **predictions)

    with open(history_path + '/cnn-history.json', 'w') as f:
        json.dump(history.history, f)

    ### Load up the mirror circuits ###
    o_circuits = np.load(proc_path + '/processed_mirrored_circuits.npz')['circuits']
    o_infidelities = np.load(proc_path + '/processed_mirrored_infidelities.npz')['infidelities']

    o_circuits = o_circuits[:, :, :len_encoding]

    o_c_t = np.transpose(o_circuits, axes = (0, 2, 1))
    o_c_t_r = o_c_t.reshape(o_c_t.shape[0], num_qubits, num_channels, -1)
    mirrored_cnn_circuits = np.transpose(o_c_t_r, axes = (0, 1, 3, 2))

    np.savez_compressed(proc_path + '/processed_cnn_mirrored_circuits.npz', circuits = mirrored_cnn_circuits)

    ### Make and save the mirror circuit predictions ###
    mirrored_preds = model.predict(mirrored_cnn_circuits).reshape((len(o_circuits),))
    np.savez_compressed(pred_path + '/processed_cnn_mirrored_predictions.npz', predictions = mirrored_preds)

    print(f'Done with the coherent noise tuning and training on experiment {exp_num}')
       
