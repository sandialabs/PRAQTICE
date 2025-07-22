"""Run GST on the combined FPR/GR experiment designs.
"""

# Set the environment to prevent fork-bombing as early as possible
# Must be before numpy imports, etc.
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['GOTO_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
os.environ['NUMEXPR_NUM_THREADS'] = '1'
os.environ['VECLIB_MAXIMUM_THREADS'] = '1'
os.environ['MKL_NUM_THREADS'] = '1'

import argparse
import copy
import importlib
import inspect
import multiprocessing as mp
import numpy as np
import pandas as pd
import pathlib
import pickle
import pygsti
import sys
import time

sys.path.append('..')
import utils
import model_utils


def rankprint(string, args, comm=None):
    if (comm is None or comm.Get_rank() == 0) and not args.no_verbose:
        print(string)

def get_args():
    """Helper function to build argparse.
    """
    # Include defaults in help message
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Global options
    parser.add_argument('--data_dir', type=str, required=True, help='Directory for experiment design/results serialization') # Takes place of file_label
    parser.add_argument('--modelpack', type=str, default='smq1Q_XYI',
                        help='Name of the GST modelpack to use for CombinedExperimentDesign and noisy model generation')
    parser.add_argument('--no_verbose', action='store_true', help='Turn off base script verbosity')

    # Compute options (listed here because they are handled early)
    parser.add_argument('--use_mpi', action='store_true', help='Load mpi4py and initialize a comm object')
    parser.add_argument('--mem_per_core', type=int, default=4, help='Memory limit per core in GiB')
    parser.add_argument('--use_mp_pool', action='store_true', help='Use a multiprocessing.Pool to distribute GST runs')
    parser.add_argument('--num_mp_workers', type=int, default=1, help='Number of workers available in the multiprocessing.Pool')
    
    # Experiment design generation options
    parser.add_argument('--run_edesign_gen', action='store_true', help='Flag to run experiment design generation')
    parser.add_argument('--maxL_power', type=int, default=4, help='Power of 2 to use for the largest max length')
    parser.add_argument('--maxL_power_step', type=int, default=1, help='Stepsize for maxL_power to use')
    parser.add_argument('--gr_types', type=str, nargs='+', default=['Full', 'Lite', 'Bare'], help='Types of germ reduction (GR) to include')
    parser.add_argument('--fpr_types', type=str, nargs='+', default=['Full', 'Global', 'PerGerm', 'PerGermPower'],
        help='Types of fiducial pair reduction (FPR) to include')
    parser.add_argument('--fpr_fractions', type=float, nargs='+', default=[0.25, 0.5, 0.75], help='Fractions to use for FPR random sampling')
    parser.add_argument('--fpr_samples', type=int, default=0, help='Number of samples to use for FPR random sampling')
    parser.add_argument('--fpr_seed', type=int, default=2021, help='Seed for random FPR sampling')
    parser.add_argument('--fpr_verbosity', type=int, default=1, help='FPR verbosity (also passed into FPR algorithms)')
    parser.add_argument('--skip_global_fpr', action='store_true', help='Whether to skip global FPR in edesign generation')
    parser.add_argument('--num_per_germ_power_randomizations', type=int, nargs='+', default=[100], help='Number of random fiducial pairs, per fiducial pair set size, to try during the per-germ power FPR algorithm')
    parser.add_argument('--per_germ_power_min_iterations', type=int, default=100, help='Minimum number of random fiducial pairs (per set size) to try during the per-germ power FPR algorithm before allowing the algorithm to exit early if it has found a sufficient candidate solution.')
    parser.add_argument('--per_germ_power_eigval_tol', type=float, nargs='+', default= [.5], help='(for per-germ power FPR) Gives the multiplicative reduction in the magnitude of the minimum        eigenvalue relative to the value for the full fiducial set the user is willing to tolerate.')
    parser.add_argument('--per_germ_power_cond_tol', type=float, nargs='+', default= [2], help='(for per-germ power FPR) Gives the multiplicative increase in the magnitude of the condition number relative to the value for the full fiducial set the user is willing to tolerate.')
    parser.add_argument('--num_soln_returned', type=int, default=5, help='Number of initial candidate solutions to return before refinement during the per-germ and per-germ power FPR algorithms')
    parser.add_argument('--type_soln_returned', type=str, default='best', help='Type of initial candidate solutions to return before refinement during the per-germ and per-germ power FPR algorithms')
    parser.add_argument('--retry_for_smaller', action='store_true', help='Whether to try resampling from the initial candidate solutions to try and find a smaller design during per-germ and per-germ power FPR algorithms')
    
    # Noisy model generation options
    parser.add_argument('--run_model_gen', action='store_true', help='Flag to run noisy model generation')
    parser.add_argument('--model_file', type=str, help='Filename for serializing the model with format <name>_<index>.pkl.' +
                                                       '\nMust be specified when --run_model_gen or --run_gst is active')
    parser.add_argument('--model_json', type=str, nargs='+', help='JSON files use to generate Lindblad error generator dictionaries')
    parser.add_argument('--model_gen_seed', type=int, default=2021, help='Seed for noisy model generation')
    parser.add_argument('--model_gen_verbose', action='store_true', help='Whether to print verbose messages during model generation')

    # Simulation/GST options
    parser.add_argument('--run_gst', action='store_true', help='Flag to run GST')
    # --model_file is also required here
    parser.add_argument('--gst_verbosity', type=int, default=2, help='GST protocol verbosity')
    parser.add_argument('--gst_modes', type=str, default='full TP,CPTP,Target,Truth', help='GST modes (parameterizations) to run. Truth is the noisy model')
    parser.add_argument('--use_wildcard', action='store_true', help='Whether to include wildcard error as a GST bad fit option')
    parser.add_argument('--num_clicks', type=int, default=1000, help='Number of simulated data samples to generate')
    parser.add_argument('--datagen_seed', type=int, default=2021, help='Seed for dataset generation')
    parser.add_argument('--param_blk_size', type=int, help='param_blk_sizes to use with blocked GST (if None, run unblocked)')

    # Analysis options
    parser.add_argument('--run_analysis', action='store_true', help='Flag to run analysis')
    parser.add_argument('--sep_edesign_mode', action='store_true', help='Flag to expect separate folders for each edesign.')
    parser.add_argument('--df_file', help='Name of file for Pandas DataFrame (either ending in .csv or .h5).\nMust be specified when --run_analysis is active')

    # FIM options
    parser.add_argument('--run_fim', action='store_true', help='Flag to run FIM')
    parser.add_argument('--fim_verbosity', type=int, default=3, help='FIM calculation verbosity')
    parser.add_argument('--fim_num_counts', type=int, default=1000, help='Number of counts in FIM calculation')
    parser.add_argument('--fim_outfile', type=str, default='fim', help='Name of NumPy file for FIM (either no ending or .npy).')

    return parser.parse_args()

def create_exp_design(args, comm=None):
    max_lengths = [2**i for i in range(0, args.maxL_power+1, args.maxL_power_step)]
    gr_types = args.gr_types
    fpr_types = args.fpr_types
    fpr_fractions = args.fpr_fractions
    num_fpr_samples = args.fpr_samples
    fpr_seed = args.fpr_seed
    fpr_verbosity = args.fpr_verbosity
    save_dir = args.data_dir
    mem_per_core = args.mem_per_core
    mp_pool = args.use_mp_pool
    mp_workers = args.num_mp_workers
    numpergermpower=args.num_per_germ_power_randomizations 
    numpergermpoweriters=args.per_germ_power_min_iterations
    param_blk_size=args.param_blk_size
    eigvaltol= args.per_germ_power_eigval_tol 
    condtol= args.per_germ_power_cond_tol
    num_soln_returned=args.num_soln_returned
    type_soln_returned=args.type_soln_returned
    if args.retry_for_smaller:
        retry_for_smaller=True
    else:
        retry_for_smaller=False
    
    # Load modelpack
    try:
        modelpack = importlib.import_module(f'pygsti.modelpacks.{args.modelpack}')
        rankprint(f'Successfully loaded modelpack {args.modelpack}', args, comm)
    except ModuleNotFoundError as err:
        raise Exception(f'Could not load "{args.modelpack}" from pygsti.modelpacks, please select a valid modelpack') from err

    if (comm is None or comm.Get_rank() == 0) and not args.no_verbose:
        print('\n>> Experiment design generation options <<')
        print(f'  Max lengths: {max_lengths}')
        print(f'  Germ reduction (GR) types: {gr_types}')
        print(f'  Fiducial pair reduction (FPR) types: {fpr_types}')
        print(f'  Random FPR fractions: {fpr_fractions}')
        print(f'  Number of random FPR samples: {num_fpr_samples}')
        print(f'  Seed for random FPR sampling: {fpr_seed}')
        print(f'  FPR verbosity: {fpr_verbosity}')
        print(f'  FPR mem limit: {mem_per_core}')
        print(f'  Parameter block size: {param_blk_size}')
        print(f'  Save directory: {save_dir}\n')

    # Generate CombinedExperimentDesign
    exp_design = None
    if comm is None or comm.Get_rank() == 0:
        start = time.time()
        exp_design = utils.create_combined_exp_design_from_modelpack(modelpack, max_lengths, fpr_fractions,
                                                                     num_fpr_samples, fpr_seed, fpr_verbosity,
                                                                     gr_types=gr_types, fpr_types=fpr_types,
                                                                     mem_limit=mem_per_core, use_mp_pool= mp_pool, num_mp_workers=mp_workers, 
                                                                     numpergermpower=numpergermpower, numpergermpoweriters=numpergermpoweriters, eigvaltol=eigvaltol, condtol=condtol,
                                                                     param_blk_size=param_blk_size, num_soln_returned=num_soln_returned, type_soln_returned=type_soln_returned, 
                                                                     retry_for_smaller=retry_for_smaller)
        utils.print_fid_germ_summary(exp_design)
        end=time.time()
        print('Elapsed Time Edesign Gen: ', end-start)
        exp_design.write(save_dir) # Write exp design to file
    
    if comm is not None:
        exp_design = comm.bcast(exp_design, root = 0)
    
    return exp_design

def create_noisy_model(args, comm=None):
    assert(args.model_file is not None), "Must specify --model_file with --run_model_gen"
    save_file = args.data_dir + '/' + args.model_file
    jsonfiles = args.model_json
    seed = args.model_gen_seed
    verbose = args.model_gen_verbose

    # Load modelpack
    try:
        modelpack = importlib.import_module(f'pygsti.modelpacks.{args.modelpack}')
        rankprint(f'Successfully loaded modelpack {args.modelpack}', args, comm)
    except ModuleNotFoundError as err:
        raise Exception(f'Could not load "{args.modelpack}" from pygsti.modelpacks, please select a valid modelpack') from err

    # All generation is done on one rank
    noisy_model = None
    if comm is None or comm.Get_rank() == 0:
        if not args.no_verbose:
            print('\n>> Noisy model generation options <<')
            print(f'  JSON inputfiles: {jsonfiles}')
            print(f'  Seed for RNG sampling: {seed}')
            print(f'  Model generation verbose: {verbose}')
            print(f'  Save file: {save_file}\n')

        target_model = modelpack.target_model(gate_type='static')

        rand_state = np.random.RandomState(seed)

        err_dicts = []
        for jsonfile in jsonfiles:
            noise_dict = model_utils.sample_errgen_json(jsonfile, rand_state, verbose, target_model.state_space.num_qubits)
            err_dicts.append(noise_dict)
        
        # Combine dicts if needed
        full_err_dict = model_utils.combine_dicts(*err_dicts)

        # Current noisy_model assignment
        noisy_model = model_utils.create_noisy_model_from_lindblad_dict(target_model, full_err_dict, verbose)
        
        # TODO: GAHHH I really don't like pickling, but there's a bug in the load_model somewhere...
        # Probably solvable, but trying to finish this up first
        with open(save_file, 'wb') as f:
            pickle.dump(noisy_model, f)
    
    if comm is not None:
        noisy_model = comm.bcast(noisy_model, root=0)
    
    return noisy_model

def simulate_data_and_run_gst(args, comm=None):
    assert(args.model_file is not None), "Must specify --model_file with --run_gst"

    try:    
        #exp_design = pygsti.protocols.StandardGSTDesign.from_dir(args.data_dir)
        exp_design= pygsti.io.read_edesign_from_dir(args.data_dir)
        #print(exp_design)
        #exp_design= exp_design['Full', 'Full']
        rankprint(f'Successfully loaded experiment design from {args.data_dir}', args, comm)
    except FileNotFoundError as err:
        raise Exception(f'Could not find an edesign to load in {args.data_dir}. Specify --run_edesign_gen to generate an edesign.') from err
    except AttributeError as err:
        raise Exception('AttributeErrors while doing loads are often due to mismatched pickles. ' + \
                            'Please run with the version of PyGSTi that generated the experiment design (or generate a new design).') from err
    
    #Load modelpack
    try:
        modelpack = importlib.import_module(f'pygsti.modelpacks.{args.modelpack}')
        rankprint(f'Successfully loaded modelpack {args.modelpack}', args, comm)
    except ModuleNotFoundError as err:
        raise Exception(f'Could not load "{args.modelpack}" from pygsti.modelpacks, please select a valid modelpack') from err
    
    #try:
    #    exp_design.processor_spec= modelpack.processor_spec([0,1])
    #except:
    #    pass
        
    #print(exp_design)
    try:
        with open(f'{args.data_dir}/{args.model_file}', 'rb') as f:
            noisy_model = pickle.load(f)
        rankprint(f'Successfully loaded noisy model from {args.data_dir}/{args.model_file}', args, comm)
    except:
        try:
            rankprint('Maybe it is located in the base folder instead? Trying to fall back.', args, comm)
            with open(f'{args.model_file}', 'rb') as f:
                noisy_model = pickle.load(f)
            rankprint(f'Successfully loaded noisy model from {args.data_dir}/{args.model_file}', args, comm)
        except  FileNotFoundError as err:
            raise Exception(f'Could not find the model in {args.data_dir}/{args.model_file} or in pwd. Specify --run_model_gen to generate a noisy model.') from err

    save_dir = args.data_dir
    memlimit = args.mem_per_core*(1024)**3
    param_blk_size = args.param_blk_size
    mp_pool = args.use_mp_pool
    mp_workers = args.num_mp_workers

    modes = args.gst_modes
    verbosity = args.gst_verbosity
    badfit_options = None
    if args.use_wildcard:
        badfit_options = {'actions': ['wildcard'], 'wildcard_methods': ['barrier']}
    num_clicks = args.num_clicks
    datagen_seed = args.datagen_seed

    # Autogenerate results name to enable easy recordkeeping
    model_name = args.model_file.replace('.pkl','')
    results_name = f'model_{model_name}_datagen_{datagen_seed}'

    rankprint('Simulating data', args, comm)
    start = time.time()

    if comm is None or comm.Get_rank() == 0:
        if not args.no_verbose:
            print('\n>> GST options <<')
            print(f'  GST modes (Truth = noisy model): {modes}')
            print(f'  Parameter block size: {param_blk_size}')
            print(f'  Verbosity: {verbosity}')
            print(f'  Badfit options: {badfit_options}')
            print(f'  Number of datagen samples: {num_clicks}')
            print(f'  Seed for datagen samples: {datagen_seed}')
            print(f'  Protocol name in results: {results_name}')
            print(f'  Save directory: {save_dir}\n')

    dataset = pygsti.data.simulate_data(noisy_model, exp_design.all_circuits_needing_data, num_clicks, seed=datagen_seed, comm=comm)
    pygsti.io.write_dataset(save_dir + '/' + results_name + '.txt', dataset)
    data = pygsti.protocols.ProtocolData(copy.deepcopy(exp_design),dataset)
    end = time.time()
    rankprint(f'Simulation time: {end - start}', args, comm)

    gaugeopt = pygsti.protocols.GSTGaugeOptSuite(gaugeopt_target=noisy_model, gaugeopt_suite_names='stdgaugeopt')
    gst_proto = pygsti.protocols.StandardGST(modes=modes, target_model= modelpack.target_model(),
                                            badfit_options=badfit_options,
                                            verbosity=verbosity,
                                            gaugeopt_suite=gaugeopt,
                                            models_to_test={'Truth':noisy_model},
                                            name=results_name)

    if param_blk_size is not None:

        # Hack in blocked target model creation
        key = list(data.edesign.keys())[0]
        blocked_target_model = data.edesign[key].create_target_model()
        blocked_target_model.sim = pygsti.forwardsims.MatrixForwardSimulator(param_blk_sizes=(param_blk_size, param_blk_size))

        def create_blocked_target_model():
            return blocked_target_model

        for k in data.edesign.keys():
            data.edesign[k].create_target_model = create_blocked_target_model

        noisy_model.sim = pygsti.forwardsims.MatrixForwardSimulator(param_blk_sizes=(param_blk_size, param_blk_size))
   

    # Switch runner for parallelism mode
    if mp_pool:
        gst_runner = utils.MultiprocessPoolRunner(gst_proto, False, pygsti.protocols.StandardGSTDesign, mp_workers)
    else:
        gst_runner = pygsti.protocols.SimpleRunner(gst_proto, False, pygsti.protocols.StandardGSTDesign)

    rankprint('Running GST', args, comm)
    start = time.time()
    results = gst_proto.run(data, memlimit=memlimit, comm=comm)
    end = time.time()
    rankprint(f'GST time: {end - start}', args, comm)
        
    # Write results to file
    if comm is None or comm.Get_rank() == 0:
        # Remove blocked matrix function from result so we can serialize
        results.data = pygsti.protocols.ProtocolData(exp_design, dataset)

        results.write(save_dir) 

def fisher_information_matrix(args, comm=None):
    try:    
        exp_design = pygsti.protocols.StandardGSTDesign.from_dir(args.data_dir)
        rankprint(f'Successfully loaded experiment design from {args.data_dir}', args, comm)
    except FileNotFoundError as err:
        raise Exception(f'Could not find an edesign to load in {args.data_dir}. Specify --run_edesign_gen to generate an edesign.') from err
    except AttributeError as err:
        raise Exception('AttributeErrors while doing loads are often due to mismatched pickles. ' + \
                            'Please run with the version of PyGSTi that generated the experiment design (or generate a new design).') from err

    # Load modelpack
    try:
        modelpack = importlib.import_module(f'pygsti.modelpacks.{args.modelpack}')
        rankprint(f'Successfully loaded modelpack {args.modelpack}', args, comm)
    except ModuleNotFoundError as err:
        raise Exception(f'Could not load "{args.modelpack}" from pygsti.modelpacks, please select a valid modelpack') from err

    save_dir = args.data_dir
    memlimit = args.mem_per_core*(1024)**3
    param_blk_size = args.param_blk_size

    verbosity = args.fim_verbosity
    num_clicks = args.fim_num_counts

    results_name = args.fim_outfile

    # Generate SPAM-depolarized model
    depol_model = modelpack.target_model('TP').depolarize(spam_noise=1e-3)
    depol_model.sim = pygsti.forwardsims.MatrixForwardSimulator(param_blk_sizes=(param_blk_size, param_blk_size))

    if comm is None or comm.Get_rank() == 0:
        if not args.no_verbose:
            print('\n>> FIM options <<')
            print(f'  Edesign directory: {args.data_dir}')
            print(f'  Modelpack: {args.modelpack}')
            print(f'  Memory per core (GB): {args.mem_per_core}')
            print(f'  Parameter block size: {param_blk_size}')
            print(f'  Verbosity: {verbosity}')
            print(f'  Results name: {results_name}')

    rankprint('Calculating Fisher information matrix', args, comm)

    start = time.time()
    ralloc = pygsti.baseobjs.ResourceAllocation(comm=comm, mem_limit=memlimit)
    fobj = utils.FisherObjective(depol_model, exp_design.all_circuits_needing_data, resource_alloc=ralloc,
        num_counts=num_clicks, verbosity=verbosity)
    
    fim = fobj.hessian()

    end = time.time()
    rankprint(f'FIM time: {end - start}', args, comm)
        
    # Write results to file
    if comm is None or comm.Get_rank() == 0:
        np.save(f'{save_dir}/{results_name}', np.array(fim))

def run_analysis(args, comm=None):
    data_dir = args.data_dir
    df_file = args.df_file
    use_csv = False
    if df_file.endswith('.csv'):
        use_csv = True
    mp_pool = args.use_mp_pool
    mp_workers = args.num_mp_workers
    num_iters= len( [2**i for i in range(0, args.maxL_power+1, args.maxL_power_step)] )

    print('Read arguments.')
    
    assert(comm is None), "Analysis is not parallelized via MPI. Use the --use_mp_pool and --num_mp_workers options instead"

    #try:    
    #    combined_data = pygsti.io.read_results_from_dir(args.data_dir)
    #    rankprint(f'Successfully loaded results from {args.data_dir}', args, comm)
    #except FileNotFoundError as err:
    #    raise Exception(f'Could not find results to load in {args.data_dir}. Specify --run_gst to generate GST results.') from err
    #except AttributeError as err:
    #    raise Exception('AttributeErrors while doing loads are often due to mismatched pickles. ' + \
    #                        'Please run with the version of PyGSTi that generated the results (or generate new results).') from err
    
    if comm is None or comm.Get_rank() == 0:
        if not args.no_verbose:
            print('\n>> Analysis options <<')
            print(f'  Dataframe file: {df_file}')
            print(f'  Dataframe format: ' + ('CSV' if use_csv else 'HDF5'))

    path = pathlib.Path(data_dir)
    
    print('Path: ', path)
    
    if args.sep_edesign_mode:
        result_dirs= [data_dir[:-1]]
    else:
        result_dirs = sorted([d for d in path.iterdir() if d.is_dir() and '_' in d.name])
    
    print('Results Directories: ', result_dirs)
    
    
    # Build tasks, i.e. traverse to all estimates
    tasks = []
    #for keys, data in combined_data.items():
    for rdir in result_dirs:
        try:    
            data = pygsti.io.read_results_from_dir(rdir)
            rankprint(f'Successfully loaded results from {rdir}', args, comm)
        except FileNotFoundError as err:
            raise Exception(f'Could not find results to load in {rdir}. Specify --run_gst to generate GST results.') from err
        except AttributeError as err:
            raise Exception('AttributeErrors while doing loads are often due to mismatched pickles. ' + \
                                'Please run with the version of PyGSTi that generated the results (or generate new results).') from err
        
        #fpr_key, gr_key = keys
        if args.sep_edesign_mode:
            #The separated data directories for each edesign will be assumed to have the FPR and GR keys as the last two elements of the name separated by underscores.
            split_rdir= rdir.split('_')
            fpr_key = split_rdir[-2]
            gr_key = split_rdir[-1]
        else:
            fpr_key, gr_key = rdir.name.split('_')
        # Split FPR samples if needed
        if '-' in fpr_key:
            fpr_key, fpr_index = fpr_key.split('-')
        else:
            fpr_index = 0
        
        for prot_key, results in data.for_protocol.items():
            # Split protocol key into model and datagen samples
            entries = prot_key.split('_')
            split_index = entries.index('datagen')
            model = '_'.join(entries[1:split_index-1])
            model_index = entries[split_index-1]
            datagen_index = entries[split_index+1]
            
            try:
                model_file = f'{data_dir}/{model}_{model_index}.pkl'
                with open(model_file, 'rb') as f:
                    noisy_model = pickle.load(f)
                #rankprint(f'Successfully loaded noisy model from {model_file}', args, comm)
            except FileNotFoundError as err:
                raise Exception(f'Could not find the model in {model_file}. Ensure result keys match model_<name>_<index>_datagen_<index>.') from err

            for est_key, est in results.estimates.items():
                for op_label in noisy_model.operations.keys():
                    metadata = [fpr_key, fpr_index, gr_key, model, model_index, datagen_index, est_key, num_iters, op_label]
                    tasks.append([est, noisy_model, metadata])
    
    # Multiprocess over queue now that we have collected all tasks
    rows = []
    if mp_pool:
        ctx = mp.get_context('spawn')
        with ctx.Pool(mp_workers) as pool:
            # Unnecessary version with some nice progress output
            it = pool.imap(utils._run_analysis_for_df_row, tasks)
            
            for i in range(len(tasks)):
                rows.append(next(it))
                print(f'Completed row {i+1} / {len(tasks)} at {time.asctime()}')
    else:
        # Serial version
        for i, task in enumerate(tasks):
            rows.append(utils._run_analysis_for_df_row(task))
            print(f'Completed row {i+1} / {len(tasks)} at {time.asctime()}')

    # Build DataFrame
    df = pd.DataFrame(rows, columns=['FPR', 'FPR Index', 'GR', 'Model', 'Model Index', 'Datagen Index',
                                     'Parameterization', 'Num GST Iterations','Gate', 'L', 'JTD', 'Eval MAE', 'Angle MAE', 'Diamond Dist'])

    # Write results to file
    if comm is None or comm.Get_rank() == 0:
        if use_csv:
            df.to_csv(df_file)
        else:
            df.to_hdf(df_file, key='GST', mode='w')

if __name__ == '__main__':
    args = get_args()

    # Initialize MPI if needed
    if args.use_mpi and args.use_mp_pool:
        raise SyntaxError('Cannot specify both --use_mpi and --use_mp_pool')

    comm = None
    if args.use_mpi:
        from mpi4py import MPI
        comm = MPI.COMM_WORLD

    # Handle global options
    if (comm is None or comm.Get_rank() == 0) and not args.no_verbose:
        print('>> Global options <<')
        print('  MPI:', 'Disabled' if comm is None else 'Enabled')
        if comm is not None: print(f' Number of MPI workers: {comm.Get_size()}')
        print('  Multiprocessing.Pool:', 'Disabled' if not args.use_mp_pool else 'Enabled')
        if args.use_mp_pool: print(f'  Number of Pool workers: {args.num_mp_workers}')
        print(f'  Memory limit: {args.mem_per_core} GiB per core')
        print(f'  Data directory: {args.data_dir}')
        print(f'  Modelpack: {args.modelpack}')
        print(f'  Run experiment design generation: {args.run_edesign_gen}')
        print(f'  Run noisy model generation: {args.run_model_gen}')
        print(f'  Run GST: {args.run_gst}')
        print(f'  Run analysis: {args.run_analysis}')
    
    start_all = time.time()

    # Create experiment design if requested
    if args.run_edesign_gen:
        rankprint('Running experiment design generation...', args, comm)
        start = time.time()
        exp_design = create_exp_design(args, comm)
        end = time.time()
        rankprint(f'Finished experiment design generation in {end - start:.2f} s', args, comm)

    # Create noisy model if requested
    if args.run_model_gen:
        rankprint('Running noisy model generation...', args, comm)
        start = time.time()
        noisy_model = create_noisy_model(args, comm)
        end = time.time()
        rankprint(f'Finished noisy model generation in {end - start:.2f} s', args, comm)

    # Run GST if requested    
    if args.run_gst:
        rankprint('Running GST...', args, comm)
        start = time.time()
        simulate_data_and_run_gst(args, comm)
        end = time.time()
        rankprint(f'Finished GST in {end - start:.2f} s', args, comm)

    # Run result extraction if requested
    if args.run_analysis:
        rankprint('Running analysis...', args, comm)
        start = time.time()
        run_analysis(args, comm)
        end = time.time()
        rankprint(f'Finished analysis in {end - start:.2f} s', args, comm)

    # Run FIM if requested
    if args.run_fim:
        rankprint('Running FIM...', args, comm)
        start = time.time()
        fisher_information_matrix(args, comm)
        end = time.time()
        rankprint(f'Finished FIM in {end - start:.2f} s', args, comm)

    end_all = time.time()
    
    rankprint(f'Completed all requested stages in {end_all - start_all:.2f} s', args, comm)
