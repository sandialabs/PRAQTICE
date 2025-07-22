"""Utilities for FPR/GR resource testing from modelpacks
"""

from collections import namedtuple
import itertools
import multiprocessing as mp
import numpy as np
import time

import pygsti
from pygsti.modelpacks import GSTModelPack
from pygsti.circuits import Circuit
from pygsti.protocols import StandardGSTDesign as GSTDesign
from pygsti.protocols.protocol import CombinedExperimentDesign

def create_combined_exp_design_from_modelpack(smqpack, maxLs, random_fracs=None, num_random_samples=1, base_seed=2021,
                                              verbosity=1, rand_interp=False, num_random_interp_samples=None,
                                              interp_fracs=None, gr_types=None, fpr_types=None, mem_limit=2,
                                              use_mp_pool='False', num_mp_workers=None, numpergermpower=100,
                                              numpergermpoweriters=100, eigvaltol=0.5, condtol=2,
                                              param_blk_size=100, num_soln_returned=5, type_soln_returned='best', 
                                              retry_for_smaller=True):
    """Create a CombinedExperimentDesign for all FPR/GR combinations of interest from a modelpack.

    Parameters
    ----------
    smqpack: GSTModelPack
        A standard multi-qubit (smq) GST model pack containing fiducial pairs and germs
    
    maxLs: list of int
        The max lengths for all GST designs

    random_fracs: list of float or None
        Fractions to use when generating random FPR samples. If None, defaults to 1/8 to 7/8.
    
    num_random_samples: int
        Number of samples for random FPR. Defaults to 1, but probably want more than that
        for better statistics
    
    base_seed: int or None
        Seed for multiple RandomState objects created
    
    verbosity: int
        If 0, no additional output. If 1, print general logging output.
        verbosity - 1 is passed into FPR/GR selection routines,
        so verbosity > 1 will also print additional FPR/GR output.
    
    Returns
    -------
    comb_design: CombinedExperimentDesign
        A combined experiment design with keys as tuple of strs ('<FPR type>', '<GR type>')
        and entries as the relevant GST design (always a subset of the full GST design)
    """
    assert(isinstance(smqpack, GSTModelPack)), "Must pass in GSTModelPack for access to fids/germs"
    
    all_gr_types = ['Full', 'Lite', 'Bare']
    if gr_types is None:
        gr_types = all_gr_types
    for gr_type in gr_types:
        assert gr_type in all_gr_types, f"GR type must be in {all_gr_types}, not {gr_type}"

    all_fpr_types = ['Full','PerGerm', 'PerGermPower', 'PerGermPowerRand']
    if fpr_types is None:
        fpr_types = ['Full', 'PerGerm', 'PerGermPowerRand']
    for fpr_type in fpr_types:
        assert fpr_type in all_fpr_types, f"FPR type must be in {all_fpr_types}, not {fpr_type}"

    #assert 'Full' in gr_types and 'Full' in fpr_types, "Full,Full must be in the experiment designs"

    target_model = smqpack.target_model()
    #target_model.sim = 'map'
    sim = pygsti.forwardsims.MatrixForwardSimulator(param_blk_sizes=(param_blk_size,param_blk_size))
    target_model.sim = sim
    prep_fids = smqpack.prep_fiducials()
    meas_fids = smqpack.meas_fiducials()
    tpb0_labels = target_model.state_space.labels[0]
    germs = {
        'Full': smqpack.germs(lite=False),
        'Lite': smqpack.germs(lite=True),
        'Bare': [Circuit([op], line_labels=tpb0_labels) for op in target_model.operations.keys()] # Works only for ExplicitOpModel
    }
    
    #Initialize the combined design to None.
    comb_design= None
    
    # Full FPR
    if 'Full' in fpr_types:
        if verbosity: print('>> Generating exp designs for full FPR <<')
        for gr in gr_types:
            if comb_design is None:
                full_design = GSTDesign(target_model, prep_fids, meas_fids, germs[gr], maxLs)
                comb_design = CombinedExperimentDesign.from_edesign(full_design, ("Full",gr))
            else:
                comb_design['Full', gr] = GSTDesign(target_model, prep_fids, meas_fids, germs[gr], maxLs)    

    # Global FPR
    #if 'Global' in fpr_types:
    #    if verbosity: print('\n>> Generating exp designs for global FPR <<')
    #    global_fid_pairs = {
    #        'Full': smqpack.global_fidpairs,
    #        'Lite': smqpack.global_fidpairs_lite,
    #    }
    #    for gkey in gr_types:
    #        fid_pairs = global_fid_pairs.get(gkey, None)

     #       # For ones that are not provided by modelpack, do global FPR ourselves
     #       if fid_pairs is None:
     #           if verbosity: print(f'  >> Generating new global fid pairs for {gkey} germs <<')
     #           fid_pairs = pygsti.alg.find_sufficient_fiducial_pairs(
     #               target_model, prep_fids, meas_fids, germs[gkey],
     #               search_mode="random", n_random=100, seed=base_seed,
     #               verbosity=verbosity-1, mem_limit=mem_limit*1024**3, minimum_pairs=2)
     #       else:
     #           if verbosity: print(f'  >> Using global fid pairs available in modelpack for {gkey} germs <<')
     #       
     #       assert(fid_pairs is not None), "Fiducial pairs should not be None"
     #       
     #       comb_design['Global', gkey] = GSTDesign(target_model, prep_fids, meas_fids, germs[gkey], maxLs,
     #                                               fiducial_pairs=fid_pairs)
    
    # TODO: Random FPR should probably be reworked for robustness testing... Maybe Corey can help us with this
    # Random FPR
    # Three variants: global and per-germ and iterpolating
    # To make seeding easy, initialize one random state for ALL random sampling
    # This repeats a bit of code (i.e. not using keep_fraction kwarg), but allows for better RNG control
    rand_state = np.random.default_rng(base_seed)

    def sample_fid_pairs(num_prep_fids, num_meas_fids, keep_fraction, rand_state):
        all_fidpairs = list(itertools.product(range(num_prep_fids), range(num_meas_fids)))
        num_to_keep = int(round(keep_fraction * len(all_fidpairs)))
        chosen_indices = rand_state.choice(range(len(all_fidpairs)), size=num_to_keep, replace=False)
        return [all_fidpairs[idx] for idx in chosen_indices]
    
    #define a sampling procedure which interpolates between the full FPR and per-germ-power FPR
    #by treating the per-germ-power FPR as a minimal set and randomly sampling pairs of fiducials
    #to add in addition.
    
    #create a dictionary of fiducial pairs indexed by the germ label which corresponds to the "full"
    #FPR for the three different germ sets. 
    full_fpr_per_germ_power_dict= {(germ,power):list(itertools.product(range(len(prep_fids)), range(len(meas_fids)))) for germ in germs['Full'] for power in maxLs}
    full_fpr_per_germ_power_lite_dict= {(germ,power):list(itertools.product(range(len(prep_fids)), range(len(meas_fids)))) for germ in germs['Lite'] for power in maxLs}
    full_fpr_per_germ_power_bare_dict= {(germ,power):list(itertools.product(range(len(prep_fids)), range(len(meas_fids)))) for germ in germs['Bare'] for power in maxLs}
    
    full_fpr_per_germ_power_dicts= {'Full':full_fpr_per_germ_power_dict, 'Lite':full_fpr_per_germ_power_lite_dict, 'Bare':full_fpr_per_germ_power_bare_dict}
    
    #Now define a function for doing the per-germ-power random FPR:
    def sample_fid_pairs_per_germ_power(full_fid_set, bare_germs, rand_frac, rand_state):
        """
        full_fid_set: a dictionary indexed by (germ,power) tuples corresponding to a "full" FPR for any given germ set set.
        rand_frac: The fraction of the fiducial pairs in full_fid_set keep.
        rand_state: an RNG object
        """        
        outputdict={}
        for germ_power_key, germ_power_fid_list in full_fid_set.items():
            #If the germ is one of the bare germs and the germ power is 1 then add in all of the fiducial pairs.
            #This should wind up being done automatically by the experiment design constructor.
            keep_idx_list=rand_state.choice(np.arange(len(germ_power_fid_list)),replace=False, size=round(rand_frac*len(germ_power_fid_list)))
            outputdict[germ_power_key]=[germ_power_fid_list[idx] for idx in keep_idx_list]
        return outputdict
        
    
#     def sample_fid_pairs_interp(full_fid_set, per_germ_power_fid_set, addl_fraction, rand_state):
#             """
#             full_fid_set: a dictionary indexed by (germ,power) tuples corresponding to a "full" FPR for any given germ set set.
#             per_germ_power_fid_set: a dictionary indexed by (germ, power) tuples corresponding to the per-germ per-power FPR scheme.
#             addl_fraction: The fraction of the fiducial pairs in full_fid_set-per_germ_power_fid_set to add back into our fiducial set.
#             rand_state: an RNG object
#             """
#             #Start by casting the dictionaries as sets, taking the set difference and recasting the result as a dictionary.
#             #python doesn't like doing the set differences when the values are lists since they aren't hashable, so cast the values as
#             #tuples first.
#             full_fid_set_tuple={key:tuple(value) for key, value in full_fid_set.items()}
#             per_germ_power_fid_set_tuple={key:tuple(value) for key, value in per_germ_power_fid_set.items()}
#             #turns out dictionary views support set operations.
#             fid_set_diff= full_fid_set_tuple.items()-per_germ_power_fid_set_tuple.items()
#             #print(fid_set_diff)
#             fid_set_diff_dict= {key:list(value) for key, value in fid_set_diff}
#             #print(fid_set_diff_dict)
            
#             #pick a random fraction of these to add back into our set of fiducials in addition to those in
#             #per_germ_power_fid_set
#             fid_set_diff_keys= list(fid_set_diff_dict.keys())
#             fid_set_diff_values= list(fid_set_diff_dict.values())
#             extra_fid_indices=rand_state.choice(np.arange(len(fid_set_diff_keys)),replace=False, size=round(addl_fraction*len(fid_set_diff_keys)))
#             #append these fiducial pairs onto output dictionary.
#             output_dict= per_germ_power_fid_set
#             for i in extra_fid_indices:
#                 output_dict[fid_set_diff_keys[i]]= fid_set_diff_values[i]
#             return output_dict
            
    if 'PerGermPowerRand' in fpr_types:          
        if random_fracs is None:
            random_fracs = np.array(range(1,8)) / 8.0

        if verbosity:
            print(f'\n>> Generating exp designs for random FPR <<')
            print(f'   ({num_random_samples} samples with base seed {base_seed})')
        for i in range(num_random_samples):
            for keep_frac in random_fracs:
                for gkey in gr_types:
                    #Sample per-germ fidpairs
                    per_germ_power_rand_fid_pairs=sample_fid_pairs_per_germ_power(full_fpr_per_germ_power_dicts[gkey],germs['Bare'], keep_frac, rand_state)
                    
                    if comb_design is None:
                        per_germ_power_rand_design=GSTDesign(target_model, prep_fids, meas_fids, germs[gkey], maxLs,
                                                            fiducial_pairs=per_germ_power_rand_fid_pairs)
                        comb_design = CombinedExperimentDesign.from_edesign(per_germ_power_rand_design, (f'PerGermRand{keep_frac:0.2f}-{i}',                                                                                                            gkey))
                    else:
                        comb_design[f'PerGermRand{keep_frac:0.2f}-{i}', gkey] = GSTDesign(target_model, prep_fids, meas_fids, germs[gkey],                                                                                              maxLs, fiducial_pairs=per_germ_power_rand_fid_pairs)
    
    
    # Per-Germ FPR
    if 'PerGerm' in fpr_types:
        if verbosity: print('\n>> Generating exp designs for per-germ FPR <<')
        
        #different version depending on whether we are parallelizing over the germs.
        
        assert(len(eigvaltol)==len(condtol))
        
        for num_randomizations in numpergermpower:

            for tol in zip(eigvaltol, condtol):

                if use_mp_pool==False:
                    fid_pairs_per_germ={}
                    for gkey in gr_types:
                        if verbosity: print(f'  >> Generating new per-germ fid pairs for {gkey} germs <<')
                        fid_pairs_per_germ[gkey] = pygsti.alg.find_sufficient_fiducial_pairs_per_germ(target_model, prep_fids, meas_fids, 
                                                                                                      germs[gkey], maxLs, search_mode="random",
                                                                                                      constrain_to_tp=True, n_random=num_randomizations,
                                                                                                      min_iterations=numpergermpoweriters, 
                                                                                                      base_loweig_tol= tol[0], condition_number_tol=tol[1],
                                                                                                      seed=base_seed, verbosity=verbosity-1,
                                                                                                      num_soln_returned= num_soln_returned, type_soln_returned=type_soln_returned,
                                                                                                      retry_for_smaller=retry_for_smaller, mem_limit=mem_limit*1024**3)

                        assert(fid_pairs_per_germ[gkey] is not None), "Fiducial pairs should not be None"

                        if comb_design is None:
                            per_germ_design=GSTDesign(target_model, prep_fids, meas_fids, germs[gkey],
                                                                        fiducial_pairs=fid_pairs_per_germ[gkey])
                            comb_design = CombinedExperimentDesign.from_edesign(per_germ_design, (f'PerGermE{tol[0]:0.3f}C{tol[1]:0.2f}N{num_randomizations}', gkey))
                        else:
                            comb_design[f'PerGermE{tol[0]:0.3f}C{tol[1]:0.2f}N{num_randomizations}', gkey] = GSTDesign(target_model, prep_fids, meas_fids, germs[gkey], maxLs,
                                                                        fiducial_pairs=fid_pairs_per_germ[gkey])
                #otherwise parallelize over the germs in each germ set using multiprocessing.
                else:
                    fid_pairs_per_germ={}
                    for gkey in gr_types:
                        fid_pairs_per_germ[gkey]={}
                        if verbosity: print(f'  >> Generating new per-germ fid pairs for {gkey} germs <<')

                        #we want to parallelize over the germs in the given germ set and then combine the germ-power fiducial dictionaries
                        #need to start by constructing a task_queue that we can pass into a helper function to unpack.
                        fcn_args= [(target_model, prep_fids, meas_fids, germ, maxLs, "random", True, num_randomizations, numpergermpoweriters, tol[0], tol[1], base_seed, verbosity-1, num_soln_returned, type_soln_returned, retry_for_smaller, mem_limit*1024**3)  for germ in germs[gkey] ]

                        with mp.Pool(num_mp_workers) as pool:
                            it = pool.imap(_run_per_germ_fpr_on_node, fcn_args)

                            # Totally unnecessary, but give us some nice progress output
                            for i in range(len(fcn_args)):
                                try:
                                    fid_pairs_per_germ[gkey].update(next(it))
                                    print(f'Completed protocol {i+1} / {len(fcn_args)} at {time.asctime()}')
                                except StopIteration:
                                    # Should not make it here
                                    break  

                        assert(fid_pairs_per_germ[gkey] is not None), "Fiducial pairs should not be None"

                        if comb_design is None:
                            per_germ_design=GSTDesign(target_model, prep_fids, meas_fids, germs[gkey], maxLs,
                                                                        fiducial_pairs=fid_pairs_per_germ[gkey])
                            comb_design = CombinedExperimentDesign.from_edesign(per_germ_design, (f'PerGermE{tol[0]:0.3f}C{tol[1]:0.2f}N{num_randomizations}', gkey))
                        else:
                            comb_design[f'PerGermE{tol[0]:0.3f}C{tol[1]:0.2f}N{num_randomizations}', gkey] = GSTDesign(target_model, prep_fids, meas_fids, germs[gkey], maxLs,
                                                                        fiducial_pairs=fid_pairs_per_germ[gkey]) 



    # Per-Germ-Power FPR (new, varies with germ AND L)
    if 'PerGermPower' in fpr_types:
        if verbosity: print('\n>> Generating exp designs for per-germ-power FPR <<')
        
        #different version depending on whether we are parallelizing over the germs.
        if use_mp_pool==False:
            fid_pairs_per_germ_power={}
            for gkey in gr_types:
                if verbosity: print(f'  >> Generating new per-germ fid pairs for {gkey} germs <<')
                fid_pairs_per_germ_power[gkey] = pygsti.alg.find_sufficient_fiducial_pairs_per_germ_power(target_model, prep_fids, meas_fids, germs[gkey], maxLs,
                    search_mode="random", constrain_to_tp=True, n_random=numpergermpower, min_iterations=numpergermpoweriters, base_loweig_tol= tol[0], condition_number_tol=tol[1], seed=base_seed, verbosity=verbosity-1, mem_limit=mem_limit*1024**3)

                assert(fid_pairs_per_germ_power[gkey] is not None), "Fiducial pairs should not be None"

                if comb_design is None:
                    per_germ_power_design=GSTDesign(target_model, prep_fids, meas_fids, germs[gkey], maxLs,
                                                                fiducial_pairs=fid_pairs_per_germ_power[gkey])
                    comb_design = CombinedExperimentDesign.from_edesign(per_germ_power_design, (f'PerGermPowerE{tol[0]:0.3f}C{tol[1]:0.2f}', gkey))
                else:
                    comb_design[f'PerGermPowerE{tol[0]:0.3f}C{tol[1]:0.2f}', gkey] = GSTDesign(target_model, prep_fids, meas_fids, germs[gkey], maxLs,
                                                                fiducial_pairs=fid_pairs_per_germ_power[gkey])
        #otherwise parallelize over the germs in each germ set using multiprocessing.
        else:
            fid_pairs_per_germ_power={}
            for gkey in gr_types:
                fid_pairs_per_germ_power[gkey]={}
                if verbosity: print(f'  >> Generating new per-germ fid pairs for {gkey} germs <<')
                    
                #we want to parallelize over the germs in the given germ set and then combine the germ-power fiducial dictionaries
                #need to start by constructing a task_queue that we can pass into a helper function to unpack.
                fcn_args= [(target_model, prep_fids, meas_fids, germ, maxLs, "random", True, numpergermpower, numpergermpoweriters, tol[0], tol[1], base_seed, verbosity-1, mem_limit*1024**3)  for germ in germs[gkey] ]
                
                with mp.Pool(num_mp_workers) as pool:
                    it = pool.imap(_run_per_germ_power_fpr_on_node, fcn_args)

                    # Totally unnecessary, but give us some nice progress output
                    for i in range(len(fcn_args)):
                        try:
                            fid_pairs_per_germ_power[gkey].update(next(it))
                            print(f'Completed protocol {i+1} / {len(fcn_args)} at {time.asctime()}')
                        except StopIteration:
                            # Should not make it here
                            break  
                    
                assert(fid_pairs_per_germ_power[gkey] is not None), "Fiducial pairs should not be None"

                if comb_design is None:
                    per_germ_power_design=GSTDesign(target_model, prep_fids, meas_fids, germs[gkey], maxLs,
                                                                fiducial_pairs=fid_pairs_per_germ_power[gkey])
                    comb_design = CombinedExperimentDesign.from_edesign(per_germ_power_design, (f'PerGermPowerE{tol[0]:0.3f}C{tol[1]:0.2f}', gkey))
                else:
                    comb_design[f'PerGermPowerE{tol[0]:0.3f}C{tol[1]:0.2f}', gkey] = GSTDesign(target_model, prep_fids, meas_fids, germs[gkey], maxLs,
                                                                fiducial_pairs=fid_pairs_per_germ_power[gkey]) 
                   
                
    #check if we have enabled randomized interpolations
    if rand_interp:
        #check to makes sure the num_random_interp_samples argument has been set, else throw an exception.
        if num_random_interp_samples is None:
            raise ValueError('Yo, buddy, if you want me to generate interpolated random FPR samples you need to tell me how many samples to take by initializing num_random_interp_samples!')
        if interp_fracs is None:
            raise ValueError('Yo, buddy, if you want me to generate interpolated random FPR samples you need to tell me a list of fractions of additional fiducials to evaluate!')    
        if verbosity:
            print(f'\n>> Generating exp designs for interpolated random FPR <<')
            print(f'   ({num_random_interp_samples} samples with base seed {base_seed})')
        #sample num_random_interp_samples times
        for i in range(num_random_interp_samples):
            for frac in interp_fracs:
                for gkey in gr_types:
                    interp_fid_pairs=sample_fid_pairs_interp(full_fpr_per_germ_power_dicts[gkey], fid_pairs_per_germ_power[gkey], frac, rand_state)
                    #Add this interpolation to the combined experiment design.
                    # TODO SS to CIO: Should this be PerGermPowerRandInterp?
                    comb_design[f'PerGermRandInterp{frac:0.2f}-{i}', gkey] = GSTDesign(target_model, prep_fids, meas_fids, germs[gkey], maxLs,
                    fiducial_pairs=interp_fid_pairs)
            
        
    return comb_design

def print_fid_germ_summary(cdesign):
    """Helper function to print high-level CombinedExperimentDesign fid-pair/germ statistics.

    Parameters
    ----------
    cdesign: CombinedExperimentDesign
        The design to print statistics for
    """
    assert(isinstance(cdesign, CombinedExperimentDesign)), "Expected a CombinedExperimentDesign"
    print(f'Combined design with total circuits: {len(cdesign.all_circuits_needing_data)}')
    print(f'Number of subdesigns: {len(list(cdesign.keys()))}')
    
    for (fpr, gr), subdesign in cdesign.items():
        print(f'\nSubdesign with FPR {fpr} and GR {gr}')
        assert(isinstance(subdesign, GSTDesign)), "Expected subdesign to be GSTDesign"

        num_prep_fids = len(subdesign.prep_fiducials)
        num_meas_fids = len(subdesign.meas_fiducials)
        print(f'  Prep Fids: {num_prep_fids} Meas Fids: {num_meas_fids}')
        
        if isinstance(subdesign.fiducial_pairs, dict):
            num_total_fid_pairs = 0
            all_fid_pairs = set()
            for fid_pairs in subdesign.fiducial_pairs.values():
                num_total_fid_pairs += len(fid_pairs)
                all_fid_pairs.update(fid_pairs)

            per_germ_avg = num_total_fid_pairs / len(subdesign.fiducial_pairs.keys())
            num_uniq_fid_pairs = len(all_fid_pairs)
            print(f'  Total Unique Fid-Pairs: {num_uniq_fid_pairs} (Avg per Germ: {per_germ_avg:0.2f})')
        else:
            num_fid_pairs = len(subdesign.fiducial_pairs) if subdesign.fiducial_pairs is not None \
                else num_prep_fids*num_meas_fids
            print(f'  Total Fid-Pairs: {num_fid_pairs}')
        
        num_germs = len(subdesign.germs)
        print(f'  Germs: {num_germs}  Max Lengths: {subdesign.maxlengths}')
        
        num_circs = len(subdesign.all_circuits_needing_data)
        print(f'  Total Circuits: {num_circs}')

# Helper function for multiprocessing Pool
def _run_protocol_on_node(args):
    node, protocol, memlimit, path = args
    
    results = protocol.run(node.data, memlimit, None)

    return path, results

# Helper function for multiprocessing Pool
def _run_per_germ_fpr_on_node(args):
    target_model, prep_fids, meas_fids, germ, maxLs, search_mode, constrain_to_tp, numpergermpower, numpergermpoweriters, eigvaltol, condtol, base_seed, verbosity, num_soln_returned,type_soln_returned, retry_for_smaller, memlimit  = args
    
    result=pygsti.alg.find_sufficient_fiducial_pairs_per_germ(target_model, prep_fids, meas_fids, [germ],
                                                              search_mode=search_mode, constrain_to_tp=constrain_to_tp,
                                                              n_random=numpergermpower, min_iterations=numpergermpoweriters,
                                                              base_loweig_tol= eigvaltol, condition_number_tol=condtol,
                                                              seed=base_seed, verbosity=verbosity,
                                                              num_soln_returned= num_soln_returned,
                                                              type_soln_returned=type_soln_returned,
                                                              retry_for_smaller=retry_for_smaller, mem_limit=memlimit)
    return result

# Helper function for multiprocessing Pool
def _run_per_germ_power_fpr_on_node(args):
    target_model, prep_fids, meas_fids, germ, maxLs, search_mode, constrain_to_tp, numpergermpower, numpergermpoweriters, eigvaltol, condtol, base_seed, verbosity, memlimit  = args
    
    result=pygsti.alg.find_sufficient_fiducial_pairs_per_germ_power(target_model, prep_fids, meas_fids, [germ], maxLs, search_mode=search_mode, constrain_to_tp=constrain_to_tp,
                    n_random=numpergermpower, min_iterations=numpergermpoweriters, base_loweig_tol= eigvaltol, condition_number_tol=condtol, seed=base_seed, verbosity=verbosity,
                    mem_limit=memlimit)

    return result

class MultiprocessPoolRunner(pygsti.protocols.ProtocolRunner):
    """
    Runs a single protocol on every data node that has no sub-nodes (possibly separately for each pass).

    This is like SimpleRunner except it uses a multiprocessing.Pool to distribute over the designs.

    Parameters
    ----------
    protocol : Protocol
        The protocol to run.

    protocol_can_handle_multipass_data : bool, optional
        Whether `protocol` is able to process multi-pass data, or
        if :class:`MultiPassProtocol` objects should be created
        implicitly.

    edesign_type : type or 'all'
        Only run `protocol` on leaves with this type.  (If 'all', then
        no filtering is performed.)
    """
    def __init__(self, protocol, protocol_can_handle_multipass_data=False, edesign_type='all', num_workers=1):
        """
        Create a new SimpleRunner object, which runs a single protocol on every
        'leaf' of the data-tree.

        Parameters
        ----------
        protocol : Protocol
            The protocol to run.

        protocol_can_handle_multipass_data : bool, optional
            Whether `protocol` is able to process multi-pass data, or
            if :class:`MultiPassProtocol` objects should be created
            implicitly.

        edesign_type : type or 'all'
            Only run `protocol` on leaves with this type.  (If 'all', then
            no filtering is performed.)
        
        num_workers : int, optional
            Number of workers to use in the Pool

        Returns
        -------
        SimpleRunner
        """
        self.protocol = protocol
        self.edesign_type = edesign_type
        self.do_passes_separately = not protocol_can_handle_multipass_data
        self.num_pool_workers = num_workers

    def run(self, data, memlimit=None, comm=None):
        """
        Run all the protocols specified by this protocol-runner on `data`.

        Parameters
        ----------
        data : ProtocolData
            The input data.

        memlimit : int, optional
            A rough per-processor memory limit in bytes.
        
        comm : mpi4py.MPI.Comm, optional
            Not used, but kept to keep interface with other ProtocolRunners

        Returns
        -------
        ProtocolResultsDir
        """
        ret = pygsti.protocols.ProtocolResultsDir(data)  # creates entire tree of nodes

        # Not using comm, i.e. no parallelism underneath

        # Task queue
        task_queue = [] # Contains (node, protocol, memlimit, path to node)

        # Select all nodes that need protocols needed
        def visit_node(node, path):
            if len(node.data) > 0:
                for subname, subnode in node.items():
                    new_path = path.copy()
                    new_path.append(subname)
                    visit_node(subnode, new_path)
            elif node.data.is_multipass() and self.do_passes_separately:
                implicit_multipassprotocol = pygsti.protocols.MultiPassProtocol(self.protocol, name=self.protocol.name)
                task_queue.append((node, implicit_multipassprotocol, memlimit,  path))
            elif self.edesign_type == 'all' or isinstance(node.data.edesign, self.edesign_type):
                task_queue.append((node, self.protocol, memlimit, path))
            else:
                pass  # don't add this node as a task, since the experiment design has the wrong type
        visit_node(ret, [])

        # A place to hold all results as (path to node, results)
        completed_queue = []

        # Multiprocess over queue now that we have collected all tasks
        with mp.Pool(self.num_pool_workers) as pool:
            it = pool.imap(_run_protocol_on_node, task_queue)

            # Totally unnecessary, but give us some nice progress output
            for i in range(len(task_queue)):
                try:
                    completed_queue.append(next(it))
                    print(f'Completed protocol {i+1} / {len(task_queue)} at {time.asctime()}')
                except StopIteration:
                    # Should not make it here
                    break
            
        # Slot all results into their proper place in the full tree structure
        def update_node(node, path, results):
            if len(path):
                update_node(node[path[0]], path[1:], results) # Traverse down to proper node
                return
            
            node.for_protocol[self.protocol.name] = results
        
        for path, results in completed_queue:
            update_node(ret, path, results)
        
        return ret

# Helper function for post-processing analysis
def _run_analysis_for_df_row(args):
    estimate, noisy_model, metadata = args

    # Gauge optimize
    go_iters = []
    #metadata[-2] should be the number of iterations.
    #metadata[6] should be the estimate key. We'll skip gauge optimizing the cptp estimates for now.
    current_estimate_key= metadata[6]
    #When we encounter the CPTP estimate we'll just return the un-gauge optimized estimate
    for k in range(metadata[-2]):
        if current_estimate_key== 'CPTP':
            gaugeopt_iter = estimate.models['iteration ' + str(k) + ' estimate']
        else:
            gaugeopt_iter = pygsti.algorithms.gaugeopt_to_target(estimate.models['iteration ' + str(k) + ' estimate'], noisy_model, item_weights={'gates': 1.0, 'spam': 0.0})
        go_iters.append(gaugeopt_iter)
    
    # Pull out operation matrices
    op_label = metadata[-1] # Last entry should be gate label
    iter_ops = [goi.operations[op_label].to_dense() for goi in go_iters]
    true_op = noisy_model.operations[op_label].to_dense()

    # Jamiolkowski trace distance
    jtd = [pygsti.tools.jtracedist(op, true_op) for op in iter_ops]

    evals = np.array([sorted(np.linalg.eigvals(op)) for op in iter_ops])
    true_evals = np.array(sorted(np.linalg.eigvals(true_op)))

    N = len(true_evals)

    # Eigenvalues MAEs
    eval_MAE = [sum(abs(iter_evals - true_evals)) / N for iter_evals in evals] 
    
    # Angle MAEs, based on the eigenvalues
    angle_MAE = [sum(abs(np.angle(iter_evals) - np.angle(true_evals))) / N for iter_evals in evals]

    max_lengths = [2**i for i in range(len(jtd))]
    
    #diamond distance:
    diamond_dist= [pygsti.tools.diamonddist(op, true_op) for op in iter_ops]
    return metadata[:-1]  + [str(op_label), max_lengths, jtd, eval_MAE, angle_MAE, diamond_dist]


class FisherObjective(pygsti.objectivefns.objectivefns.TimeIndependentMDCObjectiveFunction):
    def __init__(self, model, circuits, num_counts=1e5, resource_alloc=None, verbosity=3):
        dsDummy = pygsti.data.DataSet(static=False)
        outcomes = list(model.povms['Mdefault'].keys()); num_outcomes = len(outcomes)
        for c in circuits: dsDummy.add_count_dict(c, {outcome_lbl: num_counts / num_outcomes for outcome_lbl in outcomes})
        dsDummy.done_adding_data()
        objDummy = namedtuple('DummyRawObjFn', 'printer')(pygsti.baseobjs.VerbosityPrinter(verbosity, comm=resource_alloc))

        mdc_store = self._create_mdc_store(model, dsDummy, circuits, resource_alloc,
                                           method_names=('fn', 'hessian'), array_types=(), verbosity=verbosity)
        super().__init__(objDummy, mdc_store, None, verbosity)

    def _hessian_from_block(self, hprobs, dprobs12, probs, counts, total_counts, freqs, resource_alloc):
        # fisher_info = N * (1/probs * dprobs12 - hprobs) (re-using dprobs12 and hprobs mem)
        if resource_alloc.is_host_leader:  # hprobs, dprobs12, and probs are shared among resource_alloc procs 
            dprobs12 *= (1.0 / probs)[:, None, None]
            fisher_info = dprobs12; fisher_info -= hprobs
            fisher_info *= total_counts[:, None, None]
        else:
            fisher_info = dprobs12
        resource_alloc.host_comm_barrier()  # let root procs finish their work before moving fwd                                                                                                                                        
        return np.sum(fisher_info, axis=0)


