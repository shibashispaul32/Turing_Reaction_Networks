import re
import numpy as np
import pandas as pd
from te_bifurcation import model2te, run_bf
import sobol_seq
from rrplugins import Plugin
auto = Plugin("tel_auto2000")
import pickle
import reaction_nets
import pathlib

def run_cont(nodes, nreps=10000, if_only_seeds=True, pert_id=0):
    '''
    Perform bifurcation analysis (numerical continuation)
    nodes: a list of Node objects, each for a characteristic complex
    nreps: number of parameter sets for each reaction network
    if_only_seeds: only perform continuation for the 10 reactions known to generate bifurcations
    pert_id: perturbation id (0-11). 0 is unperturbed. The id determines the type of perturbation
    '''
    outputpath = pathlib.Path('./run_'+str(pert_id).zfill(2)+'/data_hb_lp')
    outputpath.mkdir(parents=True, exist_ok=True)
    print('Bifurcation output path: '+str(outputpath))
    seed_paths = [
            ['AAB', 1], ['ABC',0], ['AABB', 2], ['AABB', 3], ['AAAB', 0],
            ['AAAB', 2], ['AAAB', 3], ['AABC', 2], ['AABC', 3],
            ['ABCD', 1]
        ]
    for n in nodes:
        n.num_bifs = []
        for ip, p in enumerate(n.paths[::]):
            n.num_bifs.append({'LP':0, 'HB':0})
            if if_only_seeds:
                if not [n.pattern, ip] in seed_paths:
                    continue
            vs, pars = reaction_nets.make_odes(p)
            model = {'vars':vs, 'pars':dict(zip(pars, [1]*len(pars))), 'fns': {}, 'aux': [],
                          'name':n.pattern+'_'+str(ip)}
            model['pars']['koff'] = 1
            for p in model['pars']:
                if p.startswith('K_'):
                    model['pars'][p] =  700
            n.models.append(model)
            r = model2te(model, ics=dict(zip(vs.keys(), [0.1]*len(vs))))
            degps = [p for p in r.ps() if re.match(r'[a-h]_', p)]
            control_k = 'k_A'
            if 'B' in r.fs():
                control_k = 'k_B'
            else:
                control_k = 'k_A'
            ks = [p for p in r.ps() if (re.match(r'k_', p) and p != control_k)]
            Ks = [p for p in r.ps() if (re.match(r'K_', p))]
            ss = sobol_seq.i4_sobol_generate(len(degps)+len(ks)+len(Ks), int(nreps)+0)[0:,:]
            l = 10**(ss[:,:len(degps)]*(np.log10(10)-np.log10(0.1)) + np.log10(0.1))
            if pert_id == 9:
                l = 10**(ss[:,:len(degps)]*(np.log10(1)-np.log10(0.1)) + np.log10(0.1))
            if pert_id == 10:
                l = 10**(ss[:,:len(degps)]*(np.log10(10)-np.log10(1)) + np.log10(1))
            for i in range(int(nreps)):
                if i % 10 == 0:
                    print('Network: '+n.pattern+'-'+str(ip), 'Set ID: '+str(i))
                else:
                    print('Set ID: '+str(i))
                for j, p in enumerate(degps):
                    r[p] = l[i,j]
                    if pert_id == 11:
                        r[p] = l[i,0]
                    model['pars'][p] = r[p]
                for k, p in enumerate(ks):
                    r[p] = 10**(ss[i,k+len(degps)]*(np.log10(10)-np.log10(0.1)) + np.log10(0.1))
                    if pert_id == 7:
                        r[p] = 1
                    elif pert_id == 8:
                        r[p] = 10**(ss[i,0+len(degps)]*(np.log10(10)-np.log10(0.1)) + np.log10(0.1))
                    model['pars'][p] = r[p]
                for K, p in enumerate(Ks):
                    r[p] = 10**(ss[i,K+len(degps)+len(ks)]*(np.log10(1000)-np.log10(10)) + np.log10(10))
                    if pert_id == 1:
                        r[p] = 10
                    elif pert_id == 2:
                        r[p] = 100
                    elif pert_id == 3:
                        r[p] = 1000
                    elif pert_id == 4:
                        if len(p) == 5:
                            r[p] = 1000
                        elif len(p) == 4:
                            r[p] = 100
                        elif len(p) == 6:
                            r[p] = 10000
                    elif pert_id == 5:
                        if len(p) == 5:
                            r[p] = 100
                        elif len(p) == 4:
                            r[p] = 1000
                        elif len(p) == 6:
                            r[p] = 10
                    elif pert_id == 6:
                        if len(p) == 5:
                            r[p] = 10
                        elif len(p) == 4:
                            r[p] = 100
                        elif len(p) == 6:
                            r[p] = 1
                    model['pars'][p] = r[p]
                if 'B' in r.fs():
                    control_s = 'sB'
                else:
                    control_s = 'sA'
                r[control_s] = 0
                model['pars'][control_s] = r[control_s]
                uplim = 120
                data, bounds, boundsh = run_bf(r, auto, dirc="+", par=control_s, lims=[0,uplim],
                    ds=1E-2, dsmin=1E-5, dsmax=0.1)
                if len(boundsh) > 0 or len(bounds) > 3:
                    if len(boundsh) > 0:
                        print('HB point found')
                        n.num_bifs[-1]['HB'] += 1
                        model['HB_index'] = boundsh
                    if len(bounds) > 3:
                        #print('LP point found')
                        n.num_bifs[-1]['LP'] += 1
                    nameset = str(outputpath)+'/'+n.pattern+'_P'+str(ip)+'_'+str(i)
                    data.to_csv(nameset+'.csv')
                    a_file = open(nameset+"_antimony.txt", "wt")
                    a_file.write(r.getCurrentAntimony())
                    a_file.close()
                    m_file = open(nameset+"_modeldict", "wb")
                    pickle.dump(model, m_file)
                    m_file.close()
    for n in nodes[:]:
        for ip, path in enumerate(n.paths):
            print(n.pattern, ip)
            try:
                print(n.num_bifs[ip])
            except IndexError:
                print({'LP':0, 'HB':0})

perturbs = {0:'Basal range',
            1: 'Making all association constants = 10',
            2: 'Making all association constants = 100',
            3: 'Making all association constants 1000',
            4: 'Positive binding cooperativity 100, 1000, 10000',
            5: 'Negative binding cooperativity, 1000, 100, 10',
            6: 'Negative binding cooperativity 100, 10, 1',
            7: 'Making every degradation rate constant of free species 1',
            8: 'Making degradation rate constants of A, C, D equal',
            9: 'Making all regulated degradation factors <= 1',
            10: 'Making all regulated degradation factors >= 1',
            11: 'Making regulated degradation factors for A, C, D equal'
            }


if __name__ == '__main__':
    nodes = reaction_nets.make_nodes()
    for pert_id in range(len(perturbs)):
        print(pert_id, perturbs[pert_id])
        run_cont(nodes, nreps=100, if_only_seeds=True, pert_id=pert_id) # A quick test run
        #run_cont(nodes, nreps=10000, if_only_seeds=False, pert_id=pert_id) # The paper's full setup

