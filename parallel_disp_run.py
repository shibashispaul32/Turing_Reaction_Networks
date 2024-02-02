# ReadMe 'run.py'
#--------------------
# data_hb_lp contains the results from numerical continuation using tellurium and includes parameter sets exhibiting Hopf(hb) and Saddle node(lp) bifurcations
# with this run.py python script we run the dispersion codes using the Hopf bifurcation parameter (hb) sets residing inside 'data_hb_lp' directory in parallel
#-----------------------------
#-----------------------------
#-----------------------------
n_proc=15 # total number of processors used for the parallel run
#-----------------------------
#-----------------------------
#-----------------------------
import numpy as np
import pandas as pd
import os
import re
#-----------------------------
file_list=list(set([re.match(r'[A-D]+_P\d_\d+',s).group()\
                         for s in os.listdir('data_hb_lp')])) # lists the ODE-ID of parameter sets inside data_hb_lp directory
if 'data_hb_all_paths' not in os.listdir():
    os.mkdir('data_hb_all_paths')
for i in file_list:
    df=pd.read_csv(f'data_hb_lp/{i}.csv')
    if 3 in df['TY'].to_list():
        os.system(f'cp data_hb_lp/{i}.* data_hb_lp/{i}_* data_hb_all_paths/.')# filtering the hb parameter sets and copying inside data_hb_all_paths directory
#-----------------------------
file_list=list(set([re.match(r'[A-D]+_P\d_\d+',s).group()\
                         for s in os.listdir('data_hb_all_paths')]))# lists the ODE-ID of hb parameter sets
paths=sorted(sorted(list(set([re.match(r'[A-D]+_P\d',s).group()\
                         for s in os.listdir('data_hb_all_paths')]))),key=len)# lists the unique biochemical networks
if 'file_list.txt' in os.listdir():
    os.system('rm file_list.txt')
for i0 in file_list:
    print(i0,file=open('file_list.txt','a'))# write the ODE-ID s in file_list.txt file for future reference
#-----------------------------
a=np.char.array(file_list)# creating numpy string array from file_list python list
chunks=np.array_split(a,n_proc)# splitting the array into n_proc number of chunks
if 'parallel_runs' not in os.listdir():
    os.system('mkdir parallel_runs')
#-----------------------------
# Executing the codes in n_proc different processors with n_proc different chunks by creating and copying required files to n_proc different directories inside a runs directory
#-----------------------------
for i in range(len(chunks)):
    if 'run-'+str(i) not in os.listdir('parallel_runs'):
        os.system(f'mkdir parallel_runs/run-{i}')
        os.system(f'mkdir parallel_runs/run-{i}/data_hb_all_paths')
        os.system(f'mkdir parallel_runs/run-{i}/data_disp_all_path')
    for k in chunks[i]:
        os.system(f'cp data_hb_all_paths/{k}* parallel_runs/run-{i}/data_hb_all_paths/.')
    os.system(f'cp script_scan_disp_all.py parallel_runs/run-{i}/.')
    path_mem=os.getcwd()
    os.chdir(f'parallel_runs/run-{i}/')
    os.system('python script_scan_disp_all.py &') # executing the python script 'script_scan_disp_all.py'
    os.chdir(path_mem)

