#--------------------------------------------------------------
#--------------------------------------------------------------
n_proc=20
compiler='ifortran'# 'gfortran' will be considerably slower
#--------------------------------------------------------------
#--------------------------------------------------------------
import pandas as pd
import numpy as np
import os
import re
import pickle
import time
import matplotlib.pyplot as plt
#--------------------------------------
pde_ids=['AAB_P1_693_D_10482',
 'ABC_P0_782_D_926',
 'AAAB_P0_4_D_224',
 'AAAB_P2_942_D_5661',
 'AAAB_P3_737_D_17277',
 'AABB_P2_641_D_13420',
 'AABB_P3_121_D_64',
 'AABC_P2_417_D_5724',
 'AABC_P3_809_D_5018',
 'ABCD_P1_04843_D_16501']
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
def simulate_pde_id(pde_id,\
    #--------------------------------------------------------------
    #                  simulation parameters
    #--------------------------------------------------------------
    d_t=0.001,
    n_time=100000,
    d_x=0.5,n_x=100,
    d_y=0.5,n_y=100,
    n_stim=20,stim_amp=0.01,spread_stim=0.5,
    # boundary_cond='periodic'
    boundary_cond='zero flux',
    grid_size=None,run_sim=1,compiler='ifortran'
             ):
    '''The function simulate_pde_id gets input of pde-id and simulation parameters
    and generate a fortran code for simulation based on explicit forward Euler algorithm and executes the code and saves all the output files insise a directory named as the given pde-id.
    important:for this function to work properly, intel fortran should be installed in the working environment
    returns: simulated pde_id (string) and time gap between snapshots (float)'''
    if grid_size != None:
        d_x=grid_size
        d_y=grid_size
    n_wri=int(n_time/10)
    ode_id_current=re.search(r'[A-D]+_P\d_\d+',pde_id).group()
    model=pickle.load(open('data_hb_all_paths/'+ode_id_current+'_modeldict','rb'))
    try:os.system('rm -f test.* fort.*')
    except:pass
    #--------------------------------------------------------------
    code_string=''
    code_string+=f'implicit none\n!--------------------------------------\n'
    #--------------------------------------------------------------
    #                      variable declaration
    #--------------------------------------------------------------
    free_pars={}
    df=pd.read_csv("disp_data_compiled/disp_"\
                   +re.search(r'[A-D]+_P\d',ode_id_current).group()+".csv",index_col=0)
    vall=df[df['pde_id']==pde_id].values[0]
    vall=vall[2:]
    for i0,i1 in enumerate(df.columns[2:]):
        if i0>0 and i0<1+(df.shape[-1]-3)/2:
            free_pars[i1+'_0']=vall[i0]
        elif i0>=1+(df.shape[-1]-3)/2:
            free_pars['D'+i1]=vall[i0]
        else:
            free_pars[i1]=vall[i0]
    #--------------------------------------------------------------
    for k in model['vars'].keys():
        code_string+=f'double precision,allocatable,dimension(:,:)::{k},{k}_m\n'
        code_string+=f'double precision::lap_{k},r_{k}\n'
    code_string+=f'double precision,allocatable,dimension(:,:)::r\n'    
    for k,v in model['pars'].items():
        if k!='sB':
            code_string+=f'double precision::{k}={v}\n'
    for k,v in free_pars.items():
        code_string+=f'double precision::{k}={v}\n'
    code_string+=f'double precision::d_t={d_t},d_x={d_x},d_y={d_y}\n'    
    code_string+=f'double precision::stim_amp={stim_amp},spread_stim={spread_stim},r_x,r_y\n'    
    code_string+=f'integer::i_stim,n_stim={n_stim}\n'    
    code_string+=f'integer::n_time={n_time},n_x={n_x},n_y={n_y}\n'    
    code_string+=f'integer::i,j,k_time,l,n_wri={n_wri}\n'    
    code_string+=f'!--------------------------------------\n'    
    #--------------------------------------------------------------
    #                    memory allocation
    #--------------------------------------------------------------
    for k in model['vars'].keys():
        code_string+=f'allocate({k}(n_x,n_y),{k}_m(n_x,n_y))\n'
    code_string+=f'allocate(r(n_x,n_y))! random array of size (n_x,n_y)\n'
    code_string+=f'!--------------------------------------\n'    
    #--------------------------------------------------------------
    #                   assign initial conditions
    #--------------------------------------------------------------
    code_string+=f'call random_seed\n'
    for k in model['vars'].keys():
        code_string+=f'{k}={k}_0\n'
    code_string+=f'do i_stim=1,n_stim\n'
    code_string+=f'call random_number(harvest=r_x)\n'
    code_string+=f'call random_number(harvest=r_y)\n'
    code_string+=f'r_x=(0.1*n_x*d_x)+(r_x*(d_x*(n_x-(n_x*0.2))))\n'
    code_string+=f'r_y=(0.1*n_y*d_y)+(r_y*(d_y*(n_y-(n_y*0.2))))\n'
    
    code_string+=f'do i=1,n_x\ndo j=1,n_y\n'
    for k in model['vars'].keys():
        code_string+=f'{k}(i,j)={k}(i,j)+stim_amp*exp(-(((i*d_x)-r_x)**2+((j*d_y)-r_y)**2)/(2*spread_stim))\n'
    code_string+=f'end do\nend do\nend do\n'
    code_string+=f'!--------------------------------------\n'
    #--------------------------------------------------------------
    #                   assign boundary condition
    #--------------------------------------------------------------
    if (boundary_cond[0].upper())=='Z':
        for k in model['vars'].keys():
            code_string+=f'{k}(1,:)={k}(2,:);{k}(n_x,:)={k}((n_x-1),:)\n'
            code_string+=f'{k}(:,1)={k}(:,2);{k}(:,n_y)={k}(:,n_y-1)\n'
    if (boundary_cond[0].upper())=='P':
        for k in model['vars'].keys():
            code_string+=f'{k}(1,:)={k}((n_x-1),:);{k}(n_x,:)={k}(2,:)\n'
            code_string+=f'{k}(:,1)={k}(:,(n_y-1));{k}(:,n_y)={k}(:,2)\n'
    code_string+=f'!--------------------------------------\n'
    #--------------------------------------------------------------
    #                        pde integration
    #--------------------------------------------------------------
    code_string+=f'do k_time=0,n_time\n'# time loop
    code_string+=f'!--------------------------------------\n'
    #--------------------------------------------------------------
    #                       writing files
    #--------------------------------------------------------------
    for k in model['vars'].keys():
        code_string+=f'{k}_m={k}\n'
    #-----------
    code_string+=f'if(mod(k_time,n_wri).eq.0)then\n'
    multi=2+len(model['vars'].keys())
    code_string+=f"do i=1,n_x\ndo j=1,n_x\
        \nwrite(200+(k_time/n_wri),'({multi}(2x,f12.5))')(i*d_x),(j*d_y),"
    dummy_string=''
    for i,k in enumerate(sorted(list(model['vars'].keys()),key=len)):
        if i<len(model['vars'].keys())-1:
            code_string+=f"{k}_m(i,j),"
            dummy_string+=f"{k}_m(i,j),"
        else:
            code_string+=f"{k}_m(i,j)"
            dummy_string+=f"{k}_m(i,j)"
    code_string+=f'\nend do\n'
    code_string+=f"write(200+(k_time/n_wri),*)\n"
    code_string+=f'end do\n'
    code_string+=f'end if\n'
    code_string+=f'!--------------------------------------\n'
#--------------------------------------------------------------
    code_string+=f'do i=2,n_x-1\n'
    code_string+=f'do j=2,n_y-1\n'
    for k in model['vars'].keys():
        code_string+=f'lap_{k}=({k}_m(i+1,j)+{k}_m(i-1,j)-(2*{k}_m(i,j)))/(d_x**2)\n'
        code_string+=f'lap_{k}=lap_{k}+({k}_m(i,j+1)+{k}_m(i,j-1)-(2*{k}_m(i,j)))/(d_y**2)\n'
        code_string+=f'lap_{k}=lap_{k}*DD_{k}\n'
        code_string+=f'\n'
    for k,v in model['vars'].items(): 
        code_string+=f"r_{k}="+re.sub(r'(\*[A-Z]+)',r'\1_m(i,j)',v)+"\n"
        code_string+=f'\n'
    for k in model['vars'].keys():
        code_string+=f'{k}(i,j)={k}_m(i,j)+(d_t*(r_{k}+lap_{k}))\n'
        code_string+=f'\n'
    code_string+=f'\n'
    code_string+=f'end do\n'
    code_string+=f'end do\n'
    #--------------------------------------------------------------
    #                   assign boundary condition
    #--------------------------------------------------------------
    if (boundary_cond[0].upper())=='Z':
        for k in model['vars'].keys():
            code_string+=f'{k}(1,:)={k}(2,:);{k}(n_x,:)={k}((n_x-1),:)\n'
            code_string+=f'{k}(:,1)={k}(:,2);{k}(:,n_y)={k}(:,n_y-1)\n'
    if (boundary_cond[0].upper())=='P':
        for k in model['vars'].keys():
            code_string+=f'{k}(1,:)={k}((n_x-1),:);{k}(n_x,:)={k}(2,:)\n'
            code_string+=f'{k}(:,1)={k}(:,n_y-1);{k}(:,n_y)={k}(:,2)\n'
    code_string+=f'!--------------------------------------\n'
    code_string+=f'end do\n'# time-loop
    #--------------------------------------------------------------
    code_string+=f'!--------------------------------------\nend\n'
    lines=(code_string.split('\n'))
    len_lines=[len(i) for i in lines]
    n=130
    for i,line in enumerate(lines):
        if len(line) >n :
            for kk in range(int(len(line)/n)):
                line=line[:(((kk+1)*n)+(2*kk))]+'&\n&'+line[(((kk+1)*n)+(2*kk)):]
                lines[i]=line
    code_string="\n".join(lines)
    #--------------------------------------------------------------
    with open("variables.txt","a") as file:
        file.write(dummy_string)
    with open("test.f90","a") as file:
        file.write(code_string)
    if run_sim:  
        if d_x<=1:
            res="high_res"
        else:
            res="low_res"
        os.system(f'rm -rf {pde_id}_{res}')
        os.mkdir(f"{pde_id}_{res}")
        if compiler[0].lower()=='i':
            os.system('rm -f fort.* *.out && ifort test.f90 -o test.out')
        elif compiler[0].lower()=='g':
            os.system('rm -f fort.* *.out && gfortran test.f90 -o test.out')
        os.system(f'mv test.f90 test.out variables.txt {pde_id}_{res}/.')
        os.chdir(f'{pde_id}_{res}')
        os.system('./test.out & \n echo $! > pid.log')
        os.chdir('../')
    return pde_id,(d_t*n_time)/10
#-----------------------
#-----------------------
#-----------------------
def plot_pattern(input_file,pde_id,t_wri,res,set_title=0,ticks_labels=1,\
                 n_cont=20,z_col=3,c_map='RdGy_r',dpi=300,v_lim=False):
    if res[0].lower()=='h':
        res="_high_res"
    elif res[0].lower()=='l':
        res="_low_res"
    else:
        print('ildefined resolution')
        return
    if input_file not in os.listdir(pde_id+res):
        print('file does not exist yet!!')
        return
    df=pd.read_csv(f'{pde_id+res}/{input_file}',sep=r'\s+',header=None)
    dim=int(np.sqrt(len(df)))
    x=df.iloc[:,0].values.reshape(dim,dim)
    y=df.iloc[:,1].values.reshape(dim,dim)
    z=df.iloc[:,z_col].values.reshape(dim,dim)
    fig0,fig=plt.subplots()
    if v_lim==False:
        cs=fig.contourf(x,y,z,n_cont,cmap=c_map)
        cbar=fig0.colorbar(cs)
    else:
        cs=fig.contourf(x,y,z,n_cont,cmap=c_map,levels=np.linspace(v_lim[0],v_lim[1],400))
        cbar=fig0.colorbar(cs,ticks=np.linspace(v_lim[0],v_lim[1],6))
    if ticks_labels:
        fig.set_xlabel('X',size='11')
        fig.set_ylabel('Y',size='11')
    else:
        fig.set_xticks([])
        fig.set_yticks([])
    if set_title:
        fig.set_title(f'{pde_id} at t={t_wri*(int(input_file[-3:])-200):.1f} time units')
    os.chdir(pde_id+res)
    plt.savefig(pde_id+res+'.png',dpi=dpi)
    plt.close()
    os.chdir('../')
#--------------------------------------------------------------
#--------------------------------------------------------------
pde_ids_sim=[]
for grid_size in [1.0,5.0]:
    for pde_id in pde_ids[:]:
        pde_id_sim,t_wri=simulate_pde_id(pde_id=pde_id#'pde_id_goes_here'
                                ,grid_size=grid_size,\
                                d_t=10e-5,n_time=5000000,\
                                spread_stim=30.0,n_stim=10,run_sim=1,compiler=compiler)
        if grid_size <=1:
            res='_high_res'
        else:
            res='_low_res'
        pde_ids_sim.append(pde_id_sim+res)
        if np.mod(len(pde_ids_sim),n_proc) == 0:
            while 'fort.210' not in os.listdir(pde_ids_sim[-1]):
                time.sleep(0.1)
#-------------------------------------------------------                
while 'fort.210' not in os.listdir(pde_ids_sim[-1]):
    time.sleep(0.1)
#-------------------------------------------------------                
for i in pde_ids_sim:
    pde_id=re.search(r'[ABCD]+_P\d_\d+_D_\d+',i).group()
    res=re.search(r'high|low',i).group()
    plot_pattern(input_file='fort.210',pde_id=pde_id,res=res,t_wri=t_wri)