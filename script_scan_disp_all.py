import time;t0=time.time()
import tellurium as te
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from rrplugins import *
auto = Plugin("tel_auto2000")
import sobol
import sympy as sym
import os
import pickle
import re
#----------------------------------------------------------------------------------------------
file_list=[]
for i in os.listdir('data_hb_all_paths'):
    x=re.findall(r'[A-Z]+_P\d_\d+',i)
    if len(x)>0:
        file_list.append(x[0])
        del x
file_list=sorted(list(set(file_list)))
os.system('rm file_list.txt data_disp_all_path/* temp_scan_result.csv -f')
for i0 in file_list:
    print(i0,file=open('file_list.txt','a'))
#----------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------
del_files=False
#----------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------
if del_files==True:
    try:
       os.remove('scan_result.csv')
    except:
       pass
    for i in os.listdir('data_disp_all_path'):
       os.remove('data_disp_all_path/'+i)
#----------------------------------------------------------------------------------------------
model_info=[];path_info=[];file_index=[];turing=[];turing_accept=[];turing_count=[];turing_a_count=[]
#----------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------
for file_id in file_list[:]:
#----------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------
    model_info.append(re.search(r'[A-D]+',file_id).group())
    path_info.append(re.search(r'P\d',file_id).group())
    file_index.append(re.search(r'_\d+',file_id).group()[1:])
#----------------------------------------------------------------------------------------------
    s='data_hb_all_paths/'+file_id
    model=pickle.load(open(s+'_modeldict','rb'))
    with open(s+'_antimony.txt', 'r') as f:
         antstr = f.read()
    r=te.loada(antstr)  
    ics=list(r.getFloatingSpeciesConcentrations())  
    data=pd.read_csv(s+'.csv', index_col=0) 
    data_0=data.copy()
    #--------------------------------------
    if (2 in set(data['TY'].values)):
         data['PT'][data['TY']==2]=abs(data['PT'][data['TY']==2].values[0])
         model['LP_index']=list(data.index[data['TY']==2])
        #  print('lp:\n',model['LP_index'])
    if (3 in set(data['TY'].values)):
         data['PT'][data['TY']==3]=abs(data['PT'][data['TY']==3].values[0])
         model['HB_index']=list(data.index[data['TY']==3])
        #  print('hb:\n',model['HB_index'])
    model['BIF_index']=list(data[(data.TY==2) | (data.TY==3)].index)
    model['BIF_PAR']=list(data[(data.TY==2) | (data.TY==3)].PAR)
    model['BIF_PAR'], model['BIF_index'] = zip(*sorted(zip(model['BIF_PAR'], model['BIF_index'])))
    #--------------------------------------
    if (2 not in set(data['TY'].values)) & (3 in set(data['TY'].values)):
         if (data['TY'].value_counts()[3]==1):
              model['HB_index']=[data[(data['PT']>0) & (data['TY']!=9)].index[0],data[(data['PT']>0) & (data['TY']!=9)].index[-1]]
    #--------------------------------------          

    if (2 in set(data['TY'].values)) & (3 in set(data['TY'].values)):
         if (data['TY'].value_counts()[3]>=1):
              bounds=[model['BIF_index'][0],model['BIF_index'][-1]]
              drop_indx=data[(data.PAR>data.PAR[bounds[0]]) & (data.PAR<data.PAR[bounds[-1]]) & (data.PT<0)].index
              data.drop(drop_indx , inplace=True)
              data.sort_values('PAR',inplace=True)
              data.reset_index(drop=True,inplace=True)
    #----------------------------------------------------------------------------------------------
    #---------------------------------
    #     Symbolic Calculation
    #---------------------------------
    def Jacobian(v_str, f_list):
        vars = sym.symbols(v_str)
        f = sym.sympify(f_list)
        f1=sym.Matrix(f_list)
        J = sym.zeros(len(f),len(vars))
        for i, fi in enumerate(f):
            for j, s in enumerate(vars):
                J[i,j] = sym.diff(fi, s)
        return J,f1
    #---------------------------------
    var_str=''
    funcs=[]
    vars=[]
    for k,v in model['vars'].items():
        var_str=var_str+str(k)+' '
        funcs.append(v)
        vars.append(k)
    #---------------------------------
    pars=[];par_val=[]
    for k,v in model['pars'].items():
        if k not in ['sB']:
            pars.append(k)
            par_val.append(v)
    #---------------------------------
    #    lambdify Jacobian
    #---------------------------------
    j_sp,simul_f_0=Jacobian(var_str, funcs)
    f_sp=sym.lambdify([tuple(pars)],j_sp,'sympy')
    j_sp_1=f_sp(tuple(par_val))
    f_np=sym.lambdify([tuple(['sB']+vars)],j_sp_1,'numpy')
    #----------------------------------------------------------------------------------------------
    sol_dict={}
    sol_dict['s_B']=data['PAR']
    for i in vars:
        sol_dict[i]=data[str(i)]
    data_df=pd.DataFrame(sol_dict)
    #----------------------------------------------------------------------------------------------
    def get_prediction(input_var,input_val,sol_dict
                          ,sobol_seq=True,sample_size=200                           # scanning parameters
                          ,k_max=5,delta_k=0.01                                     # dispersion parameters
                          ,scan_min=1.0,scan_max=10.0                               # scanning parameters
                          ,aa=0,if_plot=False,print_out=False
                          ):
            '''
            returns sigma_R (nearest available to the input),
                    prediction (0= no, 1= yes) '''
            #----------------------------------------------------------------------------------
            #  define steady states using tellurium results loaded in sol_dict dictionary
            #----------------------------------------------------------------------------------  
            disp_ics=[]                  
            index=np.where(abs(sol_dict[input_var]-input_val)==min(abs(sol_dict[input_var]-input_val)))[0][0]
            for key in sol_dict:
                    disp_ics.append(sol_dict[key][index])
            #--------------------------------------------------------------
            #   Defining free diffusive parameters with sobol seq
            #--------------------------------------------------------------
            if sobol_seq:
                    D_free=sobol.sample(dimension=len(disp_ics)-1, n_points=sample_size)
            else:
                    D_free=np.random.rand(sample_size,len(disp_ics)-1)
            D_free=scan_min+((scan_max-scan_min)*D_free)
            #--------------------------------------------------------------
            ord_pm=[];pattern_par=[]
            for i in range(len(D_free)):
                    D0=np.diag(D_free[i])
                    max_re_eig=[]
                    ks=np.arange(0,k_max,delta_k)
                    for k in ks:
                            D=D0*k**2
                    #---------------------------
                    # Defining Jacobian
                    #---------------------------
                            jac=f_np(tuple(disp_ics))
                    #---------------------------------------
                    #   scanning the dispersion curve
                    #---------------------------------------
                            max_re_eig.append(max(np.real(np.linalg.eig(jac-D)[0])))
                    k_rising=0;k_lowering=0
                    if max_re_eig[0]>0:
                            a_0=max_re_eig[0]
                            c_0=(max(max_re_eig)-max_re_eig[0])
                            max_re_eig=max_re_eig-a_0-(0.5*c_0)
                    for j in range(len(max_re_eig)-1):
                            if (max_re_eig[j]*max_re_eig[j+1]) < 0:
                                    if (max_re_eig[j+1]-max_re_eig[j]) > 0:
                                            k_rising=(ks[j]+ks[j+1])/2
                                    if (max_re_eig[j+1]-max_re_eig[j]) < 0:
                                            k_lowering=(ks[j]+ks[j+1])/2
                    if(k_rising>0):
                            if k_rising<k_lowering:
                                    ind=1
                                    pattern_par.append(disp_ics+list(D_free[i]))
                            else:
                                    ind=2
                    if (k_rising==0):
                            ind=3
                    ord_pm.append(ind)
                    if max_re_eig[0]>0:
                            max_re_eig=max_re_eig+a_0+(0.5*c_0)
                    if ind==1 and if_plot==True:
                            fig=plt.subplot()
                            if aa==0:
                                    fig.plot(ks,np.zeros(len(ks)),ls='dotted')
                                    aa=1
                            fig.plot(ks,max_re_eig)
                            # print(D_free[i])
            if len(set(ord_pm))==1:
                   pred_t=0
            else:
                   pred_t=1
            if print_out:
                    print(disp_ics[0:2],pred_t)
            return disp_ics[0],disp_ics[1],pred_t,pattern_par   
    #----------------------------------------------------------------------------------------------
    ratio=0.2  #changed the value from 0.5 to 0.2 to decrease the range
    p1=model['BIF_PAR'][0]-(ratio*model['BIF_PAR'][0])
    p2=model['BIF_PAR'][-1]+(ratio*model['BIF_PAR'][0])
    sbs=np.linspace(p1,p2,20)
    predictions=[]
    sb_map=[]
    par_map=[]
    par_tur=[]
    for sb_in in sbs:
        sb_out,par_out,pred,tpars=get_prediction('s_B',input_val=sb_in,sol_dict=sol_dict,k_max=5.0,delta_k=0.05,scan_max=20.0,sobol_seq=True)
        predictions.append(pred)
        sb_map.append(sb_out)
        par_map.append(par_out)
        if len(tpars)>0:
            par_tur=par_tur+tpars
    #----------------------------------------------------------------------------------------------
    fig1,axs=plt.subplots(1,1,figsize=(6,5))
    axs.set_ylabel('A',size='20')
    axs.set_xlabel('$\sigma_B$',size='20')
    axs.plot(data_0["PAR"],data_0["A"],c='k',linewidth='0.8',ls='solid')
    axs.autoscale()
    r2=1
    axs.set_xlim(p1-(r2*p1),p2+(r2*p1))
    axs.scatter(sb_map,par_map,c=predictions,cmap='bwr',s=20.0)
    if (2 in set(data_0['TY'].values)):
        axs.scatter(data_0["PAR"].iloc[data_0.index[data_0['TY']==2]],data_0["A"].iloc[data_0.index[data_0['TY']==2]],\
                    c='green',marker="D",label='LP',s=50,edgecolors='k',alpha=0.5)
    if (3 in set(data_0['TY'].values)):
        axs.scatter(data_0["PAR"].iloc[data_0.index[data_0['TY']==3]],data_0["A"].iloc[data_0.index[data_0['TY']==3]],\
                    c='yellow',edgecolors='black',marker="D",label='HB',s=50,alpha=0.5)
    fig1.suptitle(s[18:]+'\n red= Turing+                          blue= Turing-')
    axs.legend()
    #-----------------------------
    # axs.set_xscale('log')
    # axs.set_yscale('log')
    #-----------------------------
    fig1.tight_layout()
    plt.savefig('data_disp_all_path/'+s[18:]+'.png')
    plt.close()
    #----------------------------------------------------------------------------------------------
    col=['sB']+list(model['vars'].keys())
    for i in model['vars'].keys():
        col.append('D_'+i)
    if len(par_tur)>0:
    #-----------------------------  
        turing.append(1)
        turing_count.append(len(par_tur))
    #-----------------------------  
        df_par_0=pd.DataFrame(par_tur,columns=col)
        lll=(df_par_0.columns.tolist())[1:]
        lll=sorted(lll)
        lll.insert(0,'sB')    
        lll_new=['sB']+sorted(lll[1:int(1+(len(lll)-1)/2)],key=len)+sorted(lll[int(1+(len(lll)-1)/2):len(lll)],key=len)
        df_par=df_par_0.reindex(lll_new,axis=1)
        df_diff=df_par.iloc[:,int(1+((df_par.shape[1]-1)/2)):].copy()
    #-----------------------------  
        def list_avg(a):
            if len(a)==0:
                return 0
            else:
                return np.average(a)
    #-----------------------------  
        accept=[]
        for i in df_diff.index:
            one_comp=[]
            two_comp=[]
            three_comp=[]
            four_comp=[]
            for j in df_diff.columns:
                if len(j)==3:
                    one_comp.append(df_diff[j][i])
                if len(j)==4:
                    two_comp.append(df_diff[j][i])
                if len(j)==5:
                    three_comp.append(df_diff[j][i])
                if len(j)==6:
                    four_comp.append(df_diff[j][i])
            tempo=[list_avg(one_comp),list_avg(two_comp),list_avg(three_comp),list_avg(four_comp)]
            tempo = [i for i in tempo if i != 0]  
            difference=np.array(tempo[:-1])-np.array(tempo[1:])
            if np.all(difference>0):
                accept.append(1)
            else:
                accept.append(0)
    #-----------------------------  
        df_diff['accept']=accept
        df_par['accept']=accept
        df_par.to_csv('data_disp_all_path/'+s[18:]+'.csv')
        if 1 in set(df_diff['accept']):
            turing_accept.append(1)
            turing_a_count.append(df_diff['accept'].sum())
        else:
            turing_accept.append(0)
            turing_a_count.append(0)
    else:
        turing.append(0)
        turing_accept.append(0)
        turing_count.append(0)
        turing_a_count.append(0)
    print(f"{model_info[-1]}, {path_info[-1]}, {file_index[-1]}, \
          {turing[-1]:d}, {turing_count[-1]:d}, {turing_accept[-1]:d},\
           {turing_a_count[-1]:d}",file=open('temp_scan_result.csv','a'))
#----------------------------------------------------------------------------------------------
df_scan=pd.DataFrame(list(zip(file_list,model_info,path_info,file_index,turing,turing_count,turing_accept,turing_a_count)), columns=['ode_id','model_info','path_info','file_index','turing','turing_count','turing_accept','turing_a_count'])
df_scan.to_csv('scan_result.csv')
