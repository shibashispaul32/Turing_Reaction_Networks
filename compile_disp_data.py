# Readme 'script-combine.py'
# after successful execution of the 'run.py' script run this script to combine all the fragments of n_proc runs
# the combined results are stored inside 'data_disp_all_path' directory
#-------------------------------------------------------------------
import os
import pandas as pd
import numpy as np
import os
import re
import pickle
import pprint
import matplotlib.pyplot as plt
plt.style.use('seaborn-v0_8-dark')
#-------------------------------------------------------------------
if 'data_disp_all_path' not in os.listdir():
    os.mkdir('data_disp_all_path')
    for i in os.listdir('parallel_runs'):
        os.system(f'cp parallel_runs/{i}/data_disp_all_path/* data_disp_all_path/.')
#-------------------------------------------------------------------
df=pd.DataFrame()
for i in os.listdir('parallel_runs'):
    df_temp=pd.read_csv(f'parallel_runs/{i}/scan_result.csv',index_col=0)
    df=pd.concat([df,df_temp],axis=0,ignore_index=1)
df=df.drop_duplicates()
df=df.reset_index()
df.drop(columns='index',inplace=True)
df.to_csv('scan_result.csv')
#-------------------------------------------------------------------
if ('scan_result_hb.csv' not in os.listdir()) & ('scan_result_lp.csv' not in os.listdir()):
    df_hb=df.copy()
    df_lp=df.copy()
    for i in sorted(df['ode_id'].to_list()):
        df_temp=pd.read_csv(f'data_hb_all_paths/{i}.csv',index_col=0)
        if 3 not in df_temp['TY'].to_list():
            drop_index=df[df['ode_id']==i].index
            df_hb.drop(drop_index,inplace=True)
        if 2 not in df_temp['TY'].to_list():
            drop_index=df[df['ode_id']==i].index
            df_lp.drop(drop_index,inplace=True)
    df_hb.reset_index(inplace=True)
    df_hb.drop(columns='index',inplace=True)
    df_lp.reset_index(inplace=True)
    df_lp.drop(columns='index',inplace=True)
    df_hb.to_csv('scan_result_hb.csv')
    df_lp.to_csv('scan_result_lp.csv')
#---------------------------------------------------------    
def get_stat():
    paths=sorted(sorted(list(set([re.match(r'[A-D]+_P\d',s).group()\
                             for s in os.listdir('data_hb_all_paths')]))),key=len)
    n_path=len(paths)
    df=pd.read_csv('scan_result_hb.csv',index_col=0)
    p_disp=((df['turing_count'].sum())/(df.shape[0]*200*20))*100
    s=f"percentage of turing enabling pattern for  data set4 : {p_disp:.4f} %\n"
    p_par=(df['turing'].sum()/df.shape[0])*100
    s+=f"--------\namongst {df.shape[0]} parameter sets that can produce Turing patterns:  {p_par} %\n--------\n"
    s+="overall pattern enabling dispersion sets (%): "+str(((df['turing_count'].sum())/(n_path*10000*4000))*100)
    f=open("turing_stat/disp_stat.txt","w")
    f.write(s)
#---------------------------------------------------------
def disp_fil_stat(s_filter=False):  
    if 'turing_stat' not in os.listdir():
        os.mkdir('turing_stat')
    scan_data=pd.read_csv('scan_result_hb.csv',index_col=0)
    scan_data['path']=scan_data['model_info']+'_'+scan_data['path_info']
    del scan_data['model_info'],scan_data['path_info']
    s="path,sampled_dispersion_pars,turing_pos_pars,turing_acceptable_pars\n"
    file_list=list(set([re.match(r'[A-D]+_P\d_\d+',s).group()\
                         for s in os.listdir('data_hb_all_paths')]))
    paths=sorted(sorted(list(set([re.match(r'[A-D]+_P\d',s).group()\
                         for s in os.listdir('data_hb_all_paths')]))),key=len)
    for i in paths:
        if s_filter==True:
            f=open('turing_stat/disp_fil_stat.csv','w')
            for j in file_list:
                if re.match(i,j):
                    df=pd.read_csv(f'data_hb_all_paths/{j}.csv',index_col=0)
                    if (df[df['TY']==3]['PAR'].min() > 30) | (df[df['TY']==3]['PAR'].max() <1):
                        scan_data.drop(scan_data[scan_data['ode_id']==j].index,inplace=True)
        else:
            f=open('turing_stat/disp_stat.csv','w')
        a=scan_data[scan_data['path'].str.contains(i)].shape[0]*4000
        b=scan_data[scan_data['path']==i]['turing_count'].sum()
        c=scan_data[scan_data['path']==i]['turing_a_count'].sum()
        s+=f"{i},{a},{b},{c}\n"
        f.write(s)
df1=disp_fil_stat(s_filter=False)
#---------------------------------------------------------
def bar_plot():
    plt.style.use('seaborn-v0_8-notebook')
    df_disp=pd.read_csv('turing_stat/disp_stat.csv')
    df_disp['path']=df_disp['path'].str.replace(r'_P',r'-',regex=True)
    df_disp['path']=df_disp['path'].str.replace(r'(\d)',lambda x: f"{str(int(x.group())+1)}",regex=True)
    #--------------------------------------
    fig,axs=plt.subplots(nrows=2,ncols=1)
    fig.suptitle('Dispersion analysis statistics')
    fig.set_figheight(10)
    fig.set_figwidth(7)
    #--------------------------------------
    df_disp.plot.bar(x='path',y=['sampled_dispersion_pars','turing_pos_pars'],rot=30,ax=axs[0])
    axs[0].legend(['Sampled Dispersion Parameters','Pattern enabling dispersion parameter sets'],loc=[0,1])
    axs[0].set_yscale('log')
    axs[0].set_ylabel('count')
    axs[0].set_xlabel('networks')
    axs[0].figure.tight_layout()
    #--------------------------------------
    df_disp.plot.bar(x='path',y=['turing_acceptable_pars'],color='indianred',rot=30,ax=axs[1],width=0.25)
    axs[1].legend(['Pattern-enabling Dispersion parameter sets with ordered diffusion coefficients'],loc=[0,1])
    axs[1].set_yscale('log')
    axs[1].set_ylabel('count')
    axs[1].set_xlabel('networks')
    axs[1].figure.tight_layout()
    plt.close()
    fig.savefig('turing_stat/dispersion-stat.png',dpi=300)
bar_plot()
get_stat()
#------------------------------------------------    
#------------------------------------------------    
## Compile Dispersion DATA
#This section compiles all the dispersion parameter sets according to the network and assign each parameter set an unique ID, i.e., PDE_ID for future reference, the result is saved inside 
#------------------------------------------------ 
# Creating ODE_ID list and list of unique paths 
#------------------------------------------------  
file_list=pd.read_csv('file_list.txt',header=None).iloc[:,0].to_list()
paths=[]
for i in file_list:
    temp=re.search(r'[A-D]+_P\d',i).group()
    if temp not in paths:
        paths.append(temp)
    del temp
paths=sorted(paths,key=lambda x:(len(re.search(r'[A-D]+',x).group()),re.search(r'[A-D]+_P\d',x).group()))
file_list=sorted(file_list,key=lambda x:(len(re.search(r'[A-D]+',x).group()),re.search(r'[A-D]+_P\d',x).group(),int(re.search(r'_\d+',x).group()[1:])))  
#------------------------------------------------ 
#------------------------------------------------ 
scan_data=pd.read_csv('scan_result.csv',index_col=0)
try:
    scan_data.drop(['turing_accept','turing_a_count'],axis=1,inplace=True)
except:
    pass
#------------------------------------------------ 
#------------------------------------------------ 
def max_val(a):
    if len(a)==0:
        return 0
    else:
        return np.max(a)
#----------------------------- 
def min_val(a):
    if len(a)==0:
        return 0
    else:
        return np.min(a)
#-----------------------------
def perf_avg(a):
    if len(a)==0:
        return 0
    else:
        return np.mean(a) 
#------------------------------------------------  
def filter_disp_scan():
    for dec_variable in ['average','stringent']: # set 'average' for using average as the condition and 'stringent' for the stringent condition
        turing_accept=[];turing_accept_count=[]
        if dec_variable not in os.listdir('data_disp_all_path'):
            os.mkdir('data_disp_all_path/'+dec_variable)
    #----------------------------- 
        for k in scan_data.index:
            index=(scan_data.loc[k,'ode_id'])
            if (scan_data.loc[k,'turing'])==0:
                turing_accept.append(0)
                turing_accept_count.append(0)
            else:
                df_temp=pd.read_csv('data_disp_all_path/'+index+'.csv',index_col=0)
                df_temp.drop('accept',axis=1,inplace=True)       
                df_diff=df_temp.iloc[:,1+int((df_temp.shape[-1]-1)/2):].copy()
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
                    tempo_max_0=[max_val(one_comp),max_val(two_comp),max_val(three_comp),max_val(four_comp)]
                    tempo_min_0=[min_val(one_comp),min_val(two_comp),min_val(three_comp),min_val(four_comp)]
                    avg_list_0=[perf_avg(one_comp),perf_avg(two_comp),perf_avg(three_comp),perf_avg(four_comp)]
                    tempo_max = [i for i in tempo_max_0 if i != 0]  
                    tempo_min = [i for i in tempo_min_0 if i != 0]  
                    avg_list = [i for i in avg_list_0 if i != 0]
                    decision=[]
                    if dec_variable=='average':
                        for low,hi in zip(avg_list[:-1],avg_list[1:]):
                            if low>=hi:decision.append(1)
                            else:decision.append(0)
                    elif dec_variable=='stringent':
                        for low,hi in zip(tempo_min[:-1],tempo_max[1:]):
                            if low>=hi:decision.append(1)
                            else:decision.append(0)
                    if 0 not in decision:
                        accept.append(1)
                    else:
                        accept.append(0)
            #----------------------------- 
                df_temp.insert(loc=1,value=accept,column='accept')
                df_temp.to_csv('data_disp_all_path/'+dec_variable+'/'+index+'.csv')
                if 1 in accept:
                    turing_accept.append(1)
                    turing_accept_count.append(np.sum(accept))
                else:
                    turing_accept.append(0)
                    turing_accept_count.append(0)
        scan_data['turing_accept']=turing_accept
        scan_data['turing_accept_count']=turing_accept_count
        scan_data.to_csv('scan_result_'+dec_variable+'.csv')
filter_disp_scan()
#------------------------------------------------  
#------------------------------------------------  
def data_compile():
    file_list=pd.read_csv('file_list.txt',header=None).iloc[:,0].to_list()
    file_list=sorted(file_list,key=lambda x:(len(re.search(r'[A-D]+',x).group()),re.search(r'[A-D]+_P\d',x).group(),int(re.search(r'_\d+',x).group()[1:])))
    if 'disp_data_compiled' in os.listdir():
        os.system('rm disp_data_compiled/*')
    for filter_var in ['average','stringent']:
        paths=sorted(sorted(list(set([re.match(r'[A-D]+_P\d',s).group()\
                         for s in os.listdir('data_hb_all_paths')]))),key=len)
        scan_data=pd.read_csv('scan_result.csv',index_col=0)
        scan_data.drop(['turing_accept','turing_a_count'],axis=1,inplace=True)
        for i in paths:
            ode_id=[];pde_id=[];l=1
            for k in file_list:
                if re.search(i,k) and scan_data[scan_data['ode_id']==k]['turing'].values[0]==1:
                        df_read=pd.read_csv(f'data_disp_all_path/{filter_var}/'+k+'.csv')
                        #--------------------------------
                        if l==1:
                            df_cat=df_read.copy()
                        else:
                            df_cat=pd.concat([df_mem,df_read],ignore_index=1)
                        df_mem=df_cat.copy()
                        for m in df_read.index:
                            ode_id.append(k)
                            pde_id.append(k+'_D_'+str(l))
                            l+=1
            df_cat.drop(df_cat.columns[[0]], axis=1, inplace=True)
            df_cat.rename(columns={'accept':'accept_'+filter_var},inplace=True)
            df_cat.insert(loc=0,value=ode_id,column='ode_id')
            df_cat.insert(loc=1,value=pde_id,column='pde_id')
            if 'disp_data_compiled' not in os.listdir():os.mkdir('disp_data_compiled')
            df_cat.to_csv('disp_data_compiled/dispersion_'+filter_var+'_'+i+'.csv')
    for i in paths:
        df1=pd.read_csv(f'disp_data_compiled/dispersion_average_{i}.csv',index_col=0)
        df2=pd.read_csv(f'disp_data_compiled/dispersion_stringent_{i}.csv',index_col=0)
        df0=df2[df2.columns[~df2.columns.str.contains(r'accept')]]
        df1=df1[df1.columns[df1.columns.str.contains(r'accept_')]]
        df2=df2[df2.columns[df2.columns.str.contains(r'accept')]]
        df_comp=pd.concat([df0,df1,df2],axis=1)
        df_comp.to_csv(f'disp_data_compiled/disp_{i}.csv')
    os.system('rm disp_data_compiled/dispersion*')
data_compile()
