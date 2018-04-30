'''
Created on Jan 6, 2017
@author: javi
'''
import os
from subprocess import call
import time
import Bio.Cluster
import util
import plot
from sklearn import preprocessing
import numpy as np 
import glob
import msa
#from sklearn.datasets.california_housing import TARGET_FILENAME
import constants as cons
import pandas 
import logging
from scipy.cluster.hierarchy import  linkage
from _ast import Num

def coevolution_analisys(method, top_df,index, zmip_natural, zmip_evol, outputh_path,contact_map_path,window, pdb_name, contact_threashold=1):
    contact_threashold_str = str(contact_threashold)
    #util.save_list_to_csv(zmip_natural_result_path+"_order.csv", zmip_natural, ['Position 1',' Position 2','ZMIP'])
    contact_map=util.load_contact_map(contact_map_path)
    #Agrego los numeros de contactos de la matriz
    contacts_count = np.count_nonzero(contact_map>=contact_threashold)
    #df.set_value(index, 'contact_threashold', contact_threashold)
    #df.set_value(index, 'contacts_count', contacts_count)
    #sincronizo las senales de coevolucion para calcular spearman rank correlation
    m,m2=util.sincronice_mi(zmip_natural, zmip_evol)
    m_=[row [2] for row in m]
    m2_=[row[2] for row in m2]
    
    #Agrego spearman entre natural y evolucionado
    value_spearman = spearman(m_,m2_)
    #df.set_value(index, 'spearman_evol_nat', value_spearman) 
    #0.100578206022
    #0.102156238538
    
    #par positions 
    #df.set_value(index, 'par_positions_count', len(zmip_natural)) 
    
    m1_norm,m2_norm=normalice_(m_,m2_)
    m_np = np.c_[ np.asarray(m), np.asarray(m1_norm) ]
    m2_np = np.c_[ np.asarray(m2), np.asarray(m2_norm) ]   
    #array for contact mi matrix comparission
    x_nat_t = []
    y_evol_t = []
    x_nat_f = []
    y_evol_f = []
    for x, y in map(None, m_np, m2_np):
        #se le resta a uno porque los arrays comienzan en la posicion 0
        pos1 = int(x[0]-1)
        pos2 = int(x[1]-1)
        v = contact_map[pos1][pos2]
        if(  v < contact_threashold):
            x_nat_f.append(x[3])
            y_evol_f.append(y[3])
        else:
            x_nat_t.append(x[3])
            y_evol_t.append(y[3])
            
    plot.contacts_with_mi(x_nat_t,y_evol_t,x_nat_f,y_evol_f,outputh_path+'best_'+method+'_contacts.png',outputh_path, 'PDB: '+pdb_name+' - Coevolution ' + method+ ' pairs Evol vs Natural')
            
    
    #ordeno zmip evolucionado sincronizado 
    util.order(m2)
    #util.save_list_to_csv(mi_result_file_path+"_order.csv", m2, ['Position 1',' Position 2','ZMIP'])
    #top_df = pandas.DataFrame()        
    print "TOTAL PAR POSITIONS " + str(len(zmip_natural))
            
    result_file = open(outputh_path+method+'.txt','w')
    result_file.write(outputh_path+method+ '\n')
    result_file.write(" SPEARMAN RANK CORRELATION " + str(value_spearman)+ '\n')
    top_coevolution(zmip_natural,zmip_evol,method,0.5,contact_map,outputh_path+method+'top_0.5percent_withcon'+contact_threashold_str+'.png',outputh_path,result_file,top_df,1, pdb_name,contact_threashold)
    top_coevolution(zmip_natural,zmip_evol,method,1,contact_map,outputh_path+method+'top_1percent_withcon'+contact_threashold_str+'.png',outputh_path,result_file,top_df,2, pdb_name,contact_threashold)
    top_coevolution(zmip_natural,zmip_evol,method,2,contact_map,outputh_path+method+'top_2percent_withcon'+contact_threashold_str+'.png',outputh_path,result_file,top_df,3, pdb_name,contact_threashold)
    top_coevolution(zmip_natural,zmip_evol,method,3,contact_map,outputh_path+method+'top_3percent_withcon'+contact_threashold_str+'.png',outputh_path,result_file,top_df,4, pdb_name,contact_threashold)
    top_coevolution(zmip_natural,zmip_evol,method,4,contact_map,outputh_path+method+'top_4percent_withcon'+contact_threashold_str+'.png',outputh_path,result_file,top_df,5, pdb_name,contact_threashold)
    top_coevolution(zmip_natural,zmip_evol,method,5,contact_map,outputh_path+method+'top_5percent_withcon'+contact_threashold_str+'.png',outputh_path,result_file,top_df,6, pdb_name,contact_threashold)
    
    #top_rank(zmip_natural,m2,3,contact_map,mi_result_file_path+'top_3percent_withcon.png',mi_result_file_path,result_file)
    #top_rank(zmip_natural,m2,4,contact_map,mi_result_file_path+'top_4percent_withcon.png',mi_result_file_path,result_file)
    #top_rank(zmip_natural,m2,5,contact_map,mi_result_file_path+'top_5percent_withcon.png',mi_result_file_path,result_file)
    result_file.close()
    
    print '************************************************************************'
    

'''
Plots information about the top_rank.
For example give information about the top_rank matches
Plot a matrix with the contact and the high values (top_rank) of the evolution and the natural msa
'''
def top_coevolution(natural_coevolution,evolutionated_coevolution,method,top,contact_map,contact_map_top_coev_path,filename,result_file, top_df,index,pdb_name,contact_threashold=1):
    num = len(natural_coevolution)*top//100
    a=natural_coevolution[0:int(num)]
    b=evolutionated_coevolution[0:int(num)]
    data=matches_coevolved_positions(a,b)
    a=remove_column(a, 2)
    b=remove_column(b, 2)
    x_nat=[]
    y_nat=[]
    x_evol=[]
    y_evol=[]
    evol_contact_pairs=[]
    nat_contact = 0
    evol_contact = 0
    for x, y in map(None, a, b):
        pos1 = int(x[0]-1)
        pos2 = int(x[1]-1)
        v = contact_map[pos1][pos2]
        if(v >= contact_threashold):
            nat_contact=nat_contact+1
        pos1 = int(y[0]-1)
        pos2 = int(y[1]-1)
        v = contact_map[pos1][pos2]
        if(v >= contact_threashold):
            evol_contact=evol_contact+1 
            evol_contact_pairs.append(y)       
        
        #se le resta a uno porque los arrays comienzan en la posicion 0
        x_nat.append(int(x[0]-1))
        y_nat.append(int(x[1]-1))
        x_evol.append(int(y[0]-1))
        y_evol.append(int(y[1]-1))
    
    plot.contact_map_with_top_rank_mi(contact_map,  x_nat, y_nat, x_evol,y_evol,contact_map_top_coev_path, filename,'PDB:' + pdb_name + ' - Top ' + str(top) +'% for '+ method + ' method' ) 
    
    #find information about secondary structure.
    #find information functional information about position.
    
    result_file.write("TOP : "  + str(top) + "% PAR POSITIONS : " + str(num)+ '\n')
    result_file.write("NATURAL CONTACTS QUANTITY : " + str(nat_contact) + " - %"+ str(nat_contact*100/num)+ '\n')
    result_file.write("EVOL CONTACTS QUANTITY : " + str(evol_contact) + " - %"+ str(evol_contact*100/num)+ '\n')
    
    result_file.write("MATCH POSITIONS BETWEEN NAT AND EVOL (NO WINDOW) : " + str(len(data))+ '\n')
    result_file.write("MATCH POSITIONS  : " + str(data) + '\n')
    
    top_df.set_value(index,'nat_'+method,str(nat_contact*100/num))
    top_df.set_value(index,'evol_'+method,str(evol_contact*100/num))
    
    data_contact=[]
    for d in data:
        pos1 = int(d[0]-1)
        pos2 = int(d[1]-1)
        v=contact_map[pos1][pos2]
        if(v >= contact_threashold):
            data_contact.append(d)
    result_file.write("MATCH POSITIONS CONTACTS BETWEEN NAT AND EVOL (NO WINDOW) : " + str(len(data_contact))+ '\n')
    top_df.set_value(index,'match_'+method,str(len(data_contact)))
    result_file.write("MATCH POSITIONS CONTACTS  : " + str(data_contact) + '\n')
    result_file.write("************************************************************************" + '\n')
    
    
    print "TOP : "  + str(top) + "% PAR POSITIONS : " + str(num)
    print "MATCH POSITIONS BETWEEN NAT AND EVOL (NO WINDOW) : " + str(len(data))
    print "NATURAL CONTACTS QUANTITY : " + str(nat_contact) + " - %"+ str(nat_contact*100/num)
    print "EVOL CONTACTS QUANTITY : " + str(evol_contact) + " - %"+ str(evol_contact*100/num)
    #print data
    return nat_contact, nat_contact*100/num, evol_contact, evol_contact*100/num, len(data),len(data_contact)

def top_coevolution_analysis(method, score_coev_conformers, top_score, contact_map_path, structures,coevolution_results, natural_coevolution,coevolution_analisys_df,index_df):
    contact_map = util.load_contact_map(contact_map_path)
    fields=["Position1","Position2","Count"]
    df_total = pandas.DataFrame([],columns=fields)
    cant = 0
    #for coevolution score dendogram
    Y=[]
    for score_coev in score_coev_conformers:
        cant = cant + 1
        num = len(score_coev)*top_score//100
        zmip_evol_top=score_coev[0:int(num)]
        df = pandas.DataFrame(zmip_evol_top,columns=fields)
        df_total=df_total.append(df)
        
        #
        mi_map=np.zeros(contact_map.shape)
        for x in map(None, zmip_evol_top):
            pos1 = int(x[0]-1)
            pos2 = int(x[1]-1)
            mi_map[pos2][pos1]=1
        mi_map=mi_map.ravel()
        Y.append(mi_map)
    Z = linkage(Y, 'single')
    plot.dendogram_matrix(Z,coevolution_results + method + "_" + str(top_score) + "_dendogram.png",' Distances score top ' + str(top_score) + ' for ' + method ,structures)   
        
    #After append all the evolution MI TOP do this
    #counts = dfres.groupby(['Position1','Position2']).size()
    #print counts
    #Count cantidad de veces que aprece en los tops
    counts_df = pandas.DataFrame(df_total.groupby(['Position1','Position2']).size().rename('Count'))
    #
    sorted_df=counts_df.sort_values(by=['Count'],ascending=[False])
    #print sorted_df
    
    sorted_df['ProbTop']=pandas.Series(0.0, index=sorted_df.index)
    sorted_df['Contacts']=pandas.Series(0.0, index=sorted_df.index)
    sorted_df['ProbContact']=pandas.Series(0.0, index=sorted_df.index)
    #print prob_contact_map
    for index,mi_par in sorted_df.iterrows():
        #probablidad de que aprezca en top
        prob_top = mi_par['Count'] * 100 / cant 
        sorted_df.set_value(index, 'ProbTop' , prob_top/100)
        #cantidad de contactos 
        pos1 = int(index[0]-1)
        pos2 = int(index[1]-1)
        v=contact_map[pos1][pos2]
        sorted_df.set_value(index, 'Contacts' , v)
        #probablidad de contactos
        prob_contact = float(v) * 100 / cant
        sorted_df.set_value(index, 'ProbContact' , float(prob_contact)/100)
        
    sorted_df=sorted_df.sort_values(by=['Count','Contacts'],ascending=[False,False])
    sorted_df.to_csv(coevolution_results + method + "_" + str(top_score) + "_coevolution.csv", sep='\t', encoding='utf-8')
    
    
    
    df = pandas.read_csv(coevolution_results + method + "_" + str(top_score) + "_coevolution.csv",delim_whitespace=True,header=0,usecols=[0,1,2,4])
    #df=df.sort(['Count', 'Contacts'], ascending=[False, False])
    df=df.sort_values(by=['Count', 'Contacts'],ascending=[False, False])
    #contact map with sum score
    cmap_triu=np.triu(contact_map, -1)
    top_score_list=df.values.tolist()
    #add mi values to map 
    for x in map(None, top_score_list):
        pos1 = int(x[0]-1)
        pos2 = int(x[1]-1)
        v = cmap_triu[pos2][pos1]
        if(v!=0):
            print ("error")
        #set the value of the smi    
        cmap_triu[pos2][pos1]=x[2]
    util.save_contact_map(cmap_triu, coevolution_results + "contact_map_top_" + str(top_score)+ "_"+method+".dat")    
    plot.contact_map_sum(cmap_triu,coevolution_results + "contact_map_top_" + str(top_score)+ "_"+method+".png",'Sum Contact Map and Sum Top '+str(top_score)+' for ' + method)   
    
    
    
    
    
    #Review how to show the conjunction of top SCORES
    result_file = open(coevolution_results+method+'.txt','w')
    result_file.write(coevolution_results+method+ '\n')
    
    
    coevolution_analisys_df.set_value(index_df,'top',str(top_score))
    coevolution_analisys_df.set_value(index_df,'contact_threashold',1)
    top_coevolution(natural_coevolution,top_score_list,method,top_score,contact_map,coevolution_results+ method+"_top_"+ str(top_score)+ 'percent_contact_threashold_1.png',coevolution_results,result_file, coevolution_analisys_df,index_df,'Thio Ecoli Conformers',1)
    index_df=index_df + 1
    coevolution_analisys_df.set_value(index_df,'top',str(top_score))
    coevolution_analisys_df.set_value(index_df,'contact_threashold',4)
    top_coevolution(natural_coevolution,top_score_list,method,top_score,contact_map,coevolution_results+ method+"_top_"+ str(top_score)+ 'percent_contact_threashold_4.png',coevolution_results,result_file, coevolution_analisys_df,index_df,'Thio Ecoli Conformers',4)
    index_df=index_df + 1
    coevolution_analisys_df.set_value(index_df,'top',str(top_score))
    coevolution_analisys_df.set_value(index_df,'contact_threashold',8)
    top_coevolution(natural_coevolution,top_score_list,method,top_score,contact_map,coevolution_results+ method+"_top_"+ str(top_score)+ 'percent_contact_threashold_8.png',coevolution_results,result_file, coevolution_analisys_df,index_df,'Thio Ecoli Conformers',8)
    
    result_file.close()
    
    
    
    '''
    correlation_p = sorted_df['Count'].corr(sorted_df['Contacts'], method='pearson')
    correlation_k = sorted_df['Count'].corr(sorted_df['Contacts'], method='kendall')
    correlation_s = sorted_df['Count'].corr(sorted_df['Contacts'], method='spearman')

    mean = sorted_df["ProbTop"].mean()
    median = sorted_df["ProbTop"].median()
    var = sorted_df["ProbTop"].var()
    mode = sorted_df["ProbTop"].mode()
    
    corr_df = pandas.DataFrame()
    corr_df.set_value(1, 'pearson', correlation_p)
    corr_df.set_value(1, 'kendall', correlation_k)
    corr_df.set_value(1, 'spearman', correlation_s)
    corr_df.to_csv(result_top  + "_corr.csv", sep='\t', encoding='utf-8')
    
    import matplotlib.pyplot as plt
    sorted_df.plot.scatter(x='Count', y='Contacts');
    plt.xlabel('Numero de veces que aparece en el top en la MI Conjunta')
    plt.ylabel('Numero de contactos de la matriz conjunta')
    plt.xticks([0,1,2,3,4,5,6,7,8])
    plt.yticks([0,1,2,3,4,5,6,7,8])
    plt.savefig(result_top  + "_corr.png");
    plt.gcf().clear()
    '''
    
def dendogram_matrix(contact_maps_paths,output_path,title,structures,clustering_type):
    Y=[]
    for cmap_pth in contact_maps_paths:
        contact_map = util.load_contact_map(cmap_pth)
        contact_map=contact_map.ravel()
        Y.append(contact_map)
    Z = linkage(Y, clustering_type)
    plot.dendogram_matrix(Z,output_path,title,structures)
    
def dendogram_top_mi(top, zmip_nat,zmip_paths, output_path,title ,structures, clustering_type, shape ,window):
    Y=[]
    zmip_natural = util.load_zmip(zmip_nat,window)
    util.order(zmip_natural)
    for zmip_path in zmip_paths:
        zmip = util.load_zmip(zmip_path,window)
        zmip_natural,zmip=util.sincronice_mi(zmip_natural, zmip)
        util.order(zmip)
        #pares de posiciones 
        num = len(zmip)*top//100
        num = int(num)
        zmip=zmip[0:num]   
        cmap_with_mi=np.zeros(shape)
        for x in map(None, zmip):
            pos1 = int(x[0]-1)
            pos2 = int(x[1]-1)
            cmap_with_mi[pos2][pos1]=1
        print zmip
        print "-------"
        print cmap_with_mi
        cmap_with_mi=cmap_with_mi.ravel()
        Y.append(cmap_with_mi)
    Z = linkage(Y, clustering_type)
    plot.dendogram_matrix(Z,output_path,title,structures)        
    
def plot_optimization(optimization_results):
    df = pandas.read_csv(optimization_results, header=0, usecols=['auc_mi','auc_di','auc_frob','auc_psicov', 'beta', 'nsus', 'run'])
    df_to_plot = df[(df['run'] == 20000)]
    df_to_plot = df_to_plot[(df_to_plot['beta'] == 0.1) | (df_to_plot['beta'] == 0.25) | (df_to_plot['beta'] == 0.5) | (df_to_plot['beta'] == 1) | (df_to_plot['beta'] == 5) | (df_to_plot['beta'] == 7) | (df_to_plot['beta'] == 10) | (df_to_plot['beta'] == 15) | (df_to_plot['beta'] == 20)]
    colors = ['tan', 'chocolate' ,'yellow', 'yellowgreen', 'red', 'orange', 'coral', 'olive', 'green']
    # df_to_plot=df[(df['beta']==0.1) | (df['beta']==0.5) | (df['beta']==1) | (df['beta']==2) | (df['beta']==3) | (df['beta']==5)]
    plot.auc_optimization(df_to_plot, 'auc_mi' , colors, optimization_results + '_mi.png','MI Optimization Values')
    plot.auc_optimization(df_to_plot, 'auc_di' , colors, optimization_results + '_di.png','DI Optimization Values')
    plot.auc_optimization(df_to_plot, 'auc_frob' , colors, optimization_results + '_frob.png','FROB Optimization Values')
    plot.auc_optimization(df_to_plot, 'auc_psicov' , colors, optimization_results + '_psicov.png','PSICOV Optimization Values')

'''
Calculate the AUC.  
For the protein family (fasta_path) and the contact_map calculates the AUC. 
'''
#interpolate false
def auc(fasta_path,contact_map):
    start_time = time.time()
    print "auc"
    #call(["julia", "auc_script.jl" ])
    call([cons.julia_exe_path, "mitos/auc.jl",fasta_path,contact_map])
    print "auc"
    print("--- %s seconds ---" % (time.time() - start_time)) 
'''
Calculate the AUC.  
For all the families in clustered_sequences_path and the contact_map_path
The results are:
a zmip for each of the families (result_zmip_path); and
a text file: result_auc_file_name wich contains the AUC for every family. 
'''
def auc_job(pdb_name,model_name,chain_name,contact_map_path,clustered_sequences_path,result_auc_file_name, result_zmip_path):
    start_time = time.time()
    print "auc_process_all"
    #call(["julia", "auc_script.jl" ])
    call([cons.julia_exe_path, "mitos/auc_process_all.jl",pdb_name,model_name,chain_name,contact_map_path,clustered_sequences_path,result_auc_file_name,result_zmip_path])
    print "auc_process_all"
    print("--- %s seconds ---" % (time.time() - start_time)) 




def evol_analisys(msa_file, mi_data_output_path, msa_conservation_path,msa_name):
    msa.buslje09(msa_file,mi_data_output_path)
    msa.msa_information(msa_file, msa_conservation_path,msa_name)

def run_analisys(df,index, zmip_natural_result_path, mi_results_path, pattern_array,contact_map_path,outputpath,window):
    #levanto el zmip natural
    zmip_natural = util.load_zmip(zmip_natural_result_path,window)
    util.order(zmip_natural)
    util.save_list_to_csv(zmip_natural_result_path+"_order.csv", zmip_natural, ['Position 1',' Position 2','ZMIP'])
    contact_map=util.load_contact_map(contact_map_path)
    
    for filename in os.listdir(mi_results_path):
        if filename.endswith(".dat") & any(r in filename for r in pattern_array):
            print " Calculation of : " + filename + " with contact map " + contact_map_path
            #levanto el zmip evolucionado 
            zmip_evol = util.load_zmip(mi_results_path + filename,window)
            #sincronizo las senales de coevolucion para calcular spearman rank correlation
            m,m2=util.sincronice_mi(zmip_natural, zmip_evol)
            m_=[row [2] for row in m]
            m2_=[row[2] for row in m2]
            value_spearman = spearman(m_,m2_)
            df.set_value(index, 'spearman_zmip_evol_nat', value_spearman) 
            #MI PLOT
            #contact_map=util.load_contact_map(contact_map_path)
            #plot.contact_map_(contact_map, outputpath)
            #Test load_contact_map
            #v=contact_map[4][2]
            #v=contact_map[4][1]
            #v=contact_map[4][3]
            m1_norm,m2_norm=normalice_(m_,m2_)
            m_np = np.c_[ np.asarray(m), np.asarray(m1_norm) ]
            m2_np = np.c_[ np.asarray(m2), np.asarray(m2_norm) ]   
            #array for contact mi matrix comparission
            x_nat_t = []
            y_evol_t = []
            x_nat_f = []
            y_evol_f = []
            #scores mi, for roc_curve 
            scores_nat = []
            scores_evol = []
            #1 contact
            #0 no contact
            y_true = []
            for x, y in map(None, m_np, m2_np):
                #se le resta a uno porque los arrays comienzan en la posicion 0
                pos1 = int(x[0]-1)
                pos2 = int(x[1]-1)
                v = contact_map[pos1][pos2]
                #x[0] = posicion 1
                #x[1] = posicion 2
                #x[2] = mi sin normalizar
                #x[3] = mi normalizado
                scores_nat.append(x[3])
                scores_evol.append(y[3])
                if(  v == 0):
                    x_nat_f.append(x[3])
                    y_evol_f.append(y[3])
                    y_true.append(0)
                else:
                    x_nat_t.append(x[3])
                    y_evol_t.append(y[3])
                    y_true.append(1)
            
            #y_true = np.array([0, 1, 1, 1])
            #y_scores = np.array([0.1, 0.4, 0.35, 0.8])
            #y_true = np.array([0, 0, 0, 1])
            #target = np.array([0.1, 0.4, 0.35, 0.8])
           
            labels=['Natural', 'Evol']
            scores=[]
            scores.append(scores_nat)
            scores.append(scores_evol)
            #plot.roc_curve(y_true,scores_nat,scores_evol)
            colors = ['blue', 'red']
            plot.roc_curve(df,index,y_true,scores,labels,colors, outputpath+filename+'_roc_curve.png')
            
            plot.contacts_with_mi(x_nat_t,y_evol_t,x_nat_f,y_evol_f,outputpath+filename+'contacts_with_mi.png',filename)
            
            #ordeno zmip evolucionado sincronizado 
            util.order(m2)
            
            util.save_list_to_csv(mi_results_path+filename+"_order.csv", m2, ['Position 1',' Position 2','ZMIP'])
            
            print "TOTAL PAR POSITIONS " + str(len(zmip_natural))
            
            result_file = open(outputpath+filename+".txt","w")
            result_file.write(filename+ '\n')
            result_file.write(" SPEARMAN RANK CORRELATION " + str(value_spearman)+ '\n')
            top_rank(zmip_natural,m2,0.5,contact_map,outputpath+filename+'top_0.5percent_withcon.png',filename,result_file)
            top_rank(zmip_natural,m2,1,contact_map,outputpath+filename+'top_1percent_withcon.png',filename,result_file)
            top_rank(zmip_natural,m2,2,contact_map,outputpath+filename+'top_2percent_withcon.png',filename,result_file)
            top_rank(zmip_natural,m2,3,contact_map,outputpath+filename+'top_3percent_withcon.png',filename,result_file)
            top_rank(zmip_natural,m2,4,contact_map,outputpath+filename+'top_4percent_withcon.png',filename,result_file)
            top_rank(zmip_natural,m2,5,contact_map,outputpath+filename+'top_5percent_withcon.png',filename,result_file)
            result_file.close()
            print '************************************************************************'

def getTargetScores(mi_file_path,contact_map, mi_data_path_natural=None,sincronized_with_natural=False,window=1,threshold=0,split_char=r'\t+',score_position=2):
    cmap=util.load_contact_map(contact_map,np.float64)
    
    cmap[cmap > threshold] = 1
    cmap[cmap <= threshold] = 0
    
    zmip_evol = util.load_zmip(mi_file_path, window,split_char)
    if(sincronized_with_natural==True):
        zmip_natural = util.load_zmip(mi_data_path_natural,window)
        m,zmip_evol=util.sincronice_mi(zmip_natural, zmip_evol)
    
    scores = []
    target = []
    for x in zmip_evol:
        v = cmap[int(x[0]-1)][int(x[1]-1)] 
        scores.append(x[score_position])
        target.append(v)
    return target,scores    
def run_analisys_singular(df,index, zmip_natural_result_path, mi_result_file_path, contact_map_path,outputpath,window, pdb_name, contact_threashold=1):
    contact_threashold_str = str(contact_threashold)
    #levanto el zmip natural
    zmip_natural = util.load_zmip(zmip_natural_result_path,window)
    util.order(zmip_natural)
    #util.save_list_to_csv(zmip_natural_result_path+"_order.csv", zmip_natural, ['Position 1',' Position 2','ZMIP'])
    contact_map=util.load_contact_map(contact_map_path)
    #Agrego los numeros de contactos de la matriz
    contacts_count = np.count_nonzero(contact_map>=contact_threashold)
    df.set_value(index, 'contact_threashold', contact_threashold)
    df.set_value(index, 'contacts_count', contacts_count)
    
    print " Calculation of : " + mi_result_file_path + " with contact map " + contact_map_path
    #levanto el zmip evolucionado 
    zmip_evol = util.load_zmip(mi_result_file_path ,window)
    #sincronizo las senales de coevolucion para calcular spearman rank correlation
    m,m2=util.sincronice_mi(zmip_natural, zmip_evol)
    m_=[row [2] for row in m]
    m2_=[row[2] for row in m2]
    
    #Agrego spearman entre natural y evolucionado
    value_spearman = spearman(m_,m2_)
    df.set_value(index, 'spearman_zmip_evol_nat', value_spearman) 
    #0.100578206022
    #0.102156238538
    
    #par positions 
    df.set_value(index, 'par_positions_count', len(zmip_natural)) 
    
    m1_norm,m2_norm=normalice_(m_,m2_)
    m_np = np.c_[ np.asarray(m), np.asarray(m1_norm) ]
    m2_np = np.c_[ np.asarray(m2), np.asarray(m2_norm) ]   
    #array for contact mi matrix comparission
    x_nat_t = []
    y_evol_t = []
    x_nat_f = []
    y_evol_f = []
    #scores mi, for roc_curve 
    scores_nat = []
    scores_evol = []
    #1 contact
    #0 no contact
    y_true = []
    for x, y in map(None, m_np, m2_np):
        #se le resta a uno porque los arrays comienzan en la posicion 0
        pos1 = int(x[0]-1)
        pos2 = int(x[1]-1)
        v = contact_map[pos1][pos2]
        scores_nat.append(x[3])
        scores_evol.append(y[3])
        if(  v < contact_threashold):
            x_nat_f.append(x[3])
            y_evol_f.append(y[3])
            y_true.append(0)
        else:
            x_nat_t.append(x[3])
            y_evol_t.append(y[3])
            y_true.append(1)
            
    labels=['Natural', 'Evol']
    scores=[]
    scores.append(scores_nat)
    scores.append(scores_evol)
    #plot.roc_curve(y_true,scores_nat,scores_evol)
    colors = ['blue', 'red']
    plot.roc_curve(df,index,y_true,scores,labels,colors, mi_result_file_path+'_roc_curve_'+contact_threashold_str+'.png')
            
    plot.contacts_with_mi(x_nat_t,y_evol_t,x_nat_f,y_evol_f,mi_result_file_path+'contacts_with_mi'+contact_threashold_str+'.png',mi_result_file_path, pdb_name)
            
    #ordeno zmip evolucionado sincronizado 
    util.order(m2)
            
    #util.save_list_to_csv(mi_result_file_path+"_order.csv", m2, ['Position 1',' Position 2','ZMIP'])
    top_df = pandas.DataFrame()        
    print "TOTAL PAR POSITIONS " + str(len(zmip_natural))
            
    result_file = open(mi_result_file_path+contact_threashold_str+'.txt','w')
    result_file.write(mi_result_file_path+ '\n')
    result_file.write(" SPEARMAN RANK CORRELATION " + str(value_spearman)+ '\n')
    top_rank_result = top_rank(zmip_natural,m2,0.5,contact_map,mi_result_file_path+'top_0.5percent_withcon'+contact_threashold_str+'.png',mi_result_file_path,result_file,top_df,1, pdb_name,contact_threashold)
    top_rank_result = top_rank(zmip_natural,m2,1,contact_map,mi_result_file_path+'top_1percent_withcon'+contact_threashold_str+'.png',mi_result_file_path,result_file,top_df,2, pdb_name,contact_threashold)
    top_rank_result = top_rank(zmip_natural,m2,2,contact_map,mi_result_file_path+'top_2percent_withcon'+contact_threashold_str+'.png',mi_result_file_path,result_file,top_df,3, pdb_name,contact_threashold)
    top_rank_result = top_rank(zmip_natural,m2,3,contact_map,mi_result_file_path+'top_3percent_withcon'+contact_threashold_str+'.png',mi_result_file_path,result_file,top_df,4, pdb_name,contact_threashold)
    top_rank_result = top_rank(zmip_natural,m2,4,contact_map,mi_result_file_path+'top_4percent_withcon'+contact_threashold_str+'.png',mi_result_file_path,result_file,top_df,5, pdb_name,contact_threashold)
    top_rank_result = top_rank(zmip_natural,m2,5,contact_map,mi_result_file_path+'top_5percent_withcon'+contact_threashold_str+'.png',mi_result_file_path,result_file,top_df,6, pdb_name,contact_threashold)
    
    #top_rank(zmip_natural,m2,3,contact_map,mi_result_file_path+'top_3percent_withcon.png',mi_result_file_path,result_file)
    #top_rank(zmip_natural,m2,4,contact_map,mi_result_file_path+'top_4percent_withcon.png',mi_result_file_path,result_file)
    #top_rank(zmip_natural,m2,5,contact_map,mi_result_file_path+'top_5percent_withcon.png',mi_result_file_path,result_file)
    result_file.close()
    top_df.to_csv(mi_result_file_path + '_top'+contact_threashold_str+'.csv')
    print '************************************************************************'

'''
Normalize information m,m2
'''
def normalice_desarrollo(m,m2,m3):
    #import numpy as np   
    #X_train = np.array([[ 1., -1.,  2.],[ 2.,  0.,  0.],[ 0.,  1., -1.]])
    min_max_scaler = preprocessing.MinMaxScaler()
    m_train_minmax = min_max_scaler.fit_transform(m)
    m2_train_minmax = min_max_scaler.fit_transform(m2) 
    m3_train_minmax = min_max_scaler.fit_transform(m3) 
    return m_train_minmax,m2_train_minmax,m3_train_minmax           
           
'''
Normalize information m,m2
'''
def normalice_(m,m2):
    #import numpy as np   
    #X_train = np.array([[ 1., -1.,  2.],[ 2.,  0.,  0.],[ 0.,  1., -1.]])
    min_max_scaler = preprocessing.MinMaxScaler()
    m_train_minmax = min_max_scaler.fit_transform( np.asarray(m).reshape(-1, 1))
    m2_train_minmax = min_max_scaler.fit_transform( np.asarray(m2).reshape(-1, 1)) 
    return m_train_minmax,m2_train_minmax   
def spearman(x,y):
    return 1 - Bio.Cluster.distancematrix((x,y), dist="s")[1][0]

def matches_coevolved_positions(matrix_ref,matrix_evol, intersection_evol=False):
    data = []
    for j in matrix_ref:
        pos  = j[0]     
        pos2 = j[1]   
        for e in matrix_evol:
            pos_  = e[0]     
            pos2_ = e[1]
            if (pos_==pos and pos2_==pos2):
                #j.append(e[2])
                aux=j[:]
                if (intersection_evol==False):
                    aux.append(e[2])   
                data.append(aux)
                break    
    return  data
def remove_column(matrix, column):
    return [row[:column] + row[column+1:] for row in matrix]

'''
Plots information about the top_rank.
For example give information about the top_rank matches
Plot a matrix with the contact and the high values (top_rank) of the evolution and the natural msa
'''
def top_rank(x,y,top,contact_map,outputpath,filename,result_file, top_df,index,pdb_name,contact_threashold=1):
    num = len(x)*top//100
    a=x[0:int(num)]
    b=y[0:int(num)]
    print a
    print b
    data=matches_coevolved_positions(a,b)
    a=remove_column(a, 2)
    b=remove_column(b, 2)
    x_nat=[]
    y_nat=[]
    x_evol=[]
    y_evol=[]
    evol_contact_pairs=[]
    nat_contact = 0
    evol_contact = 0
    for x, y in map(None, a, b):
        pos1 = int(x[0]-1)
        pos2 = int(x[1]-1)
        v = contact_map[pos1][pos2]
        if(v >= contact_threashold):
            nat_contact=nat_contact+1
        pos1 = int(y[0]-1)
        pos2 = int(y[1]-1)
        v = contact_map[pos1][pos2]
        if(v >= contact_threashold):
            evol_contact=evol_contact+1 
            evol_contact_pairs.append(y)       
        
        #se le resta a uno porque los arrays comienzan en la posicion 0
        x_nat.append(int(x[0]-1))
        y_nat.append(int(x[1]-1))
        x_evol.append(int(y[0]-1))
        y_evol.append(int(y[1]-1))
    
    plot.contact_map_with_top_rank_mi(contact_map,  x_nat, y_nat, x_evol,y_evol,outputpath,filename,pdb_name + ' Top ' + str(top) +'%')
    
    #find information about secondary structure.
    #find information functional information about position.
    
    result_file.write("TOP : "  + str(top) + "% PAR POSITIONS : " + str(num)+ '\n')
    result_file.write("NATURAL CONTACTS QUANTITY : " + str(nat_contact) + " - %"+ str(nat_contact*100/num)+ '\n')
    result_file.write("EVOL CONTACTS QUANTITY : " + str(evol_contact) + " - %"+ str(evol_contact*100/num)+ '\n')
    
    result_file.write("MATCH POSITIONS BETWEEN NAT AND EVOL (NO WINDOW) : " + str(len(data))+ '\n')
    result_file.write("MATCH POSITIONS  : " + str(data) + '\n')
    
    top_df.set_value(index,'top',str(top))
    top_df.set_value(index,'par_positions',str(num))
    top_df.set_value(index,'nat_contact',str(nat_contact))
    top_df.set_value(index,'nat_contact_%',str(nat_contact*100/num))
    top_df.set_value(index,'evol_contact',str(evol_contact))
    top_df.set_value(index,'evol_contact_%',str(evol_contact*100/num))
    
    data_contact=[]
    for d in data:
        pos1 = int(d[0]-1)
        pos2 = int(d[1]-1)
        v=contact_map[pos1][pos2]
        if(v >= contact_threashold):
            data_contact.append(d)
    result_file.write("MATCH POSITIONS CONTACTS BETWEEN NAT AND EVOL (NO WINDOW) : " + str(len(data_contact))+ '\n')
    top_df.set_value(index,'match_positions',str(len(data_contact)))
    result_file.write("MATCH POSITIONS CONTACTS  : " + str(data_contact) + '\n')
    result_file.write("************************************************************************" + '\n')
    
    
    print "TOP : "  + str(top) + "% PAR POSITIONS : " + str(num)
    print "MATCH POSITIONS BETWEEN NAT AND EVOL (NO WINDOW) : " + str(len(data))
    print "NATURAL CONTACTS QUANTITY : " + str(nat_contact) + " - %"+ str(nat_contact*100/num)
    print "EVOL CONTACTS QUANTITY : " + str(evol_contact) + " - %"+ str(evol_contact*100/num)
    #print data
    return nat_contact, nat_contact*100/num, evol_contact, evol_contact*100/num, len(data),len(data_contact)

def top_rank_intersection(execution_folder,generic_top,  contact_map_path,zmip_natural_result_path, zmip_evol_intersect_result_path, top_df, zmip_reference_result_path, zmip_prom_result_path, index, window, contact_threshold, top_threshold, sinchronize_with_natural=True):
    
    df = pandas.read_csv(zmip_evol_intersect_result_path,delim_whitespace=True,header=0,usecols=[0,1,2])
    df=df.loc[df['Count']>=top_threshold]
    num=len(df.index)
    
    contact_map= util.load_contact_map(contact_map_path)
    zmip_natural = util.load_zmip(zmip_natural_result_path,window)
    util.order(zmip_natural)
    natural=zmip_natural[0:num]
    
    evol=df[0:num]
    evol=evol.values.tolist()
    
    data=matches_coevolved_positions(natural,evol,intersection_evol=True)
    
    zmip_reference = util.load_zmip(zmip_reference_result_path,window)
    zmip_prom = util.load_zmip(zmip_prom_result_path,window)
    if(sinchronize_with_natural==True):
        zmip_natural,zmip_reference=util.sincronice_mi(zmip_natural, zmip_reference)
        zmip_natural,zmip_prom=util.sincronice_mi(zmip_natural, zmip_prom)
    
    util.order(zmip_reference)
    util.order(zmip_prom)
    
    reference=zmip_reference[0:num]
    prom=zmip_prom[0:num]
    data_reference=matches_coevolved_positions(natural,reference,intersection_evol=True)
    data_prom=matches_coevolved_positions(natural,prom,intersection_evol=True)
    
    x_nat=[]
    y_nat=[]
    x_evol=[]
    y_evol=[]
    x_ref=[]
    y_ref=[]
    x_prom=[]
    y_prom=[]
    
    
    
    evol_contact_pairs=[]
    ref_contact_pairs=[]
    prom_contact_pairs=[]
    nat_contact = 0
    evol_contact = 0
    ref_contact = 0
    prom_contact = 0
    for x, y, z, p in map(None, natural, evol,reference,prom):
        x_nat.append(int(x[0]-1))
        y_nat.append(int(x[1]-1))
        x_evol.append(int(y[0]-1))
        y_evol.append(int(y[1]-1))
        x_ref.append(int(z[0]-1))
        y_ref.append(int(z[1]-1))
        
        x_prom.append(int(p[0]-1))
        y_prom.append(int(p[1]-1))
        
        
        
        pos1 = int(x[0]-1)
        pos2 = int(x[1]-1)
        v = contact_map[pos1][pos2]
        if(v >= contact_threshold):
            nat_contact=nat_contact+1
        
        pos1 = int(y[0]-1)
        pos2 = int(y[1]-1)
        v = contact_map[pos1][pos2]
        if(v >= contact_threshold):
            evol_contact=evol_contact+1 
            evol_contact_pairs.append(y)
        
        pos1 = int(z[0]-1)
        pos2 = int(z[1]-1)
        v = contact_map[pos1][pos2]
        if(v >= contact_threshold):
            ref_contact=ref_contact+1 
            ref_contact_pairs.append(z)
        
        pos1 = int(p[0]-1)
        pos2 = int(p[1]-1)
        v = contact_map[pos1][pos2]
        if(v >= contact_threshold):
            prom_contact=prom_contact+1 
            prom_contact_pairs.append(p)
        
        
    plot.contact_map_with_top_rank_mi_sum(contact_map,  x_nat, y_nat, x_evol,y_evol,execution_folder + 'top_'+ generic_top +'contact_map_thresold_'+str(contact_threshold)+'with_intersection_top_threshold_'+str(top_threshold)+'.png','', 'top_'+ generic_top +'THIO_ECOLI Top Thresold ' + str(top_threshold) +' Contact Threshold ' + str(contact_threshold))
    plot.contact_map_with_top_rank_mi_sum(contact_map,  x_nat, y_nat, x_ref,y_ref,execution_folder + 'top_'+ generic_top +'contact_map_thresold_'+str(contact_threshold)+'with_reference_top_threshold_'+str(top_threshold)+'.png','', 'top_'+ generic_top +'THIO_ECOLI 2TRX Thresold ' + str(top_threshold) +' Contact Threshold ' + str(contact_threshold))
    plot.contact_map_with_top_rank_mi_sum(contact_map,  x_nat, y_nat, x_prom,y_prom,execution_folder + 'top_'+ generic_top +'contact_map_thresold_'+str(contact_threshold)+'with_prom_top_threshold_'+str(top_threshold)+'.png','', 'top_'+ generic_top +'THIO_ECOLI Prom Thresold ' + str(top_threshold) +' Contact Threshold ' + str(contact_threshold))
    
    top_df.set_value(index,'top_generic',generic_top)
    top_df.set_value(index,'top_threshold',str(top_threshold))
    top_df.set_value(index,'contact_threshold',str(contact_threshold))
    top_df.set_value(index,'par_positions',str(num))
    top_df.set_value(index,'nat_contact',str(nat_contact))
    top_df.set_value(index,'nat_contact_%',str(nat_contact*100/num))
    top_df.set_value(index,'evol_contact',str(evol_contact))
    top_df.set_value(index,'evol_contact_%',str(evol_contact*100/num))
    
    data_contact=[]
    for d in data:
        pos1 = int(d[0]-1)
        pos2 = int(d[1]-1)
        v=contact_map[pos1][pos2]
        if(v>=contact_threshold):
            data_contact.append(d)
    top_df.set_value(index,'match_positions',str(len(data_contact)))
    
    top_df.set_value(index,'prom_contact',str(prom_contact))
    top_df.set_value(index,'prom_contact_%',str(prom_contact*100/num))
    data_contact_prom=[]
    for d in data_prom:
        pos1 = int(d[0]-1)
        pos2 = int(d[1]-1)
        v=contact_map[pos1][pos2]
        if(v>=contact_threshold):
            data_contact_prom.append(d)
    top_df.set_value(index,'match_positions_prom',str(len(data_contact_prom)))
    
    top_df.set_value(index,'ref_contact',str(ref_contact))
    top_df.set_value(index,'ref_contact_%',str(ref_contact*100/num))
    
    data_contact_ref=[]
    for d in data_reference:
        pos1 = int(d[0]-1)
        pos2 = int(d[1]-1)
        v=contact_map[pos1][pos2]
        if(v>=contact_threshold):
            data_contact_ref.append(d)
    top_df.set_value(index,'match_positions_reference',str(len(data_contact_ref)))
    
    
def call_top_rank(x, y, top, contact_map, outputpath, filename, result_file, top_df, index, pdb_name):
    num = len(x)*top//100
    a=x[0:int(num)]
    b=y[0:int(num)]
    top_rank_generic(num, x, y, top, contact_map, outputpath, filename, result_file, top_df, index, pdb_name)

def top_rank_generic(num, x, y, top_title,contact_map,outputpath,filename,result_file, top_df,index,pdb_name, contact_threashold=1):
    a=x[0:int(num)]
    b=y[0:int(num)]
    data=matches_coevolved_positions(a,b)
    a=remove_column(a, 2)
    b=remove_column(b, 2)
    x_nat=[]
    y_nat=[]
    x_evol=[]
    y_evol=[]
    evol_contact_pairs=[]
    nat_contact = 0
    evol_contact = 0
    for x, y in map(None, a, b):
        pos1 = int(x[0]-1)
        pos2 = int(x[1]-1)
        v = contact_map[pos1][pos2]
        if(v >= contact_threashold):
            nat_contact=nat_contact+1
        pos1 = int(y[0]-1)
        pos2 = int(y[1]-1)
        v = contact_map[pos1][pos2]
        if(v >= contact_threashold):
            evol_contact=evol_contact+1 
            evol_contact_pairs.append(y)       
        
        #se le resta a uno porque los arrays comienzan en la posicion 0
        x_nat.append(int(x[0]-1))
        y_nat.append(int(x[1]-1))
        x_evol.append(int(y[0]-1))
        y_evol.append(int(y[1]-1))
    
    plot.contact_map_with_top_rank_mi(contact_map,  x_nat, y_nat, x_evol,y_evol,outputpath,filename,pdb_name + ' Top ' + str(top_title) +'%')
    
    #find information about secondary structure.
    #find information functional information about position.
    
    result_file.write("TOP : "  + str(top_title) + "% PAR POSITIONS : " + str(num)+ '\n')
    result_file.write("NATURAL CONTACTS QUANTITY : " + str(nat_contact) + " - %"+ str(nat_contact*100/num)+ '\n')
    result_file.write("EVOL CONTACTS QUANTITY : " + str(evol_contact) + " - %"+ str(evol_contact*100/num)+ '\n')
    
    result_file.write("MATCH POSITIONS BETWEEN NAT AND EVOL (NO WINDOW) : " + str(len(data))+ '\n')
    result_file.write("MATCH POSITIONS  : " + str(data) + '\n')
    
    top_df.set_value(index,'top',str(top_title))
    top_df.set_value(index,'par_positions',str(num))
    top_df.set_value(index,'nat_contact',str(nat_contact))
    top_df.set_value(index,'nat_contact_%',str(nat_contact*100/num))
    top_df.set_value(index,'evol_contact',str(evol_contact))
    top_df.set_value(index,'evol_contact_%',str(evol_contact*100/num))
    
    data_contact=[]
    for d in data:
        pos1 = int(d[0]-1)
        pos2 = int(d[1]-1)
        v=contact_map[pos1][pos2]
        if(v >= contact_threashold):
            data_contact.append(d)
    result_file.write("MATCH POSITIONS CONTACTS BETWEEN NAT AND EVOL (NO WINDOW) : " + str(len(data_contact))+ '\n')
    top_df.set_value(index,'match_positions',str(len(data_contact)))
    result_file.write("MATCH POSITIONS CONTACTS  : " + str(data_contact) + '\n')
    result_file.write("************************************************************************" + '\n')
    
    
    print "TOP : "  + str(top_title) + "% PAR POSITIONS : " + str(num)
    print "MATCH POSITIONS BETWEEN NAT AND EVOL (NO WINDOW) : " + str(len(data))
    print "NATURAL CONTACTS QUANTITY : " + str(nat_contact) + " - %"+ str(nat_contact*100/num)
    print "EVOL CONTACTS QUANTITY : " + str(evol_contact) + " - %"+ str(evol_contact*100/num)
    #print data
    return nat_contact, nat_contact*100/num, evol_contact, evol_contact*100/num, len(data),len(data_contact)

    

'''
Plots information about the top_rank.
For example give information about the top_rank matches
Plot a matrix with the contact and the high values (top_rank) of the evolution and the natural msa
'''
def top_rank_desa(x,evol1,evol2,top,contact_map,outputpath,filename,result_file):
    num = len(x)*top/100
    a=x[0:int(num)]
    b=evol1[0:int(num)]
    c=evol2[0:int(num)]
    
    data=matches_coevolved_positions(b,c)
    a=remove_column(a, 2)
    b=remove_column(b, 2)
    c=remove_column(c, 2)
    x_nat=[]
    y_nat=[]
    x_evol1=[]
    y_evol1=[]
    x_evol2=[]
    y_evol2=[]
    evol_contact_pairs=[]
    evol_contact_pairs2=[]
    nat_contact = 0
    evol_contact = 0
    evol_contact2 = 0
    for x, y, f in map(None, a, b, c):
        pos1 = int(x[0]-1)
        pos2 = int(x[1]-1)
        v = contact_map[pos1][pos2]
        if(v == 1):
            nat_contact=nat_contact+1
            
        pos1 = int(y[0]-1)
        pos2 = int(y[1]-1)
        v = contact_map[pos1][pos2]
        if(v == 1):
            evol_contact=evol_contact+1 
            evol_contact_pairs.append(y)
        
        pos1 = int(f[0]-1)
        pos2 = int(f[1]-1)
        v = contact_map[pos1][pos2]
        if(v == 1):
            evol_contact2=evol_contact2+1 
            evol_contact_pairs2.append(y)           
        
        #se le resta a uno porque los arrays comienzan en la posicion 0
        x_nat.append(int(x[0]-1))
        y_nat.append(int(x[1]-1))
        
        x_evol1.append(int(y[0]-1))
        y_evol1.append(int(y[1]-1))
        
        x_evol2.append(int(f[0]-1))
        y_evol2.append(int(f[1]-1))
    
    plot.contact_map_with_top_rank_mi_desarrollo(contact_map,  x_nat, y_nat, x_evol1,y_evol1,x_evol2,y_evol2,outputpath,filename)
    
    #find information about secondary structure.
    #find information functional information about position.
    
    result_file.write("TOP : "  + str(top) + "% PAR POSITIONS : " + str(num)+ '\n')
    result_file.write("NATURAL CONTACTS QUANTITY : " + str(nat_contact) + " - %"+ str(nat_contact*100/num)+ '\n')
    result_file.write("EVOL CONTACTS QUANTITY : " + str(evol_contact) + " - %"+ str(evol_contact*100/num)+ '\n')
    
    result_file.write("MATCH POSITIONS BETWEEN NAT AND EVOL (NO WINDOW) : " + str(len(data))+ '\n')
    result_file.write("MATCH POSITIONS  : " + str(data) + '\n')
    
    data_contact=[]
    for d in data:
        pos1 = int(d[0]-1)
        pos2 = int(d[1]-1)
        v=contact_map[pos1][pos2]
        if(v==1):
            data_contact.append(d)
    result_file.write("MATCH POSITIONS CONTACTS BETWEEN NAT AND EVOL (NO WINDOW) : " + str(len(data_contact))+ '\n')
    result_file.write("MATCH POSITIONS CONTACTS  : " + str(data_contact) + '\n')
    result_file.write("************************************************************************" + '\n')
    
    
    print "TOP : "  + str(top) + "% PAR POSITIONS : " + str(num)
    print "MATCH POSITIONS BETWEEN NAT AND EVOL (NO WINDOW) : " + str(len(data))
    print "NATURAL CONTACTS QUANTITY : " + str(nat_contact) + " - %"+ str(nat_contact*100/num)
    print "EVOL CONTACTS QUANTITY : " + str(evol_contact) + " - %"+ str(evol_contact*100/num)
    #print data    

'''
Toma todas las matrices de contactos de todas la proteinas evolucionadas,  retorna una matriz con la sumas y otra matriz con las probabilidad de los contactos.
'''
def sum_contact_map(family_folder,pdb_to_compare):
    family_folder_pdb = family_folder+"/PDB/"
    cmap_sum = None
    cant = 0 
    for index,pdb_protein_to_evolve in pdb_to_compare.iterrows():
        pdb_folder = family_folder_pdb + pdb_protein_to_evolve['pdb_folder_name']
        if(os.path.isdir(pdb_folder)): 
            contact_map = pdb_folder + "/contact_map_sync.dat"
            cmap = util.load_contact_map(contact_map)
            print contact_map
            if(cmap_sum==None):
                cmap_sum = cmap
                cant=cant+1
            else:
                if(cmap.shape==cmap_sum.shape):
                    cmap_sum = cmap_sum + cmap
                    cant=cant+1
                else:
                    print " diferent size natural : " + str(cmap_sum.shape)  + " AND " + contact_map + " : " + str(cmap.shape)
                    #pdb_to_evol_df=pdb_to_evol_df.drop(index)
    #print cmap_sum 
    util.save_contact_map(cmap_sum, family_folder + "/sum_contact_map.dat")
    cmap_sum = cmap_sum.astype(float)
    camp_prob = cmap_sum / cant
    print camp_prob   
    util.save_contact_map(camp_prob, family_folder + "/prob_contact_map.dat") 
    plot.contact_map(camp_prob,family_folder + "/prob_contact_map.png")
    conserved_contacts = np.count_nonzero(camp_prob == 1.0)  
    print conserved_contacts
    
def contact_map_sum_prob(execution_folder,contact_maps_paths):    
    cmap_sum = None
    cant = 0
    for contact_map_path in contact_maps_paths:
        cmap = util.load_contact_map(contact_map_path)
        #if(cmap_sum==None):
        if(cant==0):
            cmap_sum = cmap
            cant=cant+1
        else:
            if(cmap.shape==cmap_sum.shape):
                cmap_sum = cmap_sum + cmap
                cant=cant+1
            else:
                print " diferent size natural : " + str(cmap_sum.shape)  + " AND " + contact_map_path + " : " + str(cmap.shape)
                #pdb_to_evol_df=pdb_to_evol_df.drop(index)
    #print cmap_sum 
    util.save_contact_map(cmap_sum, execution_folder + "/sum_contact_map.dat")
    cmap_sum = cmap_sum.astype(float)
    camp_prob = cmap_sum / cant
    print camp_prob   
    util.save_contact_map(camp_prob, execution_folder + "/prob_contact_map.dat") 
    plot.contact_map(camp_prob,execution_folder + "/prob_contact_map.png")
    plot.contact_map_sum(cmap_sum,execution_folder + "/contact_map.png",'Sum Contact Map')
    
    df=pandas.DataFrame()
    '''conserved_contacts_100 = np.count_nonzero(camp_prob == 1.0)
    conserved_contacts_75_100 = np.count_nonzero(camp_prob > 0.75 and camp_prob < 1.0 )
    conserved_contacts_50_75 = np.count_nonzero(camp_prob > 0.5 and camp_prob < 0.75 )
    conserved_contacts_25_50 = np.count_nonzero(camp_prob > 0.25 and camp_prob < 0.5 )
    conserved_contacts_0_25 = np.count_nonzero(camp_prob > 0.0 and camp_prob < 0.25 )  
    '''
    
    total_contacts = np.count_nonzero(cmap_sum != 0)
    conserved_contacts_8 = np.count_nonzero(cmap_sum == 8)
    conserved_contacts_7 = np.count_nonzero(cmap_sum == 7)
    conserved_contacts_6 = np.count_nonzero(cmap_sum == 6)
    conserved_contacts_5 = np.count_nonzero(cmap_sum == 5)
    conserved_contacts_4 = np.count_nonzero(cmap_sum == 4)
    conserved_contacts_3 = np.count_nonzero(cmap_sum == 3)
    conserved_contacts_2 = np.count_nonzero(cmap_sum == 2)
    conserved_contacts_1 = np.count_nonzero(cmap_sum == 1)
    
    df.set_value(1, '#proteins', 1)
    df.set_value(1, '#contacts', conserved_contacts_1)
    df.set_value(1, '%contacts', conserved_contacts_1*100/float(total_contacts))
    df.set_value(2, '#proteins', 2)
    df.set_value(2, '#contacts', conserved_contacts_2)
    df.set_value(2, '%contacts', conserved_contacts_2*100/float(total_contacts))
    df.set_value(3, '#proteins', 3)
    df.set_value(3, '#contacts', conserved_contacts_3)
    df.set_value(3, '%contacts', conserved_contacts_3*100/float(total_contacts))
    df.set_value(4, '#proteins', 4)
    df.set_value(4, '#contacts', conserved_contacts_4)
    df.set_value(4, '%contacts', conserved_contacts_4*100/float(total_contacts))
    df.set_value(5, '#proteins', 5)
    df.set_value(5, '#contacts', conserved_contacts_5)
    df.set_value(5, '%contacts', conserved_contacts_5*100/float(total_contacts))
    df.set_value(6, '#proteins', 6)
    df.set_value(6, '#contacts', conserved_contacts_6)
    df.set_value(6, '%contacts', conserved_contacts_6*100/float(total_contacts))
    df.set_value(7, '#proteins', 7)
    df.set_value(7, '#contacts', conserved_contacts_7)
    df.set_value(7, '%contacts', conserved_contacts_7*100/float(total_contacts))
    df.set_value(8, '#proteins', 8)
    df.set_value(8, '#contacts', conserved_contacts_8)
    df.set_value(8, '%contacts', conserved_contacts_8*100/float(total_contacts))
    
    df.set_value(9, '#total_contacts', total_contacts)
    df.to_csv(execution_folder+'contact_distribution.csv')
        
"""
Lee la informacion sobre consevacion (KL) por columna de cada una de las proteinas evolucionadas. 
No se esta aplicando ningun entrecruzamiento de la informacion.
Solamente se esta ploteando la conservacion por columna para cada una de las proteinas (el grafico no puede apreciar resultados concretos) 
Luego si se realiza la media de conservacion resultado en conservation_media.png 
"""    
def comparative_conservation(family_folder, family_name, pdb_to_compare):
    natural_msa_conservation= family_folder + "/"+family_name+".fasta_data_kl.csv"
    family_folder_pdb = family_folder+"/PDB/"
    #[ name for name in os.listdir(thedir) if os.path.isdir(os.path.join(thedir, name)) ]
    msas_entropy=[]
    msa_entropy_media=[]
    df=msa.read_conservation(natural_msa_conservation)
    df = df.dropna()
    msa_entropy = [df['Entropy'].tolist(),family_name + "_NATURAL"]
    msas_entropy.append(msa_entropy)
    cant=0
    for index,pdb_protein_to_evolve in pdb_to_compare.iterrows():
        pdb_folder = family_folder_pdb + pdb_protein_to_evolve['pdb_folder_name']
        if(os.path.isdir(pdb_folder)):
            beta=str(pdb_protein_to_evolve['beta'])
            nsus=str(pdb_protein_to_evolve['nsus'])
            runs=str(int(pdb_protein_to_evolve['runs']))
            sufix = "sequences-beta"+beta+".0-nsus"+nsus+".0-runs"+runs+"_data_kl.csv"
            conservation_file = pdb_folder + "/sincronized_evol_path/conservation/"+sufix
            if(os.path.isfile(conservation_file)):
                df=msa.read_conservation(conservation_file)
                df = df.dropna()
                msa_entropy = [df['Entropy'].tolist(),pdb_folder]
                msas_entropy.append(msa_entropy)
                if(cant==0):
                    msa_entropy_media = df['Entropy'].tolist()
                else:
                    msa_entropy_media = [x + y for x, y in zip(msa_entropy_media , df['Entropy'].tolist() )]
                cant=cant+1
    msa_entropy_media =   [x / cant  for x in msa_entropy_media]     
    plot.conservation_between_msas(msas_entropy,family_folder + "/conservation.png")  
    #Media Graphic with natural
    msas_entropy=[]
    df=msa.read_conservation(natural_msa_conservation)
    df = df.dropna()
    msa_entropy = [df['Entropy'].tolist(),family_name + "_NATURAL"]
    msas_entropy.append(msa_entropy)
    msas_entropy.append([msa_entropy_media,"MEDIA"])
    plot.conservation_between_msas(msas_entropy,family_folder + "/conservation_media.png") 
    #TODO persistir conservacion media y natura por posicion en un archivo aparte
    

'''
Compute joined msas
'''
def compute_joined_msas(family_folder,pdb_to_compare,contact_map_path):
    logging.info('compute_joined_msas :: Begin ')
    start_time = time.time()
    joined_path =  family_folder + "/joined/"
    name = "joined_evol_msa"
    joined_fasta_path = joined_path + name +  ".fasta"
    if not os.path.exists(joined_path):
        os.makedirs(joined_path)
    fasta_files = []
    '''for index,pdb_protein_to_evolve in pdb_to_compare.iterrows():
        pdb_folder = family_folder + "/PDB/" +  pdb_protein_to_evolve['pdb_folder_name']
        if(os.path.isdir(pdb_folder)):
            beta=str(pdb_protein_to_evolve['beta'])
            nsus=str(pdb_protein_to_evolve['nsus'])
            runs=str(int(pdb_protein_to_evolve['runs']))
            sufix = "sequences-beta"+beta+".0-nsus"+nsus+".0-runs"+runs+".fasta"
            msa_file = pdb_folder + "/sincronized_evol_path/"+sufix
            if(os.path.isfile(msa_file)):
                fasta_files.append(msa_file)
    with open(joined_fasta_path, "w") as joined_fasta:
        for fasta in fasta_files:
            logging.info('Attach MSA  ' + fasta )
            with open(fasta) as infile:
                for line in infile:
                    if('>' not in line):
                        joined_fasta.write(line)
            infile.close() 
    joined_fasta.close()
    logging.info('End of Attach evolutionated MSAs  ' )
    '''
    #msa.clustering_singular("0.70",joined_fasta_path, joined_fasta_clustering_path)
    
    columns=["run_id","auc_threshold_0","auc_01_threshold_0","auc_threshold_25","auc_01_threshold_25","auc_threshold_50","auc_01_threshold_50","auc_threshold_75","auc_01_threshold_75"]
    df = pandas.DataFrame(columns=columns)
    for i  in xrange(1,10):
        random_fasta = joined_fasta_path + "_"+str(i)
        msa.random_seq(joined_fasta_path, random_fasta,10000)
        mi_data_output_path = random_fasta+ ".csv"
        evol_analisys(random_fasta, mi_data_output_path, random_fasta+"_conservation" , name)
        df.set_value(i, 'run_id', i)
        target,scores=getTargetScores(mi_data_output_path,contact_map_path,1,0)
        auc,auc01 = util.getAUC(target,scores)
        df.set_value(i, 'auc_threshold_0', auc)
        df.set_value(i, 'auc_01_threshold_0', auc01)
        target,scores=getTargetScores(mi_data_output_path,contact_map_path,1,0.25)
        auc,auc01 = util.getAUC(target,scores)
        df.set_value(i, 'auc_threshold_25', auc)
        df.set_value(i, 'auc_01_threshold_25', auc01)
        target,scores=getTargetScores(mi_data_output_path,contact_map_path,1,0.50)
        auc,auc01 = util.getAUC(target,scores)
        df.set_value(i, 'auc_threshold_50', auc)
        df.set_value(i, 'auc_01_threshold_50', auc01)
        target,scores=getTargetScores(mi_data_output_path,contact_map_path,1,0.75)
        auc,auc01 = util.getAUC(target,scores)
        df.set_value(i, 'auc_threshold_75', auc)
        df.set_value(i, 'auc_01_threshold_75', auc01)
    df.to_csv(joined_fasta_path+"_poblation.csv")        
    '''
        buslje09(clustered_sequences_path,mi_data_path)
        target,scores=dataanalisys.getTargetScores(mi_data_path,contact_map_path,window)
        auc,auc01 = util.getAUC(target,scores)
        '''
    #target,scores=dataanalisys.getTargetScores(mi_data_path,contact_map_sync,window)
    logging.info('End of the execution process compute_joined_msas')


def analisys_mi_with_contact_map(execution_folder, mi_paths,contact_map_path, zmip_natural_result_path,top_mi , window, sinchronize_with_natural=True):
    contact_map = util.load_contact_map(contact_map_path)
    
    zmip_natural = util.load_zmip(zmip_natural_result_path,window)
    util.order(zmip_natural)
    
    fields=["Position1","Position2","Count"]
    df_total = pandas.DataFrame([],columns=fields)
    cant = 0
    for mi_file_path in mi_paths:
        cant = cant + 1
        zmip_evol = util.load_zmip(mi_file_path,window)
        
        if(sinchronize_with_natural==True):
            zmip_natural,zmip_evol=util.sincronice_mi(zmip_natural, zmip_evol)
        
        util.order(zmip_evol)
        num = len(zmip_evol)*top_mi//100
        zmip_evol_top=zmip_evol[0:int(num)]
        df = pandas.DataFrame(zmip_evol_top,columns=fields)
        df_total=df_total.append(df)
    
    #After append all the evolution MI TOP do this
    #counts = dfres.groupby(['Position1','Position2']).size()
    #print counts
    #Count cantidad de veces que aprece en los tops
    counts_df = pandas.DataFrame(df_total.groupby(['Position1','Position2']).size().rename('Count'))
    #
    sorted_df=counts_df.sort_values(by=['Count'],ascending=[False])
    #print sorted_df
    
    sorted_df['ProbTop']=pandas.Series(0.0, index=sorted_df.index)
    sorted_df['Contacts']=pandas.Series(0.0, index=sorted_df.index)
    sorted_df['ProbContact']=pandas.Series(0.0, index=sorted_df.index)
    #print prob_contact_map
    for index,mi_par in sorted_df.iterrows():
        #probablidad de que aprezca en top
        prob_top = mi_par['Count'] * 100 / cant 
        sorted_df.set_value(index, 'ProbTop' , prob_top/100)
        #cantidad de contactos 
        pos1 = int(index[0]-1)
        pos2 = int(index[1]-1)
        v=contact_map[pos1][pos2]
        sorted_df.set_value(index, 'Contacts' , v)
        #probablidad de contactos
        prob_contact = float(v) * 100 / cant
        sorted_df.set_value(index, 'ProbContact' , float(prob_contact)/100)
        
    
    sorted_df.to_csv(execution_folder + '/top_'+str(top_mi)+'_mi.csv', sep='\t', encoding='utf-8')
    
    
    
    correlation_p = sorted_df['Count'].corr(sorted_df['Contacts'], method='pearson')
    correlation_k = sorted_df['Count'].corr(sorted_df['Contacts'], method='kendall')
    correlation_s = sorted_df['Count'].corr(sorted_df['Contacts'], method='spearman')
    
    '''correlation_sp = df.corr(method='pearson')
    correlation_sp = df.corr(method='kendall')
    correlation_sp = df.corr(method='spearman')
    '''
    
    
    mean = sorted_df["ProbTop"].mean()
    median = sorted_df["ProbTop"].median()
    var = sorted_df["ProbTop"].var()
    mode = sorted_df["ProbTop"].mode()
    
    corr_df = pandas.DataFrame()
    corr_df.set_value(1, 'pearson', correlation_p)
    corr_df.set_value(1, 'kendall', correlation_k)
    corr_df.set_value(1, 'spearman', correlation_s)
    corr_df.to_csv(execution_folder + '/top_'+str(top_mi)+'_mi_correlation.csv', sep='\t', encoding='utf-8')
    
    import matplotlib.pyplot as plt
    sorted_df.plot.scatter(x='Count', y='Contacts');
    plt.xlabel('Numero de veces que aparece en el top en la MI Conjunta')
    plt.ylabel('Numero de contactos de la matriz conjunta')
    plt.xticks([0,1,2,3,4,5,6,7,8])
    plt.yticks([0,1,2,3,4,5,6,7,8])
    plt.savefig(execution_folder + '/top_'+str(top_mi)+'_mi_correlation.png');
    plt.gcf().clear()
    
    
            
"""
Esta funcion toma el top de MI de todas las proteinas evolucionadas y luego realiza una agrupacion indicando la cantidad de veces que aparecen los pares.
Ordena los pares de forma descendente, osea los pares que mas aparecen en el top quedan arriba. 
Ademas se agrega la columna indicando la probabilidad de contacto que existen entre ellos.
"""
def comparative_mi_information(family_folder, family_name,top, window, pdb_to_compare):     
    import matplotlib.pyplot as plt      
    logging.info('Begin of the execution process family MI information')
    family_folder_pdb = family_folder+"/PDB/"
    #for protein_pdb in os.listdir(family_folder_pdb):
    fields=["Position1","Position2","Count"]
    df_total = pandas.DataFrame([],columns=fields)
    cant = 0
    for index,pdb_protein_to_evolve in pdb_to_compare.iterrows():
        pdb_folder = family_folder_pdb + pdb_protein_to_evolve['pdb_folder_name']
        if(os.path.isdir(pdb_folder)): 
            #buscar un archivo en particular directamente
            beta=str(pdb_protein_to_evolve['beta'])
            nsus=str(pdb_protein_to_evolve['nsus'])
            runs=str(int(pdb_protein_to_evolve['runs']))
            sufix = "sequences-beta"+beta+".0-nsus"+nsus+".0-runs"+runs
            zmip_file = pdb_folder + "/mi_data/zmip_"+sufix+".csv"
            if(os.path.isfile(zmip_file)):
                cant = cant + 1
                zmip_evol = util.load_zmip(zmip_file,window)
                util.order(zmip_evol)
                num = len(zmip_evol)*top/100
                zmip_evol_top=zmip_evol[0:int(num)]
                df = pandas.DataFrame(zmip_evol_top,columns=fields)
                df_total=df_total.append(df)
            #print df_total
    print df_total
    
    """a=[14.0,76.0, 45.2345]
    b=[23.0,54.0, 34.5]
    c=[[14.0,76.0, 45.2345],[23.0,54.0, 34.5],[45.0,90.0, 34.5]]
    #c.append(a)
    #c.append(b)
    df = pandas.DataFrame(c,columns=fields)
    d=[[11.0,76.0, 45.2345],[23.0,54.0, 34.5],[45.0,90.0, 34.5]]
    #d.append(a)
    #d.append(b)
    df2 = pandas.DataFrame(d,columns=fields)
    dfres=df.append(df2)
    e=[[10.0,76.0, 45.2345],[23.0,54.0, 34.5],[45.0,90.0, 34.5]]
    dfe = pandas.DataFrame(e,columns=fields)
    dfres=dfres.append(dfe)
   
    """
    #After append all the evolution MI TOP do this
    #counts = dfres.groupby(['Position1','Position2']).size()
    #print counts
    counts_df = pandas.DataFrame(df_total.groupby(['Position1','Position2']).size().rename('Count'))
    #print counts_df
    sorted_df=counts_df.sort_values(by=['Count'],ascending=[False])
    #print sorted_df
    sorted_df['ProbContact']=pandas.Series(0.0, index=sorted_df.index)
    sorted_df['ProbTop']=pandas.Series(0.0, index=sorted_df.index)
    prob_contact_map = util.load_contact_map(family_folder + "/prob_contact_map.dat",np.float64)
    #print prob_contact_map
    for index,mi_par in sorted_df.iterrows():
        #por bug arreglar
        #print mi_par
        pos1 = int(index[0]-1)
        pos2 = int(index[1]-1)
        v=prob_contact_map[pos1][pos2]
        sorted_df.set_value(index, 'ProbContact' , v)
        prob_top = mi_par['Count'] * 100 / cant 
        sorted_df.set_value(index, 'ProbTop' , prob_top/100)
        #sorted_df[index]['ProbContact']=v
    
    
    sorted_df.to_csv(family_folder + "/top_family_mi.csv", sep='\t', encoding='utf-8')
    
    
    correlation_p = sorted_df['ProbContact'].corr(sorted_df['ProbTop'], method='pearson')
    correlation_k = sorted_df['ProbContact'].corr(sorted_df['ProbTop'], method='kendall')
    correlation_s = sorted_df['ProbContact'].corr(sorted_df['ProbTop'], method='spearman')
    
    '''correlation_sp = df.corr(method='pearson')
    correlation_sp = df.corr(method='kendall')
    correlation_sp = df.corr(method='spearman')
    '''
    mean = sorted_df["ProbTop"].mean()
    median = sorted_df["ProbTop"].median()
    var = sorted_df["ProbTop"].var()
    mode = sorted_df["ProbTop"].mode()
    sorted_df.plot.scatter(x='ProbTop', y='ProbContact');
    plt.savefig(family_folder + "/top_family_mi.png");
    plt.show()
    plt.gcf().clear()
    print sorted_df
    
def prom_zmip(mi_paths, mi_prom_result,window):
    fields=["position1","position2","zmip"]
    df_total = pandas.DataFrame([],columns=fields)
    cant = 0
    for mi_file_path in mi_paths:
        cant = cant + 1
        zmip_evol = util.load_zmip(mi_file_path,window)
        df = pandas.DataFrame(zmip_evol,columns=fields)
        df_total=df_total.append(df)
    
    sum_df = pandas.DataFrame(df_total.groupby(['position1', 'position2'])['zmip'].sum().rename('zmip'))
    sum_df['zmip'] = sum_df['zmip'] / cant
    sum_df.to_csv(mi_prom_result, sep='\t', encoding='utf-8', header=False)
               
'''
Not in use.
'''    
def kendall(x,y):
    print 1 - Bio.Cluster.distancematrix((x,y), dist="k")[1][0]
    
    
    
def top_rank_comparation(execution_folder, top,contact_map_path,zmip_natural_result_path, zmip_evol_intersect_result_path, top_df, zmip_reference_result_path, zmip_prom_result_path, index, window, contact_threshold, top_threshold, sinchronize_with_natural=True):
    #matriz de contacto
    contact_map= util.load_contact_map(contact_map_path)
    
    zmip_natural = util.load_zmip(zmip_natural_result_path,window)
    util.order(zmip_natural)
    
    #pares de posiciones 
    num = len(zmip_natural)*top//100
    num = int(num)
    natural=zmip_natural[0:num]
    
    df = pandas.read_csv(zmip_evol_intersect_result_path,delim_whitespace=True,header=0,usecols=[0,1,2,4])
    df=df.sort(['Count', 'Contacts'], ascending=[False, False])
    evol=df.head(n=num)
    evol=evol.values.tolist()
    
    evol_num=num
    if(num>len(df.index)):
        evol_num=len(df.index)
    
    data=matches_coevolved_positions(natural,evol,intersection_evol=True)
    
    zmip_reference = util.load_zmip(zmip_reference_result_path,window)
    zmip_prom = util.load_zmip(zmip_prom_result_path,window)
    if(sinchronize_with_natural==True):
        zmip_natural,zmip_reference=util.sincronice_mi(zmip_natural, zmip_reference)
        zmip_natural,zmip_prom=util.sincronice_mi(zmip_natural, zmip_prom)
    
    util.order(zmip_reference)
    util.order(zmip_prom)
    
    reference=zmip_reference[0:num]
    prom=zmip_prom[0:num]
    data_reference=matches_coevolved_positions(natural,reference,intersection_evol=True)
    data_prom=matches_coevolved_positions(natural,prom,intersection_evol=True)
    
    x_nat=[]
    y_nat=[]
    x_evol=[]
    y_evol=[]
    x_ref=[]
    y_ref=[]
    x_prom=[]
    y_prom=[]
    
    
    
    evol_contact_pairs=[]
    ref_contact_pairs=[]
    prom_contact_pairs=[]
    nat_contact = 0
    evol_contact = 0
    ref_contact = 0
    prom_contact = 0
    for x, y, z, p in map(None, natural, evol,reference,prom):
        '''x_nat.append(int(x[0]-1))
        y_nat.append(int(x[1]-1))
        x_evol.append(int(y[0]-1))
        y_evol.append(int(y[1]-1))
        x_ref.append(int(z[0]-1))
        y_ref.append(int(z[1]-1))
        
        x_prom.append(int(p[0]-1))
        y_prom.append(int(p[1]-1))
        '''
        
        
        pos1 = int(x[0]-1)
        pos2 = int(x[1]-1)
        v = contact_map[pos1][pos2]
        if(v >= contact_threshold):
            nat_contact=nat_contact+1
        
        
        if (y != None):
            pos1 = int(y[0]-1)
            pos2 = int(y[1]-1)
            v = contact_map[pos1][pos2]
            if(v >= contact_threshold):
                evol_contact=evol_contact+1 
                evol_contact_pairs.append(y)
        
        pos1 = int(z[0]-1)
        pos2 = int(z[1]-1)
        v = contact_map[pos1][pos2]
        if(v >= contact_threshold):
            ref_contact=ref_contact+1 
            ref_contact_pairs.append(z)
        
        pos1 = int(p[0]-1)
        pos2 = int(p[1]-1)
        v = contact_map[pos1][pos2]
        if(v >= contact_threshold):
            prom_contact=prom_contact+1 
            prom_contact_pairs.append(p)
        
        
    #plot.contact_map_with_top_rank_mi_sum(contact_map,  x_nat, y_nat, x_evol,y_evol,execution_folder + 'top_'+ top +'contact_map_thresold_'+str(contact_threshold)+'with_intersection_top_threshold_'+str(top_threshold)+'.png','', 'top_'+ top +'THIO_ECOLI Top Thresold ' + str(top_threshold) +' Contact Threshold ' + str(contact_threshold))
    #plot.contact_map_with_top_rank_mi_sum(contact_map,  x_nat, y_nat, x_ref,y_ref,execution_folder + 'top_'+ top +'contact_map_thresold_'+str(contact_threshold)+'with_reference_top_threshold_'+str(top_threshold)+'.png','', 'top_'+ top +'THIO_ECOLI 2TRX Thresold ' + str(top_threshold) +' Contact Threshold ' + str(contact_threshold))
    #plot.contact_map_with_top_rank_mi_sum(contact_map,  x_nat, y_nat, x_prom,y_prom,execution_folder + 'top_'+ top +'contact_map_thresold_'+str(contact_threshold)+'with_prom_top_threshold_'+str(top_threshold)+'.png','', 'top_'+ top +'THIO_ECOLI Prom Thresold ' + str(top_threshold) +' Contact Threshold ' + str(contact_threshold))
    
    top_df.set_value(index,'top',top)
    top_df.set_value(index,'contact_threshold',str(contact_threshold))
    top_df.set_value(index,'par_positions',str(num))
    top_df.set_value(index,'nat_contact',str(nat_contact))
    top_df.set_value(index,'nat_contact_%',str(nat_contact*100/num))
    top_df.set_value(index,'evol_contact',str(evol_contact))
    top_df.set_value(index,'evol_contact_%',str(evol_contact*100/num))
    
    data_contact=[]
    for d in data:
        pos1 = int(d[0]-1)
        pos2 = int(d[1]-1)
        v=contact_map[pos1][pos2]
        if(v>=contact_threshold):
            data_contact.append(d)
    top_df.set_value(index,'match_positions',str(len(data_contact)))
    
    top_df.set_value(index,'prom_contact',str(prom_contact))
    top_df.set_value(index,'prom_contact_%',str(prom_contact*100/num))
    data_contact_prom=[]
    for d in data_prom:
        pos1 = int(d[0]-1)
        pos2 = int(d[1]-1)
        v=contact_map[pos1][pos2]
        if(v>=contact_threshold):
            data_contact_prom.append(d)
    top_df.set_value(index,'match_positions_prom',str(len(data_contact_prom)))
    
    top_df.set_value(index,'ref_contact',str(ref_contact))
    top_df.set_value(index,'ref_contact_%',str(ref_contact*100/num))
    
    data_contact_ref=[]
    for d in data_reference:
        pos1 = int(d[0]-1)
        pos2 = int(d[1]-1)
        v=contact_map[pos1][pos2]
        if(v>=contact_threshold):
            data_contact_ref.append(d)
    top_df.set_value(index,'match_positions_reference',str(len(data_contact_ref)))
    
    
    
    
def generate_contact_map_with_top_mi_matrix(top,sum_contact_map, top_smi_path, out_matrix):
    cmap = util.load_contact_map(sum_contact_map) 
    cmap_with_mi=np.triu(cmap, -1)
    
    df = pandas.read_csv(top_smi_path,delim_whitespace=True,header=0,usecols=[0,1,2,4])
    #df=df.sort(['Count', 'Contacts'], ascending=[False, False])
    df=df.sort_values(by=['Count', 'Contacts'],ascending=[False, False])
    
    top_smi=df.values.tolist()
    
    #add mi values to map 
    for x in map(None, top_smi):
        pos1 = int(x[0]-1)
        pos2 = int(x[1]-1)
        v = cmap_with_mi[pos2][pos1]
        if(v!=0):
            print ("error")
        #set the value of the smi    
        cmap_with_mi[pos2][pos1]=x[2]
    util.save_contact_map(cmap_with_mi, out_matrix)    
    plot.contact_map_sum(cmap_with_mi,'../THIO_ECOLI_4_107/contact_map_with_mi_top_'+top+'.png','Sum Contact Map and Sum Top '+top+' MI')    
    
    
    
def generate_contact_map_with_top_mi_two_separated_matrix(top,sum_contact_map, top_smi_path, out_matrix):
    cmap = util.load_contact_map(sum_contact_map) 
    cmap =np.triu(cmap, -1)
    mat_mi = np.zeros(cmap.shape)
    
    df = pandas.read_csv(top_smi_path,delim_whitespace=True,header=0,usecols=[0,1,2,4])
    #df=df.sort(['Count', 'Contacts'], ascending=[False, False])
    df=df.sort_values(by=['Count', 'Contacts'],ascending=[False, False])
    top_smi=df.values.tolist()
    
    #add mi values to map 
    for x in map(None, top_smi):
        pos1 = int(x[0]-1)
        pos2 = int(x[1]-1)
        v = mat_mi[pos2][pos1]
        if(v!=0):
            print ("error")
        #set the value of the smi    
        mat_mi[pos2][pos1]=x[2]
    
    print  mat_mi
    
    #util.save_contact_map(cmap_with_mi, out_matrix)    
    plot.contact_map_sum_top_mi_matrix(cmap, mat_mi ,'../THIO_ECOLI_4_107/contact_map_with_mi_top_'+top+'_two_color_bars.png','Sum Contact Map and Sum Top 1 MI') 
    

def dendogram_matrix(contact_maps_paths,output_path,title,structures,clustering_type):
    Y=[]
    for cmap_pth in contact_maps_paths:
        contact_map = util.load_contact_map(cmap_pth)
        contact_map=contact_map.ravel()
        Y.append(contact_map)
    Z = linkage(Y, clustering_type)
    plot.dendogram_matrix(Z,output_path,title,structures)
    
def dendogram_top_mi(top, zmip_nat,zmip_paths, output_path,title ,structures, clustering_type, shape ):
    Y=[]
    for zmip_path in zmip_paths:
        zmip = util.load_zmip(zmip_path,window)
        zmip_natural,zmip=util.sincronice_mi(zmip_natural, zmip)
        util.order(zmip)
        #pares de posiciones 
        num = len(zmip)*top//100
        num = int(num)
        zmip=zmip[0:num]   
        cmap_with_mi=np.zeros(shape)
        for x in map(None, zmip):
            pos1 = int(x[0]-1)
            pos2 = int(x[1]-1)
            cmap_with_mi[pos2][pos1]=1
        print zmip
        print "-------"
        print cmap_with_mi
        cmap_with_mi=cmap_with_mi.ravel()
        Y.append(cmap_with_mi)
    Z = linkage(Y, clustering_type)
    plot.dendogram_matrix(Z,output_path,title,structures)   
        
            