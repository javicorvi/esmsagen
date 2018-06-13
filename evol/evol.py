'''
Evolution module
'''
import glob
import gzip
import scpe
import urllib
import msa_analysis 
import msa 
import plot 
import util 
import os
import time
import constants as const
import pdb
import pandas
import logging 


import constants
from evol import msa_analysis

logging.basicConfig(filename=const.log_file_path, level=logging.DEBUG, format="%(asctime)s - %(levelname)s - %(filename)s - %(funcName)s - %(message)s ")
consoleHandler = logging.StreamHandler()
rootLogger = logging.getLogger()
rootLogger.addHandler(consoleHandler)

'''
Family Evolution
'''
execute_family_evol = False

'''
Calculate the MI for the natural MSA putting the protein as the reference
'''
execute_natural_mi_msa = False
'''
Calculate the Conservation of the families natural MSA
'''
execute_msa_natural_information = False

execute_download_pdbs = False

optimized_scpe_variables = False

'''
Execute the evolution of the protein with SCPE.
Generate several families; taking account de parameters beta, nsus and runs 
'''
execute_scpe = False

'''
Execute the clustering for the families generated with SCPE to avoid redundant sequences with high identity
'''
execute_clustering = False
'''
Recorta el msa
'''
execute_cut_msa = False
'''
Execute the analisys of the MSA: Seq Logo.
'''
execute_msa_information = False
'''
Execute the MI calcultation busjle09
'''
execute_mi = False
'''
Execute the analisys of the information
'''
execute_msa_analysis = False
'''
Execute the analisys of the information between all the PDBS and MSA generated. All together
'''
execute_joined_pdb_analisys = True

'''
Pattern to execute process0
'''
pattern = ["sequences"]

'''
Iterates over the structures, pdbs and execute the scpe and the clusterization 
'''        
input_families_folder = "../FAMILY_PF00085/"



'''
Main method for optimize the evolution of the structure, returns the optimization information.
Run the evolution with diferents parameters.
'''          
def optimize_evolution(execution_folder, protein, pdb_name, chain, synchronized=False):
    logging.info('Run Optimization For ' + pdb_name)
    start_time_total = time.time()
    columns = ["protein", "pdb", "chain", "beta", "nsus", "run", "auc_mi", "auc_01_mi", "auc_di", "auc_01_di", "auc_frob", "auc_01_frob","auc_psicov","auc_01_psicov","execution_time", "start_residue", "end_residue","message"]
    df = pandas.DataFrame(columns=columns)
    optimization_folder = execution_folder + "optimization/"
    optimization_file_path = optimization_folder + "optimization.csv"
    scpe_sequences = optimization_folder + "scpe_sequences/"
    clustered_sequences_path = optimization_folder + "clustered_sequences/"
    curated_sequences_path = optimization_folder + "curated_sequences/"
    mi_data_path = clustered_sequences_path + "mi/"
    di_data_path = clustered_sequences_path + "di/"
    frob_data_path = clustered_sequences_path + "frob/"
    psicov_path = clustered_sequences_path + "psicov/"
    if(synchronized==True):
        mi_data_path = curated_sequences_path + "mi/"
        di_data_path = curated_sequences_path + "di/"
        frob_data_path = curated_sequences_path + "frob/"
        psicov_path = curated_sequences_path + "psicov/"
        
    contact_map = optimization_folder + "contact_map.dat"
    contact_map_sync = optimization_folder + "contact_map_sync.dat" 
    pdb_path = execution_folder + pdb_name + ".pdb"
    pdb_clean_path = execution_folder +  pdb_name + "_chain_" + chain +".pdb"
    #create all the structure for the optimization
    if not os.path.exists(scpe_sequences):
        os.makedirs(scpe_sequences)
    if not os.path.exists(clustered_sequences_path):
        os.makedirs(clustered_sequences_path)
    if not os.path.exists(curated_sequences_path):
        os.makedirs(curated_sequences_path)
    if not os.path.exists(mi_data_path):
        os.makedirs(mi_data_path) 
    if not os.path.exists(di_data_path):
        os.makedirs(di_data_path)
    if not os.path.exists(frob_data_path):
        os.makedirs(frob_data_path)         
    if not os.path.exists(psicov_path):
        os.makedirs(psicov_path)      
    #hacer condicion para download
    #pdb.download(pdb_name, pdb_path)
    #pdb.clean_pdb(pdb_path, pdb_clean_path, chain)
    
    beta = ["0.001","0.01","0.1","0.5", "1.0", "2.0", "3.0", "5.0", "7.0", "10.0", "15.0", "20.0"]
    runs = ["1000", "5000", "10000", "20000"]
    nsus = ["1.0", "2.0", "3.0", "5.0", "7.0", "10.0", "15.0","20.0"]
    index = 1
    for b in beta:
        for sus in nsus:
            for r in runs:
                start_time = time.time()
                logging.info('Calculation of beta ' + b + ' nsus ' + sus + ' run ' + r)
                df.at[index,'protein']=protein
                df.at[index,'pdb']=pdb_name
                df.at[index,'chain']=chain
                #df.at[index,'start_residue']=start_residue
                #df.at[index,'end_residue']=end_residue
                df.at[index,'beta']=b
                df.at[index,'nsus']=sus
                df.at[index,'run']=r
                try: 
                    auc, auc01, auc_di, auc01_di, auc_frob, auc01_frob, auc_psicov, auc01_psicov, count_seq_scpe, count_seq_cluster, seq_lenght = evol(pdb_clean_path, b, r, sus, chain, scpe_sequences, clustered_sequences_path, curated_sequences_path, mi_data_path, di_data_path, frob_data_path, psicov_path,contact_map, contact_map_sync,True,synchronized)
                    df.at[index,'auc_mi']=auc
                    df.at[index,'auc_01_mi']=auc01
                    df.at[index,'auc_di']=auc_di
                    df.at[index,'auc_01_di']=auc01_di
                    df.at[index,'auc_frob']=auc_frob
                    df.at[index,'auc_01_frob']=auc01_frob
                    df.at[index,'auc_psicov']=auc_psicov
                    df.at[index,'auc_01_psicov']=auc01_psicov
                    df.at[index,'count_seq_scpe']=count_seq_scpe
                    df.at[index,'count_seq_cluster']=count_seq_cluster
                    df.at[index,'seq_lenght']=seq_lenght
                except Exception as inst:
                    print inst
                    x = inst.args
                    print x
                    df.at[index,'message']="error_running_optimization"
                    logging.error('Error with beta ' + b + ' nsus ' + sus + ' run ' + r)        
                df.at[index,'execution_time']= time.time() - start_time
                index = index + 1
                df.to_csv(optimization_file_path)    
    df['execution_time_optimization_total'] = time.time() - start_time_total            
    df.to_csv(optimization_file_path)                
   

'''
Metodo de evolucion para una estructura con unos parametros de SCPE indicados, luego se realiza todas los calculos referidos a la coevolucion con los diferentes metodos: mi, di, frob. Los ultimos dos 
Gaussianos.

'''
def evol(pdb_file_complete_filename_to_evolve, beta, runs, nsus, chain, scpe_sequences, clustered_sequences_folder_path, curated_sequences_path, mi_data_path, di_path, frob_path, psicov_path,contact_map_path, contact_map_sync,coevolution_calc=True,syncronized=False):    
    sufix = "sequences-beta" + beta + "-nsus" + nsus + "-runs" + runs
    file_name = sufix + ".fasta"
    output_msa_path = scpe_sequences + file_name
    clustered_sequences_path = clustered_sequences_folder_path + file_name + ".cluster"
    clustered_tmp_sequences_path = clustered_sequences_path + ".clstr"
    curated_sequences_path = curated_sequences_path + file_name 
    mi_data_path = mi_data_path + "mi_" + sufix + ".csv"
    dca_di = di_path + "di_" + sufix + ".csv" 
    dca_Frob = frob_path + "frob_" + sufix + ".csv"
    psicov_path = psicov_path + "psicov_" + sufix + ".csv"
    scpe.run(pdb_file_complete_filename_to_evolve, beta, runs, nsus, chain, output_msa_path, contact_map_path)
    count_seq_scpe = msa.count_sequences(output_msa_path)
    msa.clustering_singular("0.62", output_msa_path, clustered_sequences_path)
    count_seq_cluster = msa.count_sequences(clustered_sequences_path)
    seq_lenght = msa.count_aminoacids(clustered_sequences_path)
    util.delete_files(output_msa_path)
    util.delete_files(clustered_tmp_sequences_path)
    if(coevolution_calc):
        if(syncronized):
            util.sincronize_natural_evol_msa(clustered_sequences_path, curated_sequences_path, pattern, 2, -3)
            msa.buslje09(curated_sequences_path, mi_data_path)
            msa.gaussDcaFrobenius(curated_sequences_path, dca_Frob)
            msa.gaussDcaDirectInformation(curated_sequences_path, dca_di)
            msa.create_msa_without_id(curated_sequences_path, curated_sequences_path+"_noid.aln")
            msa.psicov(curated_sequences_path+"_noid.aln",psicov_path)
            util.sincronize_contact_map(contact_map_path, contact_map_sync, 2, 106)
            contact_map_path = contact_map_sync
        else:
            msa.buslje09(clustered_sequences_path, mi_data_path)
            msa.gaussDcaFrobenius(clustered_sequences_path, dca_Frob)
            msa.gaussDcaDirectInformation(clustered_sequences_path, dca_di)
            msa.create_msa_without_id(clustered_sequences_path, clustered_sequences_path+"_noid.aln")
            msa.psicov(clustered_sequences_path+"_noid.aln",psicov_path)
        # target,scores=msa_analysis.getTargetScores(mi_data_path,contact_map_sync,constants.neighbour_window)
        
        target, scores = msa_analysis.getTargetScores(mi_data_path, contact_map_path, constants.neighbour_window)
        auc, auc01 = util.getAUC(target, scores)
        
        
        target_di, scores_di = msa_analysis.getTargetScores(dca_di, contact_map_path, constants.neighbour_window, split_char=' ')
        auc_di, auc01_di = util.getAUC(target_di, scores_di)
        
        target_frob, scores_frob = msa_analysis.getTargetScores(dca_Frob, contact_map_path, constants.neighbour_window, split_char=' ')
        auc_frob, auc01_frob = util.getAUC(target_frob, scores_frob)
       
        target_psicov, scores_psicov = msa_analysis.getTargetScores(psicov_path, contact_map_path, constants.neighbour_window, split_char=' ',score_position=4)
        auc_psicov, auc01_psicov = util.getAUC(target_psicov, scores_psicov)
        
        return auc, auc01, auc_di, auc01_di, auc_frob, auc01_frob, auc_psicov, auc01_psicov, count_seq_scpe, count_seq_cluster, seq_lenght


def natural_information(execution_folder):
    optimization_folder = execution_folder + "optimization/"
    
    contact_map = optimization_folder + "contact_map.dat"
    contact_map_sync = optimization_folder + "contact_map_sync.dat" 
    util.sincronize_contact_map(contact_map, contact_map_sync, 2, 106)   
     
    #Natural information TODO set reference protein to msa.
    natural_result_path = execution_folder + "/natural/"
    natural_results = natural_result_path + "resuts.csv"
    columns = ["protein", "pdb", "chain", "auc_mi", "auc_01_mi", "auc_di", "auc_01_di", "auc_frob", "auc_01_frob","auc_psicov","auc_01_psicov","message"]
    
    natural_results_df = pandas.DataFrame(columns=columns)
    if not os.path.exists(natural_result_path):
        os.makedirs(natural_result_path)
    natural_msa = execution_folder + "PF00085_THIO_ECOLI_reference.fasta"
    natural_mi_result_path = natural_result_path + "mi.csv"
    natural_frob_result_path = natural_result_path + "frob.csv"
    natural_di_result_path= natural_result_path + "di.csv"
    natural_psicov_result_path = natural_result_path + "psicov.csv"
    
    
    msa.buslje09(natural_msa, natural_mi_result_path)
    msa.gaussDcaFrobenius(natural_msa, natural_frob_result_path)
    msa.gaussDcaDirectInformation(natural_msa, natural_di_result_path)
    msa.create_msa_without_id(natural_msa, natural_msa+"_noid.aln")
    msa.psicov(natural_msa+"_noid.aln",natural_psicov_result_path)
    
    
    
    
    target, scores = msa_analysis.getTargetScores(natural_mi_result_path, contact_map_sync, constants.neighbour_window)
    auc, auc01 = util.getAUC(target, scores)
    target_di, scores_di = msa_analysis.getTargetScores(natural_di_result_path, contact_map_sync, constants.neighbour_window, split_char=' ')
    auc_di, auc01_di = util.getAUC(target_di, scores_di)
    target_frob, scores_frob = msa_analysis.getTargetScores(natural_frob_result_path, contact_map_sync, constants.neighbour_window, split_char=' ')
    auc_frob, auc01_frob = util.getAUC(target_frob, scores_frob)
    target_psicov, scores_psicov = msa_analysis.getTargetScores(natural_psicov_result_path, contact_map_sync, constants.neighbour_window, split_char=' ',score_position=4)
    auc_psicov, auc01_psicov = util.getAUC(target_psicov, scores_psicov)
    
    natural_results_df.at[1,'auc_mi']=auc
    natural_results_df.at[1,'auc_01_mi']=auc01
    natural_results_df.at[1,'auc_di']=auc_di
    natural_results_df.at[1,'auc_01_di']=auc01_di
    natural_results_df.at[1,'auc_frob']=auc_frob
    natural_results_df.at[1,'auc_01_frob']=auc01_frob
    natural_results_df.at[1,'auc_psicov']=auc_psicov
    natural_results_df.at[1,'auc_01_psicov']=auc01_psicov
    
    natural_results_df.to_csv(natural_results)
    
    #Information to plot natural rocs
    colors = ['green','blue','yellow','red']
    labels = ['MI', 'DI','FROB','PSICOV']
    sc = []
    targets=[]
    sc.append(scores)
    targets.append(target)
    sc.append(scores_di)
    targets.append(target_di)
    sc.append(scores_frob)
    targets.append(target_frob)
    sc.append(scores_psicov)
    targets.append(target_psicov)
    
    #plot.roc_curve_(targets, sc, labels, 'MSA Natural ROC curves' ,colors, '',natural_result_path + 'rocs.png')
    
    msa.summary(natural_msa, natural_result_path + "conservation.csv", 'PF00085_THIO_ECOLI_reference')
    
def compare_evolution_with_natural(execution_folder, protein, pdb_name, chain, top_result=0.80):
    optimization_folder = execution_folder + "optimization/"
    optimization_file_path = optimization_folder + "optimization.csv"
    clustered_sequences_path = optimization_folder + "clustered_sequences/"
    curated_sequences_path = optimization_folder + "curated_sequences_path/"
    mi_data_path = curated_sequences_path + "mi/"
    di_data_path = curated_sequences_path + "di/"
    frob_data_path = curated_sequences_path + "frob/"
    psicov_data_path = curated_sequences_path + "psicov/"
    conservation_path = curated_sequences_path + "conservation/"
    curated_sequences_result = curated_sequences_path + "resuts.csv"
    if not os.path.exists(curated_sequences_path):
        os.makedirs(curated_sequences_path) 
    if not os.path.exists(conservation_path):
        os.makedirs(conservation_path)  
    if not os.path.exists(mi_data_path):
        os.makedirs(mi_data_path)  
    if not os.path.exists(di_data_path):
        os.makedirs(di_data_path)
    if not os.path.exists(frob_data_path):
        os.makedirs(frob_data_path)
    if not os.path.exists(psicov_data_path):
        os.makedirs(psicov_data_path)   
    
    contact_map_sync = optimization_folder + "contact_map_sync.dat" 
    
    #natural information
    #natural_information(execution_folder)
    #sincronized hardcode for 2trx thio ecoli.
    util.sincronize_natural_evol_msas(clustered_sequences_path, curated_sequences_path, pattern, 2, -3)
     
    optimization_df = pandas.read_csv(optimization_file_path, header=0, index_col=0)
    
    columns = ["protein", "pdb", "chain", "beta", "nsus", "run", "auc_mi", "auc_01_mi", "auc_di", "auc_01_di", "auc_frob", "auc_01_frob","auc_psicov","auc_01_psicov","execution_time", "start_residue", "end_residue","message"]
    sincronized_top_results = pandas.DataFrame(columns=columns)
    index = 1
    for index_opti, row_optimization in optimization_df.iterrows():
        beta = str(row_optimization['beta'])
        nsus = str(row_optimization['nsus'])
        runs = str(int(row_optimization['run']))
        sufix = "sequences-beta" + beta + "-nsus" + nsus + "-runs" + runs
        curated_seq = curated_sequences_path + sufix + ".fasta.cluster"
        if((row_optimization['auc_mi'] >= top_result) | (row_optimization['auc_di'] >= top_result) | (row_optimization['auc_frob'] >= top_result) | (row_optimization['auc_psicov'] >= top_result)):
            try: 
                sincronized_top_results.at[index,'protein']=row_optimization['protein']
                sincronized_top_results.at[index,'pdb']=row_optimization['pdb']
                sincronized_top_results.at[index,'chain']=row_optimization['chain']
                sincronized_top_results.at[index,'beta']=row_optimization['beta']
                sincronized_top_results.at[index,'nsus']=row_optimization['nsus']
                sincronized_top_results.at[index,'run']=row_optimization['run']
                
                mi_path = mi_data_path + "mi_" + sufix + ".csv"
                dca_di = di_data_path + "di_" + sufix + ".csv" 
                dca_Frob = frob_data_path + "frob_" + sufix + ".csv"
                psicov_path = psicov_data_path + "psicov_" + sufix + ".csv"
    
                msa.buslje09(curated_seq, mi_path)
                msa.gaussDcaFrobenius(curated_seq, dca_Frob)
                msa.gaussDcaDirectInformation(curated_seq, dca_di)
                msa.create_msa_without_id(curated_seq, curated_seq+"_noid.aln")
                msa.psicov(curated_seq+"_noid.aln",psicov_path)
                
                target, scores = msa_analysis.getTargetScores(mi_path, contact_map_sync, constants.neighbour_window)
                auc, auc01 = util.getAUC(target, scores)
                target_di, scores_di = msa_analysis.getTargetScores(dca_di, contact_map_sync, constants.neighbour_window, split_char=' ')
                auc_di, auc01_di = util.getAUC(target_di, scores_di)
                target_frob, scores_frob = msa_analysis.getTargetScores(dca_Frob, contact_map_sync, constants.neighbour_window, split_char=' ')
                auc_frob, auc01_frob = util.getAUC(target_frob, scores_frob)
                target_psicov, scores_psicov = msa_analysis.getTargetScores(psicov_path, contact_map_sync, constants.neighbour_window, split_char=' ',score_position=4)
                auc_psicov, auc01_psicov = util.getAUC(target_psicov, scores_psicov)
                
                sincronized_top_results.at[index,'auc_mi']=auc
                sincronized_top_results.at[index,'auc_01_mi']=auc01
                sincronized_top_results.at[index,'auc_di']=auc_di
                sincronized_top_results.at[index,'auc_01_di']=auc01_di
                sincronized_top_results.at[index,'auc_frob']=auc_frob
                sincronized_top_results.at[index,'auc_01_frob']=auc01_frob
                sincronized_top_results.at[index,'auc_psicov']=auc_psicov
                sincronized_top_results.at[index,'auc_01_psicov']=auc01_psicov
                
                #msa.msa_information(msa_file, msa_conservation_path,msa_name)
                #msa_analysis.run_analisys_singular(optimization_df, index, natural_mi_result_path, mi_data_file, contact_map_sync, mi_data_path, constants.neighbour_window)
                sincronized_top_results.at[index,'message']="success"
            except Exception as inst:
                print inst
                x = inst.args
                print x
                sincronized_top_results.at[index,'message']="error"
            index = index + 1    
            sincronized_top_results.to_csv(curated_sequences_result)    
#puede ser que vaya en analisis
def analyse_optimus_msa(execution_folder, pdb_name, natural_result_path):
    optimization_folder = execution_folder + "optimization/"
    curated_sequences_path = optimization_folder + "curated_sequences/"
    mi_data_path = curated_sequences_path + "mi/"
    di_data_path = curated_sequences_path + "di/"
    frob_data_path = curated_sequences_path + "frob/"
    psicov_data_path = curated_sequences_path + "psicov/"
    
    
    curated_sequences_result = curated_sequences_path + "resuts.csv"
    optimization_result = optimization_folder + "optimization.csv"
    
    curated_sequences_path_best_results = curated_sequences_path + "best_coevolution_methods/"
    contact_map_sync = optimization_folder + "contact_map_sync.dat" 
    best_coev_values_path= curated_sequences_path + "best_coev_values.csv"
    best_01_coev_values_path= curated_sequences_path + "best_01_coev_values.csv"
    
    conservation_path = curated_sequences_path_best_results + "conservation/"
    
    
    if not os.path.exists(curated_sequences_path_best_results):
        os.makedirs(curated_sequences_path_best_results)   
    
    if not os.path.exists(conservation_path):
        os.makedirs(conservation_path)   
    
    
    columns = ["method", "auc_01", "beta", "nsus", "runs"]
    best_01_coev_values = pandas.DataFrame(columns=columns)
    best_coev_values = pandas.DataFrame(columns=columns)
    
    curated_result_df = pandas.read_csv(optimization_result, header=0, index_col=0)
    curated_result_df=curated_result_df.sort_values(["auc_mi"],ascending=False)
    top_auc_mi=curated_result_df.iloc[0]
    beta = str(top_auc_mi['beta'])
    nsus = str(top_auc_mi['nsus'])
    runs = str(int(top_auc_mi['run']))
    sufix = "sequences-beta" + beta + "-nsus" + nsus + "-runs" + runs
    mi_path = mi_data_path + "mi_" + sufix + ".csv"
    
    best_coev_values.at[1,'method']='mi'
    best_coev_values.at[1,'beta']=beta 
    best_coev_values.at[1,'nsus']=nsus 
    best_coev_values.at[1,'runs']=runs
    
    
    
    
    curated_result_df = pandas.read_csv(optimization_result, header=0, index_col=0)
    curated_result_df=curated_result_df.sort_values(["auc_01_mi"],ascending=False)
    top_auc_mi=curated_result_df.iloc[0]
    beta = str(top_auc_mi['beta'])
    nsus = str(top_auc_mi['nsus'])
    runs = str(int(top_auc_mi['run']))
    sufix = "sequences-beta" + beta + "-nsus" + nsus + "-runs" + runs
    mi_01_path = mi_data_path + "mi_" + sufix + ".csv"
    
    #conservation
    mi_msa_opt= curated_sequences_path + sufix + ".fasta.cluster"
    mi_msa_conservation= conservation_path + sufix + ".csv"
    msa.summary(mi_msa_opt, mi_msa_conservation, 'Conservation ' + sufix + ".fasta.cluster")
    
    msas_entropy=[]
    df_c = pandas.read_csv(mi_msa_conservation, usecols=['Entropy'])
    l=df_c['Entropy'].tolist()
    l.insert(0, None)
    msas_entropy.append([l, 'Evol', 'red'])
    df_c = pandas.read_csv(natural_result_path + "conservation.csv", usecols=['Entropy'])
    l=df_c['Entropy'].tolist()
    l.insert(0, None)
    msas_entropy.append([l, 'Natural', 'blue'])
    plot.conservation_comparation(msas_entropy, conservation_path + sufix + ".png", 'Conservation Natural MSA and Evol MSA')
    
    msa.seq_to_logo(mi_msa_opt, 'MI Best MSA. Beta ' + beta + ' Nsus ' + nsus + ' Runs '+runs)
    
    best_01_coev_values.at[1,'method']='mi'
    best_01_coev_values.at[1,'beta']=beta 
    best_01_coev_values.at[1,'nsus']=nsus 
    best_01_coev_values.at[1,'runs']=runs 
    
    curated_result_df = pandas.read_csv(optimization_result, header=0, index_col=0)
    curated_result_df=curated_result_df.sort_values(["auc_di"],ascending=False)
    top_auc_mi=curated_result_df.iloc[0]
    beta = str(top_auc_mi['beta'])
    nsus = str(top_auc_mi['nsus'])
    runs = str(int(top_auc_mi['run']))
    sufix = "sequences-beta" + beta + "-nsus" + nsus + "-runs" + runs
    di_path = di_data_path + "di_" + sufix + ".csv"
    
    best_coev_values.at[2,'method']='di'
    best_coev_values.at[2,'beta']=beta 
    best_coev_values.at[2,'nsus']=nsus 
    best_coev_values.at[2,'runs']=runs 
    
    curated_result_df = pandas.read_csv(optimization_result, header=0, index_col=0)
    curated_result_df=curated_result_df.sort_values(["auc_01_di"],ascending=False)
    top_auc_mi=curated_result_df.iloc[0]
    beta = str(top_auc_mi['beta'])
    nsus = str(top_auc_mi['nsus'])
    runs = str(int(top_auc_mi['run']))
    sufix = "sequences-beta" + beta + "-nsus" + nsus + "-runs" + runs
    di_01_path = di_data_path + "di_" + sufix + ".csv"
    
    best_01_coev_values.at[2,'method']='di'
    best_01_coev_values.at[2,'beta']=beta 
    best_01_coev_values.at[2,'nsus']=nsus 
    best_01_coev_values.at[2,'runs']=runs 
    
    #conservation
    di_msa_opt= curated_sequences_path + sufix + ".fasta.cluster"
    di_msa_conservation= conservation_path + sufix + ".csv"
    msa.summary(di_msa_opt, di_msa_conservation, 'Conservation ' + sufix + ".fasta.cluster")
    
    msas_entropy=[]
    df_c = pandas.read_csv(di_msa_conservation, usecols=['Entropy'])
    l=df_c['Entropy'].tolist()
    l.insert(0, None)
    msas_entropy.append([l, 'Evol', 'red'])
    df_c = pandas.read_csv(natural_result_path + "conservation.csv", usecols=['Entropy'])
    l=df_c['Entropy'].tolist()
    l.insert(0, None)
    msas_entropy.append([l, 'Natural', 'blue'])
    plot.conservation_comparation(msas_entropy, conservation_path + sufix + ".png", 'Conservation Natural MSA and Evol MSA')
    msa.seq_to_logo(di_msa_opt, 'DI Best MSA. Beta ' + beta + ' Nsus ' + nsus + ' Runs '+runs)
    curated_result_df = pandas.read_csv(optimization_result, header=0, index_col=0)
    curated_result_df=curated_result_df.sort_values(["auc_frob"],ascending=False)
    top_auc_mi=curated_result_df.iloc[0]
    beta = str(top_auc_mi['beta'])
    nsus = str(top_auc_mi['nsus'])
    runs = str(int(top_auc_mi['run']))
    sufix = "sequences-beta" + beta + "-nsus" + nsus + "-runs" + runs
    frob_path = frob_data_path + "frob_" + sufix + ".csv"
    
    best_coev_values.at[3,'method']='frob'
    best_coev_values.at[3,'beta']=beta 
    best_coev_values.at[3,'nsus']=nsus 
    best_coev_values.at[3,'runs']=runs 
    
    curated_result_df = pandas.read_csv(optimization_result, header=0, index_col=0)
    curated_result_df=curated_result_df.sort_values(["auc_01_frob"],ascending=False)
    top_auc_mi=curated_result_df.iloc[0]
    beta = str(top_auc_mi['beta'])
    nsus = str(top_auc_mi['nsus'])
    runs = str(int(top_auc_mi['run']))
    sufix = "sequences-beta" + beta + "-nsus" + nsus + "-runs" + runs
    frob_01_path = frob_data_path + "frob_" + sufix + ".csv"
    
    best_01_coev_values.at[3,'method']='frob'
    best_01_coev_values.at[3,'beta']=beta 
    best_01_coev_values.at[3,'nsus']=nsus 
    best_01_coev_values.at[3,'runs']=runs 
    
    #conservation
    frob_msa_opt= curated_sequences_path + sufix + ".fasta.cluster"
    frob_msa_conservation= conservation_path + sufix + ".csv"
    msa.summary(frob_msa_opt, frob_msa_conservation, 'Conservation ' + sufix + ".fasta.cluster")
    
    msas_entropy=[]
    df_c = pandas.read_csv(frob_msa_conservation, usecols=['Entropy'])
    l=df_c['Entropy'].tolist()
    l.insert(0, None)
    msas_entropy.append([l, 'Evol', 'red'])
    df_c = pandas.read_csv(natural_result_path + "conservation.csv", usecols=['Entropy'])
    l=df_c['Entropy'].tolist()
    l.insert(0, None)
    msas_entropy.append([l, 'Natural', 'blue'])
    plot.conservation_comparation(msas_entropy, conservation_path + sufix + ".png", 'Conservation Natural MSA and Evol MSA')
    msa.seq_to_logo(frob_msa_opt, 'FROB Best MSA. Beta ' + beta + ' Nsus ' + nsus + ' Runs '+runs)
    
    curated_result_df = pandas.read_csv(optimization_result, header=0, index_col=0)
    curated_result_df=curated_result_df.sort_values(["auc_psicov"],ascending=False)
    top_auc_mi=curated_result_df.iloc[0]
    beta = str(top_auc_mi['beta'])
    nsus = str(top_auc_mi['nsus'])
    runs = str(int(top_auc_mi['run']))
    sufix = "sequences-beta" + beta + "-nsus" + nsus + "-runs" + runs
    psicov_path = psicov_data_path + "psicov_" + sufix + ".csv"
    
    best_coev_values.at[4,'method']='psicov'
    best_coev_values.at[4,'beta']=beta 
    best_coev_values.at[4,'nsus']=nsus 
    best_coev_values.at[4,'runs']=runs 
    
    curated_result_df = pandas.read_csv(optimization_result, header=0, index_col=0)
    curated_result_df=curated_result_df.sort_values(["auc_01_psicov"],ascending=False)
    top_auc_mi=curated_result_df.iloc[0]
    beta = str(top_auc_mi['beta'])
    nsus = str(top_auc_mi['nsus'])
    runs = str(int(top_auc_mi['run']))
    sufix = "sequences-beta" + beta + "-nsus" + nsus + "-runs" + runs
    psicov_01_path = psicov_data_path + "psicov_" + sufix + ".csv"
    
    best_01_coev_values.at[4,'method']='psicov'
    best_01_coev_values.at[4,'beta']=beta 
    best_01_coev_values.at[4,'nsus']=nsus 
    best_01_coev_values.at[4,'runs']=runs 
    
    #conservation
    psicov_msa_opt= curated_sequences_path + sufix + ".fasta.cluster"
    psicov_msa_conservation= conservation_path + sufix + ".csv"
    msa.summary(psicov_msa_opt, psicov_msa_conservation, 'Conservation ' + sufix + ".fasta.cluster")
    
    msas_entropy=[]
    df_c = pandas.read_csv(psicov_msa_conservation, usecols=['Entropy'])
    l=df_c['Entropy'].tolist()
    l.insert(0, None)
    msas_entropy.append([l, 'Evol', 'red'])
    df_c = pandas.read_csv(natural_result_path + "conservation.csv", usecols=['Entropy'])
    l=df_c['Entropy'].tolist()
    l.insert(0, None)
    msas_entropy.append([l, 'Natural', 'blue'])
    plot.conservation_comparation(msas_entropy, conservation_path + sufix + ".png", 'Conservation Natural MSA and Evol MSA')
    
    msa.seq_to_logo(psicov_msa_opt, 'PSICOV Best MSA. Beta ' + beta + ' Nsus ' + nsus + ' Runs '+runs)
    
    #PLOT TODO JUNTO
    
    msas_entropy=[]
    df_c = pandas.read_csv(mi_msa_conservation, usecols=['Entropy'])
    l=df_c['Entropy'].tolist()
    l.insert(0, None)
    msas_entropy.append([l, 'Evol Best MI and FROB (Beta=5 Nsus=20)', 'green'])
    
    df_c = pandas.read_csv(di_msa_conservation, usecols=['Entropy'])
    l=df_c['Entropy'].tolist()
    l.insert(0, None)
    msas_entropy.append([l, 'Evol Best DI (Beta=0.5 Nsus=5)', 'yellow'])
    
    df_c = pandas.read_csv(psicov_msa_conservation, usecols=['Entropy'])
    l=df_c['Entropy'].tolist()
    l.insert(0, None)
    msas_entropy.append([l, 'Evol Best PSICOV (Beta=7 Nsus=20)', 'red'])
    
    
    df_c = pandas.read_csv(natural_result_path + "conservation.csv", usecols=['Entropy'])
    l=df_c['Entropy'].tolist()
    l.insert(0, None)
    msas_entropy.append([l, 'Natural', 'blue'])
    plot.conservation_comparation(msas_entropy, conservation_path +  "conservation.png", 'Conservation Natural MSA and Evol MSA')
    
    
    
    target, scores = msa_analysis.getTargetScores(mi_path, contact_map_sync, constants.neighbour_window)
    auc, auc01 = util.getAUC(target, scores)
    target_di, scores_di = msa_analysis.getTargetScores(di_path, contact_map_sync, constants.neighbour_window, split_char=' ')
    auc_di, auc01_di = util.getAUC(target_di, scores_di)
    target_frob, scores_frob = msa_analysis.getTargetScores(frob_path, contact_map_sync, constants.neighbour_window, split_char=' ')
    auc_frob, auc01_frob = util.getAUC(target_frob, scores_frob)
    target_psicov, scores_psicov = msa_analysis.getTargetScores(psicov_path, contact_map_sync, constants.neighbour_window, split_char=' ',score_position=4)
    auc_psicov, auc01_psicov = util.getAUC(target_psicov, scores_psicov)
    
    best_coev_values.at[1,'auc']=auc
    best_coev_values.at[2,'auc']=auc_di
    best_coev_values.at[3,'auc']=auc_frob
    best_coev_values.at[4,'auc']=auc_psicov
    
    best_coev_values.to_csv(best_coev_values_path) 
    
    #Plot optimos AUC
    colors = ['green','blue','yellow','red']
    labels = ['MI', 'DI','FROB','PSICOV']
    sc = []
    targets=[]
    sc.append(scores)
    targets.append(target)
    sc.append(scores_di)
    targets.append(target_di)
    sc.append(scores_frob)
    targets.append(target_frob)
    sc.append(scores_psicov)
    targets.append(target_psicov)
    
    plot.roc_curve_(targets, sc, labels, 'Best AUC ROC curves for coevolution methods' ,colors, '',curated_sequences_path + 'best_rocs_auc.png')
    
    
    target, scores = msa_analysis.getTargetScores(mi_01_path, contact_map_sync, constants.neighbour_window)
    auc, auc01 = util.getAUC(target, scores)
    target_di, scores_di = msa_analysis.getTargetScores(di_01_path, contact_map_sync, constants.neighbour_window, split_char=' ')
    auc_di, auc01_di = util.getAUC(target_di, scores_di)
    target_frob, scores_frob = msa_analysis.getTargetScores(frob_01_path, contact_map_sync, constants.neighbour_window, split_char=' ')
    auc_frob, auc01_frob = util.getAUC(target_frob, scores_frob)
    target_psicov, scores_psicov = msa_analysis.getTargetScores(psicov_01_path, contact_map_sync, constants.neighbour_window, split_char=' ',score_position=4)
    auc_psicov, auc01_psicov = util.getAUC(target_psicov, scores_psicov)
    
    
    best_01_coev_values.at[1,'auc_01']=auc01
    best_01_coev_values.at[2,'auc_01']=auc01_di
    best_01_coev_values.at[3,'auc_01']=auc01_frob
    best_01_coev_values.at[4,'auc_01']=auc01_psicov
     
    best_01_coev_values.to_csv(best_01_coev_values_path)    
    
    #Plot optimos AUC_01
    colors = ['green','blue','yellow','red']
    labels = ['MI', 'DI','FROB','PSICOV']
    sc = []
    targets=[]
    sc.append(scores)
    targets.append(target)
    sc.append(scores_di)
    targets.append(target_di)
    sc.append(scores_frob)
    targets.append(target_frob)
    sc.append(scores_psicov)
    targets.append(target_psicov)
    
    plot.roc_curve_(targets, sc, labels, 'Best AUC_01 ROC curves for coevolution methods' ,colors, '',curated_sequences_path + 'best_rocs_auc_01.png')
    
    #top_rank_result = top_rank(zmip_natural,m2,1,contact_map,mi_result_file_path+'top_1percent_withcon'+contact_threashold_str+'.png',mi_result_file_path,result_file,top_df,2, pdb_name,contact_threashold)
    natural_mi_result_path = natural_result_path + "mi.csv"
    natural_frob_result_path = natural_result_path + "frob.csv"
    natural_di_result_path= natural_result_path + "di.csv"
    natural_psicov_result_path = natural_result_path + "psicov.csv"
    
    
    target_nat, scores_nat = msa_analysis.getTargetScores(natural_mi_result_path, contact_map_sync, constants.neighbour_window)
    #auc, auc01 = util.getAUC(target, scores)
    sc = []
    targets=[]
    sc.append(scores_nat)
    targets.append(target_nat)
    sc.append(scores)
    targets.append(target)
    colors = ['blue','red']
    labels = ['Natural', 'Evol']
    plot.roc_curve_(targets,sc,labels,'Best Roc Curve Evol vs Natural for MI',colors, '' ,curated_sequences_path_best_results+'best_mi_roc.png')
    
    target_di_nat, scores_di_nat = msa_analysis.getTargetScores(natural_di_result_path, contact_map_sync, constants.neighbour_window, split_char=' ')
    #auc_di, auc01_di = util.getAUC(target_di, scores_di)
    sc = []
    targets=[]
    sc.append(scores_di_nat)
    targets.append(target_di_nat)
    sc.append(scores_di)
    targets.append(target_di)
    
    plot.roc_curve_(targets,sc,labels,'Best Roc Curve Evol vs Natural for DI',colors, '' ,curated_sequences_path_best_results+'best_di_roc.png')
    
    
    target_frob_nat, scores_frob_nat = msa_analysis.getTargetScores(natural_frob_result_path, contact_map_sync, constants.neighbour_window, split_char=' ')
    #auc_frob, auc01_frob = util.getAUC(target_frob, scores_frob)
    sc = []
    targets=[]
    sc.append(scores_frob_nat)
    targets.append(target_frob_nat)
    sc.append(scores_frob)
    targets.append(target_frob)
    plot.roc_curve_(targets,sc,labels,'Best Roc Curve Evol vs Natural for FROB',colors, '' ,curated_sequences_path_best_results+'best_frob_roc.png')
    
    
    target_psicov_nat, scores_psicov_nat = msa_analysis.getTargetScores(natural_psicov_result_path, contact_map_sync, constants.neighbour_window, split_char=' ',score_position=4)
    #auc_psicov, auc01_psicov = util.getAUC(target_psicov, scores_psicov)
    sc = []
    targets=[]
    sc.append(scores_psicov_nat)
    targets.append(target_psicov_nat)
    sc.append(scores_psicov)
    targets.append(target_psicov)
    plot.roc_curve_(targets,sc,labels,'Best Roc Curve Evol vs Natural for PSICOV',colors, '' ,curated_sequences_path_best_results+'best_psicov_roc.png')
    
    
    
    #MI 
    mi_natural = util.load_zmip(natural_mi_result_path, constants.neighbour_window)
    util.order(mi_natural)
    mi_evol = util.load_zmip(mi_01_path , constants.neighbour_window)
    util.order(mi_evol)
    
    di_natural = util.load_zmip(natural_di_result_path, constants.neighbour_window,split_char=' ')
    di_evol = util.load_zmip(di_01_path, constants.neighbour_window,split_char=' ')
    
    frob_natural = util.load_zmip(natural_frob_result_path, constants.neighbour_window,split_char=' ')
    frob_evol = util.load_zmip(frob_01_path, constants.neighbour_window,split_char=' ')
    
    psicov_natural = util.load_zmip(natural_psicov_result_path, constants.neighbour_window,split_char=' ')
    psicov_evol = util.load_zmip(psicov_01_path, constants.neighbour_window,split_char=' ')
    
    [r.pop(2) for r in psicov_evol]
    [r.pop(2) for r in psicov_evol]
    [r.pop(2) for r in psicov_natural]
    [r.pop(2) for r in psicov_natural]
    
    coevolution_analisys_df = pandas.DataFrame()
    coevolution_analisys_df.set_value(1,'top',0.5)
    coevolution_analisys_df.set_value(2,'top',1)
    coevolution_analisys_df.set_value(3,'top',2)
    coevolution_analisys_df.set_value(4,'top',3)
    coevolution_analisys_df.set_value(5,'top',4)
    coevolution_analisys_df.set_value(6,'top',5)
    msa_analysis.coevolution_analisys('MI',coevolution_analisys_df, 1, mi_natural, mi_evol, curated_sequences_path_best_results, contact_map_sync,  constants.neighbour_window, pdb_name)
    msa_analysis.coevolution_analisys('DI',coevolution_analisys_df, 1, di_natural, di_evol, curated_sequences_path_best_results, contact_map_sync,  constants.neighbour_window, pdb_name)
    msa_analysis.coevolution_analisys('FROB',coevolution_analisys_df, 1, frob_natural, frob_evol, curated_sequences_path_best_results, contact_map_sync,  constants.neighbour_window, pdb_name)
    msa_analysis.coevolution_analisys('PSICOV',coevolution_analisys_df, 1, psicov_natural, psicov_evol, curated_sequences_path_best_results, contact_map_sync,  constants.neighbour_window, pdb_name)
    
    coevolution_analisys_df.to_csv(curated_sequences_path_best_results +'tops_contact_threashold_1.csv')
    
    #msa_analysis.top_coevolution(natural_coevolution,evolutionated_coevolution,top,contact_map,contact_map_top_coev_path,filename,result_file, top_df,index,pdb_name,contact_threashold=1)
    
    
    
    
def optimize_evolution_conformers(execution_folder,structures,protein, chain):
    for pdb_name in structures:
        execution_folder_conf = execution_folder + pdb_name + "/"
        optimize_evolution(execution_folder_conf, protein, pdb_name, chain, synchronized=True)

def analyse_optimus_msa_conformers(execution_folder,structures, natural_result_path):
    for pdb_name in structures:
        execution_folder_conf = execution_folder + pdb_name + "/"
        analyse_optimus_msa(execution_folder_conf, pdb_name, natural_result_path)


def analyse_top_coevolution_conformers(execution_folder,structures,sum_contact_map_path,natural_result_path):
    mi_coev = []
    di_coev = []
    frob_coev = []
    psicov_coev = []
    coevolution_results = execution_folder + "coevolution_results/"
    if not os.path.exists(coevolution_results):
        os.makedirs(coevolution_results)
    for pdb_name in structures:
        execution_folder_conf = execution_folder + pdb_name + "/"
        optimization_folder = execution_folder_conf + "optimization/"
        curated_sequences_path = optimization_folder + "curated_sequences/"
        mi_data_path = curated_sequences_path + "mi/"
        di_data_path = curated_sequences_path + "di/"
        frob_data_path = curated_sequences_path + "frob/"
        psicov_data_path = curated_sequences_path + "psicov/"
        best_01_coev_values_path= curated_sequences_path + "best_01_coev_values.csv"
        best_01_coev_df = pandas.read_csv(best_01_coev_values_path, header=0, index_col=0)
        for index_opti, best_conf in best_01_coev_df.iterrows(): 
            beta = str(best_conf['beta'])
            nsus = str(best_conf['nsus'])
            runs = str(int(best_conf['runs']))
            sufix = "sequences-beta" + beta + "-nsus" + nsus + "-runs" + runs + ".csv"
            if(best_conf['method']=='mi'):
                mi_evol = util.load_zmip(mi_data_path + "mi_"+sufix , constants.neighbour_window)
                util.order(mi_evol)
                mi_coev.append(mi_evol)
            if(best_conf['method']=='di'):
                di_evol = util.load_zmip(di_data_path+ "di_"+sufix, constants.neighbour_window,split_char=' ') 
                di_coev.append(di_evol)
            if(best_conf['method']=='frob'):
                frob_evol = util.load_zmip(frob_data_path + "frob_"+ sufix, constants.neighbour_window,split_char=' ')
                frob_coev.append(frob_evol)
            if(best_conf['method']=='psicov'):
                psicov_evol = util.load_zmip(psicov_data_path + "psicov_"+ sufix, constants.neighbour_window,split_char=' ')
                [r.pop(2) for r in psicov_evol]
                [r.pop(2) for r in psicov_evol]
                psicov_coev.append(psicov_evol)
    
    
    
    
    
    natural_mi_result_path = natural_result_path + "mi.csv"
    natural_frob_result_path = natural_result_path + "frob.csv"
    natural_di_result_path= natural_result_path + "di.csv"
    natural_psicov_result_path = natural_result_path + "psicov.csv"
    
    #MI 
    mi_natural = util.load_zmip(natural_mi_result_path, constants.neighbour_window)
    util.order(mi_natural)
    di_natural = util.load_zmip(natural_di_result_path, constants.neighbour_window,split_char=' ')
    frob_natural = util.load_zmip(natural_frob_result_path, constants.neighbour_window,split_char=' ')
    psicov_natural = util.load_zmip(natural_psicov_result_path, constants.neighbour_window,split_char=' ')
    [r.pop(2) for r in psicov_natural]
    [r.pop(2) for r in psicov_natural]
    
    
    coevolution_analisys_df = pandas.DataFrame()
    index=1
    
    corr_df = pandas.DataFrame()
    
    msa_analysis.top_coevolution_analysis('mi', mi_coev,0.5, sum_contact_map_path, structures,coevolution_results,mi_natural,coevolution_analisys_df,index,corr_df)
    coevolution_analisys_df.to_csv(coevolution_results +'top_thio_ecoli_conformers.csv')
    msa_analysis.top_coevolution_analysis('di', di_coev,0.5, sum_contact_map_path, structures,coevolution_results ,di_natural,coevolution_analisys_df,index,corr_df)     
    coevolution_analisys_df.to_csv(coevolution_results +'top_thio_ecoli_conformers.csv')
    msa_analysis.top_coevolution_analysis('frob', frob_coev,0.5, sum_contact_map_path, structures,coevolution_results,frob_natural,coevolution_analisys_df,index,corr_df)            
    coevolution_analisys_df.to_csv(coevolution_results +'top_thio_ecoli_conformers.csv')
    msa_analysis.top_coevolution_analysis('psicov', psicov_coev,0.5, sum_contact_map_path, structures,coevolution_results,psicov_natural,coevolution_analisys_df,index,corr_df)            
    coevolution_analisys_df.to_csv(coevolution_results +'top_thio_ecoli_conformers.csv')
    index=index*3 +1
    msa_analysis.top_coevolution_analysis('mi', mi_coev,1, sum_contact_map_path, structures,coevolution_results,mi_natural,coevolution_analisys_df,index,corr_df)  
    msa_analysis.top_coevolution_analysis('di', di_coev,1, sum_contact_map_path, structures,coevolution_results ,di_natural,coevolution_analisys_df,index,corr_df)       
    msa_analysis.top_coevolution_analysis('frob', frob_coev,1, sum_contact_map_path, structures,coevolution_results,frob_natural,coevolution_analisys_df,index,corr_df)             
    msa_analysis.top_coevolution_analysis('psicov', psicov_coev,1, sum_contact_map_path, structures,coevolution_results,psicov_natural,coevolution_analisys_df,index,corr_df)  
    index=index*3 + 1
    msa_analysis.top_coevolution_analysis('mi', mi_coev,2, sum_contact_map_path, structures,coevolution_results,mi_natural,coevolution_analisys_df,index,corr_df)  
    msa_analysis.top_coevolution_analysis('di', di_coev,2, sum_contact_map_path, structures,coevolution_results ,di_natural,coevolution_analisys_df,index,corr_df)   
    msa_analysis.top_coevolution_analysis('frob', frob_coev,2, sum_contact_map_path, structures,coevolution_results,frob_natural,coevolution_analisys_df,index,corr_df)              
    msa_analysis.top_coevolution_analysis('psicov', psicov_coev,2, sum_contact_map_path, structures,coevolution_results,psicov_natural,coevolution_analisys_df,index,corr_df)  
    index=index*3 +1
    msa_analysis.top_coevolution_analysis('mi', mi_coev,3, sum_contact_map_path, structures,coevolution_results,mi_natural,coevolution_analisys_df,index,corr_df)  
    msa_analysis.top_coevolution_analysis('di', di_coev,3, sum_contact_map_path, structures,coevolution_results ,di_natural,coevolution_analisys_df,index,corr_df)       
    msa_analysis.top_coevolution_analysis('frob', frob_coev,3, sum_contact_map_path, structures,coevolution_results,frob_natural,coevolution_analisys_df,index,corr_df)             
    msa_analysis.top_coevolution_analysis('psicov', psicov_coev,3, sum_contact_map_path, structures,coevolution_results,psicov_natural,coevolution_analisys_df,index,corr_df)              
    index=index*3 +1
    msa_analysis.top_coevolution_analysis('mi', mi_coev,4, sum_contact_map_path, structures,coevolution_results,mi_natural,coevolution_analisys_df,index,corr_df)  
    msa_analysis.top_coevolution_analysis('di', di_coev,4, sum_contact_map_path, structures,coevolution_results ,di_natural,coevolution_analisys_df,index,corr_df)       
    msa_analysis.top_coevolution_analysis('frob', frob_coev,4, sum_contact_map_path, structures,coevolution_results,frob_natural,coevolution_analisys_df,index,corr_df)             
    msa_analysis.top_coevolution_analysis('psicov', psicov_coev,4, sum_contact_map_path, structures,coevolution_results,psicov_natural,coevolution_analisys_df,index,corr_df)              
    index=index*3 + 1
    msa_analysis.top_coevolution_analysis('mi', mi_coev,5, sum_contact_map_path, structures,coevolution_results,mi_natural,coevolution_analisys_df,index,corr_df)  
    msa_analysis.top_coevolution_analysis('di', di_coev,5, sum_contact_map_path, structures,coevolution_results ,di_natural,coevolution_analisys_df,index,corr_df)       
    msa_analysis.top_coevolution_analysis('frob', frob_coev,5, sum_contact_map_path, structures,coevolution_results,frob_natural,coevolution_analisys_df,index,corr_df)             
    msa_analysis.top_coevolution_analysis('psicov', psicov_coev,5, sum_contact_map_path, structures,coevolution_results,psicov_natural,coevolution_analisys_df,index,corr_df)              

    coevolution_analisys_df.to_csv(coevolution_results +'top_thio_ecoli_conformers.csv')
    corr_df.to_csv(coevolution_results +'top_thio_ecoli_conformers_corr.csv')
   
def generate_combined_msa(execution_folder,structures,num):
    mi_msa=[]
    di_msa=[]
    frob_msa=[]
    psicov_msa=[]
    mi_msa_noid=[]
    di_msa_noid=[]
    frob_msa_noid=[]
    psicov_msa_noid=[]
    for pdb_name in structures:
        execution_folder_conf = execution_folder + pdb_name + "/"
        optimization_folder = execution_folder_conf + "optimization/"
        curated_sequences_path = optimization_folder + "curated_sequences/"
        best_01_coev_values_path= curated_sequences_path + "best_01_coev_values.csv"
        best_01_coev_df = pandas.read_csv(best_01_coev_values_path, header=0, index_col=0)
        for index_opti, best_conf in best_01_coev_df.iterrows(): 
            beta = str(best_conf['beta'])
            nsus = str(best_conf['nsus'])
            runs = str(int(best_conf['runs']))
            sufix = "sequences-beta" + beta + "-nsus" + nsus + "-runs" + runs + ".fasta"
            msa_path=curated_sequences_path+sufix
            msa_path_no_id=curated_sequences_path+sufix+ '_no_seq_ids.fasta'
            if(best_conf['method']=='mi'):
                mi_msa.append(msa_path)
                mi_msa_noid.append(msa_path_no_id)
            if(best_conf['method']=='di'):
                di_msa.append(msa_path)
                di_msa_noid.append(msa_path_no_id)
            if(best_conf['method']=='frob'):
                frob_msa.append(msa_path)
                frob_msa_noid.append(msa_path_no_id)
            if(best_conf['method']=='psicov'):
                psicov_msa.append(msa_path)
                psicov_msa_noid.append(msa_path_no_id)
    for msa_ in mi_msa:
        msa.create_msa_without_id(msa_, str(msa_) + '_no_seq_ids.fasta')
    for msa_ in di_msa:
        msa.create_msa_without_id(msa_, str(msa_) + '_no_seq_ids.fasta')
    for msa_ in frob_msa:
        msa.create_msa_without_id(msa_, str(msa_) + '_no_seq_ids.fasta')
    for msa_ in psicov_msa:
        msa.create_msa_without_id(msa_, str(msa_) + '_no_seq_ids.fasta')
    
    bootstraping_folder = execution_folder + 'msa_bootstraping/'
    if not os.path.exists(bootstraping_folder):
        os.makedirs(bootstraping_folder)
    
    msa_conjuntion_bootstrap_path = bootstraping_folder + 'msa_conjuntion_' + str(num) + '/'
    if not os.path.exists(msa_conjuntion_bootstrap_path):
        os.makedirs(msa_conjuntion_bootstrap_path)
        
    msa.create_msa_bootstrap(mi_msa_noid, msa_conjuntion_bootstrap_path + 'mi_msa_bootstrap_' + str(num) + '.fasta', num)
    msa.create_msa_bootstrap(di_msa_noid, msa_conjuntion_bootstrap_path + 'di_msa_bootstrap_' + str(num) + '.fasta', num)
    msa.create_msa_bootstrap(frob_msa_noid, msa_conjuntion_bootstrap_path + 'frob_msa_bootstrap_' + str(num) + '.fasta', num)
    msa.create_msa_bootstrap(psicov_msa_noid, msa_conjuntion_bootstrap_path + 'psicov_msa_bootstrap_' + str(num) + '.fasta', num)
    

def conjunction_analysis(execution_folder, structures,contact_map_path,natural_result_path):
    num_=[20000]
    df = pandas.DataFrame()
    for num in num_:
        generate_combined_msa(execution_folder, structures, num)
        analisys_msa_conjuntion_thio_ecoli_conformeros(df,execution_folder, structures,num,contact_map_path,natural_result_path)
        
        
def analisys_msa_conjuntion_thio_ecoli_conformeros(df,execution_folder,structures,num, contact_map_path,natural_result_path):
    bootstraping_folder = execution_folder + 'msa_bootstraping/'
    msa_conjuntion_bootstrap_path = bootstraping_folder + 'msa_conjuntion_' + str(num) + '/'
    msa.buslje09(msa_conjuntion_bootstrap_path + 'mi_msa_bootstrap_' + str(num) + '.fasta', msa_conjuntion_bootstrap_path + 'mi_msa_bootstrap_' + str(num) + '.csv')
    msa.gaussDcaFrobenius(msa_conjuntion_bootstrap_path + 'frob_msa_bootstrap_' + str(num) + '.fasta', msa_conjuntion_bootstrap_path + 'frob_msa_bootstrap_' + str(num) + '.csv')
    msa.gaussDcaDirectInformation(msa_conjuntion_bootstrap_path + 'di_msa_bootstrap_' + str(num) + '.fasta', msa_conjuntion_bootstrap_path + 'di_msa_bootstrap_' + str(num) + '.csv')
    msa.create_msa_without_id(msa_conjuntion_bootstrap_path + 'psicov_msa_bootstrap_' + str(num) + '.fasta', msa_conjuntion_bootstrap_path + 'psicov_msa_bootstrap_' + str(num) + '.fasta'+"_noid.aln")
    msa.psicov(msa_conjuntion_bootstrap_path + 'psicov_msa_bootstrap_' + str(num) + '.fasta'+"_noid.aln",msa_conjuntion_bootstrap_path + 'psicov_msa_bootstrap_' + str(num) + '.csv')
    target, scores = msa_analysis.getTargetScores(msa_conjuntion_bootstrap_path + 'mi_msa_bootstrap_' + str(num) + '.csv', contact_map_path, constants.neighbour_window)
    auc, auc01 = util.getAUC(target, scores)
    target_di, scores_di = msa_analysis.getTargetScores(msa_conjuntion_bootstrap_path + 'di_msa_bootstrap_' + str(num) + '.csv', contact_map_path, constants.neighbour_window, split_char=' ')
    auc_di, auc01_di = util.getAUC(target_di, scores_di)
    target_frob, scores_frob = msa_analysis.getTargetScores(msa_conjuntion_bootstrap_path + 'frob_msa_bootstrap_' + str(num) + '.csv', contact_map_path, constants.neighbour_window, split_char=' ')
    auc_frob, auc01_frob = util.getAUC(target_frob, scores_frob)
    target_psicov, scores_psicov = msa_analysis.getTargetScores(msa_conjuntion_bootstrap_path + 'psicov_msa_bootstrap_' + str(num) + '.csv', contact_map_path, constants.neighbour_window, split_char=' ',score_position=4)
    auc_psicov, auc01_psicov = util.getAUC(target_psicov, scores_psicov)
    
    index = len(df.index) + 1
    df.at[index,'num']=num
    df.at[index,'auc_mi_01']=auc01
    df.at[index,'auc_mi']=auc
    df.at[index,'auc_di_01']=auc01_di
    df.at[index,'auc_di']=auc_di
    df.at[index,'auc_frob_01']=auc01_frob
    df.at[index,'auc_frob']=auc_frob
    df.at[index,'auc_psicov_01']=auc01_psicov
    df.at[index,'auc_psicov']=auc_psicov
    df.to_csv(msa_conjuntion_bootstrap_path + 'result_conjunction.csv')
    
    #Plot optimos AUC_01
    colors = ['green','blue','yellow','red']
    labels = ['MI', 'DI','FROB','PSICOV']
    sc = []
    targets=[]
    sc.append(scores)
    targets.append(target)
    sc.append(scores_di)
    targets.append(target_di)
    sc.append(scores_frob)
    targets.append(target_frob)
    sc.append(scores_psicov)
    targets.append(target_psicov)
    
    plot.roc_curve_(targets, sc, labels, 'Bootstrap ' + str(num) +  'MSA Conjunction ROCs ' ,colors, '',msa_conjuntion_bootstrap_path + 'rocs_auc.png')
    
    natural_mi_result_path = natural_result_path + "mi.csv"
    natural_frob_result_path = natural_result_path + "frob.csv"
    natural_di_result_path= natural_result_path + "di.csv"
    natural_psicov_result_path = natural_result_path + "psicov.csv"
    
    #MI 
    mi_natural = util.load_zmip(natural_mi_result_path, constants.neighbour_window)
    util.order(mi_natural)
    mi_evol = util.load_zmip(msa_conjuntion_bootstrap_path + 'mi_msa_bootstrap_' + str(num) + '.csv' , constants.neighbour_window)
    util.order(mi_evol)
    
    di_natural = util.load_zmip(natural_di_result_path, constants.neighbour_window,split_char=' ')
    di_evol = util.load_zmip(msa_conjuntion_bootstrap_path + 'di_msa_bootstrap_' + str(num) + '.csv', constants.neighbour_window,split_char=' ')
    
    frob_natural = util.load_zmip(natural_frob_result_path, constants.neighbour_window,split_char=' ')
    frob_evol = util.load_zmip(msa_conjuntion_bootstrap_path + 'frob_msa_bootstrap_' + str(num) + '.csv', constants.neighbour_window,split_char=' ')
    
    psicov_natural = util.load_zmip(natural_psicov_result_path, constants.neighbour_window,split_char=' ')
    psicov_evol = util.load_zmip(msa_conjuntion_bootstrap_path + 'psicov_msa_bootstrap_' + str(num) + '.csv', constants.neighbour_window,split_char=' ')
    
    [r.pop(2) for r in psicov_evol]
    [r.pop(2) for r in psicov_evol]
    [r.pop(2) for r in psicov_natural]
    [r.pop(2) for r in psicov_natural]
    
    coevolution_analisys_df = pandas.DataFrame()
    coevolution_analisys_df.set_value(1,'top',0.5)
    coevolution_analisys_df.set_value(2,'top',1)
    coevolution_analisys_df.set_value(3,'top',2)
    coevolution_analisys_df.set_value(4,'top',3)
    coevolution_analisys_df.set_value(5,'top',4)
    coevolution_analisys_df.set_value(6,'top',5)
    msa_analysis.coevolution_analisys('MI',coevolution_analisys_df, 1, mi_natural, mi_evol, msa_conjuntion_bootstrap_path, contact_map_path,  constants.neighbour_window, 'Conjuntion MSA')
    msa_analysis.coevolution_analisys('DI',coevolution_analisys_df, 1, di_natural, di_evol, msa_conjuntion_bootstrap_path, contact_map_path,  constants.neighbour_window, 'Conjuntion MSA')
    msa_analysis.coevolution_analisys('FROB',coevolution_analisys_df, 1, frob_natural, frob_evol, msa_conjuntion_bootstrap_path, contact_map_path,  constants.neighbour_window, 'Conjuntion MSA')
    msa_analysis.coevolution_analisys('PSICOV',coevolution_analisys_df, 1, psicov_natural, psicov_evol, msa_conjuntion_bootstrap_path, contact_map_path,  constants.neighbour_window, 'Conjuntion MSA')
    
    coevolution_analisys_df.to_csv(msa_conjuntion_bootstrap_path +'tops_contact_threashold_1.csv')
    
    
 
def analisys_singular_conjunction_thio_ecoli_conformeros(execution_folder,structures, num=1500):
    mi_data_output_path = execution_folder + 'mi_data_path/'
    mi_data_path_file = mi_data_output_path + "zmip_bootstrap.csv"
    zmip_natural_result_file = "../THIO_ECOLI_4_107_2TRX_A/natural/zmip_PF00085_THIO_ECOLI_reference.csv"
    contact_map_path = '../THIO_ECOLI_4_107/sum_contact_map.dat'
    # zmip_reference_result_path = '../THIO_ECOLI_4_107/2TRX/mi_data_path/zmip_sequences-beta5.0-nsus15.0-runs20000.csv'
    top_df = pandas.DataFrame()
    pdb_name = 'THIO_ECOLI CONJUNCTION CONFORMEROS'
    dataanalisys.run_analisys_singular(top_df, 1, zmip_natural_result_file, mi_data_path_file, contact_map_path, mi_data_output_path, window, pdb_name) 
    # dataanalisys.top_rank_intersection(execution_folder, contact_map_path,zmip_natural_result_path, mi_data_path_file, top_df, zmip_reference_result_path, index=1, window, contact_threshold=1, top_threshold=1, sinchronize_with_natural=True)
    
    top_df.to_csv(execution_folder + 'result_conjunction.csv')
    
    
   
    
   
   
   
          