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
from paramiko.util import constant_time_bytes_eq

import constants

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
def optimize_evolution(execution_folder, protein, pdb_name, chain):
    logging.info('Run Optimization For ' + pdb_name)
    start_time_total = time.time()
    columns = ["protein", "pdb", "chain", "beta", "nsus", "run", "auc_mi", "auc_01_mi", "auc_di", "auc_01_di", "auc_frob", "auc_01_frob","auc_psicov","auc_01_psicov","execution_time", "start_residue", "end_residue","message"]
    df = pandas.DataFrame(columns=columns)
    optimization_folder = execution_folder + "optimization/"
    optimization_file_path = optimization_folder + "optimization.csv"
    scpe_sequences = optimization_folder + "scpe_sequences/"
    clustered_sequences_path = optimization_folder + "clustered_sequences/"
    sincronized_evol_path = optimization_folder + "sincronized_sequences/"
    mi_data_path = clustered_sequences_path + "mi/"
    di_data_path = clustered_sequences_path + "di/"
    frob_data_path = clustered_sequences_path + "frob/"
    psicov_path = clustered_sequences_path + "psicov/"
    contact_map = optimization_folder + "contact_map.dat"
    contact_map_sync = optimization_folder + "contact_map_sync.dat" 
    pdb_path = execution_folder + pdb_name + ".pdb"
    pdb_clean_path = execution_folder +  pdb_name + "_chain_" + chain +".pdb"
    #create all the structure for the optimization
    if not os.path.exists(scpe_sequences):
        os.makedirs(scpe_sequences)
    if not os.path.exists(clustered_sequences_path):
        os.makedirs(clustered_sequences_path)
    if not os.path.exists(sincronized_evol_path):
        os.makedirs(sincronized_evol_path)
    if not os.path.exists(mi_data_path):
        os.makedirs(mi_data_path) 
    if not os.path.exists(di_data_path):
        os.makedirs(di_data_path)
    if not os.path.exists(frob_data_path):
        os.makedirs(frob_data_path)         
    if not os.path.exists(psicov_path):
        os.makedirs(psicov_path)      
    pdb.download(pdb_name, pdb_path)
    pdb.clean_pdb(pdb_path, pdb_clean_path, chain)
    
    beta = ["0.001","0.01","0.1","0.5", "1.0", "2.0", "3.0", "5.0", "7.0", "10.0", "15.0", "20.0"]
    runs = ["1000", "5000", "10000", "20000"]
    nsus = ["1.0", "2.0", "3.0", "5.0", "7.0", "10.0", "15.0","20.0"]
    auc_max = 0
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
                    auc, auc01, auc_di, auc01_di, auc_frob, auc01_frob, auc_psicov, auc01_psicov, count_seq_scpe, count_seq_cluster, seq_lenght = evol(pdb_clean_path, b, r, sus, chain, scpe_sequences, clustered_sequences_path, sincronized_evol_path, mi_data_path, di_data_path, frob_data_path, psicov_path,contact_map, contact_map_sync)
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
                    if(auc > auc_max):
                        parameters = (b, sus, r) 
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
    return parameters

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
        msa.buslje09(curated_sequences_path, mi_data_path)
        msa.gaussDcaFrobenius(curated_sequences_path, dca_Frob)
        msa.gaussDcaDirectInformation(curated_sequences_path, dca_di)
        msa.create_msa_without_id(curated_sequences_path, curated_sequences_path+"_noid.aln")
        msa.psicov(curated_sequences_path+"_noid.aln",psicov_path)
        
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
    
    plot.roc_curve_(targets, sc, labels, 'MSA Natural ROC curves' ,colors, '',natural_result_path + 'rocs.png')
    
    
    
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
def analyse_optimus_msa(execution_folder, pdb_name):
    optimization_folder = execution_folder + "optimization/"
    curated_sequences_path = optimization_folder + "curated_sequences_path/"
    mi_data_path = curated_sequences_path + "mi/"
    di_data_path = curated_sequences_path + "di/"
    frob_data_path = curated_sequences_path + "frob/"
    psicov_data_path = curated_sequences_path + "psicov/"
    conservation_path = curated_sequences_path + "conservation/"
    curated_sequences_result = curated_sequences_path + "resuts.csv"
    curated_sequences_path_best_results = curated_sequences_path + "best_coevolution_methods/"
    contact_map_sync = optimization_folder + "contact_map_sync.dat" 
    best_coev_values_path= curated_sequences_path + "best_coev_values.csv"
    best_01_coev_values_path= curated_sequences_path + "best_01_coev_values.csv"
    
    if not os.path.exists(curated_sequences_path_best_results):
        os.makedirs(curated_sequences_path_best_results)   
    
    columns = ["method", "auc_01", "beta", "nsus", "runs"]
    best_01_coev_values = pandas.DataFrame(columns=columns)
    best_coev_values = pandas.DataFrame(columns=columns)
    
    curated_result_df = pandas.read_csv(curated_sequences_result, header=0, index_col=0)
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
    
    
    curated_result_df = pandas.read_csv(curated_sequences_result, header=0, index_col=0)
    curated_result_df=curated_result_df.sort_values(["auc_01_mi"],ascending=False)
    top_auc_mi=curated_result_df.iloc[0]
    beta = str(top_auc_mi['beta'])
    nsus = str(top_auc_mi['nsus'])
    runs = str(int(top_auc_mi['run']))
    sufix = "sequences-beta" + beta + "-nsus" + nsus + "-runs" + runs
    mi_01_path = mi_data_path + "mi_" + sufix + ".csv"
    
    best_01_coev_values.at[1,'method']='mi'
    best_01_coev_values.at[1,'beta']=beta 
    best_01_coev_values.at[1,'nsus']=nsus 
    best_01_coev_values.at[1,'runs']=runs 
    
    curated_result_df = pandas.read_csv(curated_sequences_result, header=0, index_col=0)
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
    
    curated_result_df = pandas.read_csv(curated_sequences_result, header=0, index_col=0)
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
    
    
    curated_result_df = pandas.read_csv(curated_sequences_result, header=0, index_col=0)
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
    
    curated_result_df = pandas.read_csv(curated_sequences_result, header=0, index_col=0)
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
    
    curated_result_df = pandas.read_csv(curated_sequences_result, header=0, index_col=0)
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
    
    curated_result_df = pandas.read_csv(curated_sequences_result, header=0, index_col=0)
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
    
    #TODO Conservation MSA
    #conservation_msa
    
    #
    
    #top_rank_result = top_rank(zmip_natural,m2,1,contact_map,mi_result_file_path+'top_1percent_withcon'+contact_threashold_str+'.png',mi_result_file_path,result_file,top_df,2, pdb_name,contact_threashold)
    natural_result_path = execution_folder + "/natural/"
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
    
    coevolution_analisys_df.to_csv(curated_sequences_path_best_results +'tops_contact_threashold_1csv')
    
    #msa_analysis.top_coevolution(natural_coevolution,evolutionated_coevolution,top,contact_map,contact_map_top_coev_path,filename,result_file, top_df,index,pdb_name,contact_threashold=1)
    
def evol_conformers(execution_folder, best_coev_results_paths, structures):
    chain = 'A'
    
    protein = 'THIO_ECOLI_4_107'
    columns = ["protein", "pdb", "chain", "beta", "nsus", "run", "auc_mi", "auc_01_mi", "auc_di", "auc_01_di", "auc_frob", "auc_01_frob","auc_psicov","auc_01_psicov","execution_time", "start_residue", "end_residue","message"]
    
    for pdb_name in structures:
        index = 1
        df = pandas.DataFrame(columns=columns)
        conf_path = execution_folder + pdb_name + "/"
        scpe_sequences = conf_path + "scpe_sequences/"
        clustered_sequences_path = conf_path + "clustered_sequences/"
        curated_sequences_path = conf_path + "curated_sequences/"
        mi_data_path = curated_sequences_path + "mi/"
        di_data_path = curated_sequences_path + "di/"
        frob_data_path = curated_sequences_path + "frob/"
        psicov_data_path = curated_sequences_path + "psicov/"
        contact_map = conf_path + "contact_map.dat"
        contact_map_sync = conf_path + "contact_map_sync.dat" 
        pdb_clean_path = conf_path +  pdb_name + "_clean.pdb"
        
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
        if not os.path.exists(psicov_data_path):
            os.makedirs(psicov_data_path) 
            
        columns = ["protein", "pdb", "chain", "beta", "nsus", "run", "auc_mi", "auc_01_mi", "auc_di", "auc_01_di", "auc_frob", "auc_01_frob","auc_psicov","auc_01_psicov","execution_time", "start_residue", "end_residue","message"]    
        result_path =  curated_sequences_path + "results.csv"        
        best_coev_results = pandas.read_csv(best_coev_results_paths, header=0, index_col=0)
        for index_opti, method_coev in best_coev_results.iterrows():
            method=method_coev['method']
            beta = str(method_coev['beta'])
            nsus = str(method_coev['nsus'])
            runs = str(method_coev['runs'])
            sufix = "sequences-beta" + beta + "-nsus" + nsus + "-runs" + runs
            curated_seq = curated_sequences_path + sufix + ".fasta.cluster"
            start_time = time.time()
            logging.info('Calculation of beta ' + beta + ' nsus ' + nsus + ' run ' + runs)
            df.at[index,'protein']=protein
            df.at[index,'pdb']=pdb_name
            df.at[index,'chain']=chain
            df.at[index,'method']=method
            df.at[index,'beta']=beta
            df.at[index,'nsus']=nsus
            df.at[index,'run']=runs
            try: 
                if(method=='mi'):
                    evol(pdb_clean_path, beta, runs, nsus, chain, scpe_sequences, clustered_sequences_path, curated_sequences_path, mi_data_path, di_data_path, frob_data_path, psicov_data_path,contact_map, contact_map_sync, coevolution_calc=False, syncronized=True)
                    util.sincronize_natural_evol_msas(clustered_sequences_path, curated_sequences_path, pattern, 2, -3)
                    util.sincronize_contact_map(contact_map, contact_map_sync, 2, 106)
                    mi_path = mi_data_path + "mi_" + sufix + ".csv"
                    msa.buslje09(curated_seq, mi_path)
                    target, scores = msa_analysis.getTargetScores(mi_path, contact_map_sync, constants.neighbour_window)
                    auc, auc01 = util.getAUC(target, scores)
                    df.at[index,'score']=auc
                    df.at[index,'score_01']=auc01
                if(method=='psicov'):
                    evol(pdb_clean_path, beta, runs, nsus, chain, scpe_sequences, clustered_sequences_path, curated_sequences_path, mi_data_path, di_data_path, frob_data_path, psicov_data_path,contact_map, contact_map_sync, coevolution_calc=False, syncronized=True)
                    util.sincronize_natural_evol_msas(clustered_sequences_path, curated_sequences_path, pattern, 2, -3)
                    util.sincronize_contact_map(contact_map, contact_map_sync, 2, 106)
                    psicov_path = psicov_data_path + "psicov_" + sufix + ".csv"
                    msa.create_msa_without_id(curated_seq, curated_seq+"_noid.aln")
                    msa.psicov(curated_seq+"_noid.aln",psicov_path) 
                    target_psicov, scores_psicov = msa_analysis.getTargetScores(psicov_path, contact_map_sync, constants.neighbour_window, split_char=' ',score_position=4)
                    auc_psicov, auc01_psicov = util.getAUC(target_psicov, scores_psicov)
                    df.at[index,'score']=auc_psicov
                    df.at[index,'score_01']=auc01_psicov
                #dca_di = di_data_path + "di_" + sufix + ".csv" 
                #dca_Frob = frob_data_path + "frob_" + sufix + ".csv"
                #msa.gaussDcaFrobenius(curated_seq, dca_Frob)
                #msa.gaussDcaDirectInformation(curated_seq, dca_di)
                #target_di, scores_di = msa_analysis.getTargetScores(dca_di, contact_map_sync, constants.neighbour_window, split_char=' ')
                #auc_di, auc01_di = util.getAUC(target_di, scores_di)
                #target_frob, scores_frob = msa_analysis.getTargetScores(dca_Frob, contact_map_sync, constants.neighbour_window, split_char=' ')
                #auc_frob, auc01_frob = util.getAUC(target_frob, scores_frob)
            except Exception as inst:
                print inst
                x = inst.args
                print x
                df.at[index,'message']="error_running_optimization"
                logging.error('Error with beta ' + beta + ' nsus ' + nsus + ' run ' + runs)        
            df.at[index,'execution_time']= time.time() - start_time
            index = index + 1
            df.to_csv(result_path)    
        df.to_csv(result_path)     

def analyse_confomers(execution_folder, structures, natural_result_path):
    for pdb_name in structures:
        conf_path = execution_folder + pdb_name + "/"
        curated_sequences_path = conf_path + "curated_sequences/"
        mi_data_path = curated_sequences_path + "mi/"
        #di_data_path = curated_sequences_path + "di/"
        #frob_data_path = curated_sequences_path + "frob/"
        psicov_data_path = curated_sequences_path + "psicov/"
        contact_map_sync = conf_path + "contact_map_sync.dat" 
        
        natural_mi_result_path = natural_result_path + "mi.csv"
        #natural_frob_result_path = natural_result_path + "frob.csv"
        #natural_di_result_path= natural_result_path + "di.csv"
        natural_psicov_result_path = natural_result_path + "psicov.csv"
        
        target_nat, scores_nat = msa_analysis.getTargetScores(natural_mi_result_path, contact_map_sync, constants.neighbour_window)
        target_psicov_nat, scores_psicov_nat = msa_analysis.getTargetScores(natural_psicov_result_path, contact_map_sync, constants.neighbour_window, split_char=' ',score_position=4)
        
        result_path =  curated_sequences_path + "results.csv"
        results_conf = pandas.read_csv(result_path, header=0, index_col=0)
        
        coevolution_analisys_df = pandas.DataFrame()
        coevolution_analisys_df.set_value(1,'top',0.5)
        coevolution_analisys_df.set_value(2,'top',1)
        coevolution_analisys_df.set_value(3,'top',2)
        coevolution_analisys_df.set_value(4,'top',3)
        coevolution_analisys_df.set_value(5,'top',4)
        coevolution_analisys_df.set_value(6,'top',5)
        
        
        sc = []
        targets=[]
        colors = ['green','orange','blue','red']
        labels = ['Natural MI', 'Evol MI', 'Natural PSICOV', 'Evol PSICOV']
        for index_opti, conf_result in results_conf.iterrows():
            beta = str(conf_result['beta'])
            nsus = str(conf_result['nsus'])
            runs = str(int(conf_result['run']))
            sufix = "sequences-beta" + beta + "-nsus" + nsus + "-runs" + runs
            if(conf_result['method']=='mi'):
                mi_path = mi_data_path + "mi_" + sufix + ".csv"
                target, scores = msa_analysis.getTargetScores(mi_path, contact_map_sync, constants.neighbour_window)
                sc.append(scores_nat)
                targets.append(target_nat)
                sc.append(scores)
                targets.append(target)
                
                #MI 
                mi_natural = util.load_zmip(natural_mi_result_path, constants.neighbour_window)
                util.order(mi_natural)
                mi_evol = util.load_zmip(mi_path , constants.neighbour_window)
                util.order(mi_evol)
                
                msa_analysis.coevolution_analisys('MI',coevolution_analisys_df, 1, mi_natural, mi_evol, curated_sequences_path, contact_map_sync,  constants.neighbour_window, pdb_name)
            if(conf_result['method']=='psicov'):
                psicov_path = psicov_data_path + "psicov_" + sufix + ".csv"
                target_psicov, scores_psicov = msa_analysis.getTargetScores(psicov_path, contact_map_sync, constants.neighbour_window, split_char=' ',score_position=4)
                sc.append(scores_psicov_nat)
                targets.append(target_psicov_nat)
                sc.append(scores_psicov)
                targets.append(target_psicov)
                
                plot.roc_curve_(targets,sc,labels,'PDB:' + pdb_name +  ' - Best Roc Curve Evol vs Natural for MI and PSICOV',colors, '',curated_sequences_path+'best_psicov_roc.png')
                
                psicov_natural = util.load_zmip(natural_psicov_result_path, constants.neighbour_window,split_char=' ')
                psicov_evol = util.load_zmip(psicov_path, constants.neighbour_window,split_char=' ')
                
                [r.pop(2) for r in psicov_evol]
                [r.pop(2) for r in psicov_evol]
                [r.pop(2) for r in psicov_natural]
                [r.pop(2) for r in psicov_natural]
                
                msa_analysis.coevolution_analisys('PSICOV',coevolution_analisys_df, 1, psicov_natural, psicov_evol, curated_sequences_path, contact_map_sync,  constants.neighbour_window, pdb_name)
                
        coevolution_analisys_df.to_csv(curated_sequences_path +'tops_contact_threashold_1csv')    
                
                
def run_families_evol():
    logging.info('Begin of the execution process')
    start_time = time.time()
    for family_folder in os.listdir(input_families_folder):
        logging.info(family_folder)
        msa_gz_path = glob.glob(input_families_folder + family_folder + "/*.final.gz")
        if(len(msa_gz_path) == 0):
            logging.error('No existe alineamiento de la familia ' + family_folder)
            return
        if(len(msa_gz_path) > 1):
            logging.error('Existe mas de un alineamiento de la familia ' + family_folder)
            return
        msa_gz_path = msa_gz_path[0]
        aux_path = msa_gz_path.split('/')
        family_pdb_information = aux_path[0] + "/" + aux_path[1] + "/" + aux_path[2] + "/" + aux_path[2] + "_pdb_level.csv" 
        family_pdb_evol_info_path = aux_path[0] + "/" + aux_path[1] + "/" + aux_path[2] + "/" + aux_path[2] + "_evol_info.csv"
        
        if(not os.path.isfile(family_pdb_evol_info_path)):
            pdb_to_evol_df = util.find_pdb_to_evolve(family_pdb_information)
            pdb_to_evol_df.to_csv(family_pdb_evol_info_path)
        else:
            pdb_to_evol_df = pandas.read_csv(family_pdb_evol_info_path, header=0, index_col='cluster')    
        
        if(execute_download_pdbs):
            pdb.download_pdbs_df(input_families_folder, family_folder, pdb_to_evol_df)
        
        if(execute_family_evol):
            family_evol(input_families_folder, family_folder, pdb_to_evol_df, family_pdb_evol_info_path)
            pdb_to_evol_df.to_csv(family_pdb_evol_info_path)

        if(execute_joined_pdb_analisys):
            pdb_to_compare = pdb_to_evol_df.loc[pdb_to_evol_df['status'] == 'okey']
            
            # msa_analysis.comparative_conservation(input_families_folder + family_folder, family_folder, pdb_to_compare)
            
            # msa_analysis.sum_contact_map(input_families_folder + family_folder, pdb_to_compare)
            
            # msa_analysis.comparative_mi_information(input_families_folder + family_folder, family_folder , 1 ,constants.neighbour_window, pdb_to_compare)  
            
            msa_analysis.compute_joined_msas(input_families_folder + family_folder , pdb_to_compare, input_families_folder + family_folder + "/prob_contact_map.dat")
            
    logging.info('End of the execution process')
    logging.info('Time of Execution --- %s seconds ---' % (time.time() - start_time))

        
def family_evol(input_families_folder, family_folder, pdb_to_evol_df, family_pdb_evol_info_path):
    start_time = time.time()
    try:
        logging.info('Family Evol ' + family_folder)
        # todo manejar errores
        msa_gz_path = glob.glob(input_families_folder + family_folder + "/*.final.gz")
        if(len(msa_gz_path) == 0):
            logging.error('No existe alineamiento de la familia ' + family_folder)
            return
        if(len(msa_gz_path) > 1):
            logging.error('Existe mas de un alineamiento de la familia ' + family_folder)
            return
        
        msa_gz_path = msa_gz_path[0]
        aux_path = msa_gz_path.split('/')
        msa_file_name_fasta = aux_path[0] + "/" + aux_path[1] + "/" + aux_path[2] + "/" + aux_path[2] + ".fasta"    
        zmip_natural_path = msa_file_name_fasta + "_zmip.dat"
        
        family_pdb_information = aux_path[0] + "/" + aux_path[1] + "/" + aux_path[2] + "/" + aux_path[2] + "_pdb_level.csv"    
        
        with gzip.open(msa_gz_path, 'rb') as f:
            aux_path = f.filename.split('/')
            msa_filename = os.path.basename(f.filename)
            msa_complete_filename_stock = aux_path[0] + "/" + aux_path[1] + "/" + aux_path[2] + "/" + msa_filename[:-3]
            msa_file = open(msa_complete_filename_stock , "w")
            file_content = f.read()
            msa_file.write(file_content)
            msa_file.flush()
            msa_file.close()
        
        msa.convertMSAToFasta(msa_complete_filename_stock, msa_file_name_fasta)
        
        if(execute_natural_mi_msa):
            msa.natural_msa_mi(msa_file_name_fasta, zmip_natural_path)
        # Natural Conservation Information
        if(execute_msa_natural_information):
            msa.msa_information(msa_file_name_fasta, msa_file_name_fasta, aux_path[2])
        
        # web_logo.create_web_logo(msa_file_name_fasta, msa_file_name_fasta + "_logo_sh.png",msa_file_name_fasta + "_data_sh.csv", 'png', aux_path[2], logo_type='SHANNON')
        # web_logo.create_web_logo(msa_file_name_fasta, msa_file_name_fasta + "_logo_kl.png",msa_file_name_fasta + "_data_kl.csv", 'png', aux_path[2], logo_type='KL')
        
        # pdb_paths_files = input_families_folder +  family_folder  +"/PDB/*.pdb.gz"
        pdb_paths_files = input_families_folder + family_folder + "/PDB/"
        
        # optimizacion
        optimized_family = False
        optimization_folder = input_families_folder + family_folder + "/optimization_folder/"
        optimization_file_path = optimization_folder + "optimization.csv"
        # si existe el arhivo de optimizacion entonces no hay que optimizar
        if(os.path.isfile(optimization_file_path)):
            optimized_family = True
        
        sufix = "_superimposed.pdb.gz"
        # pdb_to_evol_df = util.find_pdb_to_evolve(family_pdb_information)
        logging.info('Begin of the PDBs Evolution ')
        for index, pdb_protein_to_evolve in pdb_to_evol_df.iterrows():
            pdb_name = pdb_protein_to_evolve['pdb']
            file_name_pdb = pdb_protein_to_evolve['seq'].replace("/", "_").replace("-", "_") + "_" + pdb_protein_to_evolve['pdb'] + "_" + pdb_protein_to_evolve['chain'] + sufix
            complete_file_name_pdb = pdb_paths_files + file_name_pdb
            logging.info('Begin of the PDB ' + file_name_pdb)
            pdb_file_complete = pdb_paths_files + pdb_protein_to_evolve['pdb_folder_name'] + "/" + pdb_protein_to_evolve['pdb_folder_name'] + "_complete.pdb"
            pdb_file_complete_filename_to_evolve = pdb_paths_files + pdb_protein_to_evolve['pdb_folder_name'] + "/" + pdb_protein_to_evolve['pdb_folder_name'] + "_clean.pdb"
            # util.remove_header(pdb_file_complete_filename_to_evolve)
            protein = pdb_protein_to_evolve['seq']
            
            # for pdb_gz in glob.glob(pdb_paths_files):
            # aca arrancar otro try
            # unzip pdb and move to pdb folder
            try:
                with gzip.open(complete_file_name_pdb, 'rb') as f:
                    aux_path = f.filename.split('/')
                    pdb_file_name = os.path.basename(f.filename[:-3])
                    pdb_folder = aux_path[0] + "/" + aux_path[1] + "/" + aux_path[2] + "/" + aux_path[3] + "/" + pdb_file_name[:-17]
                    if not os.path.exists(pdb_folder):
                        os.makedirs(pdb_folder)
                    cutted_pdb_path = pdb_folder + "/" + pdb_file_name    
                    pdb_file = open(cutted_pdb_path , "w")
                    file_content = f.read()
                    pdb_file.write(file_content)
                    pdb_file.close()
                # chain name to evol
                chain_name = pdb_file_name[-18:-17]
                
                # the contact map will be created by scpe 
                contact_map = pdb_folder + "/contact_map.dat"
                
                contact_map_syncronized = pdb_folder + "/contact_map_sync.dat"
                # the folder to put de evol scpe sequences
                scpe_sequences = pdb_folder + "/scpe_sequences/"
                # the folder to put de evol clustered sequences
                clustered_sequences_path = pdb_folder + "/clustered_sequences/"
                # the folder to put de evol clustered sequences sincronized
                sincronized_evol_path = pdb_folder + "/sincronized_evol_path/"
                # MSA information
                msa_information_path = sincronized_evol_path + "conservation/"
                mi_data = pdb_folder + "/mi_data/"
                mi_data_analisys = mi_data + "info/"
                if not os.path.exists(scpe_sequences):
                    os.makedirs(scpe_sequences)
                if not os.path.exists(clustered_sequences_path):
                    os.makedirs(clustered_sequences_path)
                if not os.path.exists(sincronized_evol_path):
                    os.makedirs(sincronized_evol_path)
                if not os.path.exists(msa_information_path):
                    os.makedirs(msa_information_path)    
                if not os.path.exists(mi_data):
                    os.makedirs(mi_data)
                if not os.path.exists(mi_data_analisys):
                    os.makedirs(mi_data_analisys)
                if not os.path.exists(optimization_folder):
                    os.makedirs(optimization_folder)
                 
                cd_secuence = util.getSequence(msa_file_name_fasta, "CRBB1_HUMAN/150-232")
                # start_residue = int(pdb_protein_to_evolve['seq'][pdb_protein_to_evolve['seq'].index("/")+1:pdb_protein_to_evolve['seq'].index("-")]) 
                # end_residue = int(pdb_protein_to_evolve['seq'][pdb_protein_to_evolve['seq'].index("-")+1:]) 
                if(pdb_protein_to_evolve['status'] != 'okey'):
                    residue_position, residue_name = util.getPDBSequence(pdb_name, pdb_file_complete, chain_name)
                    print residue_position
                    print residue_name
                    util.MusclePairAlign("pdb", ''.join(residue_name), "cd_sequence", cd_secuence)
                    start_residue, end_residue = util.find_pdb_start_end_for_protein(msa_complete_filename_stock, pdb_protein_to_evolve['seq'], pdb_name, chain_name)
                # pdb_to_evol_df.set_value(index,"start_residue",start_residue)
                # pdb_to_evol_df.set_value(index,"end_residue",end_residue)
                    
                util.clean_pdb(pdb_file_complete, pdb_file_complete_filename_to_evolve, chain_name, start_residue, end_residue)   
                
                if(optimized_scpe_variables and not optimized_family):
                    optimize_evolution(protein, pdb_name, optimization_folder, optimization_file_path, pdb_file_complete_filename_to_evolve, cutted_pdb_path, chain_name, start_residue, end_residue)
                    optimized_family = True
                    
                optimization_df = pandas.read_csv(optimization_file_path)
                max_scpe_params = optimization_df.ix[optimization_df['auc'].idxmax()]
                beta = max_scpe_params['beta']
                nsus = max_scpe_params['nsus']
                runs = max_scpe_params['run']
                
                if(pdb_protein_to_evolve['status'] != 'okey'):
                    pdb_to_evol_df.set_value(index, "beta", beta)
                    pdb_to_evol_df.set_value(index, "nsus", nsus)
                    pdb_to_evol_df.set_value(index, "runs", runs)
                    pdb_to_evol_df.set_value(index, "start_pdb_residue", start_residue)
                    pdb_to_evol_df.set_value(index, "end_pdb_residue", end_residue)
                    evol_protein(pdb_to_evol_df, index, pdb_file_complete_filename_to_evolve, cutted_pdb_path, beta, runs, nsus, chain_name, scpe_sequences, clustered_sequences_path, sincronized_evol_path, zmip_natural_path, mi_data, mi_data_analisys , msa_information_path, contact_map, contact_map_syncronized)
                
                pdb_to_evol_df.set_value(index, "status", "okey")
                logging.info('End of the PDB ' + file_name_pdb)
            except Exception as inst:
                print inst
                x = inst.args
                print x
                logging.error('The PDB was not evolutionated ' + file_name_pdb)
                logging.error(inst)
                pdb_to_evol_df.set_value(index, "status", "error")
            pdb_to_evol_df.to_csv(family_pdb_evol_info_path)    
    except Exception as inst:
        print inst
        x = inst.args
        print x
        logging.error('The family was not evolutionated  ' + family_folder)
    logging.info('End family Evol  ' + family_folder)
    logging.info('--- %s seconds ---' % (time.time() - start_time))   



  



def evol_protein(data_frame_evol, index, pdb_file_complete_filename_to_evolve, cutted_pdb_path, beta, runs, nsus, chain, scpe_sequences, clustered_sequences_folder_path, sincronized_evol_path, zmip_natural_result_path, mi_data_path, mi_data_analisys, msa_conservation_path, contact_map_path, contact_map_sync):    
    start_time = time.time()
    beta = str(beta)
    nsus = str(nsus)
    runs = str(runs)
    sufix = "sequences-beta" + beta + "-nsus" + nsus + "-runs" + runs
    file_name = sufix + ".fasta"
    output_msa_path = scpe_sequences + file_name
    clustered_sequences_path = clustered_sequences_folder_path + file_name + ".cluster"
    clustered_tmp_sequences_path = clustered_sequences_path + ".clstr"
    sincronized_evol_file_path = sincronized_evol_path + file_name 
    mi_data_path = mi_data_path + "zmip_" + sufix + ".csv"
    scpe.run(pdb_file_complete_filename_to_evolve, beta, runs, nsus, chain, output_msa_path, contact_map_path)
    count_seq_scpe = msa.count_sequences(output_msa_path)
    data_frame_evol.set_value(index, 'count_seq_scpe', count_seq_scpe)
    msa.clustering_singular("0.62", output_msa_path, clustered_sequences_path)
    count_seq_cluster = msa.count_sequences(clustered_sequences_path)
    data_frame_evol.set_value(index, 'count_seq_cluster', count_seq_cluster)
    seq_lenght = msa.count_aminoacids(clustered_sequences_path)
    data_frame_evol.set_value(index, 'seq_lenght', seq_lenght)
    util.delete_files(output_msa_path)
    util.delete_files(clustered_tmp_sequences_path)
    # se realiza la optimizacion sobre el msa ya recortado
    util.synchronize_evol_with_cutted_pdb_singular(pdb_file_complete_filename_to_evolve, cutted_pdb_path, clustered_sequences_path, sincronized_evol_file_path, contact_map_path, contact_map_sync)
    seq_cutted_lenght = msa.count_aminoacids(sincronized_evol_file_path)
    data_frame_evol.set_value(index, 'seq_cutted_lenght', seq_cutted_lenght)
    # util.delete_files(clustered_sequences_path)
    msa_analysis.evol_analisys(sincronized_evol_file_path, mi_data_path, msa_conservation_path + sufix, file_name)
    msa_analysis.run_analisys_singular(data_frame_evol, index, zmip_natural_result_path, mi_data_path, contact_map_sync, mi_data_analisys, constants.neighbour_window)
    data_frame_evol.set_value(index, 'execution_time', time.time() - start_time)




'''Evoluciona una estructura de la thio_ecoli'''


def evol_protein_for_thio_ecoli_conf(pdb_file_complete_filename_to_evolve, beta, runs, nsus, chain, scpe_sequences, clustered_sequences_folder_path, curated_sequences_folder_path, mi_data_path, contact_map_path, contact_map_sync):    
    sufix = "sequences-beta" + beta + "-nsus" + nsus + "-runs" + runs
    file_name = sufix + ".fasta"
    output_msa_path = scpe_sequences + file_name
    clustered_sequences_path = clustered_sequences_folder_path + file_name + ".cluster"
    clustered_tmp_sequences_path = clustered_sequences_path + ".clstr"
    curated_sequences_path = curated_sequences_folder_path + file_name + ".cluster"
    mi_data_path = mi_data_path + "zmip_" + sufix + ".csv"
    scpe.run(pdb_file_complete_filename_to_evolve, beta, runs, nsus, chain, output_msa_path, contact_map_path)
    count_seq_scpe = msa.count_sequences(output_msa_path)
    msa.clustering_singular("0.62", output_msa_path, clustered_sequences_path)
    count_seq_cluster = msa.count_sequences(clustered_sequences_path)
    seq_lenght = msa.count_aminoacids(clustered_sequences_path)
    util.delete_files(clustered_tmp_sequences_path)
    util.sincronize_natural_evol_msas(clustered_sequences_folder_path, curated_sequences_folder_path, pattern, 2, -3)
    util.sincronize_contact_map(contact_map_path, contact_map_sync, 2, 106)
    msa.buslje09(curated_sequences_path, mi_data_path)
    target, scores = msa_analysis.getTargetScores(mi_data_path, contact_map_sync, constants.neighbour_window)
    auc, auc01 = util.getAUC(target, scores)
    return auc, auc01, count_seq_scpe, count_seq_cluster, seq_lenght
    
    


###############################
##    Methods to be moved
#################################    
   



    
#plot_auc_optimization_dca()    
# plot_rocs()
# plot_roc_natural_2trx()

# plot_rocs()
# run_families_evol()


def evol_thio_ecoli_conformeros():
    execution_folder = '../THIO_ECOLI_4_107/'
    structures = ['1THO', '2H6X', '2TRX']
    structures = ['1SRX', '1XOA', '1XOB', '2H74', '1KEB', '2H6Z', '2H76']
    # structures = ['1SRX']#solo tiene carbonos alfa ??
    structures = ['1KEB']
    structures = ['2H6Z', '2H76']
    structures = ['1XOA', '1THO']
    structures = ['1XOB', '2H73', '2TIR']
    structures = ['2TRX', '1XOA', '1XOB', '2H74','1XOB_M5','1XOB_M7','1XOB_M16','1XOA_M17']
    chain = 'A'
    beta = '5.0'
    nsus = '15.0'
    runs = '20000'
    index = 1
    protein = 'THIO_ECOLI_4_107'
    columns = ["protein", "pdb", "chain", "beta", "nsus", "run", "auc", "auc_01", "execution_time"]
    df = pandas.DataFrame(columns=columns)
    for pdb_name in structures:
        pdb_folder = execution_folder + pdb_name + "/"
        pdb_complete_path = pdb_folder + pdb_name + "_complete.pdb"
        pdb_file_complete_filename_to_evolve = pdb_folder + pdb_name + "_clean.pdb"
        if not os.path.exists(pdb_folder):
            os.makedirs(pdb_folder)
            pdb_data = urllib.urlopen('http://files.rcsb.org/download/' + pdb_name + '.pdb').read()
            pdb_file = open(pdb_complete_path, "w")
            pdb_file.write(pdb_data)
            pdb_file.close()
            util.clean_pdb(pdb_complete_path, pdb_file_complete_filename_to_evolve, chain)    
        scpe_sequences = pdb_folder + "scpe_sequences/"
        clustered_sequences_path = pdb_folder + "clustered_sequences_path/"
        curated_sequences_path = pdb_folder + "curated_sequences_path/"
        mi_data_path = pdb_folder + "mi_data_path/"
        contact_map = pdb_folder + "contact_map.dat"
        contact_map_sync = pdb_folder + "contact_map_sync.dat"
        if not os.path.exists(scpe_sequences):
            os.makedirs(scpe_sequences)
        if not os.path.exists(clustered_sequences_path):
            os.makedirs(clustered_sequences_path)
        if not os.path.exists(curated_sequences_path):
            os.makedirs(curated_sequences_path)
        if not os.path.exists(mi_data_path):
            os.makedirs(mi_data_path)     
            
        start_time = time.time()
        logging.info('Calculation of beta ' + beta + ' nsus ' + nsus + ' run ' + runs)
        df.set_value(index, 'protein', protein)
        df.set_value(index, 'pdb', pdb_name)
        df.set_value(index, 'chain', chain)
        # df.set_value(index, 'start_residue', start_residue)
        # df.set_value(index, 'end_residue', end_residue)
        df.set_value(index, 'beta', beta)
        df.set_value(index, 'nsus', nsus)
        df.set_value(index, 'run', runs)
        try: 
            auc, auc01, count_seq_scpe, count_seq_cluster, seq_lenght = evol_protein_for_thio_ecoli_conf(pdb_file_complete_filename_to_evolve, beta, runs, nsus, chain, scpe_sequences, clustered_sequences_path, curated_sequences_path, mi_data_path, contact_map, contact_map_sync)
            df.set_value(index, 'auc', auc)
            df.set_value(index, 'auc_01', auc01)
            df.set_value(index, 'count_seq_scpe', count_seq_scpe)
            df.set_value(index, 'count_seq_cluster', count_seq_cluster)
            df.set_value(index, 'seq_lenght', seq_lenght)
            
        except Exception as inst:
            print inst
            x = inst.args
            print x
            df.set_value(index, 'auc', 'error')
            df.set_value(index, 'auc_01', 'error')
            logging.error('Error with beta ' + beta + ' nsus ' + nsus + ' run ' + runs)        
        df.set_value(index, 'execution_time', time.time() - start_time)    
        index = index + 1
        df.to_csv(pdb_folder + 'results.csv')      


#evol_thio_ecoli_conformeros()
            
def analisys_thio_ecoli_conformeros():
    execution_folder = '../THIO_ECOLI_4_107/'
    #structures = ['2TRX', '1XOA', '1XOB', '2H74', '1KEB', '2H6Z', '2H6X', '2H76', '1THO']
    structures = ['2TRX', '1XOA', '1XOB', '2H74','1XOB_M5','1XOB_M7','1XOB_M16','1XOA_M17']
    chain = 'A'
    beta = '5.0'
    nsus = '15.0'
    runs = '20000'
    df_result = pandas.DataFrame()
    msa_file_name_fasta = "../THIO_ECOLI_4_107_2TRX_A/natural/PF00085_THIO_ECOLI_reference.fasta"
    zmip_natural_result_file = "../THIO_ECOLI_4_107_2TRX_A/natural/zmip_PF00085_THIO_ECOLI_reference.csv"
    index = 1
    for pdb_name in structures:
        pdb_folder = execution_folder + pdb_name + "/"
        # conservation de la evolucion
        msa_conservation_path = pdb_folder + "conservation/"
        if not os.path.exists(msa_conservation_path):
            os.makedirs(msa_conservation_path)
        
        sufix = "sequences-beta" + beta + "-nsus" + nsus + "-runs" + runs
        msa_file = pdb_folder + "curated_sequences_path/" + sufix + ".fasta.cluster"
        mi_data_path = pdb_folder + "mi_data_path/"
        mi_data_path_file = mi_data_path + "zmip_" + sufix + ".csv"
        contact_map_sync = pdb_folder + "contact_map_sync.dat"
        msa_conservation_path = pdb_folder + "conservation/"
        df_result.set_value(index, 'protein', 'THIO_ECOLI/4_107')
        df_result.set_value(index, 'pdb', pdb_name)
        df_result.set_value(index, 'chain', chain)
        df_result.set_value(index, 'beta', beta)
        df_result.set_value(index, 'nsus', nsus)
        df_result.set_value(index, 'run', runs)
        count_seq_scpe = msa.count_sequences(pdb_folder + "scpe_sequences/" + sufix + ".fasta")
        count_seq_cluster = msa.count_sequences(msa_file)
        df_result.set_value(index, 'count_scpe', count_seq_scpe)
        df_result.set_value(index, 'count_seq_cluster', count_seq_cluster)
        try:
            msa.msa_information(msa_file, msa_conservation_path, sufix)  
            msa_analysis.run_analisys_singular(df_result, index, zmip_natural_result_file, mi_data_path_file, contact_map_sync, mi_data_path, constants.neighbour_window, pdb_name)     
        except Exception as inst:
            print inst
            x = inst.args
            print x
            df_result.set_value(index, 'auc', 'error')
            df_result.set_value(index, 'auc_01', 'error')
            logging.error('Error with beta ' + beta + ' nsus ' + nsus + ' run ' + runs)
        index = index + 1
        df_result.to_csv(pdb_folder + 'results.csv')
    df_result.to_csv(execution_folder + 'results_conformeros.csv')       


'''Analisis de los conformeros de forma separada.'''
#analisys_thio_ecoli_conformeros()                                     


'''Analisis de la matriz de contacto conjunta sumada'''


def analisys_contact_map_thio_ecoli_conformeros():
    execution_folder = '../THIO_ECOLI_4_107/'
    structures = ['2TRX', '1XOA', '1XOB', '2H74','1XOB_M5','1XOB_M7','1XOB_M16','1XOA_M17']
    contact_maps_paths = [execution_folder + pdb + '/contact_map_sync.dat' for pdb in structures]
    msa_analysis.contact_map_sum_prob(execution_folder, contact_maps_paths) 

#analisys_contact_map_thio_ecoli_conformeros()


'''Agrupamiento a traves de top mi de los conformeros y analisis de resultados '''


def analisys_top_mi_thio_ecoli_conformeros():
    execution_folder = '../THIO_ECOLI_4_107/'
    execution_folder_agrup = execution_folder + 'conjunction_mi/'
    if not os.path.exists(execution_folder_agrup):
        os.makedirs(execution_folder_agrup)
     
    structures = ['2TRX', '1XOA', '1XOB', '2H74','1XOB_M5','1XOB_M7','1XOB_M16','1XOA_M17']
    mi_paths = [execution_folder + pdb + '/mi_data_path/zmip_sequences-beta5.0-nsus15.0-runs20000.csv' for pdb in structures]
    contact_map_path = execution_folder + 'sum_contact_map.dat'
    zmip_natural_result_path = "../THIO_ECOLI_4_107_2TRX_A/natural/zmip_PF00085_THIO_ECOLI_reference.csv"
    
    msa_analysis.analisys_mi_with_contact_map(execution_folder_agrup, mi_paths, contact_map_path, zmip_natural_result_path, 0.5, constants.neighbour_window, sinchronize_with_natural=True)
    msa_analysis.analisys_mi_with_contact_map(execution_folder_agrup, mi_paths, contact_map_path, zmip_natural_result_path, 1, constants.neighbour_window, sinchronize_with_natural=True)
    msa_analysis.analisys_mi_with_contact_map(execution_folder_agrup, mi_paths, contact_map_path, zmip_natural_result_path, 2, constants.neighbour_window, sinchronize_with_natural=True)
    msa_analysis.analisys_mi_with_contact_map(execution_folder_agrup, mi_paths, contact_map_path, zmip_natural_result_path, 3, constants.neighbour_window, sinchronize_with_natural=True)
    msa_analysis.analisys_mi_with_contact_map(execution_folder_agrup, mi_paths, contact_map_path, zmip_natural_result_path, 4, constants.neighbour_window, sinchronize_with_natural=True)
    msa_analysis.analisys_mi_with_contact_map(execution_folder_agrup, mi_paths, contact_map_path, zmip_natural_result_path, 5, constants.neighbour_window, sinchronize_with_natural=True)
    # hacer el desarrollo para mostras la matriz de contactos con los mi top, que se hayan encontrado en al menos 4 msa de proteinas
    
    zmip_reference_result_path = execution_folder + '2TRX/mi_data_path/zmip_sequences-beta5.0-nsus15.0-runs20000.csv'
    
    zmip_prom_result_path = execution_folder + '/prom/zmip_prom.csv'
    
    zmip_evol_intersect_result_path = execution_folder_agrup + 'top_1_mi.csv'
    # zmip_evol_intersect_result_path = execution_folder + '/top_2_mi.csv'
    top_df = pandas.DataFrame()
    generic_top = '1'
    msa_analysis.top_rank_intersection(execution_folder_agrup, generic_top, contact_map_path, zmip_natural_result_path, zmip_evol_intersect_result_path, top_df, zmip_reference_result_path, zmip_prom_result_path, 1, constants.neighbour_window, contact_threshold=8, top_threshold=8, sinchronize_with_natural=True)
    msa_analysis.top_rank_intersection(execution_folder_agrup, generic_top, contact_map_path, zmip_natural_result_path, zmip_evol_intersect_result_path, top_df, zmip_reference_result_path, zmip_prom_result_path, 2, constants.neighbour_window, contact_threshold=7, top_threshold=7, sinchronize_with_natural=True)
    msa_analysis.top_rank_intersection(execution_folder_agrup, generic_top, contact_map_path, zmip_natural_result_path, zmip_evol_intersect_result_path, top_df, zmip_reference_result_path, zmip_prom_result_path, 3, constants.neighbour_window, contact_threshold=6, top_threshold=6, sinchronize_with_natural=True)
    msa_analysis.top_rank_intersection(execution_folder_agrup, generic_top, contact_map_path, zmip_natural_result_path, zmip_evol_intersect_result_path, top_df, zmip_reference_result_path, zmip_prom_result_path, 4, constants.neighbour_window, contact_threshold=5, top_threshold=5, sinchronize_with_natural=True)
    msa_analysis.top_rank_intersection(execution_folder_agrup, generic_top, contact_map_path, zmip_natural_result_path, zmip_evol_intersect_result_path, top_df, zmip_reference_result_path, zmip_prom_result_path, 5, constants.neighbour_window, contact_threshold=4, top_threshold=4, sinchronize_with_natural=True)
    msa_analysis.top_rank_intersection(execution_folder_agrup, generic_top, contact_map_path, zmip_natural_result_path, zmip_evol_intersect_result_path, top_df, zmip_reference_result_path, zmip_prom_result_path, 6, constants.neighbour_window, contact_threshold=3, top_threshold=3, sinchronize_with_natural=True)
    msa_analysis.top_rank_intersection(execution_folder_agrup, generic_top, contact_map_path, zmip_natural_result_path, zmip_evol_intersect_result_path, top_df, zmip_reference_result_path, zmip_prom_result_path, 7, constants.neighbour_window, contact_threshold=2, top_threshold=2, sinchronize_with_natural=True)
    msa_analysis.top_rank_intersection(execution_folder_agrup, generic_top, contact_map_path, zmip_natural_result_path, zmip_evol_intersect_result_path, top_df, zmip_reference_result_path, zmip_prom_result_path, 8, constants.neighbour_window, contact_threshold=1, top_threshold=1, sinchronize_with_natural=True)
    
    top_df.to_csv(execution_folder_agrup + 'top_1_information.csv')
    
    zmip_evol_intersect_result_path = execution_folder_agrup + 'top_2_mi.csv'
    # zmip_evol_intersect_result_path = execution_folder + '/top_2_mi.csv'
    top_df = pandas.DataFrame()
    generic_top = '2'
    msa_analysis.top_rank_intersection(execution_folder_agrup, generic_top, contact_map_path, zmip_natural_result_path, zmip_evol_intersect_result_path, top_df, zmip_reference_result_path, zmip_prom_result_path, 1, constants.neighbour_window, contact_threshold=8, top_threshold=8, sinchronize_with_natural=True)
    msa_analysis.top_rank_intersection(execution_folder_agrup, generic_top, contact_map_path, zmip_natural_result_path, zmip_evol_intersect_result_path, top_df, zmip_reference_result_path, zmip_prom_result_path, 2, constants.neighbour_window, contact_threshold=7, top_threshold=7, sinchronize_with_natural=True)
    msa_analysis.top_rank_intersection(execution_folder_agrup, generic_top, contact_map_path, zmip_natural_result_path, zmip_evol_intersect_result_path, top_df, zmip_reference_result_path, zmip_prom_result_path, 3, constants.neighbour_window, contact_threshold=6, top_threshold=6, sinchronize_with_natural=True)
    msa_analysis.top_rank_intersection(execution_folder_agrup, generic_top, contact_map_path, zmip_natural_result_path, zmip_evol_intersect_result_path, top_df, zmip_reference_result_path, zmip_prom_result_path, 4, constants.neighbour_window, contact_threshold=5, top_threshold=5, sinchronize_with_natural=True)
    msa_analysis.top_rank_intersection(execution_folder_agrup, generic_top, contact_map_path, zmip_natural_result_path, zmip_evol_intersect_result_path, top_df, zmip_reference_result_path, zmip_prom_result_path, 5, constants.neighbour_window, contact_threshold=4, top_threshold=4, sinchronize_with_natural=True)
    msa_analysis.top_rank_intersection(execution_folder_agrup, generic_top, contact_map_path, zmip_natural_result_path, zmip_evol_intersect_result_path, top_df, zmip_reference_result_path, zmip_prom_result_path, 6, constants.neighbour_window, contact_threshold=3, top_threshold=3, sinchronize_with_natural=True)
    msa_analysis.top_rank_intersection(execution_folder_agrup, generic_top, contact_map_path, zmip_natural_result_path, zmip_evol_intersect_result_path, top_df, zmip_reference_result_path, zmip_prom_result_path, 7, constants.neighbour_window, contact_threshold=2, top_threshold=2, sinchronize_with_natural=True)
    msa_analysis.top_rank_intersection(execution_folder_agrup, generic_top, contact_map_path, zmip_natural_result_path, zmip_evol_intersect_result_path, top_df, zmip_reference_result_path, zmip_prom_result_path, 8, constants.neighbour_window, contact_threshold=1, top_threshold=1, sinchronize_with_natural=True)
    
    top_df.to_csv(execution_folder_agrup + 'top_2_information.csv')
#analisys_top_mi_thio_ecoli_conformeros()    
    

'''Crea los msas sin ids para luego unirlos '''


def create_msa_without_id_thio_ecoli_conformeros():
    execution_folder = '../THIO_ECOLI_4_107/'
    structures = ['2TRX', '1XOA', '1XOB', '2H74', '1KEB', '2H6Z', '2H6X', '2H76']
    msas = [execution_folder + pdb + '/curated_sequences_path/sequences-beta5.0-nsus15.0-runs20000.fasta.cluster' for pdb in structures]
    for msa_ in msas:
        msa.create_msa_without_id(msa_, str(msa_) + '_no_seq_ids.fasta')
    
    
def create_msa_conjuntion_thio_ecoli_conformeros(num=1500):
    execution_folder = '../THIO_ECOLI_4_107/'
    structures = ['2TRX', '1XOA', '1XOB', '2H74', '1KEB', '2H6Z', '2H6X', '2H76']
    msas = [execution_folder + pdb + '/curated_sequences_path/sequences-beta5.0-nsus15.0-runs20000.fasta.cluster_no_seq_ids.fasta' for pdb in structures]
    msa_conjuntion_bootstrap_path = execution_folder + 'msa_conjuntion_' + str(num) + '/'
    if not os.path.exists(msa_conjuntion_bootstrap_path):
        os.makedirs(msa_conjuntion_bootstrap_path)
    msa.create_msa_bootstrap(msas, msa_conjuntion_bootstrap_path + 'msa_bootstrap_' + str(num) + '.fasta', num)

    
def analisys_msa_conjuntion_thio_ecoli_conformeros(num=1500):
    execution_folder = '../THIO_ECOLI_4_107/'
    msa_conjuntion_bootstrap_path = execution_folder + 'msa_conjuntion_' + str(num) + '/'
    
    msa_bootstrap = msa_conjuntion_bootstrap_path + 'msa_bootstrap_' + str(num) + '.fasta'
    mi_data_output_path = msa_conjuntion_bootstrap_path + 'mi_data_path/'
    msa_conservation_path = msa_conjuntion_bootstrap_path + 'conservation/'
    if not os.path.exists(mi_data_output_path):
        os.makedirs(mi_data_output_path)
    if not os.path.exists(msa_conservation_path):
        os.makedirs(msa_conservation_path)
    mi_data_path = mi_data_output_path + "zmip_bootstrap.csv"
    msa_conservation_path = msa_conservation_path + 'bootstrap_'
    msa_bootstrap_clustered = msa_bootstrap + '_cluster_62.cluster'
    msa.clustering_singular("0.62", msa_bootstrap, msa_bootstrap_clustered)
    msa_analysis.evol_analisys(msa_bootstrap_clustered, mi_data_path, msa_conservation_path, 'Bootstrap THIO_ECOLI Conformeros')
 
 
def analisys_singular_conjunction_thio_ecoli_conformeros(num=1500):
    execution_folder = '../THIO_ECOLI_4_107/msa_conjuntion_' + str(num) + '/'
    mi_data_output_path = execution_folder + 'mi_data_path/'
    mi_data_path_file = mi_data_output_path + "zmip_bootstrap.csv"
    zmip_natural_result_file = "../THIO_ECOLI_4_107_2TRX_A/natural/zmip_PF00085_THIO_ECOLI_reference.csv"
    contact_map_path = '../THIO_ECOLI_4_107/sum_contact_map.dat'
    # zmip_reference_result_path = '../THIO_ECOLI_4_107/2TRX/mi_data_path/zmip_sequences-beta5.0-nsus15.0-runs20000.csv'
    top_df = pandas.DataFrame()
    pdb_name = 'THIO_ECOLI CONJUNCTION CONFORMEROS'
    msa_analysis.run_analisys_singular(top_df, 1, zmip_natural_result_file, mi_data_path_file, contact_map_path, mi_data_output_path, constants.neighbour_window, pdb_name) 
    # msa_analysis.top_rank_intersection(execution_folder, contact_map_path,zmip_natural_result_path, mi_data_path_file, top_df, zmip_reference_result_path, index=1, constants.neighbour_window, contact_threshold=1, top_threshold=1, sinchronize_with_natural=True)
    
    top_df.to_csv(execution_folder + 'result_conjunction.csv')


'''Realiza los calculos promediando la informacion mutua obtenida de todas las evoluciones de los conformeros'''


def analisys_prom_zmip_thio_ecoli_conformeros():
    execution_folder = '../THIO_ECOLI_4_107/'
    structures = ['2TRX', '1XOA', '1XOB', '2H74','1XOB_M5','1XOB_M7','1XOB_M16','1XOA_M17']
    mi_paths = [execution_folder + pdb + '/mi_data_path/zmip_sequences-beta5.0-nsus15.0-runs20000.csv' for pdb in structures]
    # contact_map_path = execution_folder + 'sum_contact_map.dat'
    # zmip_natural_result_path = "../THIO_ECOLI_4_107_2TRX_A/natural/zmip_PF00085_THIO_ECOLI_reference.csv"
    
    execution_prom_folder = execution_folder + 'prom/'
    mi_prom_result = execution_prom_folder + 'zmip_prom.csv'
    if not os.path.exists(execution_prom_folder):
        os.makedirs(execution_prom_folder)
    
    msa_analysis.prom_zmip(mi_paths, mi_prom_result, constants.neighbour_window)
    
    zmip_natural_result_file = "../THIO_ECOLI_4_107_2TRX_A/natural/zmip_PF00085_THIO_ECOLI_reference.csv"
    contact_map_path = '../THIO_ECOLI_4_107/sum_contact_map.dat'
    
    top_df = pandas.DataFrame()
    pdb_name = 'THIO_ECOLI PROM CONFORMEROS'
    
    msa_analysis.run_analisys_singular(top_df, 1, zmip_natural_result_file, mi_prom_result, contact_map_path, execution_prom_folder, constants.neighbour_window, pdb_name) 
    top_df.to_csv(execution_prom_folder + 'result_prom.csv')
 
#analisys_prom_zmip_thio_ecoli_conformeros()
'''Promedios de informacion mutua analisis conjunto''' 
# promedio thio_ecoli
 
# analisys_prom_zmip_thio_ecoli_conformeros()   
   
# conjuncion 

   
def conjunction_analisys(num):
   
    create_msa_without_id_thio_ecoli_conformeros()

    create_msa_conjuntion_thio_ecoli_conformeros(num)    
    
    analisys_msa_conjuntion_thio_ecoli_conformeros(num)

    analisys_singular_conjunction_thio_ecoli_conformeros(num)
    
# conjunction_analisys(1500)    
# conjunction_analisys(3000)
# conjunction_analisys(5000)
# conjunction_analisys(10000)
# top thio_ecoli

# analisys_contact_map_thio_ecoli_conformeros()


'''Analisis agrupamiento con top de conformeros'''
# analisys_top_mi_thio_ecoli_conformeros()

# draw_contact_maps()

# pdb.rms_list(unit_prot_id='P0AA25',reference='2TRX')
# util.clean_pdb('../THIO_ECOLI_4_107/all_structures/1THO.pdb', '../THIO_ECOLI_4_107/all_structures/1THO_clean.pdb', 'A')
# r = pdb.align_pdb('../THIO_ECOLI_4_107/2TRX/2TRX_clean.pdb', '../THIO_ECOLI_4_107/2TRX/2TRX_clean.pdb')
# print r


def resultados_top_mi_comparation():
    execution_folder = '../THIO_ECOLI_4_107/'
    execution_folder_smi = execution_folder + 'conjunction_mi/'
    if not os.path.exists(execution_folder_smi):
        os.makedirs(execution_folder_smi)
     
    structures = ['2TRX', '1XOA', '1XOB', '2H74','1XOB_M5','1XOB_M7','1XOB_M16','1XOA_M17']
    mi_paths = [execution_folder + pdb + '/mi_data_path/zmip_sequences-beta5.0-nsus15.0-runs20000.csv' for pdb in structures]
    contact_map_path = execution_folder + 'sum_contact_map.dat'
    zmip_natural_result_path = "../THIO_ECOLI_4_107_2TRX_A/natural/zmip_PF00085_THIO_ECOLI_reference.csv"
    
    # hacer el desarrollo para mostras la matriz de contactos con los mi top, que se hayan encontrado en al menos 4 msa de proteinas
    
    zmip_reference_result_path = execution_folder + '2TRX/mi_data_path/zmip_sequences-beta5.0-nsus15.0-runs20000.csv'
    
    zmip_prom_result_path = execution_folder + 'prom/zmip_prom.csv'
    
    zmip_evol_intersect_result_path = execution_folder_smi + 'top_1_mi.csv'
    # zmip_evol_intersect_result_path = execution_folder + '/top_2_mi.csv'
    top_df = pandas.DataFrame()
    generic_top = '1'
    
    msa_analysis.top_rank_comparation(execution_folder, 0.5, contact_map_path, zmip_natural_result_path, zmip_evol_intersect_result_path, top_df, zmip_reference_result_path, zmip_prom_result_path, 1, constants.neighbour_window, contact_threshold=1, top_threshold=1, sinchronize_with_natural=True)
    msa_analysis.top_rank_comparation(execution_folder, 1, contact_map_path, zmip_natural_result_path, zmip_evol_intersect_result_path, top_df, zmip_reference_result_path, zmip_prom_result_path, 2, constants.neighbour_window, contact_threshold=1, top_threshold=1, sinchronize_with_natural=True)
    msa_analysis.top_rank_comparation(execution_folder, 2, contact_map_path, zmip_natural_result_path, zmip_evol_intersect_result_path, top_df, zmip_reference_result_path, zmip_prom_result_path, 3, constants.neighbour_window, contact_threshold=1, top_threshold=1, sinchronize_with_natural=True)
    msa_analysis.top_rank_comparation(execution_folder, 3, contact_map_path, zmip_natural_result_path, zmip_evol_intersect_result_path, top_df, zmip_reference_result_path, zmip_prom_result_path, 4, constants.neighbour_window, contact_threshold=1, top_threshold=1, sinchronize_with_natural=True)
    msa_analysis.top_rank_comparation(execution_folder, 4, contact_map_path, zmip_natural_result_path, zmip_evol_intersect_result_path, top_df, zmip_reference_result_path, zmip_prom_result_path, 5, constants.neighbour_window, contact_threshold=1, top_threshold=1, sinchronize_with_natural=True)
    msa_analysis.top_rank_comparation(execution_folder, 5, contact_map_path, zmip_natural_result_path, zmip_evol_intersect_result_path, top_df, zmip_reference_result_path, zmip_prom_result_path, 6, constants.neighbour_window, contact_threshold=1, top_threshold=1, sinchronize_with_natural=True)        
    
    msa_analysis.top_rank_comparation(execution_folder, 0.5, contact_map_path, zmip_natural_result_path, zmip_evol_intersect_result_path, top_df, zmip_reference_result_path, zmip_prom_result_path, 7, constants.neighbour_window, contact_threshold=4, top_threshold=1, sinchronize_with_natural=True)
    msa_analysis.top_rank_comparation(execution_folder, 1, contact_map_path, zmip_natural_result_path, zmip_evol_intersect_result_path, top_df, zmip_reference_result_path, zmip_prom_result_path, 8, constants.neighbour_window, contact_threshold=4, top_threshold=1, sinchronize_with_natural=True)
    msa_analysis.top_rank_comparation(execution_folder, 2, contact_map_path, zmip_natural_result_path, zmip_evol_intersect_result_path, top_df, zmip_reference_result_path, zmip_prom_result_path, 9, constants.neighbour_window, contact_threshold=4, top_threshold=1, sinchronize_with_natural=True)
    msa_analysis.top_rank_comparation(execution_folder, 3, contact_map_path, zmip_natural_result_path, zmip_evol_intersect_result_path, top_df, zmip_reference_result_path, zmip_prom_result_path, 10, constants.neighbour_window, contact_threshold=4, top_threshold=1, sinchronize_with_natural=True)
    msa_analysis.top_rank_comparation(execution_folder, 4, contact_map_path, zmip_natural_result_path, zmip_evol_intersect_result_path, top_df, zmip_reference_result_path, zmip_prom_result_path, 11, constants.neighbour_window, contact_threshold=4, top_threshold=1, sinchronize_with_natural=True)
    msa_analysis.top_rank_comparation(execution_folder, 5, contact_map_path, zmip_natural_result_path, zmip_evol_intersect_result_path, top_df, zmip_reference_result_path, zmip_prom_result_path, 12, constants.neighbour_window, contact_threshold=4, top_threshold=1, sinchronize_with_natural=True)        
    
    msa_analysis.top_rank_comparation(execution_folder, 0.5, contact_map_path, zmip_natural_result_path, zmip_evol_intersect_result_path, top_df, zmip_reference_result_path, zmip_prom_result_path, 13, constants.neighbour_window, contact_threshold=8, top_threshold=1, sinchronize_with_natural=True)
    msa_analysis.top_rank_comparation(execution_folder, 1, contact_map_path, zmip_natural_result_path, zmip_evol_intersect_result_path, top_df, zmip_reference_result_path, zmip_prom_result_path, 14, constants.neighbour_window, contact_threshold=8, top_threshold=1, sinchronize_with_natural=True)
    msa_analysis.top_rank_comparation(execution_folder, 2, contact_map_path, zmip_natural_result_path, zmip_evol_intersect_result_path, top_df, zmip_reference_result_path, zmip_prom_result_path, 15, constants.neighbour_window, contact_threshold=8, top_threshold=1, sinchronize_with_natural=True)
    msa_analysis.top_rank_comparation(execution_folder, 3, contact_map_path, zmip_natural_result_path, zmip_evol_intersect_result_path, top_df, zmip_reference_result_path, zmip_prom_result_path, 16, constants.neighbour_window, contact_threshold=8, top_threshold=1, sinchronize_with_natural=True)
    msa_analysis.top_rank_comparation(execution_folder, 4, contact_map_path, zmip_natural_result_path, zmip_evol_intersect_result_path, top_df, zmip_reference_result_path, zmip_prom_result_path, 17, constants.neighbour_window, contact_threshold=8, top_threshold=1, sinchronize_with_natural=True)
    msa_analysis.top_rank_comparation(execution_folder, 5, contact_map_path, zmip_natural_result_path, zmip_evol_intersect_result_path, top_df, zmip_reference_result_path, zmip_prom_result_path, 18, constants.neighbour_window, contact_threshold=8, top_threshold=1, sinchronize_with_natural=True)
    
    top_df.to_csv(execution_folder + 'result_mi_comparation.csv')
    
#resultados_top_mi_comparation()





def seq_to_logo():
    execution_folder = '../THIO_ECOLI_4_107/'
    structures = ['2TRX', '1XOA', '1XOB', '2H74','1XOB_M5','1XOB_M7','1XOB_M16','1XOA_M17']
    msas = [execution_folder + pdb + '/curated_sequences_path/sequences-beta5.0-nsus15.0-runs20000.fasta.cluster' for pdb in structures]
    msa.seq_to_logo(msas, structures)
#seq_to_logo()


def entropy_by_column():
    # alignment = AlignIO.read('fasta.txt', "fasta", alphabet=Alphabet.ProteinAlphabet())
    execution_folder = '../THIO_ECOLI_4_107/'
    structures = ['2TRX', '1XOA', '1XOB', '2H74','1XOB_M5','1XOB_M7','1XOB_M16','1XOA_M17']
    # structures = ['1KEB']
    msas = [execution_folder + pdb + '/curated_sequences_path/sequences-beta5.0-nsus15.0-runs20000.fasta.cluster' for pdb in structures]
    outputs = [execution_folder + pdb + '/conservation/information_content.csv' for pdb in structures]
    i = 0
    for m in msas:
        msa.summary(m, outputs[i], structures[i])
        i = i + 1
    # p_data= df.value_counts()/len(df) # calculates the probabilities
    # entropy=sc.stats.entropy(p_data)  # input probabilities to get the entropy 
    # return entropy
    
#entropy_by_column()


'''
Entropia del alienamiento natural pf00085 gapstrip thio_ecoli
'''


def natural_entropy_by_column():
    msa_natural = '../THIO_ECOLI_4_107_2TRX_A/PF00085_THIO_ECOLI_reference.fasta'
    msa.summary(msa_natural, '../THIO_ECOLI_4_107_2TRX_A/information_content_PF00085_THIO_ECOLI_reference.csv', 'PF00085_THIO_ECOLI_reference')

# natural_entropy_by_column()


'''
Calculo de la entropia media por columna 
'''


def entropy_by_column_media():
    execution_folder = '../THIO_ECOLI_4_107/'
    structures = ['2TRX', '1XOA', '1XOB', '2H74','1XOB_M5','1XOB_M7','1XOB_M16','1XOA_M17']
    # structures = ['1KEB']
    msas_summary = [execution_folder + pdb + '/conservation/information_content.csv' for pdb in structures]
    output = execution_folder + 'information_content_media.csv'
    msa.conservation_media(msas_summary, output)
    
#entropy_by_column_media()







'''
Entropy by Column Family 
'''
def entropy_by_column_family():
    file = '../FAMILY_PF00085/PF00085/PF00085_evol_info.csv'
    df = pandas.read_csv(file, header=0, usecols=['pdb','pdb_folder_name','auc', 'auc_nat', 'auc_01', 'auc_nat_01', 'status'])
    pdb_to_compare = df.loc[df['status'] == 'okey']
    structures = pdb_to_compare['pdb_folder_name'].tolist()
    execution_folder = '../FAMILY_PF00085/PF00085/PDB/'
    msas = [execution_folder + pdb + '/sincronized_evol_path/sequences-beta5.0-nsus15.0-runs20000.fasta' for pdb in structures]
    outputs = [execution_folder + pdb + '/sincronized_evol_path/conservation/information_content.csv' for pdb in structures]
    i = 0
    for m in msas:
        print ("BEGIN " + m)
        msa.summary(m, outputs[i], structures[i])
        i = i + 1
        print ("END " + m)
    print ("fin")    
#entropy_by_column_family()

'''
Calculo de la entropia media por columna para toda la familia
'''
def entropy_by_column_media_family():
    file = '../FAMILY_PF00085/PF00085/PF00085_evol_info.csv'
    df = pandas.read_csv(file, header=0, usecols=['pdb','pdb_folder_name','status'])
    pdb_to_compare = df.loc[df['status'] == 'okey']
    structures = pdb_to_compare['pdb_folder_name'].tolist()
    execution_folder = '../FAMILY_PF00085/PF00085/PDB/'
    
    msas_summary = [execution_folder + pdb + '/sincronized_evol_path/conservation/information_content.csv' for pdb in structures]
    output = '../FAMILY_PF00085/PF00085/information_content_media.csv'
    msa.conservation_media(msas_summary, output)
    
#entropy_by_column_media_family ()

def natural_entropy_by_column_family():
    msa_natural = '../FAMILY_PF00085/PF00085/PF00085.fasta'
    msa.summary(msa_natural, '../FAMILY_PF00085/PF00085/information_content_superimposed_PF00085.csv', 'PF00085_superimposed')

#natural_entropy_by_column_family()


 
'''Analisis de la matriz de contacto conjunta sumada'''


def analisys_contact_map_thio_ecoli_family():
    file = '../FAMILY_PF00085/PF00085/PF00085_evol_info.csv'
    df = pandas.read_csv(file, header=0, usecols=['pdb','pdb_folder_name','status'])
    pdb_to_compare = df.loc[df['status'] == 'okey']
    structures = pdb_to_compare['pdb_folder_name'].tolist()
    execution_folder = '../FAMILY_PF00085/PF00085/PDB/'
    
    contact_maps_paths = [execution_folder + pdb + '/contact_map_sync.dat' for pdb in structures]
    msa_analysis.contact_map_sum_prob('../FAMILY_PF00085/PF00085', contact_maps_paths)     

#analisys_contact_map_thio_ecoli_family()

    
'''Realiza los calculos promediando la informacion mutua obtenida de todas las evoluciones de los conformeros'''


def analisys_prom_zmip_thio_ecoli_family():
    file = '../FAMILY_PF00085/PF00085/PF00085_evol_info.csv'
    df = pandas.read_csv(file, header=0, usecols=['pdb','pdb_folder_name','status'])
    pdb_to_compare = df.loc[df['status'] == 'okey']
    structures = pdb_to_compare['pdb_folder_name'].tolist()
    execution_folder = '../FAMILY_PF00085/PF00085/PDB/'
    
    mi_paths = [execution_folder + pdb + '/mi_data/zmip_sequences-beta5.0-nsus15.0-runs20000.csv' for pdb in structures]
    # contact_map_path = execution_folder + 'sum_contact_map.dat'
    # zmip_natural_result_path = "../THIO_ECOLI_4_107_2TRX_A/natural/zmip_PF00085_THIO_ECOLI_reference.csv"
    
    execution_prom_folder = '../FAMILY_PF00085/PF00085/prom/'
    mi_prom_result = execution_prom_folder + 'zmip_prom.csv'
    if not os.path.exists(execution_prom_folder):
        os.makedirs(execution_prom_folder)
    
    msa_analysis.prom_zmip(mi_paths, mi_prom_result, 0)
    
    zmip_natural_result_file = '../FAMILY_PF00085/PF00085/PF00085.fasta_zmip.dat'
    
    contact_map_path = '../FAMILY_PF00085/PF00085/sum_contact_map.dat'
    
    top_df = pandas.DataFrame()
    pdb_name = 'THIO_ECOLI PROM FAMILY'
    
    msa_analysis.run_analisys_singular(top_df, 1, zmip_natural_result_file, mi_prom_result, contact_map_path, execution_prom_folder, 0, pdb_name,1) 
    msa_analysis.run_analisys_singular(top_df, 2, zmip_natural_result_file, mi_prom_result, contact_map_path, execution_prom_folder, 0, pdb_name,16)
    msa_analysis.run_analisys_singular(top_df, 3, zmip_natural_result_file, mi_prom_result, contact_map_path, execution_prom_folder, 0, pdb_name,32)
    msa_analysis.run_analisys_singular(top_df, 4, zmip_natural_result_file, mi_prom_result, contact_map_path, execution_prom_folder, 0, pdb_name,48)
    msa_analysis.run_analisys_singular(top_df, 5, zmip_natural_result_file, mi_prom_result, contact_map_path, execution_prom_folder, 0, pdb_name,64)
    top_df.to_csv(execution_prom_folder + 'result_prom.csv')
    
#analisys_prom_zmip_thio_ecoli_family()

   

'''

'''
def create_msa_without_id_thio_ecoli_family():
    file = '../FAMILY_PF00085/PF00085/PF00085_evol_info.csv'
    df = pandas.read_csv(file, header=0, usecols=['pdb','pdb_folder_name','status'])
    pdb_to_compare = df.loc[df['status'] == 'okey']
    structures = pdb_to_compare['pdb_folder_name'].tolist()
    execution_folder = '../FAMILY_PF00085/PF00085/PDB/'
    
    msas = [execution_folder + pdb + '/sincronized_evol_path/sequences-beta5.0-nsus15.0-runs20000.fasta' for pdb in structures]
    for msa_ in msas:
        msa.create_msa_without_id(msa_, str(msa_) + '_no_seq_ids.fasta')
        
def create_msa_conjuntion_thio_ecoli_family(num=1500):
    file = '../FAMILY_PF00085/PF00085/PF00085_evol_info.csv'
    df = pandas.read_csv(file, header=0, usecols=['pdb','pdb_folder_name','status'])
    pdb_to_compare = df.loc[df['status'] == 'okey']
    structures = pdb_to_compare['pdb_folder_name'].tolist()
    execution_folder = '../FAMILY_PF00085/PF00085/PDB/'
    
    msas = [execution_folder + pdb + '/sincronized_evol_path/sequences-beta5.0-nsus15.0-runs20000.fasta_no_seq_ids.fasta' for pdb in structures]
    msa_conjuntion_bootstrap_path = '../FAMILY_PF00085/PF00085/msa_conjuntion_' + str(num) + '/'
    if not os.path.exists(msa_conjuntion_bootstrap_path):
        os.makedirs(msa_conjuntion_bootstrap_path)
    msa.create_msa_bootstrap(msas, msa_conjuntion_bootstrap_path + 'msa_bootstrap_' + str(num) + '.fasta', num)

def analisys_msa_conjuntion_thio_ecoli_family(num=1500):
    execution_folder = '../FAMILY_PF00085/PF00085/'
    msa_conjuntion_bootstrap_path = execution_folder + 'msa_conjuntion_' + str(num) + '/'
    
    msa_bootstrap = msa_conjuntion_bootstrap_path + 'msa_bootstrap_' + str(num) + '.fasta'
    mi_data_output_path = msa_conjuntion_bootstrap_path + 'mi_data_path/'
    msa_conservation_path = msa_conjuntion_bootstrap_path + 'conservation/'
    if not os.path.exists(mi_data_output_path):
        os.makedirs(mi_data_output_path)
    if not os.path.exists(msa_conservation_path):
        os.makedirs(msa_conservation_path)
    mi_data_path = mi_data_output_path + "zmip_bootstrap.csv"
    msa_conservation_path = msa_conservation_path + 'bootstrap_'
    msa_bootstrap_clustered = msa_bootstrap + '_cluster_62.cluster'
    
    msa.clustering_singular("0.62", msa_bootstrap, msa_bootstrap_clustered)
    msa_analysis.evol_analisys(msa_bootstrap_clustered, mi_data_path, msa_conservation_path, 'Bootstrap Family PF00085')
    
    #msa_analysis.evol_analisys(msa_bootstrap, mi_data_path, msa_conservation_path, 'Bootstrap Family PF00085')

def analisys_singular_conjunction_thio_ecoli_family(num=1500):
    execution_folder = '../FAMILY_PF00085/PF00085/msa_conjuntion_' + str(num) + '/'
    mi_data_output_path = execution_folder + 'mi_data_path/'
    mi_data_path_file = mi_data_output_path + "zmip_bootstrap.csv"
    zmip_natural_result_file = '../FAMILY_PF00085/PF00085/PF00085.fasta_zmip.dat'
    contact_map_path = '../FAMILY_PF00085/PF00085/sum_contact_map.dat'
    # zmip_reference_result_path = '../THIO_ECOLI_4_107/2TRX/mi_data_path/zmip_sequences-beta5.0-nsus15.0-runs20000.csv'
    top_df = pandas.DataFrame()
    pdb_name = 'Family PF00085 conjunction'
    
    #msa_analysis.run_analisys_singular(top_df, 1, zmip_natural_result_file, mi_data_path_file, contact_map_path, mi_data_output_path, constants.neighbour_window, pdb_name) 
    
    msa_analysis.run_analisys_singular(top_df, 1, zmip_natural_result_file, mi_data_path_file, contact_map_path, mi_data_output_path, 0, pdb_name,1) 
    msa_analysis.run_analisys_singular(top_df, 2, zmip_natural_result_file, mi_data_path_file, contact_map_path, mi_data_output_path, 0, pdb_name,16)
    msa_analysis.run_analisys_singular(top_df, 3, zmip_natural_result_file, mi_data_path_file, contact_map_path, mi_data_output_path, 0, pdb_name,32)
    msa_analysis.run_analisys_singular(top_df, 4, zmip_natural_result_file, mi_data_path_file, contact_map_path, mi_data_output_path, 0, pdb_name,48)
    msa_analysis.run_analisys_singular(top_df, 5, zmip_natural_result_file, mi_data_path_file, contact_map_path, mi_data_output_path, 0, pdb_name,64)
    
    
    top_df.to_csv(execution_folder + 'result_conjunction.csv')


'''

'''
def conjunction_analisys_family(num):
   
    create_msa_without_id_thio_ecoli_family()

    create_msa_conjuntion_thio_ecoli_family(num)    
    
    analisys_msa_conjuntion_thio_ecoli_family(num)

    analisys_singular_conjunction_thio_ecoli_family(num)





def select_thio_ecoli_conf():
    structures_path = '../THIO_ECOLI_4_107/all_structures_2/'
    pdb.download_pdbs('P0AA25','2TRX', 'A', structures_path)
#select_thio_ecoli_conf()




def dca_pfam_natu_2trx_reference():
    msa_natural = '../THIO_ECOLI_4_107_2TRX_A/PF00085_THIO_ECOLI_reference.fasta'
    msa_natural_dca_Frob ='../THIO_ECOLI_4_107_2TRX_A/natural/PF00085_THIO_ECOLI_reference.fasta_dca_Frob.txt'
    msa_natural_dca_DI='../THIO_ECOLI_4_107_2TRX_A/natural/PF00085_THIO_ECOLI_reference.fasta_dca_DI.txt'
    msa.dca_gauss_frob(msa_natural, msa_natural_dca_Frob)
    msa.dca_gauss_di(msa_natural, msa_natural_dca_DI)
#dca_pfam_natu_2trx_reference()  


def dca_thio_ecoli_conformeros():
    execution_folder = '../THIO_ECOLI_4_107/'
    structures = ['2TRX', '1XOA', '1XOB', '2H74','1XOB_M5','1XOB_M7','1XOB_M16','1XOA_M17']
    # structures = ['1KEB']
    msas = [execution_folder + pdb + '/curated_sequences_path/sequences-beta5.0-nsus15.0-runs20000.fasta.cluster' for pdb in structures]
    outputs_frob = [execution_folder + pdb + '/dca/frob.txt' for pdb in structures]
    outputs_di = [execution_folder + pdb + '/dca/di.txt' for pdb in structures]
    i=0
    for m in msas:
        if not os.path.exists(execution_folder + structures[i] + '/dca/'):
            os.makedirs(execution_folder + structures[i] + '/dca/')
        msa.dca_gauss_frob(m, outputs_frob[i])
        msa.dca_gauss_di(m, outputs_di[i])
        i=i+1

'''
Calcula el DCA para los MSA evolucionados completos sin sincronizar. Mover y ademas agregar este calculo dentro de la optimizacion.
'''    
def dca_2trx_evol():
    # the contact map will be created by scpe 
    pdb_folder = "../THIO_ECOLI_4_107_2TRX_A"
    optimization_folder = pdb_folder + "/optimization_folder/"
    curated_sequences_path = optimization_folder + "clustered_sequences_path/"
    dca_path = curated_sequences_path + "dca/"
    
    optimization_folder = pdb_folder + "/optimization_folder/"
    optimization_file_path = optimization_folder + "optimization_data.csv"
    optimization_df = pandas.read_csv(optimization_file_path, header=0)
    contact_map_path = '../THIO_ECOLI_4_107_2TRX_A/optimization_folder/contact_map.dat'
    
    for index, row_optimization in optimization_df.iterrows():
        try:
            if(row_optimization['analysis']!='dca_okey'):
                beta = str(row_optimization['beta'])
                nsus = str(row_optimization['nsus'])
                runs = str(int(row_optimization['run']))
                sufix = "sequences-beta" + beta + "-nsus" + nsus + "-runs" + runs
                curated_seq = curated_sequences_path + sufix + ".fasta.cluster"
                dca_di = dca_path + "dca_di" + sufix + ".csv" 
                dca_Frob = dca_path + "dca_Frob" + sufix + ".csv"
                
                #msa.dca_gauss_frob(curated_seq, dca_Frob)
                #msa.dca_gauss_di(curated_seq, dca_di)
                
                target_di, scores_di = msa_analysis.getTargetScores(dca_di, contact_map_path, constants.neighbour_window,split_char=' ')
                auc_di_, auc01_di_ = util.getAUC(target_di, scores_di)
                
                target_frob, scores_frob = msa_analysis.getTargetScores(dca_Frob, contact_map_path, constants.neighbour_window,split_char=' ')
                auc_frob_, auc01_frob_ = util.getAUC(target_frob, scores_frob)
                
                optimization_df.set_value(index, 'auc_di_', auc_di_)
                optimization_df.set_value(index, 'auc01_di_', auc01_di_)
                optimization_df.set_value(index, 'auc_frob_', auc_frob_)
                optimization_df.set_value(index, 'auc01_frob_', auc01_frob_)
                optimization_df.set_value(index, 'analysis', 'dca_okey')
                
                optimization_df.to_csv(optimization_file_path)
        except Exception as inst:
            print inst
            x = inst.args
            print x
            optimization_df.set_value(index, 'analysis', 'error')
                
#dca_2trx_evol()  


#def auc_dca_evol_2trx():
    

'''
Calcula el DCA para los MSA evolucionados completos  sincronizado
'''    
def dca_2trx_evol_sync():
    # the contact map will be created by scpe 
    pdb_folder = "../THIO_ECOLI_4_107_2TRX_A"
    optimization_folder = pdb_folder + "/optimization_folder/"
    clustered_sequences_path = optimization_folder + "clustered_sequences_path/"
    curated_sequences_path = optimization_folder + "curated_sequences_path/"
    dca_path = curated_sequences_path + "dca/"
    
    optimization_folder = pdb_folder + "/optimization_folder/"
    optimization_file_path = optimization_folder + "optimization_data.csv"
    optimization_df = pandas.read_csv(optimization_file_path, header=0, index_col=0)
    contact_map_path = '../THIO_ECOLI_4_107_2TRX_A/optimization_folder/contact_map_sync.dat'
    
    #pattern=["sequences-beta0.1","sequences-beta0.01","sequences-beta0.001","sequences-beta0.25"]
    #util.sincronize_natural_evol_msas(clustered_sequences_path, curated_sequences_path, pattern, 2, -3)
    msa_conservation_path = dca_path + "conservation/"
    zmip_natural_result_path = "../THIO_ECOLI_4_107_2TRX_A/natural/PF00085_THIO_ECOLI_reference.fasta_dca_DI.txt"
    for index, row_optimization in optimization_df.iterrows():
        try:
            #if(row_optimization['auc_di']!=''):
            if(row_optimization['auc_di_']>=0.8):
                beta = str(row_optimization['beta'])
                nsus = str(row_optimization['nsus'])
                runs = str(int(row_optimization['run']))
                sufix = "sequences-beta" + beta + "-nsus" + nsus + "-runs" + runs
                
                curated_seq = curated_sequences_path + sufix + ".fasta.cluster"
                
                dca_di = dca_path + "dca_di" + sufix + ".csv" 
                dca_Frob = dca_path + "dca_Frob" + sufix + ".csv"
                
                #msa.dca_gauss_frob(curated_seq, dca_Frob)
                #msa.dca_gauss_di(curated_seq, dca_di)
                
                target_di, scores_di = msa_analysis.getTargetScores(dca_di, contact_map_path, zmip_natural_result_path,True,constants.neighbour_window,split_char=' ')
                auc_di, auc01_di = util.getAUC(target_di, scores_di)
                target_frob, scores_frob = msa_analysis.getTargetScores(dca_Frob, contact_map_path,zmip_natural_result_path,True, constants.neighbour_window,split_char=' ')
                auc_frob, auc01_frob = util.getAUC(target_frob, scores_frob)
                
                msa_analysis.run_analisys_singular_dca(optimization_df, index, zmip_natural_result_path, dca_di, contact_map_path, dca_path, constants.neighbour_window,'2TRX')
                
                #msa.msa_information(curated_seq, msa_conservation_path,sufix)
                
                
                 
                optimization_df.set_value(index, 'auc_di', auc_di)
                optimization_df.set_value(index, 'auc01_di', auc01_di)
                optimization_df.set_value(index, 'auc_frob', auc_frob)
                optimization_df.set_value(index, 'auc01_frob', auc01_frob)
                optimization_df.set_value(index, 'analysis', 'okey')
                
                optimization_df.to_csv(optimization_file_path)
        except Exception as inst:
            print inst
            x = inst.args
            print x
            optimization_df.set_value(index, 'analysis', 'error')
                
#dca_2trx_evol_sync()  

 
def dca_2trx_evol_analisys():
    # the contact map will be created by scpe 
    pdb_folder = "../THIO_ECOLI_4_107_2TRX_A"
    optimization_folder = pdb_folder + "/optimization_folder/"
    #clustered_sequences_path = optimization_folder + "clustered_sequences_path/"
    curated_sequences_path = optimization_folder + "curated_sequences_path/"
    dca_path = curated_sequences_path + "dca/"
    
    optimization_folder = pdb_folder + "/optimization_folder/"
    optimization_file_path = optimization_folder + "optimization_data.csv"
    optimization_df = pandas.read_csv(optimization_file_path, header=0, index_col=0)
    contact_map_sync = '../THIO_ECOLI_4_107_2TRX_A/optimization_folder/contact_map_sync.dat'
    
    #pattern=["sequences-beta0.1","sequences-beta0.01","sequences-beta0.001","sequences-beta0.25"]
    #util.sincronize_natural_evol_msas(clustered_sequences_path, curated_sequences_path, pattern, 2, -3)
    zmip_natural_result_path = "../THIO_ECOLI_4_107_2TRX_A/natural/zmip_PF00085_THIO_ECOLI_reference.csv"
    
    msa_conservation_path = dca_path + "conservation/"
    if not os.path.exists(msa_conservation_path):
        os.makedirs(msa_conservation_path)
    
    for index, row_optimization in optimization_df.iterrows():
        try:
            if(row_optimization['auc_di']>=0.8):
                beta = str(row_optimization['beta'])
                nsus = str(row_optimization['nsus'])
                runs = str(int(row_optimization['run']))
                sufix = "sequences-beta" + beta + "-nsus" + nsus + "-runs" + runs
                curated_seq = curated_sequences_path + sufix + ".fasta.cluster"
                dca_di = dca_path + "dca_di" + sufix + ".csv" 
                #dca_Frob = dca_path + "dca_Frob" + sufix + ".csv"
                #msa_analysis.run_analisys_singular_dca(optimization_df, index, zmip_natural_result_path, dca_di, contact_map_sync, dca_path, constants.neighbour_window,'2TRX')
                target_di, scores_di = msa_analysis.getTargetScores(dca_di, contact_map_sync, constants.neighbour_window,split_char=' ')
                auc_di, auc01_di = util.getAUC(target_di, scores_di)
                msa.msa_information(curated_seq, msa_conservation_path,sufix)
                optimization_df.set_value(index, 'analysis', 'okey')
                optimization_df.to_csv(optimization_file_path)
        except Exception as inst:
            print inst
            x = inst.args
            print x
            optimization_df.set_value(index, 'analysis', 'error')
#dca_2trx_evol_analisys()







    
    
    


          