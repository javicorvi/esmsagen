'''
Scripts for execute evolutions and for plot diferent information.
'''

import msa_analysis
import util
import plot
import evol
import constants

'''
2TRX evolution
'''
def evol_2trx():
    chain_name = "A"
    pdb = "2TRX"
    # the contact map will be created by scpe 
    execution_folder = constants.data_path + "THIO_ECOLI_4_107_2TRX_A/"
    evol.optimize_evolution(execution_folder,"THIO_ECOLI/4-107", pdb, chain_name)
#evol_2trx()

def plot_auc_optimization():
    execution_folder = constants.data_path + "THIO_ECOLI_4_107_2TRX_A/"
    optimization_results = execution_folder + "optimization/optimization.csv"    
    msa_analysis.plot_optimization(optimization_results)
#plot_auc_optimization()    
'''
Natural information of THIO_ECOLI
'''
def natural_information():       
    #natural information
    execution_folder = constants.data_path + "THIO_ECOLI_4_107_2TRX_A/"
    evol.natural_information(execution_folder)
#natural_information()    

def compare_evolution_2trx():
    chain_name = "A"
    pdb = "2TRX"
    # the contact map will be created by scpe 
    execution_folder = constants.data_path + "THIO_ECOLI_4_107_2TRX_A/"
    evol.compare_evolution_with_natural(execution_folder,"THIO_ECOLI/4-107", pdb, chain_name)
#compare_evolution_2trx()
def analyse_optimus_msa():       
    #natural information
    execution_folder = constants.data_path + "THIO_ECOLI_4_107_2TRX_A/"
    natural_result_path = execution_folder + "/natural/"
    #natural_information()  
    evol.analyse_optimus_msa(execution_folder,'2TRX',natural_result_path)
    
#analyze_optimus_msa()

def evol_thio_ecoli_conformers():
    chain_name = "A"
    protein = "THIO_ECOLI/4-107"
    execution_folder = constants.data_path + "THIO_ECOLI_4_107_CONFORMERS/"
    structures = ['1XOB_M16','1XOA_M17','2TRX', '1XOA', '1XOB', '2H74','1XOB_M5','1XOB_M7']
    evol.optimize_evolution_conformers(execution_folder,structures, protein,chain_name)
#evol_thio_ecoli_conformers() 

'''
def evol_thio_ecoli_conformers():
    execution_folder = constants.data_path + "THIO_ECOLI_4_107_CONFORMERS/"
    best_coev_results_paths = constants.data_path + "THIO_ECOLI_4_107_2TRX_A/optimization/curated_sequences_path/best_01_coev_values.csv"
    structures = ['2TRX', '1XOA', '1XOB', '2H74','1XOB_M5','1XOB_M7','1XOB_M16','1XOA_M17']
    evol.evol_conformers(execution_folder, best_coev_results_paths, structures)
'''    
#evol_thio_ecoli_conformeros()    

def analyse_thio_ecoli_confomers():
    execution_folder = constants.data_path + "THIO_ECOLI_4_107_CONFORMERS/"
    structures = [ '1XOA', '1XOB', '2H74','1XOB_M5','1XOB_M7','1XOB_M16','1XOA_M17','2TRX']
    natural_result_path = constants.data_path + "THIO_ECOLI_4_107_2TRX_A/natural/"
    evol.analyse_optimus_msa_conformers(execution_folder, structures, natural_result_path)

#analyse_thio_ecoli_confomers()   

def analyse_top_coevolution_conformers():
    execution_folder = constants.data_path + "THIO_ECOLI_4_107_CONFORMERS/"
    structures = [ '1XOA', '1XOB', '2H74','1XOB_M5','1XOB_M7','1XOB_M16','1XOA_M17','2TRX']
    contact_map_path = execution_folder + 'sum_contact_map.dat'
    natural_result_path = constants.data_path + "THIO_ECOLI_4_107_2TRX_A/natural/"
    evol.analyse_top_coevolution_conformers(execution_folder, structures,contact_map_path,natural_result_path)
    
    contact_map_dendogram_output = execution_folder + 'contact_map_dendogram.png'
    contact_maps_paths = [execution_folder + pdb + '/contact_map_sync.dat' for pdb in structures]
    msa_analysis.dendogram_matrix(contact_maps_paths, contact_map_dendogram_output,'Contact Map Clustering entre Conformeros',structures,'single')   
analyse_top_coevolution_conformers()

