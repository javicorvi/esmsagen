'''
Scripts for execute evolutions and for plot diferent information.
'''

import msa_analysis
import util
import plot
import evol
import constants
import msa
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

#analyse_optimus_msa()


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
#analyse_top_coevolution_conformers()

def plot_comparation_top():
    execution_folder = constants.data_path + "THIO_ECOLI_4_107_CONFORMERS/"
    top_result = execution_folder + "coevolution_results/top_thio_ecoli_conformers.csv"
    msa_analysis.plot_comparation_top(top_result, execution_folder + "coevolution_results/")
#plot_comparation_top()    

def conjunction_analisys():
    execution_folder = constants.data_path + "THIO_ECOLI_4_107_CONFORMERS/"
    structures = [ '1XOA', '1XOB', '2H74','1XOB_M5','1XOB_M7','1XOB_M16','1XOA_M17','2TRX']
    contact_map_path = execution_folder + 'sum_contact_map.dat'
    natural_result_path = constants.data_path + "THIO_ECOLI_4_107_2TRX_A/natural/"
    evol.conjunction_analysis(execution_folder, structures,contact_map_path,natural_result_path)
    
    #evol.analisys_singular_conjunction_thio_ecoli_conformeros(execution_folder, structures,num)
    
#conjunction_analisys()    
# conjunction_analisys(3000)
# conjunction_analisys(5000)
# conjunction_analisys(10000)


def conservation_conformeros():
    execution_folder = constants.data_path + "THIO_ECOLI_4_107_CONFORMERS/"
    structures = [ '1XOB']
    msas = [execution_folder + pdb + '/optimization/curated_sequences/sequences-beta5.0-nsus20.0-runs20000_copy.fasta' for pdb in structures]
    msa.seq_to_logo_msas(msas,structures)
    
    structures = [ '2TRX']
    msas = [execution_folder + pdb + '/optimization/curated_sequences/sequences-beta7.0-nsus20.0-runs20000_copy.fasta' for pdb in structures]
    msa.seq_to_logo_msas(msas,structures)
    
#conservation_conformeros()
def plot_conformers_contacts():
    execution_folder = constants.data_path + "THIO_ECOLI_4_107_CONFORMERS/"
    structures = [ '1XOA', '1XOB', '2H74','1XOB_M5','1XOB_M7','1XOB_M16','1XOA_M17','2TRX']
    contact_maps_paths = [execution_folder + pdb + '/contact_map_sync.dat' for pdb in structures]
    contact_map_path = execution_folder + 'sum_contact_map.dat'
    output = execution_folder + "contacts_shared_b_cof.png"
    cmap_sum = util.load_contact_map(contact_map_path)
    import numpy as np 
    mat = np.zeros((8, 8))
    sum_total=float(np.count_nonzero(cmap_sum>=1))
    prom = 0
    cant=0
    for i in range(len(structures)):
        for j in range((i+1), len(structures)):  
            print(structures[i])
            print(structures[j])
            total_pair=util.contacts_match(contact_maps_paths[i],contact_maps_paths[j])
            val=round(total_pair * 100 / sum_total,2)
            prom = prom + val
            mat[i,j]=val
            mat[j,i]=val
            cant=cant+1
    print prom / cant
    x=['', '1XOA', '1XOB', '2H74','1XOB_M5','1XOB_M7','1XOB_M16','1XOA_M17','2TRX']
    y=['', '1XOA', '1XOB', '2H74','1XOB_M5','1XOB_M7','1XOB_M16','1XOA_M17','2TRX']
    plot.conformers_contact(mat,x,y,output)       
#plot_conformers_contacts() 
    
    
def evol_family():
    execution_folder = constants.data_path + "FAMILIES/"
    families = [ 'PF00085']
    evol.evol_families(families, execution_folder)   
#evol_family()

def analisys_family():
    execution_folder = constants.data_path + "FAMILIES/"
    families = [ 'PF00085']
    evol.analisys_families(families, execution_folder)   
analisys_family()  
