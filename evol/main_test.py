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

def analisis_para_evolciones_espeficas():
    mi_data_path = '../THIO_ECOLI_4_107_2TRX_A/optimization_folder/curated_sequences_path/mi_data_path/'
    contact_map_sync = '../THIO_ECOLI_4_107_2TRX_A/optimization_folder/contact_map_sync.dat'
    zmip_natural_result_file = '../THIO_ECOLI_4_107_2TRX_A/natural/zmip_PF00085_THIO_ECOLI_reference.csv'
    sufix = "zmip_sequences-beta5.0-nsus15.0-runs20000"
    file_name = sufix + ".csv"
    df_result = pandas.DataFrame()
    msa_analysis.run_analisys_singular(df_result, 1, zmip_natural_result_file, mi_data_path + file_name, contact_map_sync, mi_data_path, constants.neighbour_window, 'Evol 2TRX Beta 5 Nsus 15')
    
    sufix = "zmip_sequences-beta5.0-nsus10.0-runs20000"
    file_name = sufix + ".csv"
    df_result = pandas.DataFrame()
    msa_analysis.run_analisys_singular(df_result, 2, zmip_natural_result_file, mi_data_path + file_name, contact_map_sync, mi_data_path, constants.neighbour_window, 'Evol 2TRX Beta 5 Nsus 10')
    
    sufix = "zmip_sequences-beta10.0-nsus15.0-runs20000"
    file_name = sufix + ".csv"
    df_result = pandas.DataFrame()
    msa_analysis.run_analisys_singular(df_result, 3, zmip_natural_result_file, mi_data_path + file_name, contact_map_sync, mi_data_path, constants.neighbour_window, 'Evol 2TRX Beta 10 Nsus 15')
    
    sufix = "zmip_sequences-beta2.0-nsus3.0-runs20000"
    file_name = sufix + ".csv"
    df_result = pandas.DataFrame()
    msa_analysis.run_analisys_singular(df_result, 4, zmip_natural_result_file, mi_data_path + file_name, contact_map_sync, mi_data_path, constants.neighbour_window, 'Evol 2TRX Beta 2 Nsus 3')
    
    
 

def draw_contact_maps():
    execution_folder = '../THIO_ECOLI_4_107/'
    structures = ['2TRX', '1XOA', '1XOB', '2H74','1XOB_M5','1XOB_M7','1XOB_M16','1XOA_M17']
    contact_maps_paths = [execution_folder + pdb + '/contact_map_sync.dat' for pdb in structures]
    i = 0
    for cmap_pth in contact_maps_paths:
        contact_map = util.load_contact_map(cmap_pth)
        plot.contact_map(contact_map, cmap_pth + '.png', structures[i] + ' Contact Map')
        i = i + 1
        
#draw_contact_maps()

'''
Plot roc cur familia pf00085 estructuras superpuestas.
'''
def plot_roc_natural_family_superpost():
    zmip_natural_result_file = '../FAMILY_PF00085/PF00085/PF00085.fasta_zmip.dat'
    contact_map_path = '../FAMILY_PF00085/PF00085/sum_contact_map.dat'
    
    target, scores = msa_analysis.getTargetScores(zmip_natural_result_file, contact_map_path, 0)
    sc = []
    targets = []
    sc.append(scores)
    targets.append(target)
    auc, auc01 = util.getAUC(target, scores)
    colors = ['blue']
    labels = ['Familia PF00085']
    output_file = '../FAMILY_PF00085/PF00085/natural_spuerpuesta_roc'
    plot.roc_curve_(targets, sc, labels, 'Familia PF00085 estructuras superpuestas' , colors, 'Runs', output_file)
    
#plot_roc_natural_family_superpost() 

def plot_roc_natural_2trx():
    
    contact_map_path = '../THIO_ECOLI_4_107_2TRX_A/optimization_folder/contact_map_sync.dat'
    output_file = '../THIO_ECOLI_4_107_2TRX_A/natural/natural_roc_mianddca'
    
    colors = ['green' , 'blue']
    labels = ['DCA DI', 'MI']
    sc = []
    targets=[]
    
    mi_data_path_DI = '../THIO_ECOLI_4_107_2TRX_A/natural/PF00085_THIO_ECOLI_reference.fasta_dca_DI.txt'
    target_di, scores_di = msa_analysis.getTargetScores(mi_data_path_DI, contact_map_path, constants.neighbour_window,split_char=' ')
    sc.append(scores_di)
    targets.append(target_di)
    
    '''mi_data_path_Frob = '../THIO_ECOLI_4_107_2TRX_A/natural/PF00085_THIO_ECOLI_reference.fasta_dca_Frob.txt' 
    target_Frob, scores_Frob = msa_analysis.getTargetScores(mi_data_path_Frob, contact_map_path, constants.neighbour_window,split_char=' ')
    sc.append(scores_Frob)
    targets.append(target_Frob)'''
    
    mi_data_path = '../THIO_ECOLI_4_107_2TRX_A/natural/zmip_PF00085_THIO_ECOLI_reference.csv'
    target_mi, scores_mi = msa_analysis.getTargetScores(mi_data_path, contact_map_path, constants.neighbour_window)
    sc.append(scores_mi)
    targets.append(target_mi)
    
    
    plot.roc_curve_(targets, sc, labels, 'MSA Natural ROC curves' ,colors, '',output_file)
    
    
def plot_roc_evol_2trx_dca():
    
    pdb_folder = "../THIO_ECOLI_4_107_2TRX_A"
    
    contact_map_path = '../THIO_ECOLI_4_107_2TRX_A/optimization_folder/contact_map_sync.dat'
    output_file = '../THIO_ECOLI_4_107_2TRX_A/optimization_folder/curated_sequences_path/dca/dca_sequences-beta5.0-nsus15.0-runs20000.png'
    
    colors = ['orange', 'green' , 'blue']
    labels = ['DCA DI','DCA Frob', 'MI']
    sc = []
    targets=[]
    
    mi_data_path_DI = '../THIO_ECOLI_4_107_2TRX_A/optimization_folder/curated_sequences_path/dca/dca_disequences-beta5.0-nsus15.0-runs20000.txt'
    target_di, scores_di = msa_analysis.getTargetScores(mi_data_path_DI, contact_map_path, constants.neighbour_window,split_char=' ')
    sc.append(scores_di)
    targets.append(target_di)
    
    mi_data_path_Frob = '../THIO_ECOLI_4_107_2TRX_A/optimization_folder/curated_sequences_path/dca/dca_Frobsequences-beta5.0-nsus15.0-runs20000.txt' 
    target_Frob, scores_Frob = msa_analysis.getTargetScores(mi_data_path_Frob, contact_map_path, constants.neighbour_window,split_char=' ')
    sc.append(scores_Frob)
    targets.append(target_Frob)
    
    mi_data_path = '../THIO_ECOLI_4_107_2TRX_A/optimization_folder/curated_sequences_path/mi_data_path/zmip_sequences-beta5.0-nsus15.0-runs20000.csv'
    target_mi, scores_mi = msa_analysis.getTargetScores(mi_data_path, contact_map_path, constants.neighbour_window)
    sc.append(scores_mi)
    targets.append(target_mi)
    
    plot.roc_curve_(targets, sc, labels, 'MSA  ROC curves' ,colors, '',output_file)

#plot_roc_evol_2trx_dca()


def dendogram_contact_maps():
    execution_folder = '../THIO_ECOLI_4_107/'
    output_path = execution_folder + 'contact_dendogram.png'
    structures = ['2TRX', '1XOA', '1XOB', '2H74','1XOB_M5','1XOB_M7','1XOB_M16','1XOA_M17']
    contact_maps_paths = [execution_folder + pdb + '/contact_map_sync.dat' for pdb in structures]
    msa_analysis.dendogram_matrix(contact_maps_paths,output_path,'Contact Map Clustering entre Conformeros',structures,'single')
#dendogram_contact_maps()


def dendogram_top_mi():
    execution_folder = '../THIO_ECOLI_4_107/'
    output_path = execution_folder + 'top_2_mi_dendogram.png'
    zmip_natural_result_path = "../THIO_ECOLI_4_107_2TRX_A/natural/zmip_PF00085_THIO_ECOLI_reference.csv"
    structures = ['2TRX', '1XOA', '1XOB', '2H74','1XOB_M5','1XOB_M7','1XOB_M16','1XOA_M17']
    mi_paths = [execution_folder + pdb + '/mi_data_path/zmip_sequences-beta5.0-nsus15.0-runs20000.csv' for pdb in structures]
    sum_contact_map='../THIO_ECOLI_4_107/sum_contact_map.dat'
    contact_map = util.load_contact_map(sum_contact_map)
    msa_analysis.dendogram_top_mi(2,zmip_natural_result_path,mi_paths,output_path,'Top MI Clustering entre Conformeros',structures,'single', contact_map.shape, constants.neighbour_window)
#dendogram_top_mi()    


def generate_contact_map_with_top_mi_matrix():
    sum_contact_map='../THIO_ECOLI_4_107/sum_contact_map.dat'
    #sum_contact_map='../THIO_ECOLI_4_107/sum_contact_map.dat'
    #top_smi = '../top_smi.csv'
    top_1_mi = '../THIO_ECOLI_4_107/conjunction_mi/top_1_mi.csv'
    top_2_mi = '../THIO_ECOLI_4_107/conjunction_mi/top_2_mi.csv'
    top_3_mi = '../THIO_ECOLI_4_107/conjunction_mi/top_3_mi.csv'
    top_4_mi = '../THIO_ECOLI_4_107/conjunction_mi/top_4_mi.csv'
    top_5_mi = '../THIO_ECOLI_4_107/conjunction_mi/top_5_mi.csv'
    out_matrix = '../THIO_ECOLI_4_107/sum_contact_map_and_top_mi.dat'
    msa_analysis.generate_contact_map_with_top_mi_matrix(str(1),sum_contact_map, top_1_mi, out_matrix+'1')
    msa_analysis.generate_contact_map_with_top_mi_matrix(str(2),sum_contact_map, top_2_mi, out_matrix+'2')
    msa_analysis.generate_contact_map_with_top_mi_matrix(str(3),sum_contact_map, top_3_mi, out_matrix+'3')
    msa_analysis.generate_contact_map_with_top_mi_matrix(str(4),sum_contact_map, top_4_mi, out_matrix+'4')
    msa_analysis.generate_contact_map_with_top_mi_matrix(str(5),sum_contact_map, top_5_mi, out_matrix+'5')
    
    #msa_analysis.generate_contact_map_with_top_mi_two_separated_matrix(str(1),sum_contact_map, top_1_mi, out_matrix+'1')
    
    #plot.contact_map(cmap,'THIO_ECOLI_4_107/prob_contact_map_.png')
    
#generate_contact_map_with_top_mi_matrix()


#conjunction_analisys_family(100)   
#conjunction_analisys_family(200)   
#conjunction_analisys_family(400)   
#conjunction_analisys_family(600)
#conjunction_analisys_family(800)   


# conjunction_analisys(3000)
# conjunction_analisys(5000)
# conjunction_analisys(10000)


'''
Plot conservacion media con natural, de toda la familia 
'''    
def plot_conservation_media_with_natural_family():
    execution_folder = '../FAMILY_PF00085/PF00085/'
    conservation_media = '../FAMILY_PF00085/PF00085/information_content_media.csv'
    #TODO la conservacion natural debe ser del MSA natural superpuesto y recortado.
    conservation_natural = '../FAMILY_PF00085/PF00085/information_content_superimposed_PF00085.csv'
    df_media = pandas.read_csv(conservation_media, usecols=['Entropy'])
    df_natural = pandas.read_csv(conservation_natural, usecols=['Entropy'])
    natural = df_natural['Entropy'].tolist()
    natural.insert(0, None)
    media = df_media['Entropy'].tolist()
    media.insert(0, None)
    msas_entropy = [[natural, 'Natural', 'blue'], [media, 'Media', 'red']]
    plot.conservation_comparation(msas_entropy, execution_folder + 'conservation_media.png', 'Conservacion Media y Natural')
    
#plot_conservation_media_with_natural_family()    
'''
Plot conservacion media con natural
'''
def plot_conservation_media_with_natural():
    execution_folder = '../THIO_ECOLI_4_107/'
    conservation_media = execution_folder + 'information_content_media.csv'
    conservation_natural = '../THIO_ECOLI_4_107_2TRX_A/information_content_PF00085_THIO_ECOLI_reference.csv'
    df_media = pandas.read_csv(conservation_media, usecols=['Entropy'])
    df_natural = pandas.read_csv(conservation_natural, usecols=['Entropy'])
    natural = df_natural['Entropy'].tolist()
    natural.insert(0, None)
    media = df_media['Entropy'].tolist()
    media.insert(0, None)
    msas_entropy = [[natural, 'Natural', 'blue'], [media, 'Media', 'red']]
    plot.conservation_comparation(msas_entropy, execution_folder + 'conservation_media.png', 'Conservacion Media y Natural')

#plot_conservation_media_with_natural()    

'''
Plot conservacion de conformeros
'''
def plot_conservation_conformeros():
    execution_folder = '../THIO_ECOLI_4_107/'
    structures = ['2TRX', '1XOA', '1XOB', '2H74','1XOB_M5','1XOB_M7','1XOB_M16','1XOA_M17']
    msas_summary = [execution_folder + pdb + '/conservation/information_content.csv' for pdb in structures]
    output = execution_folder + 'convervation_conformeros.png'
    msas_entropy = []
    colors = ['orange', 'green', 'grey', 'yellow', 'brown', 'goldenrod', 'magenta', 'crimson']
    i = 0
    for m in msas_summary:
        df_c = pandas.read_csv(m, usecols=['Entropy'])
        msas_entropy.append([df_c['Entropy'].tolist(), structures[i], colors[i]])
        i = i + 1
    plot.conservation_comparation(msas_entropy, output, 'Conservacion Conformeros')
# plot_conservation_conformeros()



def plot_comparation_top():
    execution_folder = '../THIO_ECOLI_4_107/'
    top_df = pandas.read_csv(execution_folder + 'result_mi_comparation.csv', header=0)

    df_1 = top_df.loc[top_df['contact_threshold'] == 1]
    plot.top_comparation(df_1, '>= 1', execution_folder + 'top_comparation_contact_1.png')
    df_4 = top_df.loc[top_df['contact_threshold'] == 4]
    plot.top_comparation(df_4, '>= 4', execution_folder + 'top_comparation_contact_4.png')
    df_8 = top_df.loc[top_df['contact_threshold'] == 8]
    plot.top_comparation(df_8, '= 8', execution_folder + 'top_comparation_contact_8.png')
#plot_comparation_top()


'''
Plot de curvas rocs para evoluciones especificas
Move this method to a test file or main file
'''
def plot_rocs_example_2trx():
    mi_data_path = '../THIO_ECOLI_4_107_2TRX_A/optimization_folder/clustered_sequences_path/mi_data_path/'
    contact_map_path = '../THIO_ECOLI_4_107_2TRX_A/optimization_folder/contact_map.dat'
    dca_data_path = '../THIO_ECOLI_4_107_2TRX_A/optimization_folder/clustered_sequences_path/dca/'
    labels = ['runs 20000', 'runs 10000', 'runs 5000','runs 1000']
    colors = ['red', 'orange', 'yellow', 'chocolate']
    sc = []
    targets = []
    mi_data_path_natural = '../THIO_ECOLI_4_107_2TRX_A/natural/zmip_PF00085_THIO_ECOLI_reference.csv'
    
    sufix = "zmip_sequences-beta0.5-nsus2.0-runs20000"
    file_name = sufix + ".csv"
    target, scores = msa_analysis.getTargetScores(mi_data_path + file_name, contact_map_path, mi_data_path_natural, False, constants.neighbour_window)
    sc.append(scores)
    targets.append(target)
    
    sufix = "zmip_sequences-beta0.5-nsus2.0-runs10000"
    file_name = sufix + ".csv"
    target, scores = msa_analysis.getTargetScores(mi_data_path + file_name, contact_map_path, mi_data_path_natural, False, constants.neighbour_window)
    sc.append(scores)
    targets.append(target)
    
    sufix = "zmip_sequences-beta0.5-nsus2.0-runs5000"
    file_name = sufix + ".csv"
    target, scores = msa_analysis.getTargetScores(mi_data_path + file_name, contact_map_path, mi_data_path_natural, False, constants.neighbour_window)
    sc.append(scores)
    targets.append(target)
    
    sufix = "zmip_sequences-beta0.5-nsus2.0-runs1000"
    file_name = sufix + ".csv"
    target, scores = msa_analysis.getTargetScores(mi_data_path + file_name, contact_map_path, mi_data_path_natural, False, constants.neighbour_window)
    sc.append(scores)
    targets.append(target)
    
    
    
    
    output_file = '../THIO_ECOLI_4_107_2TRX_A/beta0.5_nsus_2_MI.png'
    plot.roc_curve_(targets, sc, labels, 'ROC curves beta=0.5 nsus=2 metodo MI' , colors, 'Runs', output_file)

    #dca
    
    sc = []
    targets = []
    
    
    sufix = "dca_disequences-beta0.5-nsus2.0-runs20000"
    file_name = sufix + ".csv"
    target, scores = msa_analysis.getTargetScores(dca_data_path + file_name, contact_map_path, mi_data_path_natural, False, constants.neighbour_window, split_char=' ')
    sc.append(scores)
    targets.append(target)
    
    sufix = "dca_disequences-beta0.5-nsus2.0-runs10000"
    file_name = sufix + ".csv"
    target, scores = msa_analysis.getTargetScores(dca_data_path + file_name, contact_map_path, mi_data_path_natural, False, constants.neighbour_window, split_char=' ')
    sc.append(scores)
    targets.append(target)
    
    sufix = "dca_disequences-beta0.5-nsus2.0-runs5000"
    file_name = sufix + ".csv"
    target, scores = msa_analysis.getTargetScores(dca_data_path + file_name, contact_map_path, mi_data_path_natural, False, constants.neighbour_window, split_char=' ')
    sc.append(scores)
    targets.append(target)
    
    sufix = "dca_disequences-beta0.5-nsus2.0-runs1000"
    file_name = sufix + ".csv"
    target, scores = msa_analysis.getTargetScores(dca_data_path + file_name, contact_map_path, mi_data_path_natural, False, constants.neighbour_window, split_char=' ')
    sc.append(scores)
    targets.append(target)
    
    output_file = '../THIO_ECOLI_4_107_2TRX_A/beta0.5_nsus_2_DCA.png'
    plot.roc_curve_(targets, sc, labels, 'ROC curves beta=0.5 nsus=2 metodo DCA' , colors, 'Runs', output_file)


#plot_rocs_example_2trx()


def plot_rocs():
    mi_data_path = '../THIO_ECOLI_4_107_2TRX_A/optimization_folder/curated_sequences_path/mi_data_path/'
    contact_map_path = '../THIO_ECOLI_4_107_2TRX_A/optimization_folder/contact_map_sync.dat'

    labels = ['Natural', 'Evol Beta 5 Nsus 15', 'Evol Beta 5 Nsus 10', 'Evol Beta 10 Nsus 15']
    colors = ['blue', 'red', 'orange', 'yellow']
    sc = []
    targets = []
    mi_data_path_natural = '../THIO_ECOLI_4_107_2TRX_A/natural/zmip_PF00085_THIO_ECOLI_reference.csv'
    target, scores = msa_analysis.getTargetScores(mi_data_path_natural, contact_map_path, constants.neighbour_window)
    sc.append(scores)
    targets.append(target)
    
    sufix = "zmip_sequences-beta5.0-nsus15.0-runs20000"
    file_name = sufix + ".csv"
    target, scores = msa_analysis.getTargetScores(mi_data_path + file_name, contact_map_path, mi_data_path_natural, False, constants.neighbour_window)
    sc.append(scores)
    targets.append(target)
    
    sufix = "zmip_sequences-beta5.0-nsus10.0-runs20000"
    file_name = sufix + ".csv"
    target, scores = msa_analysis.getTargetScores(mi_data_path + file_name, contact_map_path, mi_data_path_natural, False, constants.neighbour_window)
    sc.append(scores)
    targets.append(target)
    
    sufix = "zmip_sequences-beta10.0-nsus15.0-runs20000"
    file_name = sufix + ".csv"
    target, scores = msa_analysis.getTargetScores(mi_data_path + file_name, contact_map_path, mi_data_path_natural, False, constants.neighbour_window)
    sc.append(scores)
    targets.append(target)
    
    output_file = '../THIO_ECOLI_4_107_2TRX_A/comparation_beta_nsus_with_natural'
    plot.roc_curve_(targets, sc, labels, 'ROC curve beta/nsus comparations with natural' , colors, 'Runs', output_file)


def plot_rocs_2trx():
    mi_data_path = '../THIO_ECOLI_4_107_2TRX_A/optimization_folder/curated_sequences_path/mi_data_path/'
    dca_data_path = '../THIO_ECOLI_4_107_2TRX_A/optimization_folder/curated_sequences_path/dca/'
    contact_map_path = '../THIO_ECOLI_4_107_2TRX_A/optimization_folder/contact_map_sync.dat'

    labels = ['Natural MI', 'Natural DCA', 'MI Beta 5 Nsus 15', 'DCA Beta 0.5 Nsus 5']
    colors = ['blue', 'green', 'red', 'orange']
    sc = []
    targets = []
    
    mi_data_path_natural = '../THIO_ECOLI_4_107_2TRX_A/natural/zmip_PF00085_THIO_ECOLI_reference.csv'
    target, scores = msa_analysis.getTargetScores(mi_data_path_natural, contact_map_path, constants.neighbour_window)
    sc.append(scores)
    targets.append(target)

    mi_data_path_DI = '../THIO_ECOLI_4_107_2TRX_A/natural/PF00085_THIO_ECOLI_reference.fasta_dca_DI.txt'
    target_di, scores_di = msa_analysis.getTargetScores(mi_data_path_DI, contact_map_path, constants.neighbour_window,split_char=' ')
    sc.append(scores_di)
    targets.append(target_di)
    
    
    sufix = "zmip_sequences-beta5.0-nsus15.0-runs20000"
    file_name = sufix + ".csv"
    target, scores = msa_analysis.getTargetScores(mi_data_path + file_name, contact_map_path, mi_data_path_natural, True, constants.neighbour_window)
    sc.append(scores)
    targets.append(target)
    
    
    
    sufix = "dca_disequences-beta0.5-nsus5.0-runs20000"
    file_name = sufix + ".csv"
    target, scores = msa_analysis.getTargetScores(dca_data_path + file_name, contact_map_path, mi_data_path_DI, True, constants.neighbour_window,split_char=' ')
    sc.append(scores)
    targets.append(target)
    
    
        
    output_file = '../THIO_ECOLI_4_107_2TRX_A/comparation_beta_nsus_with_natural'
    plot.roc_curve_(targets, sc, labels, 'Best ROC curves for MI and DCA' , colors, 'Runs', output_file)

#plot_rocs_2trx()

'''Plot de curvas rocs para evoluciones especificas. Igual pero plotea al promedio con los conformeros'''


def plot_rocs_3():
    execution_folder = '../THIO_ECOLI_4_107/'
    structures = ['2TRX', '1XOA', '1XOB', '2H74', '1KEB', '2H6Z', '2H6X', '2H76']
    mi_paths = [execution_folder + pdb + '/mi_data_path/zmip_sequences-beta5.0-nsus15.0-runs20000.csv' for pdb in structures]
    
    contact_map_paths = [execution_folder + pdb + '/contact_map_sync.dat' for pdb in structures]
    contact_map_path = execution_folder + 'sum_contact_map.dat'
    
    sc = []
    targets = []
    mi_data_path_natural = '../THIO_ECOLI_4_107_2TRX_A/natural/zmip_PF00085_THIO_ECOLI_reference.csv'
    labels = ['2TRX', '1XOA', '1XOB', '2H74', '1KEB', '2H6Z', '2H6X', '2H76', 'Natural', 'Evol Prom']
    colors = ['bisque', 'lightcyan', 'lavender', 'cornsilk', 'moccasin', 'thistle', 'lemonchiffon', 'lightyellow', 'blue', 'red']
    
    # conformeros
    i = 0
    for mi_p in mi_paths:
        target, scores = msa_analysis.getTargetScores(mi_p, contact_map_paths[i], mi_data_path_natural, True, constants.neighbour_window)
        sc.append(scores)
        targets.append(target)
        i = i + 1
    # natural
    target, scores = msa_analysis.getTargetScores(mi_data_path_natural, contact_map_path, constants.neighbour_window)
    sc.append(scores)
    targets.append(target)
    # promedio
    target, scores = msa_analysis.getTargetScores(execution_folder + 'prom/zmip_prom.csv', contact_map_path, mi_data_path_natural, True, constants.neighbour_window)
    sc.append(scores)
    targets.append(target)
    
    output_file = execution_folder + '/prom/roc_curve_prom_comparations.png'
    plot.roc_curve_(targets, sc, labels, 'ROC curve Prom comparacion con natural y conformeros' , colors, 'Runs', output_file)

# plot_rocs_3()


'''Plot de curvas rocs para evoluciones especificas. Igual pero plotea al promedio con los conformeros'''


def plot_rocs_2():
    execution_folder = '../THIO_ECOLI_4_107/'
    structures = ['2TRX', '1XOA', '1XOB', '2H74','1XOB_M5','1XOB_M7','1XOB_M16','1XOA_M17']
    mi_paths = [execution_folder + pdb + '/mi_data_path/zmip_sequences-beta5.0-nsus15.0-runs20000.csv' for pdb in structures]
    contact_map_path = execution_folder + 'sum_contact_map.dat'
    labels = ['2TRX', '1XOA', '1XOB', '2H74','1XOB_M5','1XOB_M7','1XOB_M16','1XOA_M17', 'Natural', 'Evol Prom']
    colors = ['bisque', 'lightcyan', 'lavender', 'cornsilk', 'moccasin', 'thistle', 'lemonchiffon', 'lightyellow', 'blue', 'red']
    sc = []
    targets = []
    mi_data_path_natural = '../THIO_ECOLI_4_107_2TRX_A/natural/zmip_PF00085_THIO_ECOLI_reference.csv'
    
    # conformeros
    for mi_p in mi_paths:
        target, scores = msa_analysis.getTargetScores(mi_p, contact_map_path, mi_data_path_natural, True, constants.neighbour_window)
        sc.append(scores)
        targets.append(target)
    
    # natural
   
    target, scores = msa_analysis.getTargetScores(mi_data_path_natural, contact_map_path, constants.neighbour_window)
    sc.append(scores)
    targets.append(target)
    # promedio
    target, scores = msa_analysis.getTargetScores(execution_folder + 'prom/zmip_prom.csv', contact_map_path, mi_data_path_natural, True, constants.neighbour_window)
    sc.append(scores)
    targets.append(target)
    
    output_file = execution_folder + '/prom/roc_curve_prom_comparations.png'
    plot.roc_curve_(targets, sc, labels, 'ROC curve Prom comparacion con natural y conformeros' , colors, 'Runs', output_file)

#plot_rocs_2()


def plot_auc_family():
    file = '../FAMILY_PF00085/PF00085/PF00085_evol_info.csv'
    output = '../FAMILY_PF00085/PF00085/PF00085_AUC'
    df = pandas.read_csv(file, header=0, usecols=['auc', 'auc_nat', 'auc_01', 'auc_nat_01', 'status'])
    pdb_to_compare = df.loc[df['status'] == 'okey']
    plot.auc_family(pdb_to_compare['auc'], pdb_to_compare['auc_nat'], output, 'PF00085 AUC')
    output = '../FAMILY_PF00085/PF00085/PF00085_AUC_01'
    plot.auc_family(pdb_to_compare['auc_01'], pdb_to_compare['auc_nat_01'], output, 'PF00085 AUC 0.1' ,)
    
    
def plot_auc_optimization():
    optimization_results = '../THIO_ECOLI_4_107_2TRX_A/optimization_folder/optimizacion.csv'    
    df = pandas.read_csv(optimization_results, header=0, usecols=['auc_', 'beta', 'nsus', 'run'])
    df_to_plot = df[(df['run'] == 20000)]
    df_to_plot = df_to_plot[(df_to_plot['beta'] == 0.1) | (df_to_plot['beta'] == 0.25) | (df_to_plot['beta'] == 0.5) | (df_to_plot['beta'] == 1) | (df_to_plot['beta'] == 5) | (df_to_plot['beta'] == 7) | (df_to_plot['beta'] == 10) | (df_to_plot['beta'] == 15) | (df_to_plot['beta'] == 20)]
    colors = ['tan', 'chocolate' ,'yellow', 'yellowgreen', 'red', 'orange', 'coral', 'olive', 'green']
    # df_to_plot=df[(df['beta']==0.1) | (df['beta']==0.5) | (df['beta']==1) | (df['beta']==2) | (df['beta']==3) | (df['beta']==5)]
    plot.auc_optimization(df_to_plot, 'auc_' , colors, optimization_results + '.png','Optimizacion Metodo MI')
#plot_auc_optimization()
def plot_auc_optimization_dca():
    optimization_results = '../THIO_ECOLI_4_107_2TRX_A/optimization_folder/optimization_data.csv'    
    df = pandas.read_csv(optimization_results, header=0, usecols=['auc_di_', 'beta', 'nsus', 'run'])
    df_to_plot = df[(df['run'] == 20000)]
    df_to_plot = df_to_plot[(df_to_plot['beta'] == 0.1) | (df_to_plot['beta'] == 0.25) | (df_to_plot['beta'] == 0.5) | (df_to_plot['beta'] == 1) | (df_to_plot['beta'] == 5) | (df_to_plot['beta'] == 7) | (df_to_plot['beta'] == 10) | (df_to_plot['beta'] == 15) | (df_to_plot['beta'] == 20)]
    colors = ['tan', 'chocolate' ,'yellow', 'yellowgreen', 'red', 'orange', 'coral', 'olive', 'green']
    # df_to_plot=df[(df['beta']==0.1) | (df['beta']==0.5) | (df['beta']==1) | (df['beta']==2) | (df['beta']==3) | (df['beta']==5)]
    plot.auc_optimization(df_to_plot, 'auc_di_',colors, optimization_results + '_di.png','Optimizacion Metodo DCA')

