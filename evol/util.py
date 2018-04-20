'''
Created on Jan 12, 2017

@author: javi
'''
import numpy as np 
import os
import time
import re
import csv
import glob
import zipfile
import pandas
import logging
from sklearn import metrics

def getPartialAUC(fpr, tpr, max_fpr):
    idx = np.where(fpr <= max_fpr)[0]
    # linearly interpolate the ROC curve until max_fpr
    idx_last = idx.max()
    idx_next = idx_last + 1
    xc = [fpr[idx_last], fpr[idx_next]]
    yc = [tpr[idx_last], fpr[idx_next]]
    tpr = np.r_[tpr[idx], np.interp(max_fpr, xc, yc)]
    fpr = np.r_[fpr[idx], max_fpr]
    partial_roc = metrics.auc(fpr, tpr, reorder=True)
    # standardize result to lie between 0.5 and 1
    min_area = max_fpr**2/2
    max_area = max_fpr
    return 0.5*(1+(partial_roc-min_area)/(max_area-min_area))    

def getAUC(target,scores):
    if(target==[]):
        return 0.5,0.5
    fpr, tpr, _ = metrics.roc_curve(target, scores)
    auc = metrics.auc(fpr, tpr)
    auc01 = getPartialAUC(fpr, tpr, 0.1)
    return auc,auc01
def load_zmip(zmip_result_path,window_mi_neightboards=0,split_char=r'\t+'):
    column = []
    file = open(zmip_result_path)
    for line in file:
        line = line.replace('NaN','0')
        line = line.replace('\n','')
        #column.append(map(float,re.split(r'\t+', line)))
        line_float = map(float,re.split(split_char, line))
        if (line_float[0] + window_mi_neightboards >= line_float[1]):
            line_float[2]=0
        column.append(line_float)
    file.close()
    return column
def order(zmip_natural):    
    zmip_natural.sort(key=lambda x:x[2], reverse=True)
    

def save_list_to_csv(filename, zmip, header):
    zmip.insert(0, header)
    with open(filename, "wb") as f:
        writer = csv.writer(f)
        writer.writerows(zmip)    
    zmip.pop(0)
'''
Syncronize matrix_ref with matrix_evol.
'''    
def sincronice_mi(matrix_ref, matrix_evol):
    column = []
    column2 = [] 
    for j in matrix_ref:
        pos  = j[0]     
        pos2 = j[1]   
        #value = j[2]
        for e in matrix_evol:
            pos_  = e[0]     
            pos2_ = e[1]
            if (pos_==pos and pos2_==pos2):   
                column2.append(e)
                column.append(j)
                break    
    return  column, column2
def sincronize_natural_evol_msas(input_folder,output_folder,pattern_array,reg_init,reg_end_back):
    start_time = time.time()
    print "sincronize_natural_evol_alignments"
    count=0
    for filename in os.listdir(input_folder):
        if filename.endswith(".cluster") & any(r in filename for r in pattern_array):
            with open(output_folder+"/"+filename,'w') as new_file:
                with open(input_folder+"/"+filename) as old_file:
                    for line in old_file:
                        if('>' in line):
                            line = line.replace('\n','_'+str(count)+'\n')
                            new_file.write(line)
                            count=count+1
                        else:
                            new_file.write(line[reg_init:reg_end_back]+'\n')    
            old_file.close()
            new_file.close()
    print "sincronize_natural_evol_alignments"
    print("--- %s seconds ---" % (time.time() - start_time))

'''
Sincroniza y corta los msa evolucionados teniendo en cuenta solo la escructura en comun que existe entre los pdb de la familia
Recibe el pdf recortado con las posiciones que deben mantenerse luego elimina las posiciones del msa evolucionado.
'''
def synchronize_evol_with_cutted_pdb(pdb_complete_path, pdb_cutted_path, clustered_sequences_path, sincronized_evol_path, contact_map_path, sincronized_contact_map):
    start_time = time.time()
    print "sincronize_natural_evol_alignments"  
    df = pandas.read_csv(pdb_complete_path, delim_whitespace=True,header=None)
    df=df.loc[df[4] == 'A']
    df=df.dropna()
    df[5] = df[5].astype('int32')
    df=df.groupby(5).first().reset_index()
    start = df[5].min()
    end = df[5].max()
    df_columnas = pandas.read_csv(pdb_cutted_path, delim_whitespace=True,header=None,usecols=[5])
    df_columnas=df_columnas.dropna()
    df_columnas[5] = df_columnas[5].astype('int32')
    df_columnas=df_columnas.groupby(5).first().reset_index()
    df_columnas[5] = df_columnas[5].apply(lambda x: x - start)
    count=0
    for filename in os.listdir(clustered_sequences_path):
        if filename.endswith(".cluster"):
            with open(sincronized_evol_path+"/"+filename,'w') as new_file:
                with open(clustered_sequences_path+"/"+filename) as old_file:
                    for line in old_file:
                        if('>' in line):
                            line = line.replace('\n','_'+str(count)+'\n')
                            new_file.write(line)
                            count=count+1
                        else:
                            line_array=np.array(list(line))
                            new_line = line_array[df_columnas[5]]
                            new_file.write(new_line.tostring()+'\n')
            old_file.close()
            new_file.close()
    
    
    #Y = np.arange(36).reshape(6,6)
    #test = Y[np.ix_([0,3,5],[0,3,5])]
    cmap = load_contact_map(contact_map_path)
    cmap_sync=cmap[np.ix_(df_columnas[5].tolist(),df_columnas[5].tolist())]
    save_contact_map(cmap_sync, sincronized_contact_map)        
    
    print "sincronize_natural_evol_alignments"
    print("--- %s seconds ---" % (time.time() - start_time))

def find_first_residue(pdb_complete_path):
    with open(pdb_complete_path,'r') as pdb:
        first_line = pdb.readline()
        return int(first_line[23:26])
    

'''
Sincroniza y corta los msa evolucionados teniendo en cuenta solo la escructura en comun que existe entre los pdb de la familia
Recibe el pdf recortado con las posiciones que deben mantenerse luego elimina las posiciones del msa evolucionado.
'''
def synchronize_evol_with_cutted_pdb_singular(pdb_complete_path, pdb_cutted_path, clustered_sequences_path, sincronized_evol_path, contact_map_path, sincronized_contact_map):
    start_time = time.time()
    print "sincronize_natural_evol_alignments"  
    start = find_first_residue(pdb_complete_path)
    df_columnas = pandas.read_csv(pdb_cutted_path, delim_whitespace=True,header=None,usecols=[5])
    df_columnas=df_columnas.dropna()
    df_columnas[5] = df_columnas[5].astype('int32')
    df_columnas=df_columnas.groupby(5).first().reset_index()
    df_columnas[5] = df_columnas[5].apply(lambda x: x - start)
    count=0
    with open(sincronized_evol_path,'w') as new_file:
        with open(clustered_sequences_path) as old_file:
            for line in old_file:
                if('>' in line):
                    line = line.replace('\n','_'+str(count)+'\n')
                    new_file.write(line)
                    count=count+1
                else:
                    line_array=np.array(list(line))
                    new_line = line_array[df_columnas[5].tolist()]
                    new_file.write(new_line.tostring()+'\n')
    old_file.close()
    new_file.close()
    #Y = np.arange(36).reshape(6,6)
    #test = Y[np.ix_([0,3,5],[0,3,5])]
    cmap = load_contact_map(contact_map_path)
    cmap_sync=cmap[np.ix_(df_columnas[5].tolist(),df_columnas[5].tolist())]
    save_contact_map(cmap_sync, sincronized_contact_map)        
    
    print "sincronize_natural_evol_alignments"
    print("--- %s seconds ---" % (time.time() - start_time))

import random
import string 
def random_char(y):
    return ''.join(random.choice(string.ascii_letters) for x in range(y))

def load_contact_map(contact_map_path, dtype='i4'):
    cmap = np.loadtxt(contact_map_path, dtype=dtype)
    np.set_printoptions(threshold='nan')
    return cmap
def load_contact_map_deprecated(contact_map_path):
    with open(contact_map_path) as file:
        l = [map(str,line.split(' ')) for line in file ]
        file.close()
        
        a = np.loadtxt(contact_map_path, dtype='i4')
        np.set_printoptions(threshold='nan')
        #print (cmap)
        #test=[[1,2,3,4,5,6,7,8],[1,2,3,4,5,6,7,8],[1,2,3,4,5,6,7,8],[1,2,3,4,5,6,7,8],[1,2,3,4,5,6,7,8],[1,2,3,4,5,6,7,8],[1,2,3,4,5,6,7,8],[1,2,3,4,5,6,7,8]]
        a=np.array(l, order='F')
        a[a == 'false\n']=0
        a[a == 'false']=0
        a[a == '0\n']=0
        a[a == '0']=0
        a[a == 'true\n']=1
        a[a == 'true']=1
        a[a == '1']=1
        a[a == '1\n']=1
        a[a == '2']=1
        a[a == '2\n']=1
        cmap=np.array(a, dtype='i4')
        #np.set_printoptions(threshold='nan')
        #print (cmap)
        return cmap    
def load_contact_map_(contact_map_path):
    with open(contact_map_path) as file:
        l = [map(str,line.split(' ')) for line in file ]
        file.close()
        #test=[[1,2,3,4,5,6,7,8],[1,2,3,4,5,6,7,8],[1,2,3,4,5,6,7,8],[1,2,3,4,5,6,7,8],[1,2,3,4,5,6,7,8],[1,2,3,4,5,6,7,8],[1,2,3,4,5,6,7,8],[1,2,3,4,5,6,7,8]]
        a=np.array(l, order='F')
        a[a !=0 ] = 1
        cmap=np.array(a, dtype='i4')
        return cmap
        
'''
Sincronize contact map adjusting with the information of reg_init and reg_end_back
contact_map_path: input contact map
contact_map_output: out contact map sincronized
'''        
def sincronize_contact_map(contact_map_path, contact_map_output, reg_init, reg_end_back):
    start_time = time.time()
    print "sincronize_contact_map"
    #with open(contact_map_output,'w') as new_file:
    with open(contact_map_path) as old_file:
        l = [map(str,line.split(' ')) for line in old_file ]
        old_file.close()
        #test=[[1,2,3,4,5,6,7,8],[1,2,3,4,5,6,7,8],[1,2,3,4,5,6,7,8],[1,2,3,4,5,6,7,8],[1,2,3,4,5,6,7,8],[1,2,3,4,5,6,7,8],[1,2,3,4,5,6,7,8],[1,2,3,4,5,6,7,8]]
        a=np.array(l, order='F')
        b=a[reg_init:reg_end_back, reg_init:reg_end_back]
        contacts = np.count_nonzero(b == 'true')
        print "Contacts " + str(contacts)
        np.savetxt(contact_map_output, b, delimiter=' ',fmt="%s")
        numrows = len(b)
        numcols = len(b[0])
        print "ROWS " + str(numrows)
        print "COLUMNS " +str(numcols)
    print "sincronize_contact_map"
    print("--- %s seconds ---" % (time.time() - start_time))
    

    
def add_matrix(m1,m2):
    (m, n)=m1.shape
    m3=np.zeros((m, n), dtype='i4')
    for i in range(m):
        for j in range(n):
            m3[i][j] = m1[i][j]+ m2[i][j]
    return m3       

def save_contact_map(m, path):
    np.savetxt(path, m, delimiter=' ',fmt="%s")   
    
def delete_files(files_pattern):
    files = glob.glob(files_pattern)
    for f in files:
        os.remove(f)
        
def zip_files(files_pattern):
    files = glob.glob(files_pattern)
    for f in files:
        zf = zipfile.ZipFile(f+'.zip', mode='w')
        zf.write(f)       
        
"""
Retorna los PDB/Proteinas a evolucionar.  
Se queda con un pdb por cluster para evitar redundancia.
Con el primer PDB encontrado del cluster
"""        
def find_pdb_to_evolve(family_pdb_information):
    #fields = ["pdb"]
    fields = ['seq',"pdb","chain","cluster","n_residues"]
    df = pandas.read_csv(family_pdb_information,header=0,usecols=fields)
    logging.info("Cantidad Total Proteinas/PDB: " + str(len(df.index)))
    df=df.sort(["cluster","pdb"])
    df=df.groupby("cluster").first()
    df['pdb_folder_name']=df['seq'].str.replace("/","_").str.replace("-","_") + "_" + df['pdb']+"_"+df['chain']
    df['status']= 'pending'
    df['contacts_count']= np.NaN
    df['beta']= np.NaN
    df['nsus']= np.NaN
    df['runs']= np.NaN
    df['auc']= np.NaN
    df['auc_01']= np.NaN
    df['auc_nat']= np.NaN
    df['auc_nat_01']= np.NaN
    df['spearman_zmip_evol_nat']= np.NaN
    df['par_positions_count']= np.NaN
    df['execution_time']= np.NaN
    #print df
    logging.info("Cantidad de Proteinas/PDB a evolucionar (Uno por cluster):" + str(len(df.index)))
    return df



def find_pdb_start_end_for_protein(stockholm_msa_file, protein, pdb_name, chain):
    try:
        with open(stockholm_msa_file,'r') as file:
            for line in file:
                if(protein in line and "PDB; "+ pdb_name + " " + chain in line):  
                    spl = line.split(';')
                    sub = spl[2]
                    range=sub.split('-')
                    start= int(range[0])
                    end = int(range[1])
                    file.close()
                    return start,end
               
        file.close()
        raise Exception ("Error: no existe la informacion para leer el rango de residuos de la proteina " + protein  +" pdb " + pdb_name + " cadena " + chain)    
    except Exception as inst:
        print inst
        logging.error("Error no controlado intentando leer el rango de residuos de la proteina " + protein  +" pdb " + pdb_name + " cadena " + chain)

def getSequence(fasta, seq_id):
    return "ISLFEGANFKGNTIEIQDDAPSLWVFSVGSVKVSSGTWVGYQYPGYRGYQYLLEPGDFRHWNEWGAFQPQMQSL"
def getPDBSequence(pdb_name, pdb_path, chain):
    logging.info("getPDBSequence pdb " + pdb_name + " cadena " + chain)
    from Bio.PDB.PDBParser import PDBParser
    from Bio.PDB.Polypeptide import three_to_one
    from Bio.PDB.Polypeptide import is_aa
    residue_position = []
    residue_name = list()
    try:
        parser = PDBParser(PERMISSIVE=1)
        structure = parser.get_structure(pdb_name, pdb_path)
        model = structure[0]
        chain = model[chain]
        for residue in chain:
            if is_aa(residue.get_resname(), standard=True):
                residue_name.append(three_to_one(residue.get_resname()))
                residue_position.append(residue.get_full_id()[3][1])
            #else:
                #residue_name.append("X")
                #residue_position.append(residue.get_full_id()[3][1])
                #raise Exception("Secuencia no valida, error en la posicion: " + str(residue.get_full_id()[3][1]))
     
    except Exception as inst:
        print inst
        logging.error("Error no controlado intentando leer la sequencia del pdb "  + pdb_name + " cadena " + chain + " path " + pdb_path)            
        raise Exception ("PDB Invalido pdb "  + pdb_name + " cadena " + chain + " path " + pdb_path)
    return residue_position, residue_name
    
    
    '''
    df_columnas = pandas.read_csv(pdb_path, delim_whitespace=True,header=None,usecols=[5])
    df_columnas=df_columnas.dropna()
    df_columnas[5] = df_columnas[5].astype('int32')
    df_columnas=df_columnas.groupby(5).first().reset_index()
    
    with open(pdb_path,'r') as file:
        for line in file:
            if(line.startswith("ATOM") and line[21]==chain):  
                residue_position.append(line[23:26])
                residue_name.append(line[17:19])
    file.close()
    return residue_position, residue_name
    '''
'''def mapCDtoPDB(sequence, pdb, chain):
    pdb_seq, positions = getPDBSequence()
    MusclePairAlign("protein_seq",sequence,"pdb",pdb_seq)
'''    
def MusclePairAlign(id1, seq1, id2, seq2):
    from Bio.Align.Applications import MuscleCommandline
    from Bio import SeqIO
    tmpname = "test"
    if not os.path.exists(os.path.join("tmp")):
        os.makedirs(os.path.join("tmp"))
    with open(os.path.join("tmp", tmpname+".fasta"), 'w') as o:
        o.write(">"+id1+"\n")
        o.write(seq1+"\n")
        o.write(">"+id2+"\n")
        o.write(seq2+"\n")
    muscle_line = MuscleCommandline(input=os.path.join("tmp", tmpname+".fasta"), out=os.path.join("tmp", tmpname+".fasta.out"))
    stdout, stderr = muscle_line()
    handle = open(os.path.join("tmp", tmpname+".fasta.out"), "rU")
    records = list(SeqIO.parse(handle, "fasta"))
    handle.close()
    # print str(records[0].id)
    # print str(records[0].seq)
    # print str(records[1].id)
    # print str(records[1].seq)
    return(records[0].seq,records[1].seq)


'''    
import math, string, sys, fileinput

def range_bytes (): return range(256)
def range_printable(): return (ord(c) for c in string.printable)
def H(data, iterator=range_bytes):
    if not data:
        return 0
    entropy = 0
    for x in iterator():
        p_x = float(data.count(chr(x)))/len(data)
        if p_x > 0:
            entropy += - p_x*math.log(p_x, 2)
    return entropy

def main ():
    for row in fileinput.input():
        string = row.rstrip('\n')
        print ("%s: %f" % (string, H(string, range_printable)))

for str in ['gargleblaster', 'tripleee', 'magnus', 'lkjasdlk',
               'aaaaaaaa', 'sadfasdfasdf', '7&wS/p(']:
    print ("%s: %f" % (str, H(str, range_printable)))
'''    