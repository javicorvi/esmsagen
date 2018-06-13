'''
Created on Jan 6, 2017

@author: javi
'''
import os
from subprocess import call
from Bio.Align.AlignInfo import PSSM
import glob
import gzip
import time
import re
import msa_analysis
#import web_logo
import pandas
import logging
import constants as cons
import random
from Bio import AlignIO
from Bio.Align import AlignInfo
from Bio.Alphabet import IUPAC
from Bio.SubsMat import FreqTable
import Bio.Alphabet as Alphabet
from Bio import motifs


'''
Calculates de buslje09 MI corrected.
For the protein family in falsta_path
The result is stored in zmip_result_path
'''     
def buslje09(fasta_path, zmip_result_path):
    start_time = time.time()
    print "buslje09"
    #call(["julia"])
    #call(["julia04/bin/julia", "mitos/buslje09.jl" ])
    #h=heatmap(full(zmip), yflip=true)
    #call(["julia04/bin/julia", "mitos/buslje09.jl",fasta_path, zmip_result_path])
    call([cons.julia_exe_path, "../lib/mitos/buslje09.jl",fasta_path, zmip_result_path])
    print "buslje09"
    print("--- %s seconds ---" % (time.time() - start_time))   
        
def buslje09_(input_folder, zmip_result_path,pattern_array=["sequences"]): 
    start_time = time.time()
    print "buslje09_"
    for filename in os.listdir(input_folder):
        if filename.endswith(".cluster") & any(r in filename for r in pattern_array):
            buslje09(input_folder + filename , zmip_result_path + "zmip" + filename + ".dat")
    print "buslje09_"
    print("--- %s seconds ---" % (time.time() - start_time))
 
'''
Calcula el DCA por pares de posiciones , pseudocount = Frobenius 
'''     
def gaussDcaFrobenius(fasta_path,output_path):
    start_time = time.time()
    logging.info('Begin of the method')
    #call(["julia", "auc_script.jl" ])
    call([cons.julia_exe_path, "../lib/gaussdca/frobenius.jl",fasta_path, output_path])
    logging.info('Time of Execution --- %s seconds ---' % (time.time() - start_time)) 
'''
Calcula el DCA por pares de posiciones, pseudocount = Direct Information
'''    
def gaussDcaDirectInformation(fasta_path,output_path):
    start_time = time.time()
    logging.info('Begin of the method')
    #call(["julia", "auc_script.jl" ])
    call([cons.julia_exe_path, "../lib/gaussdca/direct_information.jl",fasta_path, output_path])
    logging.info('Time of Execution --- %s seconds ---' % (time.time() - start_time))  

'''
Calcula el informacion referida a la coevolucion con sparse inverse covariance estimation 
'''    
def psicov(fasta_path,output_path):
    start_time = time.time()
    logging.info('Begin of the method')
    f = open(output_path, "w")
    call(["../lib/psicov/psicov", "-p" ,"-d","0.03","-j","1",fasta_path],stdout=f)
    logging.info('Time of Execution --- %s seconds ---' % (time.time() - start_time))  


 
  
'''
Global Variables
'''
exten=".fasta"
letters = {'A':0,'R':0,'N':0,'D':0,'B':0,'C':0,'E':0,'Q':0,'Z':0,'G':0,'H':0,'I':0,'L':0,'K':0,'M':0,'F':0,'P':0,'S':0,'T':0,'W':0,'Y':0,'V':0}


def count_sequences(msa_path):
    count = 0
    with open(msa_path,'r') as msa:
        for line in msa:
            if('>' not in line):
                count=count+1
    msa.close()
    return count
def count_aminoacids(msa_path):
    #TO OJO COUNT 
    count = 0
    with open(msa_path,'r') as msa:
        for line in msa:
            if('>' not in line):
                count = len(line.strip())
                break    
    msa.close()
    return count
    
'''
Set the protein as reference in the msa.
This process called MIToS script to acomplish this function 
'''
def setProteinReference(fasta_path):
    start_time = time.time()
    logging.info('Begin of the method')
    #call(["julia", "auc_script.jl" ])
    call([cons.julia_exe_path, "../lib/mitos/setProteinReferenceToMSA.jl",fasta_path])
    logging.info('Time of Execution --- %s seconds ---' % (time.time() - start_time)) 
    
def clustering(clust, input_folder, output_folder, pattern_array=[".fasta"]):
    start_time = time.time()
    print "clustering_begins"
    for filename in os.listdir(input_folder):
        if filename.endswith(exten) & any(r in filename for r in pattern_array):
            print(filename)
            filenameclust = filename + "_"+clust + ".cluster"
            print(filenameclust)
            try:
                call(["cdhit", "-i" , input_folder+"/"+filename ,"-o", output_folder+"/"+filenameclust,"-c",clust,"-n", "4", "-M", "6000"])
            except Exception:
                print "The clusterization  get an exception with de pdb file " + input_folder+"/"+filename
                raise Exception ("The clusterization  get an exception with de pdb file " + input_folder+"/"+filename)    

    print "clustering_ends"
    print("--- %s seconds ---" % (time.time() - start_time))

def clustering_singular(clust, input_file, output_file):
    start_time = time.time()
    print(input_file)
    try:
        call(["cdhit", "-i" , input_file ,"-o", output_file,"-c",clust,"-n", "4", "-M", "6000"])
    except Exception:
        print "The clusterization  get an exception with de pdb file " + input_file
        raise Exception ("The clusterization  get an exception with de pdb file " + input_file)    
    print "clustering_ends"
    print("--- %s seconds ---" % (time.time() - start_time))

'''
Crea los logo shannon y kl 
'''
def msa_information_process(input_folder, output_folder, pattern_array=[".cluster"]):
    start_time = time.time()
    print "msa_information_process"
    for filename in os.listdir(input_folder):
        if filename.endswith(".cluster") & any(r in filename for r in pattern_array):
            print(filename)
            try:
                web_logo.create_web_logo(input_folder + filename, output_folder + filename + "_logo_sh.png",output_folder + filename + "_data_sh.csv", 'png', filename, logo_type='SHANNON')
                #web_logo.create_web_logo(input_folder + filename, output_folder + filename + "_logo_kl.png",output_folder + filename + "_data_kl.csv", 'png', filename, logo_type='KL')
                #web_logo.create_web_logo(input_folder + filename, output_folder + filename + "_logo_rob.png",output_folder + filename + "_data_rob.csv", 'png', filename, logo_type='ROB')
            except Exception as inst:
                print inst
                print "LOGO  get an exception with the msa  " + input_folder+"/"+filename
                raise Exception ("LOGO  get an exception with the msa  " + input_folder+"/"+filename)    

    print "msa_information_process"
    print("--- %s seconds ---" % (time.time() - start_time))

def msa_information(input_msa, output_msa, msa_name):
    start_time = time.time()
    print "msa_information"
    print(msa_name)
    try:
        web_logo.create_web_logo(input_msa, output_msa + "_logo_sh.png",output_msa + "_data_sh.csv", 'png', msa_name, logo_type='SHANNON')
        #web_logo.create_web_logo(input_msa, output_msa + "_logo_kl.png",output_msa + "_data_kl.csv", 'png', msa_name, logo_type='KL')
    except Exception as inst:
        print inst
        print "LOGO  get an exception with the msa  " + input_msa
        raise Exception ("LOGO  get an exception with the msa  " + input_msa)    

    print "msa_information"
    print("--- %s seconds ---" % (time.time() - start_time))

def seq_to_logo(msa,title):
    print "seq to logo"
    import subprocess 
    #subprocess.call("../lib/seq2logo-2.0/Seq2Logo.py -f " + msa + " -o "+msa+'_logo'+" -I 2 -u Bits -t '"+title+"' --format 'JPEG,PNG,SVG,PDF'", shell=True)
    print "end"


def seq_to_logo_msas(msas,structures):
    print "seq to logo"
    import subprocess 
    i=0
    for msa in msas:
        #subprocess.call("seq2logo-2.0/Seq2Logo.py -f seq2logo-2.0/test_data/fasta.txt -o logo_example -I 2 -u Bits -t 'Titulo_Logo' --formats 'JPEG,PNG,SVG,PDF'", shell=True)
        seq_to_log(msa, structures[i])
        i=i+1
    #call(["seq2logo-2.0/Seq2Logo.py", "-f" ,"seq2logo-2.0/test_data/fasta.txt"])
    #os.system("script2.py 1")
    print "end"
'''
Read the conservation info stored in file_conservation_path  
'''
def read_conservation(file_conservation_path):
    start_time = time.time()
    print "read_conservation"
    print file_conservation_path
    fields=["#","Entropy"]
    df = pandas.read_csv(file_conservation_path, delim_whitespace=True,header=7,usecols=fields)
    #df.shape[0]
    print "msa_information"
    print("--- %s seconds ---" % (time.time() - start_time))
    return df
    
def conservation(msa_path):
    import numpy as np
    import scipy.stats as sc
    from Bio import AlignIO
    from Bio.Align import AlignInfo
    from Bio.Alphabet import IUPAC
    from Bio.SubsMat import FreqTable
    import Bio.Alphabet as Alphabet
    from Bio import motifs
    for filename in os.listdir(msa_path):
        if filename.endswith(".cluster"):
            alignment = AlignIO.read(msa_path+filename, "fasta", alphabet=Alphabet.ProteinAlphabet())
            columns_quantity = []
            columns_frequency = []
            #summary_align = AlignInfo.SummaryInfo(alignment)
            #pssm = summary_align.pos_specific_score_matrix()
            #print pssm
            for x in range(0, len(alignment[0].seq)-1):
                column = alignment[:, x]
                quantity=letters
                for f in column:
                    print(f)
                    quantity[f]+=1
                double = 20/len(alignment)
                print len(alignment)
                print (quantity)
                #frequency=list(map(lambda x: x/len(alignment), quantity)) 
                frequency = dict(map(lambda (k,v): (k, v/len(alignment)), quantity.iteritems()))
                print frequency
                columns_quantity.append(quantity)
                columns_frequency.append(frequency)
            print (columns_quantity)
               
            
            '''
            m = motifs.create(alignment,alphabet=Alphabet.ProteinAlphabet())
            print (m)
            
            alfa = summary_align.alignment._alphabet
            base_alpha = Alphabet._get_base_alphabet(alfa) 
            print(summary_align)
            print(alfa)
            print(base_alpha)
            data=summary_align.information_content(5,30)
            print(data)'''
            
    #n is the number of data points
    ''''n=10
    kld = np.zeros(n, n)
    for i in range(0, n):
        for j in range(0, n):
            if(i != j):
                kld[i, j] = sc.entropy(distributions[i, :], distributions[j, :])'''

#Conver from Stockholm MSA to Fasta MSA            
def convertMSAToFasta(msa_, new_msa):
    with open(new_msa,'w') as new_file:
        with open(msa_) as old_file:
            for line in old_file:
                if('#' not in line):
                    new_line = re.split(r'\t+', line.rstrip('\t'))
                    if(len(new_line)==2):
                        new_file.write(">"+new_line[0]+"\n")
                        new_file.write(new_line[1])
    old_file.close()
    new_file.close()

# Unzip file and the convert file to fasta format and save it.  Return the path of de msa fasta format 
def natural_msa_mi(msa_file_name_fasta, result_zmip_path):
    try:
        dataanalisys.buslje09(msa_file_name_fasta, result_zmip_path)
    except BaseException as inst:
        logging.error('Error execution MI form the natural MSA ' )
        raise Exception('Error execution MI form the natural MSA')
    
def lettercount(pos):
    return {c: pos.count(c) for c in pos}




def frequency():
    sequences = ['AATC','GCCT','ATCA']
    f = zip(*sequences)
    counts = [{letter: column.count(letter) for letter in column} for column in f]
    print(counts)
    import csv
    with open('test', "wb") as f:
        writer = csv.writer(f)
        for row in counts:
            for l in ['A','T','G','C']:
                if(row.has_key(l)):
                    print l  + str(row[l])
                else:
                    print l + "0"    
        f.close()   

def random_seq(input, ouput,count):
    lines = open(input).read().splitlines()  
    with open(ouput, "w") as f:
        for i in xrange(0,count):
            f.write(">SEQ_"+str(i)+"\n")
            line = random.choice(lines)
            f.write(line+"\n") 
        f.close()    

def create_msa_bootstrap(msas_path, out_put_msa, random_number):
    i=0
    with open(out_put_msa, "w") as output:
        for msa in msas_path:
            with open(msa,"r") as f:
                lines = random.sample(f.readlines(),random_number)
                for line in lines:
                    output.write(">SEQ_"+str(i)+"\n")
                    output.write(line)
                    i=i+1
            f.close()      
    output.close()        
            
'''Elimina los ids de las secuencias'''
def create_msa_without_id(msa_path, msa_output_path):
    with open(msa_output_path, "w") as msa_output:
        with open(msa_path,"r") as msa:
            for line in msa:
                if('>' not in line):
                    msa_output.write(line)
        msa.close()
    msa_output.close()                
    
'''Informacion relacionada con el MSA '''    
def summary(msa,output,title):
    '''
       Ala (A) 9.10   Gln (Q) 3.79   Leu (L) 9.87   Ser (S) 6.69
       Arg (R) 5.71   Glu (E) 6.16   Lys (K) 4.99   Thr (T) 5.57
       Asn (N) 3.88   Gly (G) 7.26   Met (M) 2.38   Trp (W) 1.29
       Asp (D) 5.45   His (H) 2.19   Phe (F) 3.92   Tyr (Y) 2.93
       Cys (C) 1.21   Ile (I) 5.70   Pro (P) 4.85   Val (V) 6.88
    '''
    #unit_prot freq table of aminoacids 23/11/2017   
    e_freq_dict={'A':0.091, 'R':0.0571, 'N':0.0388, 'D':0.0545,'C':0.0121,'Q':0.0379,'E':0.0616,'G':0.0726,'H':0.0219,'I':0.0570,'L':0.0987,'K':0.0499,'M':0.0238,'F':0.0392,'P':0.0485,'S':0.0669,'T':0.0557,'W':0.0129,'Y':0.0293,'V':0.0688}
    #e_freq_dict={'A': 0.175, 'B': 0.325, 'C': 0.5}
    e_freq_table = FreqTable.FreqTable(e_freq_dict, FreqTable.FREQ, alphabet=Alphabet.ProteinAlphabet())
    #e_freq_table=None
    df = pandas.DataFrame()
    alignment = AlignIO.read(msa, "fasta", alphabet=Alphabet.ProteinAlphabet())
    summary_align = AlignInfo.SummaryInfo(alignment)
    total_entropy, entropy_columns, freq_dict_columns = information_content(summary_align,e_freq_table=e_freq_table)
    '''Print File de resultados'''
    for i in range(len(entropy_columns.values())):
        freq_dict = freq_dict_columns[i]
        df_2 = pandas.DataFrame([freq_dict], columns=freq_dict.keys())
        df_2['Entropy']=entropy_columns[i]
        df=df.append(df_2,ignore_index=True)
        #df.set_value(i, 'Entropy' , entropy_columns[i])
    df.to_csv(output)  
        
        
'''Informacion Entropia y Frecuencia de Amoniocidos por columna
Metodo sobreesctrito y adaptado de BioPython'''        
def information_content(self, start=0,
                            end=None,
                            e_freq_table=None, log_base=2,
                            chars_to_ignore=None):
    
        """Calculate the information content for each residue along an alignment.

        Arguments:
            - start, end - The starting an ending points to calculate the
              information content. These points should be relative to the first
              sequence in the alignment, starting at zero (ie. even if the 'real'
              first position in the seq is 203 in the initial sequence, for
              the info content, we need to use zero). This defaults to the entire
              length of the first sequence.
            - e_freq_table - A FreqTable object specifying the expected frequencies
              for each letter in the alphabet we are using (e.g. {'G' : 0.4,
              'C' : 0.4, 'T' : 0.1, 'A' : 0.1}). Gap characters should not be
              included, since these should not have expected frequencies.
            - log_base - The base of the logathrim to use in calculating the
              information content. This defaults to 2 so the info is in bits.
            - chars_to_ignore - A listing of characters which should be ignored
              in calculating the info content. Defaults to none.

        Returns:
            - A number representing the info content for the specified region.

        Please see the Biopython manual for more information on how information
        content is calculated.
        """
        # if no end was specified, then we default to the end of the sequence
        Protein20Random = 0.05
        Nucleotide4Random = 0.25
        if end is None:
            end = len(self.alignment[0].seq)
        if chars_to_ignore is None:
            chars_to_ignore = []

        if start < 0 or end > len(self.alignment[0].seq):
            raise ValueError("Start (%s) and end (%s) are not in the \
                    range %s to %s"
                    % (start, end, 0, len(self.alignment[0].seq)))
        # determine random expected frequencies, if necessary
        random_expected = None
        if not e_freq_table:
            # TODO - What about ambiguous alphabets?
            base_alpha = Alphabet._get_base_alphabet(self.alignment._alphabet)
            if isinstance(base_alpha, Alphabet.ProteinAlphabet):
                random_expected = Protein20Random
            elif isinstance(base_alpha, Alphabet.NucleotideAlphabet):
                random_expected = Nucleotide4Random
            else:
                errstr = "Error in alphabet: not Nucleotide or Protein, "
                errstr += "supply expected frequencies"
                raise ValueError(errstr)
            del base_alpha
        elif not isinstance(e_freq_table, FreqTable.FreqTable):
            raise ValueError("e_freq_table should be a FreqTable object")

        # determine all of the letters we have to deal with
        all_letters = self._get_all_letters()
        for char in chars_to_ignore:
            all_letters = all_letters.replace(char, '')

        info_content = {}
        freq_dict_ = {}
        for residue_num in range(start, end):
            freq_dict = self._get_letter_freqs(residue_num,
                                               self.alignment,
                                               all_letters, chars_to_ignore)
            # print freq_dict,
            column_score = self._get_column_info_content(freq_dict,
                                                         e_freq_table,
                                                         log_base,
                                                         random_expected)

            info_content[residue_num] = column_score
            freq_dict_[residue_num] = freq_dict
        # sum up the score
        total_info = sum(info_content.values())
        # fill in the ic_vector member: holds IC for each column
        #for i in info_content:
        #    self.ic_vector[i] = info_content[i]
        return total_info, info_content, freq_dict_ 

'''
Calcula la conservacion media
'''    
def conservation_media(msas_summary, output):    
    entropy_media=[]
    i=0
    for summary in msas_summary:
        df = pandas.read_csv(summary, usecols=['Entropy'])
        print (df)
        if(i==0):
            entropy_media = df['Entropy'].tolist()
        else:
            entropy_media = [x + y for x, y in zip(entropy_media , df['Entropy'].tolist() )]
        i=i+1
    entropy_media =   [x / i  for x in entropy_media]
    df = pandas.DataFrame(entropy_media,columns=['Entropy'])
    df.to_csv(output)

      