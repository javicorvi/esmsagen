'''

Invoques to the scpe program 

'''

from subprocess import call


'''
Invoques to the scpe several times for the beta, nsus and runs parameters
'''
def runs(pdb_file,beta,runs,nsus,chain,result_file,contact_map_path):
    for b in beta:
        for sus in nsus:
            for r in runs:
                try: 
                    call(["../lib/scpe/scpe.exe", "-i" ,pdb_file,"-q","../lib/scpe/QijCodonesEmp" ,"-c", chain,"-b",b,"-A", sus, "-a", sus, "-R", r,"-1","1","-2", result_file,"-3","25.0","-4","0","-5",contact_map_path])
                except Exception:
                    print "The SCPE execution get an exception with de pdb file " + pdb_file
                    raise Exception ("The SCPE execution get an exception with de pdb file " + pdb_file)   
                
                

'''
Run the scpe for the pdb and the parameters
'''
def run(pdb_file,beta,run,nsus,chain,output_msa_path,contact_map_path):
    try:
        call(["../lib/scpe/scpe.exe", "-i" ,pdb_file,"-q","../lib/scpe/QijCodonesEmp" ,"-c", chain,"-b",beta,"-A", nsus, "-a", nsus, "-R", run,"-1","1","-2", output_msa_path,"-3","25.0","-4","0","-5",contact_map_path])
    except Exception as inst:
        print inst
        print "The SCPE execution get an exception with de pdb file " + pdb_file
        raise Exception ("The SCPE execution get an exception with de pdb file " + pdb_file)                   



