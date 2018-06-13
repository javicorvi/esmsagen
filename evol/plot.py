'''
Created on Jan 12, 2017

@author: javi
'''
import itertools
import matplotlib.pyplot as plt
import numpy as np
from sklearn import metrics
import matplotlib.cm as cm
from scipy.cluster.hierarchy import dendrogram
from Canvas import Line
'''
Generates a plot describe the contacts with de MI of the natural and evolution msa
'''
def contacts_with_mi(x_nat_t,y_evol_t,x_nat_f,y_evol_f,output_path,filename,pdb_name):
    #trues
    #x_nat_t = [27.0,89.0,25.0]
    #y_evol_t = [5.9,3.2,4.0]
    #falses
    #x_nat_f = [3.0,1.0,3.0]
    #y_evol_f = [1.1,0.32,1.0]
    #plt.scatter(x_nat_t, y_evol_t,color="b")
    #plt.scatter(x_nat_f, y_evol_f,color="r")
    plt.axis([0, 1, 0, 1])
    #plt.axis([0, 0.6, 0, 0.6])
    no = plt.scatter(y_evol_f, x_nat_f,color="r")
    co = plt.scatter(y_evol_t, x_nat_t,color="b")
    plt.legend((no, co),
           ('No Contact', 'Contact'),
           scatterpoints=1,
           loc='lower right',
           ncol=3,
           fontsize=8)
    #plt.set_title(filename)
    plt.ylabel('Natural Covariation Values')
    plt.xlabel('Evolution Covariation Values')
    plt.title(pdb_name)
    
    
    ''''a = plt.axes([.6, .6, .25, .25])
    plt.axis([0, 1, 0, 1])
    no = plt.scatter(y_evol_f, x_nat_f,color="r")
    co = plt.scatter(y_evol_t, x_nat_t,color="b")
    plt.axvline(x=.6, ymin=0.0, ymax=0.6,linestyle='--',color='black')
    #plt.axhline(y=.4, xmin=0.25, xmax=0.402, linewidth=2, color = 'k')
    plt.axhline(y=.6, xmin=0.0, xmax=0.6,linestyle='--',color='black')
    plt.title('Zoom Out')
    plt.xticks([0,0.6,1])
    plt.yticks([0,0.6,1])
    '''
    
    #x=[[1.0,2.4]]
    #plt.scatter(x,color="r")
    plt.savefig(output_path)
    #plt.show()
    plt.gcf().clear()

'''
Deprecated
'''
def roc_curve_dep(y_true,scores_nat,scores_evol):
    fpr, tpr, _ = metrics.roc_curve(y_true, scores_nat)
    roc_auc = metrics.auc(fpr, tpr)
    
    fpr2, tpr2, _ = metrics.roc_curve(y_true, scores_evol)
    roc_auc2 = metrics.auc(fpr2, tpr2)      
            
    plt.figure()
    lw = 2
    
    
    plt.plot(fpr, tpr, color='blue', lw=lw, label='ROC curve Nat (area = %0.2f)' % roc_auc)
    plt.plot(fpr2, tpr2, color='red', lw=lw, label='ROC curve Evol (area = %0.2f)' % roc_auc2)
    
    plt.plot([0, 1], [0, 1], color='black', lw=lw, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Receiver operating characteristic example')
    plt.legend(loc="lower right")
    #plt.show()
    plt.gcf().clear()
    
def roc_curve(df,index,y_true,scores,labels,colors,output_file,title):
    fpr = dict()
    tpr = dict()
    roc_auc = dict()
    lw = 2
    plt.figure()
    #for i, (a, b) in enumerate(zip(alist, blist)):
    for i,(label,color) in enumerate(zip(labels, colors)):
        fpr[i], tpr[i], _ = metrics.roc_curve(y_true, scores[i])
        partial_auc_value_0_1 = partial_auc(fpr[i], tpr[i], 0.1)
        auc = metrics.auc(fpr[i], tpr[i])
        roc_auc[i] = auc
        if(i==0):
            df.set_value(index,'auc_nat',auc)
            df.set_value(index,'auc_nat_01',partial_auc_value_0_1)
        else:
            df.set_value(index,'auc',auc)
            df.set_value(index,'auc_01',partial_auc_value_0_1)
        plt.plot(fpr[i], tpr[i], color=color, lw=lw,label='ROC curve  {0} (auc = {1:0.2f} | auc 0.1 = {2:0.2f})'''.format(label, roc_auc[i], partial_auc_value_0_1))
    plt.plot([0, 1], [0, 1], color='black', lw=lw, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title(title)
    plt.legend(loc="lower right")
    plt.savefig(output_file)
    #plt.show()
    plt.gcf().clear()

def roc_curve_(y_true,scores,labels,title,colors,legend_title,output_file, vertical_at_01=True,vertical_line_color='grey',xy_line=True,xy_line_color='grey'):
    fpr = dict()
    tpr = dict()
    roc_auc = dict()
    lw = 2
    plt.figure()
    #for i, (a, b) in enumerate(zip(alist, blist)):
    for i,(label,color) in enumerate(zip(labels, colors)):
        '''
        fix empty scores and targets
        if((not y_true[i]) & (not scores[i])):
            fpr[i], tpr[i] = 0.5,0.5
        else:
        '''
        fpr[i], tpr[i], _ = metrics.roc_curve(y_true[i], scores[i])
        partial_auc_value_0_1 = partial_auc(fpr[i], tpr[i], 0.1)
        auc = metrics.auc(fpr[i], tpr[i])
        roc_auc[i] = auc
        plt.plot(fpr[i], tpr[i], color=color, lw=lw,label='{0} (auc = {1:0.4f} | auc 0.1 = {2:0.4f})'''.format(label, roc_auc[i], partial_auc_value_0_1))
    
    plt.plot([0, 1], [0, 1], color='black', lw=lw, linestyle='--')
    if(vertical_at_01):
        plt.axvline(x=0.1,color=vertical_line_color)
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.0])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title(title)
    plt.legend(loc="lower right",prop={'size':10},title=legend_title)
    plt.savefig(output_file)
    #plt.show()
    plt.gcf().clear()

def partial_auc(fpr, tpr, max_fpr):
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
'''
Generate the plot for the diferents auc taking into account the beta, nsus and runs
Open the result_auc_path parse all the results and plot them into 4 subplots
one for each runs
'''
def plot_auc(result_auc_path,output_path,beta,runs,nsus):
    ys_run = np.zeros((len(runs),len(beta), len(nsus)))
    for id_r,r in enumerate(runs):  
        for id_b,b in enumerate(beta):   
            #data_beta = [] * len(nsus)
            for id_n,n in enumerate(nsus):    
                source = open(result_auc_path, 'r')
                for line in source:
                    if (('beta'+b in line) & ('nsus'+n in line) & ('runs'+r+'.' in line)):
                        str,value=line.split(",")
                        #data_beta[id_n]=value
                        ys_run[id_r,id_b,id_n]=value
                source.close()  
    print ys_run  
    
    f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2 , sharex='col', sharey='row')
    #f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2)
    colors = itertools.cycle([ "yellow","orange" ,"green", "red", "blue",  "black"])
    subplt = [ax1,ax2,ax3,ax4]
    
    axis=[0, 6, 0, 1.0]
    ax1.axis(axis)
    ax2.axis(axis)
    ax3.axis(axis)
    ax4.axis(axis)
    
    #for splt in subplt:
    for id_splt,splt in enumerate(subplt):
        i=0
        splt.set_title('Runs ' + runs[id_splt])
        
        for y, c in zip(ys_run[id_splt], colors):
            splt.scatter(nsus, y, color=c,picker=True)  
            splt.plot(nsus, y,color=c, label=beta[i])
            i=i+1    
    
    
    # Now add the legend with some customizations.
    
    f.text(0.5, 0.04, 'NSUS', ha='center')
    f.text(0.04, 0.5, 'AUC', va='center', rotation='vertical')
    #f.canvas.mpl_connect('pick_event', onpick3)
    plt.legend(loc=1,prop={'size':10},title="Betas")
    #plt.show()   
    #plt.savefig(output_path+'/auc.png')
    plt.gcf().clear()




'''

'''
def contact_map_with_top_rank_mi(contact_map, x_nat, y_nat, x_evol,y_evol,output_path,filename,title):
    cmap=contact_map
    n=plt.scatter(x_nat, y_nat,color="b",s=40,  marker=(5, 2))
    e=plt.scatter(y_evol, x_evol, color="r",s=40,  marker=(5, 2))
    plt.legend((n, e), ('Natural', 'Evol'),scatterpoints=1,loc='upper right',ncol=3,fontsize=8)
    g = np.floor(cmap)
    plt.imshow(g, cmap='Greys')
    cbar=plt.colorbar(ticks=[])
    cbar.set_label('contacts')
    plt.title(title)
    plt.savefig(output_path)
    plt.gcf().clear()

def auc_family(x,y,output,title,color="black",line_xy=True,color_line='grey',axis_label=['AUC MSA Evolucionado','AUC MSA Natural']):
    plt.scatter(x, y,color=color)
    plt.axis([0, 1, 0, 1])
    if(line_xy==True):
        plt.plot([0,1],color=color_line)
    plt.xlabel(axis_label[0])
    plt.ylabel(axis_label[1])
    plt.title(title)
    plt.savefig(output+'.pdf')
    plt.savefig(output+'.png')
    plt.savefig(output+'.eps')
    plt.gcf().clear()    
    
        
    
def contact_map_with_top_rank_mi_desarrollo(contact_map, x_nat, y_nat, x_evol,y_evol,x_evol2,y_evol2,output_path,filename):
    cmap=contact_map
    plt.scatter(x_nat, y_nat,color="b",s=40,  marker=(5, 2))
    plt.scatter(y_evol, x_evol, color="r",s=40,  marker=(5, 2))
    plt.scatter(y_evol2, x_evol2, color="g",s=40,  marker=(5, 2))
    g = np.floor(cmap)
    plt.imshow(g, cmap='Greys')
    plt.title(filename)
    plt.savefig(output_path)
    #plt.show()
    plt.gcf().clear()

'''
Plot the contact map and saved in output_file
'''
def contact_map(contact_map, output_file, title='Contact Map'):
    
    '''
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    from matplotlib.colors import LogNorm
    import numpy as np
    x, y, z = np.loadtxt('data.txt', unpack=True)
    '''
    cmap=contact_map
    
    #plt.imshow(cmap,cmap=cm.Blues)
    plt.imshow(cmap,cmap=cm.Greys,interpolation='nearest')
    #plt.imshow(cmap,cmap=cm.Set1,interpolation='nearest')
    #plt.colorbar()
    #plt.imshow(cmap,cmap=cm.hot)
    plt.title(title)
    plt.savefig(output_file)
    #plt.show()
    plt.gcf().clear()

def contact_map_sum(contact_map, output_file,title='Contact Map'):
    
    '''
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    from matplotlib.colors import LogNorm
    import numpy as np
    x, y, z = np.loadtxt('data.txt', unpack=True)
    '''
    cmap=contact_map
    
    #plt.imshow(cmap,cmap=cm.Blues)
    plt.imshow(cmap,cmap=cm.YlOrRd,interpolation='nearest')
    plt.colorbar(ticks=[0,1,2,3,4,5,6,7,8])
    #plt.imshow(cmap,cmap=cm.hot)
    plt.title(title)
    plt.savefig(output_file)
    #plt.show()
    plt.gcf().clear()

def contact_map_with_top_rank_mi_sum(contact_map, x_nat, y_nat, x_evol,y_evol,output_path,filename,title):
    cmap=contact_map
    plt.scatter(x_nat, y_nat,color="b",s=40,  marker=(5, 2))
    plt.scatter(y_evol, x_evol, color="r",s=40,  marker=(5, 2))
    g = np.floor(cmap)
    plt.imshow(cmap,cmap=cm.Greys,interpolation='nearest')
    plt.colorbar(ticks=[0,1,2,3,4,5,6,7,8])
    plt.title(title)
    plt.savefig(output_path)
    plt.gcf().clear()


'''
Generates a plot describe the contacts with de MI of the natural and evolution msa
'''
def contacts_with_mi_desarrollo(x_nat_t,y_evol_t,y_evol2_t,x_nat_f,y_evol_f,y_evol2_f,output_path,filename):
    #trues
    #x_nat_t = [27.0,89.0,25.0]
    #y_evol_t = [5.9,3.2,4.0]
    #falses
    #x_nat_f = [3.0,1.0,3.0]
    #y_evol_f = [1.1,0.32,1.0]
    #plt.scatter(x_nat_t, y_evol_t,color="b")
    #plt.scatter(x_nat_f, y_evol_f,color="r")
    plt.axis([0, 1, 0, 1])
    no = plt.scatter(y_evol_f, x_nat_f,color="r")
    no2= plt.scatter(y_evol2_f, x_nat_f,color="r")
    co = plt.scatter(y_evol_t, x_nat_t,color="b")
    co2 = plt.scatter(y_evol2_t, x_nat_t,color="g")
    plt.legend((no, co,co2),
           ('No Contact', 'Contact 2TRX', 'CONTACT 1THX'),
           scatterpoints=1,
           loc='lower left',
           ncol=3,
           fontsize=8)
    #plt.set_title(filename)
    plt.ylabel('Natural MI Values')
    plt.xlabel('Evolution MI Values')
    plt.title(filename)
    
    #x=[[1.0,2.4]]
    #plt.scatter(x,color="r")
    plt.savefig(output_path)
    #plt.show()
    plt.gcf().clear()
'''
Plot the conservation of the MSA evolutionated and the Natural MSA
'''
def conservation_between_msas(msas_entropy, output_file,natural_line_style='-'):
    for index,msa_entropy in enumerate(msas_entropy):
        if(index==0):
            plt.plot(msa_entropy[0],color=np.random.rand(3,1), label=msa_entropy[1],linewidth=5,linestyle=natural_line_style)
        else:
            plt.plot(msa_entropy[0],color=np.random.rand(3,1), label=msa_entropy[1])
    plt.ylabel('Bits Entropy')
    plt.xlabel('Position')
    #plt.legend(loc=1,prop={'size':10},title="PDB")
    plt.savefig(output_file)
    #plt.show()  
    plt.gcf().clear()


'''
Plot the conservation of the MSA evolutionated and the Natural MSA
'''
def conservation_comparation(msas_entropy, output_file, title):
    plt.figure(num=None, figsize=(10, 6), dpi=80, facecolor='w', edgecolor='k')
    for index,msa_entropy in enumerate(msas_entropy):
        plt.plot(msa_entropy[0],label=msa_entropy[1],color=msa_entropy[2],linewidth=2)
    
    #plt.plot([30,33],[6,6],color='black')
    plt.ylabel('Bits Entropy')
    plt.xlabel('Position')
    plt.legend(loc=1,prop={'size':10})
    
    #np.arange(min(x), max(x)+1, 1.0)
    
    plt.axis([0, 106, 0, 7])
    plt.xticks(np.arange(0, 106, 5.0))
    plt.title(title)
    plt.savefig(output_file)
    #plt.show()  
    plt.gcf().clear()

def top_comparation(df, contact_threshold ,output_path):
    #df_to_plot=df[(df['beta']==b)]
    plt.scatter(df['top'], df['evol_mi'],color='green',label=None)
    plt.plot(df['top'], df['evol_mi'],color='green',label='MI EVOL')
    
    plt.scatter(df['top'], df['evol_di'],color='blue',label=None)
    plt.plot(df['top'], df['evol_di'],color='blue',label='DI EVOL')
    
    plt.scatter(df['top'], df['evol_frob'],color='yellow',label=None)
    plt.plot(df['top'], df['evol_frob'],color='yellow',label='FROB EVOL')
    
    plt.scatter(df['top'], df['evol_psicov'],color='red',label=None)
    plt.plot(df['top'], df['evol_psicov'],color='red',label='PSICOV EVOL')
    
    plt.scatter(df['top'], df['nat_mi'],color='green',label=None)
    plt.plot(df['top'], df['nat_mi'],color='green',label='MI NAT',linestyle=':')
    
    plt.scatter(df['top'], df['nat_di'],color='blue',label=None)
    plt.plot(df['top'], df['nat_di'],color='blue',label='DI NAT',linestyle=':')
    
    plt.scatter(df['top'], df['nat_frob'],color='yellow',label=None)
    plt.plot(df['top'], df['nat_frob'],color='yellow',label='FROB NAT',linestyle=':')
    
    plt.scatter(df['top'], df['nat_psicov'],color='red',label=None)
    plt.plot(df['top'], df['nat_psicov'],color='red',label='PSICOV NAT',linestyle=':')
    
    
    plt.legend(loc="upper right",prop={'size':8})
    plt.xticks([0,0.5,1,2,3,4,5])
    plt.yticks([0,10,20,30,40,50,60,70,80,90,100])
    #plt.axis([0, 20, 0.5, 1])
    plt.xlabel('Tops %')
    plt.ylabel('% de contactos')
    plt.title('Procentaje de contactos por Top. Contacts ' + contact_threshold)
    #plt.show()   
    plt.savefig(output_path +  "_with_nat.png")
    plt.gcf().clear()
    
    '''
    
    plt.scatter(df['top'], df['match_positions_reference'],color='green',label=None)
    plt.plot(df['top'], df['match_positions_reference'],color='green',label='2TRX')
    
    plt.scatter(df['top'], df['match_positions'],color='orange',label=None)
    plt.plot(df['top'], df['match_positions'],color='orange',label='Sumatoria de TOP MI')
    
    plt.scatter(df['top'], df['match_positions_prom'],color='red',label=None)
    plt.plot(df['top'], df['match_positions_prom'],color='red',label='Promedio')
    
    plt.legend(loc="upper right",prop={'size':10})
    plt.xticks([0,0.5,1,2,3,4,5])
    plt.yticks([0,5,10,20,30,40,50,60,70,80,90,100])
    #plt.axis([0, 20, 0.5, 1])
    plt.xlabel('Tops %')
    plt.ylabel('# contactos compartidos con el natural')
    plt.title('Contactos compartidos con el natural. Contacts ' + contact_threshold)
    #plt.show()   
    plt.savefig(output_path + '_match.png')
    
    
    plt.gcf().clear()
    '''

def top_comparation_sub(df, contact_threshold ,output_path):
    #df_to_plot=df[(df['beta']==b)]
    plt.scatter(df['top'], df['nat_contact_%'],color='blue',label=None)
    plt.plot(df['top'], df['nat_contact_%'],color='blue',label='Natural')
    
    
    
    plt.scatter(df['top'], df['ref_contact_%'],color='green',label=None)
    plt.plot(df['top'], df['ref_contact_%'],color='green',label='2TRX')
    
    plt.scatter(df['top'], df['evol_contact_%'],color='orange',label=None)
    plt.plot(df['top'], df['evol_contact_%'],color='orange',label='Sumatoria de TOP MI')
    
    plt.scatter(df['top'], df['prom_contact_%'],color='red',label=None)
    plt.plot(df['top'], df['prom_contact_%'],color='red',label='Promedio')
    
    plt.legend(loc="lower left",prop={'size':8})
    plt.xticks([0,0.5,1,2,3,4,5])
    plt.yticks([0,10,20,30,40,50,60,70,80,90,100])
    #plt.axis([0, 20, 0.5, 1])
    plt.xlabel('Tops %')
    plt.ylabel('% de contactos')
    plt.title('Procentaje de contactos por Top. Contacts ' + contact_threshold)
    #plt.show()   
    
    a = plt.axes([.6, .6, .25, .25])
    #plt.axis([0, 1, 0, 1])
    #no = plt.scatter(y_evol_f, x_nat_f,color="r")
    #co = plt.scatter(y_evol_t, x_nat_t,color="b")
    plt.scatter(df['top'], df['match_positions_reference'],color='green',label=None)
    plt.plot(df['top'], df['match_positions_reference'],color='green',label='2TRX')
    
    plt.scatter(df['top'], df['match_positions'],color='orange',label=None)
    plt.plot(df['top'], df['match_positions'],color='orange',label='Sumatoria de TOP MI')
    
    plt.scatter(df['top'], df['match_positions_prom'],color='red',label=None)
    plt.plot(df['top'], df['match_positions_prom'],color='red',label='Promedio')
    
    plt.legend(loc="upper left",prop={'size':4})
    plt.xticks([0,0.5,1,2,3,4,5])
    #plt.yticks([0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,90,100])
    plt.yticks([0,10,20,30,40,50,60,70,80])
    #plt.axis([0, 20, 0.5, 1])
    #plt.xlabel('Tops %')
    #plt.ylabel('# contactos compartidos con el natural')
    #plt.title('Contactos compartidos con el natural. Contacts ' + contact_threshold)
    #plt.show()   
    plt.savefig(output_path + '_sub.png')
    
    
    plt.gcf().clear()

def auc_optimization(df, field_to_plot,colors,output_path, title):
    betas=df.groupby('beta').first().reset_index()
    betas_list=betas['beta'].tolist()
    i=0
    for b in betas_list:
        df_to_plot=df[(df['beta']==b)]
        plt.scatter(df_to_plot['nsus'], df_to_plot[field_to_plot],color=colors[i],label=None)
        plt.plot(df_to_plot['nsus'], df_to_plot[field_to_plot],color=colors[i],label=b)
        i=i+1
    '''for y in zip(df):
        splt.scatter(nsus, y, color=c,picker=True)  
        splt.plot(nsus, y,color=c, label=beta[i])
        i=i+1    
    '''
    plt.legend(loc="lower right",prop={'size':6},title="Beta")
    plt.xticks([1,2,3,5,7,10,15,20])
    plt.yticks([0.5,0.6,0.7,0.8,0.9,1])
    plt.axis([0, 20, 0.5, 1])
    plt.xlabel('NSUS (Sustituciones no sinonimas por sitio)')
    plt.ylabel('AUC')
    plt.title(title)
    #plt.show()   
    plt.savefig(output_path)
    plt.gcf().clear()
    
    
def contact_map_sum_top_mi_matrix(mat1, mat2,output_file,title='Contact Map'):
    mat1 =np.triu(mat1, -1)
    mat2 =np.triu(mat2, 1)
    '''
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    from matplotlib.colors import LogNorm
    import numpy as np
    x, y, z = np.loadtxt('data.txt', unpack=True)
    '''
    
    #fig = plt.figure()
    #ax1 = fig.add_subplot(111)
    
    fig,ax = plt.subplots()
    pa = ax.imshow(mat1,cmap=cm.YlOrRd,interpolation='nearest')
    cba = plt.colorbar(pa,shrink=0.25,ticks=[0,1,2,3,4,5,6,7,8])
    
    pb = ax.imshow(mat2,cmap=cm.Blues,interpolation='nearest')
    cbb = plt.colorbar(pb,shrink=0.25,ticks=[0,1,2,3,4,5,6,7,8])

    plt.show()
    
    
    #plt.imshow(cmap,cmap=cm.hot)
    plt.title(title)
    plt.savefig(output_file)
    #plt.show()
    plt.gcf().clear() 

def dendogram_matrix(Z, output_path,title,labels):
    
    plt.figure(figsize=(8, 8))
    plt.title(title)
    #plt.xlabel('Indice de entrada (1-50,51-100,101-150)')
    #plt.ylabel('Distancia')
    #max_d = 10
    dendrogram(Z,leaf_rotation=90.,leaf_font_size=10.,show_contracted=True,labels=labels, orientation='top')
    #plt.axhline(y=max_d, c='k') 
    plt.savefig(output_path)
    plt.gcf().clear()