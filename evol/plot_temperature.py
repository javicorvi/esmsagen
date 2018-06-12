'''
Show all different interpolation methods for imshow
'''

import matplotlib.pyplot as plt

import numpy as np
import pandas as pd
# from the docs:

# If interpolation is None, default to rc image.interpolation. See also
# the filternorm and filterrad parameters. If interpolation is 'none', then
# no interpolation is performed on the Agg, ps and pdf backends. Other
# backends will fall back to 'nearest'.
#
# http://matplotlib.org/api/pyplot_api.html#matplotlib.pyplot.imshow


def plot_temperature(method):
    
    fig, axes = plt.subplots(nrows=2, ncols=2)
    mul=1
    df=pd.read_csv(method+"_0.5_coevolution.csv",header=0,sep='\t',usecols=['Count','Contacts'],dtype={'Contacts': np.int})    
    counts_df = df.groupby(['Count','Contacts']).size().reset_index(name='Size')
    x=counts_df['Count'].tolist()
    y=counts_df['Contacts'].tolist()
    area=(mul*counts_df['Size'])**2
    #plt.imshow(matrix,cmap=cm.Greys,interpolation='catrom')
    #axes[0,0].rcParams['axes.facecolor'] = 'orange'
    axes[0,0].scatter(x, y, s=area, c='purple', alpha=1.0)
    axes[0,0].set_axis_bgcolor('skyblue')
    axes[0,0].set_ylabel('Contacts',fontsize=8)
    axes[0,0].set_title(method+' top 0.5',fontsize=10)
    axes[0,0].set_xticks([0,1,2,3,4,5,6,7,8])
    axes[0,0].set_yticks([0,1,2,3,4,5,6,7,8])
    axes[0,0].set_xlabel('Sum Top COV',fontsize=8)
    axes[0,0].plot([0,1,2,3,4,5,6,7,8], color='grey')
    df=pd.read_csv(method+"_1_coevolution.csv",header=0,sep='\t',usecols=['Count','Contacts'],dtype={'Contacts': np.int})    
    counts_df = df.groupby(['Count','Contacts']).size().reset_index(name='Size')
    x=counts_df['Count'].tolist()
    y=counts_df['Contacts'].tolist()
    area=(mul*counts_df['Size'])**2
    #plt.imshow(matrix,cmap=cm.Greys,interpolation='catrom')
    #axes[0,1].rcParams['axes.facecolor'] = 'orange'
    axes[0,1].scatter(x, y, s=area, c='blue', alpha=10.0)
    axes[0,1].set_title(method+' top 1',fontsize=10)
    
    
    axes[0,1].set_axis_bgcolor('skyblue')
    axes[0,1].set_xticks([0,1,2,3,4,5,6,7,8])
    axes[0,1].set_yticks([0,1,2,3,4,5,6,7,8])
    axes[0,1].set_ylabel('Contacts',fontsize=8)
    axes[0,1].set_xlabel('Sum Top COV',fontsize=8)
    axes[0,1].plot([0,1,2,3,4,5,6,7,8], color='grey')
    df=pd.read_csv(method+"_2_coevolution.csv",header=0,sep='\t',usecols=['Count','Contacts'],dtype={'Contacts': np.int})    
    counts_df = df.groupby(['Count','Contacts']).size().reset_index(name='Size')
    x=counts_df['Count'].tolist()
    y=counts_df['Contacts'].tolist()
    area=(mul*counts_df['Size'])**2
    #plt.imshow(matrix,cmap=cm.Greys,interpolation='catrom')
    #axes[1,0].rcParams['axes.facecolor'] = 'orange'
    axes[1,0].scatter(x, y, s=area, c='red', alpha=10.0)
    axes[1,0].set_axis_bgcolor('skyblue')
    axes[1,0].set_title(method+' top 2',fontsize=10)
    axes[1,0].set_xticks([0,1,2,3,4,5,6,7,8])
    axes[1,0].set_yticks([0,1,2,3,4,5,6,7,8])
    axes[1,0].set_xlabel('Sum Top COV',fontsize=8)
    axes[1,0].set_ylabel('Contacts',fontsize=8)
    axes[1,0].plot([0,1,2,3,4,5,6,7,8], color='grey')
    df=pd.read_csv(method+"_3_coevolution.csv",header=0,sep='\t',usecols=['Count','Contacts'],dtype={'Contacts': np.int})    
    counts_df = df.groupby(['Count','Contacts']).size().reset_index(name='Size')
    x=counts_df['Count'].tolist()
    y=counts_df['Contacts'].tolist()
    area=(mul*counts_df['Size'])**2
    #plt.imshow(matrix,cmap=cm.Greys,interpolation='catrom')
    #axes[1,1].rcParams['axes.facecolor'] = 'orange'
    axes[1,1].scatter(x, y, s=area, c='green', alpha=10.0)
    axes[1,1].set_axis_bgcolor('skyblue')
    axes[1,1].set_title(method+' top 3',fontsize=10)
    axes[1,1].set_xticks([0,1,2,3,4,5,6,7,8])
    axes[1,1].set_yticks([0,1,2,3,4,5,6,7,8])
    axes[1,1].set_xlabel('Sum Top COV',fontsize=8)
    axes[1,1].set_ylabel('Contacts',fontsize=8)
    axes[1,1].plot([0,1,2,3,4,5,6,7,8], color='grey')
    plt.subplots_adjust(hspace = 0.4)
    
    # plt.show()
    plt.savefig(method+"_correlation.png")
    plt.gcf().clear()
plot_temperature("mi")
plot_temperature("di")
plot_temperature("frob")
plot_temperature("psicov")