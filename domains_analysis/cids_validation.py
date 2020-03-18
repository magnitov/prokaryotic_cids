from scipy import stats
import pandas as pd
import numpy as np
import statistics
import matplotlib.pyplot as plt
import scipy.signal
import matplotlib
import os
from matplotlib.backends.backend_pdf import PdfPages

map1 = pd.read_csv('final_cids/Directionality_Index/Escherichia_rep1', sep = '\s+', header=None)
map2 = pd.read_csv('final_cids/Directionality_Index/Escherichia_rep2', sep = '\s+', header=None)
RNA_seq_data = pd.read_csv('rna_seq/Escherichia_RNAseq.bed', sep = '\s+', header=None)

def RNA_valid(RNA, contact_map1,contact_map2, N, resolution):
    s=sum(RNA[3].values[:])
    cpm=s/1000000
    RNA[3]=RNA[3].values[:]/cpm
    RNA[3]=scipy.signal.savgol_filter(RNA[3], 5, 1)
    vals = (np.arange(0, 2*N, 2) + 1) * resolution/2
    x_value = np.hstack((((-1)*vals)[::-1], vals))
    df1=[]
    bound1=[]
    df2=[]
    bound2=[]


    for i in contact_map1[2]:
        if (i-N*resolution > contact_map1[1].values[0]) and (i+N*resolution < contact_map1[2].values[len(contact_map1)-1]):
            for ind,j in enumerate (RNA[2]):
                if i == j:
                    df1.append(list(RNA[3][ind-N+1:ind+N+1]))
                    bound1.append(j)    
    for i in contact_map2[2]:
        if (i-N*resolution > contact_map2[1].values[0]) and (i+N*resolution < contact_map2[2].values[len(contact_map2)-1]):
            for ind,j in enumerate (RNA[2]):
                if i == j:
                    df2.append(list(RNA[3][ind-N+1:ind+N+1]))
                    bound2.append(j)        
                
    B1 = []
    Opposite1 = [] 
    a1 = np.array(df1)
    number1 = a1.shape[0]
    size1 = a1.shape[1]//2 #size==N    
    
    B2 = []
    Opposite2 = []
    a2 = np.array(df2)
    number2 = a2.shape[0]
    size2 = a2.shape[1]//2 #size==N


    for a_row in a1:
        Opposite1 += list(a_row[size1-10:size1-5])
        B1 += list(a_row[size1-5:size1+5])
        Opposite1 += list(a_row[size1+5:size1+10])
    
    for a_row in a2:
        Opposite2 += list(a_row[size2-10:size2-5])
        B2 += list(a_row[size2-5:size2+5])
        Opposite2 += list(a_row[size2+5:size2+10])


    variable1 = np.subtract(B1, Opposite1)
    c1=np.mean (a1, axis=0) 
    pvalue1 = scipy.stats.wilcoxon(B1, Opposite1).pvalue

    
    variable2 = np.subtract(B2, Opposite2)
    c2=np.mean (a2, axis=0) 
    pvalue2 = scipy.stats.wilcoxon(B2, Opposite2).pvalue

    
    pp = PdfPages('graph-1902/Escherichia_HiCseg.pdf')

    plt.figure(figsize = (7, 5))
    
    plt.plot(x_value,c1, color='b', linewidth=3.0, label='Mean Trend rep1');    
    plt.fill_between(x_value, c1-np.std(c1), c1+np.std(c1),alpha=.1, color='b')
    
    plt.plot(x_value,c2, color='r', linewidth=3.0, label='Mean Trend rep2');    
    plt.fill_between(x_value, c2-np.std(c2), c2+np.std(c2),alpha=.1, color='r')
    res=resolution/1000
    plt.xticks([w for w in range(-12*resolution,13*resolution,3*resolution)],
               ['%i kb'%w for w in range(-int(12*resolution/1000),int(13*resolution/1000),int(3*resolution/1000))]);


 
    if pvalue1 < 10**(-5):
        pvalue1 = '<1e$^{-5}$'
    else:
        pvalue1 = '=%s' %(round(pvalue1, 5))
        
        
    if pvalue2 < 10**(-5):
        pvalue2 = '<1e$^{-5}$'
    else:
        pvalue2 = '=%s' %(round(pvalue2, 5))
    plt.title('Escherichia_HiCseg \n p-rep1%s  p-rep2%s' %(pvalue1, pvalue2))
    plt.xlabel('kb')
    plt.ylabel('RNA-seq level (reads)')
    plt.xlim(-9.5*resolution, 9.5*resolution)
    plt.legend()
    plt.savefig('graph-1902/Escherichia_HiCseg.png')
    plt.savefig(pp, format='pdf')
    pp.close();

RNA_valid(RNA_seq_data, map1, map2, 10, 5000)