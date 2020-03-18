import os
import argparse
from scipy import stats
import pandas as pd
import numpy as np 
import math
import cooler

# Function by Axel Cournac (https://github.com/koszullab/E_coli_analysis/blob/master/python_codes/directional_indice.py)
def directional(A, nw):
    n1 = A.shape[0]  
    signal1 = np.zeros((n1, 1));
    
    for i in range(0,n1) :
        vect_left = [];
        vect_right = [];
        
        for k in range(i-1,i-nw-1,-1) :
            kp =k; 
            if k < 0 :
                kp = n1 +k ;
            if A[i,kp] > 0 :
                vect_left.append(math.log(A[i,kp]));    
            else :
                vect_left.append(0);  
                    
        for k in range(i+1,i+nw+1) : 
            kp =k;
            if k >= n1 :
                kp = k - n1;
            if A[i,kp] > 0 :
                vect_right.append(math.log(A[i,kp]));    
            else :
                vect_right.append(0);  
                           
        if sum(vect_left) != 0 and sum(vect_right) != 0 :
            signal1[i] =  stats.ttest_rel(vect_right,vect_left)[0];
        else :
            signal1[i] =  0;
                           
    return signal1

# Parsing input arguments
parser = argparse.ArgumentParser()
parser.add_argument('tab_path', type = str, help = 'File name')
parser.add_argument('resolution', type = int, help = 'Resolution')
args = parser.parse_args()

# Path to the input files
tab_path = args.tab_path
resolution = args.resolution
matrix = pd.read_table(tab_path, sep = '\t', header = None).as_matrix()

# Create folders to store output files
out_dir = tab_path.split('/')[1]
replicate = tab_path.split('/')[-1].replace('.tab', '')
if not os.path.exists('./' + out_dir + '/TADs/directionality/'):
    os.mkdir('./' + out_dir + '/TADs/directionality/')

# Call TADs
for scale in [10]:
    print('Calling CIDs at scale:', scale)
    di = directional(matrix, scale)
    di = [item for sublist in di for item in sublist]
    np.save('./' + out_dir + '/TADs/directionality/' + replicate + '.DI.npy', np.array(di))

