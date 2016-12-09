# 
# Usage: python theseus_10x.py  PATH2H5 OUTPUTPATH TPT_THRESHOLD NUMBER_OF_CELLS

import pandas as pd
import numpy as np
import scipy
from scipy import sparse, io
import h5py
import re
import time
import sys

def read_10x_h5(path2file):
    
    with h5py.File(path2file,'r') as hf:
        print('List of arrays in this file: \n', hf.keys())
        start = time.time()
        bc = np.array(hf.get('barcode'))
        bc_cor = np.array(hf.get('barcode_corrected_reads'))
        umi = np.array(hf.get('umi'))
        gene = np.array(hf.get('gene'))
        reads = np.array(hf.get('reads'))
        nonconf_reads = np.array(hf.get('nonconf_mapped_reads'))
        unmap_reads = np.array(hf.get('unmapped_reads'))
        map_pos = np.array(hf.get('conf_mapped_uniq_read_pos'))

        #TABLE = pd.DataFrame([bc,gene,np.int64(reads),nonconf_reads,unmap_reads,map_pos==1])
        stop = time.time()
        
    print 'reading file took ' + str(stop-start)
    
    TABLE=pd.DataFrame()
    start = time.time()
    TABLE['bc']=bc
    TABLE['umi']=umi
    TABLE['gene']=gene
    TABLE['bcumi']=zip(bc,umi)
    TABLE['unique']=[1]*len(TABLE)
    TABLE['map_logical']=map_pos>0
    TABLE['read_counts']=reads+nonconf_reads+unmap_reads

    stop = time.time()

    print 'pandaifying took ' + str(stop-start)
    
    return TABLE
    
#filters to top numcells
def filterT(TABLE,numcells):
    
    start = time.time()
    
    map_rate = np.divide(np.multiply(1.0,np.sum(TABLE['read_counts'][TABLE['map_logical']])),np.sum(TABLE['read_counts']))
    print('mapping rate ',map_rate)
    
    TABLE_mapped = TABLE[TABLE['map_logical']]
    TABLE_mapped_merge = TABLE_mapped.groupby('bc').sum()
    TABLE_mapped_merge=TABLE_mapped_merge.sort_values('read_counts',ascending=False).reset_index()


    TABLE_topcells = TABLE_mapped_merge.iloc[range(numcells)]
    top_bcs        = set(TABLE_topcells['bc'])

    top_bc_logical = [True if x in top_bcs else False for x in TABLE_mapped['bc']]

    TABLE_mapped_filtered = TABLE_mapped[top_bc_logical]

    top_cells_mapping = np.divide(np.double(np.sum(TABLE_mapped_filtered['read_counts'])),np.sum(TABLE_mapped_merge['read_counts']))

    print('top cells mapping rate ',top_cells_mapping)

    stop = time.time()

    total_mapping = np.multiply(map_rate,top_cells_mapping)
    print('total mapping rate ',total_mapping)

    print 'filtering took ' + str(stop-start)
    
    return TABLE_mapped_filtered

path2h5=sys.argv[1]
pathout=sys.argv[2]
T = read_10x_h5(path2h5) 

tpt_filter=float(sys.argv[3])
number_of_cells=int(sys.argv[4])


print('The input path is '+path2h5)
print('The output path is '+pathout)
print('TPT threshold is '+str(tpt_filter))
print('Number of Cells is '+str(number_of_cells))


with h5py.File(path2h5,'r') as hf:
    genenames = np.array(hf.get('gene_names'))

#calculate TPT
BCUMI_group=T.groupby('bcumi').sum()
BCUMI_group=pd.DataFrame(BCUMI_group['read_counts'])
BCUMI_group.columns=['total_reads']
T_tot=T.copy()
T_tot.index=T_tot['bcumi']
T_tot=T_tot.join(BCUMI_group)
T_tot['TPT']=1.0*(np.divide(1.0*T_tot['read_counts'],T_tot['total_reads']))

tpt_logical=T_tot['TPT']>tpt_filter
T_tot_filt=T_tot[tpt_logical]
print('Filtered '+str(np.round(100*(1.0-np.mean(tpt_logical)),4))+'%')


T_F=filterT(T_tot_filt,number_of_cells)

T_F['unique']=np.array([1.0]*len(T_F))
T_F['bcgene']=[(x,y) for x,y in zip(T_F['bc'],T_F['gene'])]
T_FF = T_F.groupby('bcgene').sum()

def str2index(strlist):
    reduced=pd.DataFrame(list(set(np.sort(strlist))))
    reduced=reduced.reset_index()
    reduced.index=reduced[0]
    dftmp=pd.DataFrame(strlist,index=strlist)
    dftmp=dftmp.merge(reduced,how='left')['index']
    return np.array(dftmp),list(reduced[0])

#convert to full expression matrix 

tfrow=[x[1] for x in T_FF.index]
tfrow,gnames=str2index(tfrow)
tfcol=[x[0] for x in T_FF.index]
tfcol,cnames=str2index(tfcol)
tfdata=np.array(T_FF['unique'])
tmpcol=pd.DataFrame(np.unique(tfcol))
tmpcol['unind']=range(len(tmpcol))

dftfcol=pd.DataFrame(tfcol)
dftfcol=dftfcol.merge(tmpcol,on=0)
tfcol=np.array(dftfcol['unind'])
EXPR_MAT=scipy.sparse.csr_matrix((tfdata,(tfrow,tfcol)),shape=(np.max(tfrow)+1,np.max(tfcol)+1)).toarray()
EXPR_MAT=pd.DataFrame(EXPR_MAT)
EXPR_MAT.index=gnames
EXPR_MAT.columns=cnames

gene_labels=pd.DataFrame(genenames)
gene_labels.index=gene_labels.index+1
EXPR_MAT.index=gene_labels.ix[EXPR_MAT.index,0]
EXPR_MAT.to_csv(pathout+'expr_TPTfilt.txt',sep='\t')