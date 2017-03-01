from __future__ import print_function
import sys
import pandas as pd
import numpy.matlib
import numpy as np
import scipy
import time
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sklearn
from sklearn import preprocessing
from sklearn import mixture
from sklearn.neighbors.kde import KernelDensity
import glob
import seaborn as sns
import collections 
sns.set_context('talk')
sns.set_style('white')
sns.set_style("ticks")
import re
from scipy import sparse, io
import os
import math
import csv
import fbpca
from matplotlib import rcParams
import numpy as np
import scipy.stats as stats
from scipy.stats import gaussian_kde
import statsmodels.api as sm
import statsmodels
from statsmodels.distributions.empirical_distribution import ECDF

#GO imports
from goatools.obo_parser import GODag
from goatools.associations import read_ncbi_gene2go
from goatools.go_enrichment import GOEnrichmentStudy

# The following packages are typically not installed by default in Python installations, but would enable some additional functionality
#import Levenshtein (edit_dist)
#import infomap (info_cluster) 
#import networkx as nx (info_cluster)



## progress bar
def update_progress(progress):
    barLength = 10 # Modify this to change the length of the progress bar
    status = ""
    if isinstance(progress, int):
        progress = float(progress)
    if not isinstance(progress, float):
        progress = 0
        status = "error: progress var must be float\r\n"
    if progress < 0:
        progress = 0
        status = "Halt...\r\n"
    if progress >= 1:
        progress = 1
        status = "Done...\r\n"
    block = int(round(barLength*progress))
    text = "\rPercent: [{0}] {1}% {2}".format( "#"*block + "-"*(barLength-block), progress*100, status)
    sys.stdout.write(text)
    sys.stdout.flush()

#######read input############
def read_10x(pathin):
    """Return Pandas Dataframe containing 10x dataset """
    
    mat=scipy.io.mmread(os.path.join(pathin, "matrix.mtx"))
    genes_path = os.path.join(pathin, "genes.tsv")
    gene_ids = [row[0] for row in csv.reader(open(genes_path), delimiter="\t")]
    gene_names = [row[1] for row in csv.reader(open(genes_path), delimiter="\t")]

    gene_final = [x+'_'+y for x,y in zip(gene_ids,gene_names)]
    
    barcodes_path = os.path.join(pathin, "barcodes.tsv")
    barcodes = [row[0][0:14] for row in csv.reader(open(barcodes_path), delimiter="\t")]
    
    DGE=pd.DataFrame(mat.toarray())
    
    DGE.index=gene_final
    DGE.columns=barcodes
    
    return DGE


def genenames_from10x(genelist):
    """Return gene names from 10x index generated with read_10x """
    genesymbol=[]
    #ensemblid=[]
    for i in range(len(genelist)):
    	curgene=genelist[i]
        starts=[]
        for x in re.finditer('_',curgene):
            starts.append(x.start()+1)
        genesymbol.append(curgene[starts[-1]:])
        
    return genesymbol#,ensemblid

def genenames_from10x_mod(genelist):
    """Return gene names from 10x index generated with read_10x """
    genesymbol=[]
    #ensemblid=[]
    for i in range(len(genelist)):
        curgene=genelist[i]
        starts=[]
        for x in re.finditer('_',curgene):
            starts.append(x.start()+1)
        genesymbol.append(curgene[starts[0]:])
        
    return genesymbol#,ensemblid

def collapse2gene(DGE):

    DGE_gene=DGE.copy()
    DGE_gene.index=genenames_from10x(DGE_gene.index)
    DGE_gene=DGE_gene.groupby(DGE_gene.index).sum()

    return DGE_gene

def guide2gene(guide):
    """get genename between underscores"""
    underscore_pos = []
    count=0
    if ('INTERGENIC' in guide):
        nameout='INTERGENIC'
    elif ('_' in guide):
        for x in re.finditer('_',guide):
            if count<2:
                underscore_pos.append(x.span()[1])

        nameout=re.sub('sg','',guide[underscore_pos[0]:underscore_pos[1]-1])
    else:
        nameout=guide
    return nameout

def get_batches(cbcs):
    """Return batch - last underscore in column names"""
    batchvec=[]
    for cell in cbcs:
        starts=[]
        for x in re.finditer('_',cell):
            starts.append(x.start()+1)
        batchvec.append(cell[starts[-2]:])
        
    return np.array(batchvec)
  
def genelevel_dict(GUIDES_DICT):
    """Collapse guide level dictionary to gene level using guide2gene"""
    genes=[guide2gene(x) for x in GUIDES_DICT.keys()]
    GUIDES_DICT_GENES={}
    for gene in genes:
        GUIDES_DICT_GENES[gene]=[]

    for key in GUIDES_DICT.keys():
        GUIDES_DICT_GENES[guide2gene(key)].extend(GUIDES_DICT[key])
    return GUIDES_DICT_GENES

####transform data#######
def tp10k_transform(DGE,norm_factor=1.0e4):
    """normalize columns of pandas dataframe to sum to a constant, by default 10,000"""
    return(norm_factor*(DGE / DGE.sum()))

def Zcells(DGE):
    """Z transformation of columns of pandas"""
    DGEZ=DGE.copy()
    DGEZ=pd.DataFrame(sklearn.preprocessing.scale(DGE,axis=0))
    DGEZ.index=DGE.index
    DGEZ.columns=DGE.columns
    return DGEZ
 
def Zgenes(DGE,batchvec=None):
    """Z transformation of rows of pandas, option for per batch normalization"""
    DGEZ=DGE.copy()
    if batchvec is None:
        DGEZ=pd.DataFrame(sklearn.preprocessing.scale(DGEZ,axis=1))
        DGEZ.columns=DGE.columns
        DGEZ.index=DGE.index
    else:
        batch=np.unique(batchvec)
        for curbatch in batch:
            DGEZ.ix[:,np.array(batchvec)==curbatch]=sklearn.preprocessing.scale(DGEZ.ix[:,np.array(batchvec)==curbatch],axis=1)
    return DGEZ

def Zgenes_floor(DGE,floor=0,batchvec=None):
    """Z transformation of rows of pandas dataframe, with flooring of std dev, option for per batch normalization"""
    DGEZ=DGE.copy()
    if batchvec is None:
        curstd=DGE.std(axis=1)+floor
        curmean=DGE.mean(axis=1)
        curZ=(DGEZ.subtract(curmean,axis=0)).divide(curstd,axis=0)
        DGEZ=curZ
        DGEZ.columns=DGE.columns
        DGEZ.index=DGE.index
    else:
        batch=np.unique(batchvec)
        for curbatch in batch:
            curDGE=DGEZ.ix[:,np.array(batchvec)==curbatch]
            curstd=curDGE.std(axis=1)+floor
            curmean=curDGE.mean(axis=1)
            curZ=(curDGE.subtract(curmean,axis=0)).divide(curstd,axis=0)

            DGEZ.ix[:,np.array(batchvec)==curbatch]=np.array(curZ)
    return DGEZ



def Centergenes(DGE,batchvec=None):
    """Median centering of rows of pandas, option for per batch normalization"""

    DGEC=DGE.copy()
    if batchvec is None:
        DGEC=DGEC.subtract(DGEC.median(axis=1),axis='rows')
    else:
        batch=np.unique(batchvec)
        for curbatch in batch:
            DGEC.ix[:,np.array(batchvec)==curbatch]=DGEC.ix[:,np.array(batchvec)==curbatch].subtract(DGEC.ix[:,np.array(batchvec)==curbatch].median(axis=1),axis='rows')
    return DGEC

def permute_matrix(DGE,bins=20,verbose=0):
    """Permute genes based on similar expression levels"""
    DGE_perm=DGE.copy()

    GSUMS=np.sum(DGE,axis=1)

    breakvec = np.linspace(1,100,bins)

    breaks=[]
    for breaker in breakvec:
        breaks.append(np.percentile(GSUMS,breaker))
    breaks=np.unique(breaks)

    for i in range(len(breaks)-1):
        if verbose==1:
            print(np.round((1.0*i)/(len(breaks)-1)))
        for j in range(len(DGE.columns)):
            curlogical=np.logical_and(GSUMS>breaks[i],GSUMS<=breaks[i+1])
            DGE_perm.ix[curlogical,j]=np.random.permutation(DGE_perm.ix[curlogical,j])
    return DGE_perm

def downsample_reads(DF,per_reads=1.0,nrpc=None):

    DF_mod=DF.copy()


    numgenes=np.shape(DF_mod)[0]
    genenames=DF_mod.index
    DF_mod.index=range(numgenes)
    cells=DF_mod.columns
    
    readspercell=np.sum(DF_mod,axis=0)
    totalreads  =np.sum(readspercell)
    newreads    =np.round(totalreads*per_reads)
    cellpercents=np.divide(1.0*readspercell,totalreads)
    if nrpc:
        newreadspercell=nrpc
    else:
        newreadspercell=[int(x) for x in np.round(np.multiply(cellpercents,newreads))]
    
    DF_out=pd.DataFrame()
 
    for i in range(len(cells)):
        vectorize=[]
        curcell=DF_mod[cells[i]]
        curcell=curcell[curcell!=0]

        for j in curcell.index:
            vectorize.extend([j]*curcell[j])

        vec_sample=np.random.choice(vectorize,size=newreadspercell[i],replace=False)
        sampled_vec=np.histogram(vec_sample,bins=range(numgenes+1))[0]
        DF_out[cells[i]]=sampled_vec
            
    DF_out.index=genenames

    return DF_out
        
def downsampler(DF,percell=1.0,perreads=1.0):
        
    if percell==1.0:
        DF_sampled=DF.copy()
    else:
        newcells=int(np.round(np.shape(DF)[1]*percell))
        DF_sampled=DF.sample(newcells,axis=1)
    
    if perreads==1.0:
        return DF_sampled
    else:
        return downsample_reads(DF_sampled,perreads)

###########generate covariates#########
def dict2X(GUIDES_DICT,cbcs):
    """convert guide cbc dictionary into covariate matrix"""
    X=pd.DataFrame()
    
    for key in GUIDES_DICT.keys():
        curkey=[]
        for cbc in cbcs:
            if cbc in GUIDES_DICT[key]:
                curkey.append(1)
            else:
                curkey.append(0)
        X[key]=np.array(curkey)
        
    X.index=cbcs
    
    return X


def clusters2X(clusters,cbcs):
    """convert cell cluster cbc dictionary into covariate matrix"""
    clusterun=clusters.columns
    X=pd.DataFrame(np.zeros((len(cbcs),len(clusterun))))
    X.index=cbcs

    clusters_intersect=clusters.loc[list(set(clusters.index).intersection(set(cbcs)))]



    X.loc[clusters_intersect.index]=clusters_intersect

    return X


def Xguides2genes(DF):

    Xgene=DF.copy()

    Xgene=Xgene.T
    Xgene.index=[guide2gene(x) for x in Xgene.index]
    Xgene_group=(Xgene.groupby(Xgene.index).sum()>0).sum()
    XgeneF=1.0*(Xgene.groupby(Xgene.index).sum()>0).T

    return XgeneF


def Y2FlatCov(Y,verbose=0):
    ngenes=np.shape(Y)[0]

    triuind=np.triu_indices(ngenes)
    curnames=Y.index
    covgenes=[curnames[x]+'-'+curnames[y] for x,y in zip(triuind[0],triuind[1])]

    triu_mask=np.triu(np.ones((ngenes,ngenes))).astype(np.bool)

    ncells=np.shape(Y)[1]
    i=0

    COVout=pd.DataFrame(np.zeros((len(triuind[0]),ncells)))
    COVout.columns=Y.columns

    for col in Y:
        update_progress(np.divide(1.0*i,ncells))
        cell=pd.DataFrame(Y[col])
        #cell=np.divide(cell,np.linalg.norm(cell))
        cellcov=cell.dot(cell.T)
        triucellcov=cellcov.where(np.triu(np.ones(cellcov.shape)).astype(np.bool)).values.flatten()
        triucellcov=triucellcov[~numpy.isnan(triucellcov)]

        COVout[col]=triucellcov
        i+=1
    COVout.index=covgenes

    return COVout

def create_interactions(DF):
    """Take covariate matrix and generate pairwise interaction matrix between covariates"""
    INTERACT=pd.DataFrame()
    dfcolumns=DF.columns
    groupthese=[]

    for i in range(len(dfcolumns)):
        for j in range(len(dfcolumns)):
            
            name1=dfcolumns[i]
            name2=dfcolumns[j] 
            
            if(i<j):
                twonames=np.sort(list(set([str(name1),str(name2)])))
                if len(twonames)==2:
                    INTERACT[str(name1)+'-'+str(name2)]=np.array(DF.ix[:,i])*np.array(DF.ix[:,j])
                    groupthese.append(str(twonames[0])+'-'+str(twonames[1]))

    #INTERACT.columns=[guide2gene(x.split('-')[0])+'-'+guide2gene(x.split('-')[1]) for x in INTERACT.columns]
    INTERACT=INTERACT.T
    INTERACT['genes']=INTERACT.index
    INTERACT=INTERACT.groupby(groupthese).sum().T
    INTERACT=INTERACT>0
    INTERACT.index=DF.index
    
    return(1.0*INTERACT)

def create_interactions_nothresh(DF):
    """Take covariate matrix and generate pairwise interaction matrix between covariates"""
    INTERACT=pd.DataFrame()
    dfcolumns=DF.columns
    groupthese=[]

    for i in range(len(dfcolumns)):
        for j in range(len(dfcolumns)):
            
            name1=dfcolumns[i]
            name2=dfcolumns[j] 
            
            if(i<j):
                twonames=np.sort(list(set([str(name1),str(name2)])))
                if len(twonames)==2:
                    INTERACT[str(name1)+'-'+str(name2)]=np.array(DF.ix[:,i])*np.array(DF.ix[:,j])
                    groupthese.append(str(twonames[0])+'-'+str(twonames[1]))

    #INTERACT.columns=[guide2gene(x.split('-')[0])+'-'+guide2gene(x.split('-')[1]) for x in INTERACT.columns]
    INTERACT=INTERACT.T
    INTERACT=INTERACT.groupby(groupthese).sum().T
    INTERACT.index=DF.index
    
    return(1.0*INTERACT)


def create_3_interactions(DF):
    """Take covariate matrix and generate three-way interaction matrix between covariates"""
    INTERACT=pd.DataFrame()
    dfcolumns=DF.columns
    groupthese=[]

    for i in range(len(dfcolumns)):
        for j in range(len(dfcolumns)):
            for k in range(len(dfcolumns)):
                
                if((i<j)&(i<k)):
                    name1=dfcolumns[i]
                    name2=dfcolumns[j] 
                    name3=dfcolumns[k] 
                    threenames=np.sort(list(set([str(name1),str(name2),str(name3)])))
                    if len(threenames)==3:
                        INTERACT[str(name1)+'-'+str(name2)+'-'+str(name3)]=np.array(DF.ix[:,i])*np.array(DF.ix[:,j])*np.array(DF.ix[:,k])
                        groupthese.append(str(threenames[0])+'-'+str(threenames[1])+'-'+str(threenames[2]))
    
    #INTERACT.columns=[guide2gene(x.split('-')[0])+'-'+guide2gene(x.split('-')[1])+'-'+guide2gene(x.split('-')[2]) for x in INTERACT.columns]
    INTERACT=INTERACT.T
    INTERACT['genes']=INTERACT.index
    INTERACT=INTERACT.groupby(groupthese).sum().T
    INTERACT=INTERACT>0
    INTERACT.index=DF.index
    
    return(1.0*INTERACT)

def create_3_interactions_nothresh(DF):
    """Take covariate matrix and generate three-way interaction matrix between covariates"""
    INTERACT=pd.DataFrame()
    dfcolumns=DF.columns
    groupthese=[]

    for i in range(len(dfcolumns)):
        for j in range(len(dfcolumns)):
            for k in range(len(dfcolumns)):
                
                if((i<j)&(i<k)):
                    name1=dfcolumns[i]
                    name2=dfcolumns[j] 
                    name3=dfcolumns[k] 
                    threenames=np.sort(list(set([str(name1),str(name2),str(name3)])))
                    if len(threenames)==3:
                        INTERACT[str(name1)+'-'+str(name2)+'-'+str(name3)]=np.array(DF.ix[:,i])*np.array(DF.ix[:,j])*np.array(DF.ix[:,k])
                        groupthese.append(str(threenames[0])+'-'+str(threenames[1])+'-'+str(threenames[2]))
    
    #INTERACT.columns=[guide2gene(x.split('-')[0])+'-'+guide2gene(x.split('-')[1])+'-'+guide2gene(x.split('-')[2]) for x in INTERACT.columns]
    INTERACT=INTERACT.T
    INTERACT['genes']=INTERACT.index
    INTERACT=INTERACT.groupby(groupthese).sum().T
    INTERACT.index=DF.index
    
    return(1.0*INTERACT)

#############Linear Model Stuff########

def cv_rsq(Y,X,k=5,per=0.8,adj=[],relcel=[]):
    Y_tmp=Y.copy()
    X_tmp=X.copy()
    rsq=[]
    for i in range(k):
        print(i)
        numsamples=int(np.round(per*len(Y_tmp)))
        train=np.random.choice(range(len(Y_tmp)),size=numsamples,replace=False)
        traincells=Y_tmp.index[train]
        testcells=list(set(Y_tmp.index)-set(traincells))
        print('1',len(testcells))
        if len(relcel)>0:
            testcells=list(set(testcells).intersection(set(relcel)))
        print('2',len(testcells))
        Y_train=Y_tmp.loc[traincells]
        Y_test=Y_tmp.loc[testcells]
        flag=0

        X_train=X_tmp.loc[traincells]
        X_test=X_tmp.loc[testcells]

        lmfit=sklearn.linear_model.ElasticNet(precompute=True,l1_ratio=0.5,alpha=0.0005,max_iter=10000)#linear_model.Ridge
        lmfit.fit(X_train,Y_train)

        if len(adj)>0:
            X_train_adj=bayes_cov_col(Y_train,X_train,adj,lmfit)
            lmfit.fit(X_train_adj,Y_train)
            X_test_adj=bayes_cov_col(Y_test,X_test,adj,lmfit)
            rsq.append(lmfit.score(X_test_adj,Y_test))
        else:
            rsq.append(lmfit.score(X_test,Y_test))
    return rsq

def marginal_covariates(y,x,k=4,percent=0.8):
    """Input is observations and list of covariates
    like guides, qc, batch, guide interactions, cell types, cell type interactions
    perform k-fold CV on xx percent of data
    for each of the 2^n combinations of covariates
    """
    if isinstance(x,list):
        numsamples=int(np.round(percent*len(y)))
        X=pd.concat(x,axis=1)
#        rsqall=[]
#        for i in range(k):
#            print(i)
#            train=np.random.choice(range(len(y)),size=numsamples,replace=False)
#            traincells=y.index[train]
#            testcells=list(set(y.index)-set(traincells))
#            X_train=X.loc[traincells]
#            Y_train=y.loc[traincells]
#            X_test=X.loc[testcells]
#            Y_test=y.loc[testcells]

#            enet=sklearn.linear_model.ElasticNet(precompute=True,l1_ratio=0.5,alpha=0.0012,max_iter=10000)
#            enet.fit(X_train,Y_train)
#            print('model has been fit')

#            rsqall.append(enet.score(X_test,Y_test))
        rsqind=[]
        big_resid=[]
        for j in range(len(x)):
            print(j)
            rsqk=[]
            for i in range(k):
                print(k)
                train=np.random.choice(range(len(y)),size=numsamples,replace=False)
                traincells=y.index[train]
                testcells=list(set(y.index)-set(traincells))
                Y_train=y.loc[traincells]
                Y_test=y.loc[testcells]
                flag=0

                if j==0:
                    X_train=x[j].loc[traincells]
                    X_test=x[j].loc[testcells]
                    lmfit=sklearn.linear_model.ElasticNet(precompute=True,l1_ratio=0.5,alpha=0.0005,max_iter=10000)
                else:
                    X=pd.concat(x[0:j],axis=1)
                    X_train=X.loc[traincells]
                    X_test=X.loc[testcells]
                
                lmfit.fit(X_train,Y_train)
                rsqk.append(lmfit.score(X_test,Y_test))
                Yhat=lmfit.predict(X_test)

                if flag==0:
                    df_resid=Yhat-Y_test
                    flag=1
                else:
                    df_resid = (df_resid + (Yhat-Y_test)) / 2.0 

            rsqind.append(rsqk)
            big_resid.append(df_resid)

    else:
        print('x is not a list')
        return

    #df_rsq=pd.concat([pd.DataFrame(rsqind)],axis=0)
    return rsqind

def crosscov_interactions(X1,X2):

    cols1=X1.columns
    cols2=X2.columns

    Xout=pd.DataFrame()

    for i in range(len(cols1)):
        for j in range(len(cols2)):
            if i>j:

                curi=cols1[i]
                curj=cols2[j]
                Xout[str(curi)+'_'+str(curj)]=X1[curi]*X2[curj]

    return Xout

def nonzeroX2dict(X):
    dict_out={}
    for col in X.columns:
        curcol=X[col]
        dict_out[col]=curcol[curcol>0].index
    return dict_out


def bayes_cov_col(Y,X,cols,lm):

    #EM iterateit
    Yhat=pd.DataFrame(lm.predict(X))
    Yhat.index=Y.index
    Yhat.columns=Y.columns
    SSE_all=np.square(Y.subtract(Yhat))
    X_adjust=X.copy()


    df_SSE   = []
    df_logit = []

    for curcov in cols:

        curcells=X[X[curcov]>0].index

        if len(curcells)>2:

            X_notcur=X.copy()
            X_notcur[curcov]=[0]*len(X_notcur)

            X_sub=X_notcur.loc[curcells]

            Y_sub=Y.loc[curcells]

            Yhat_notcur=pd.DataFrame(lm.predict(X_sub))
            Yhat_notcur.index=Y_sub.index
            Yhat_notcur.columns=Y_sub.columns

            SSE_notcur=np.square(Y_sub.subtract(Yhat_notcur))
            SSE=SSE_all.loc[curcells].subtract(SSE_notcur)
            SSE_sum=SSE.sum(axis=1)

            #SSE_transform=SSE.div(GENE_std)[vargenes].sum(axis=1)
            logitify=np.divide(1.0,1.0+np.exp(SSE_sum))

            df_SSE.append(SSE_sum)
            df_logit.append(logitify)

            X_adjust[curcov].loc[curcells]=logitify

    return X_adjust


def run_model(Y,X,EM_DICT=None,verbose=0,modalpha=0.0005,removecells=1):

    enet=sklearn.linear_model.ElasticNet(precompute=True,l1_ratio=0.5,alpha=modalpha,max_iter=10000)
    enet.fit(X,Y)
    if verbose==1:
        print(enet.score(X,Y))

    Be=pd.DataFrame(enet.coef_)
    Be.columns=X.columns
    Be.index=Y.columns

    #EM iterateit
    Yhat=pd.DataFrame(enet.predict(X))
    Yhat.index=Y.index
    Yhat.columns=Y.columns
    SSE_all=np.square(Y.subtract(Yhat))

    X_adjust=X.copy()
    X_adjust['unperturbed']=[0]*len(X)

    df_SSE   = []
    df_logit = []
    df_pf    = []

    if EM_DICT is not None:

        for curcov in EM_DICT.keys():

            curcells=EM_DICT[curcov]

            X_notcur=X.copy()
            X_notcur[curcov]=[0]*len(X_notcur)

            X_sub=X_notcur.loc[curcells]

            Y_sub=Y.loc[curcells]

            GENE_var=2.0*Y_sub.var(axis=0)
            vargenes=GENE_var[GENE_var>0].index


            Yhat_notcur=pd.DataFrame(enet.predict(X_sub))
            Yhat_notcur.index=Y_sub.index
            Yhat_notcur.columns=Y_sub.columns

            SSE_notcur=np.square(Y_sub.subtract(Yhat_notcur))
            SSE=SSE_all.loc[curcells].subtract(SSE_notcur)
            SSE_sum=SSE.sum(axis=1)

            SSE_transform=SSE.div(GENE_var+0.5)[vargenes].sum(axis=1)
            logitify=np.divide(1.0,1.0+np.exp(SSE_sum))#SSE_transform))#sum))

            df_SSE.append(SSE_sum)
            df_logit.append(logitify)
            pf=np.mean(logitify>0.99)

            if verbose==1:
                
                print(curcov,pf)
            df_pf.append([curcov,pf])
            weak_perturb=1.0*(logitify<0.1)
            X_adjust[curcov].loc[curcells]=logitify
            X_adjust['unperturbed'].loc[curcells]=weak_perturb

        print('done with EM')

        #refit model

        enet=sklearn.linear_model.ElasticNet(precompute=True,l1_ratio=0.5,alpha=0.0005,max_iter=10000)

        if removecells==1:
            goodcells=X_adjust['unperturbed']!=1
            print(np.mean(goodcells))
            Y=Y[goodcells]
            X_adjust=X[goodcells]
     
        enet.fit(X_adjust,Y)
        Yhat=pd.DataFrame(enet.predict(X_adjust))
        Yhat.index=Y.index
        Yhat.columns=Y.columns

        if verbose==1:
            print(enet.score(X_adjust,Y))

        Be=pd.DataFrame(enet.coef_)
        Be.columns=X_adjust.columns
        Be.index=Y.columns
    RES_out=Y.subtract(Yhat)  

    if EM_DICT is not None:
        return(Be,X_adjust,RES_out,df_pf)#,df_SSE,df_logit)

    return(Be,X_adjust,RES_out)#,df_SSE,df_logit)


def run_model_bycol(Y,X,EM_cols=None,modalpha=0.005,verbose=0):

    enet=sklearn.linear_model.ElasticNet(precompute=True,l1_ratio=0.5,alpha=modalpha,max_iter=10000)
    enet.fit(X,Y)
    if verbose==1:
        print(enet.score(X,Y))

    Be=pd.DataFrame(enet.coef_)
    Be.columns=X.columns
    Be.index=Y.columns

    Yhat=pd.DataFrame(enet.predict(X))
    Yhat.index=Y.index
    Yhat.columns=Y.columns

    if EM_cols is not None:

        X_adjust=bayes_cov_col(Y,X,EM_cols,enet)

        #print('done with EM')

        #refit model

        enet=sklearn.linear_model.ElasticNet(precompute=True,l1_ratio=0.5,alpha=0.0004,max_iter=10000)
     
        enet.fit(X_adjust,Y)
        Yhat=pd.DataFrame(enet.predict(X_adjust))
        Yhat.index=Y.index
        Yhat.columns=Y.columns

        if verbose==1:
            print(enet.score(X_adjust,Y))

        Be=pd.DataFrame(enet.coef_)
        Be.columns=X_adjust.columns
        Be.index=Y.columns
    else:
        X_adjust=X.copy()

    RES_out=Y.subtract(Yhat)  

    return(Be,X_adjust,RES_out)


def count_27(B1,B2,B3,thresh=0.01):
    vecs1=[B1<(-thresh),np.abs(B1)<=thresh,B1>thresh]
    vecs2=[B2<(-thresh),np.abs(B2)<=thresh,B2>thresh]
    vecs3=[B3<(-thresh),np.abs(B3)<=thresh,B3>thresh]
    COUNTER=[]
    for i in range(3):
        for j in range(3):
            for k in range(3):
                COUNTER.append(np.sum(np.logical_and(np.logical_and(vecs1[i],vecs2[j]),vecs3[k])))
    return COUNTER

def return_sorted_list(in1):
    output = [0] * len(in1)
    for i, x in enumerate(sorted(range(len(in1)), key=lambda y: in1[y])):
        output[x] = i
    return np.array(output)

def index_27(B1,B2,B3,df_order,thresh=0.01):
    vecs1=[B1<(-thresh),np.abs(B1)<=thresh,B1>thresh]
    vecs2=[B2<(-thresh),np.abs(B2)<=thresh,B2>thresh]
    vecs3=[B3<(-thresh),np.abs(B3)<=thresh,B3>thresh]

    Ball=pd.concat([B1,B2,B3],axis=1)

    iarray=pd.DataFrame(['none']*len(B1))
    iarray.index=B1.index

    for i in range(3):
        for j in range(3):
            for k in range(3):
                totsum=int(np.sum(np.logical_and(np.logical_and(vecs1[i],vecs2[j]),vecs3[k])))
                iarray[np.logical_and(np.logical_and(vecs1[i],vecs2[j]),vecs3[k])]=str(i-1)+' '+str(j-1)+' '+str(k-1)

    iarray['type']=['none']*len(B1)
    iarray['order']=[0]*len(B1)
    iarray['effect']=[0]*len(B1)

    numbering=0

    for i in range(len(df_order)):

        curgroup=df_order.index[i]
        curtype=df_order.ix[i,'type']

        matches=iarray[0]==curgroup
        nummatches=np.sum(matches)

        if nummatches>0:

            Bmatch=Ball[matches]
            intarray=[int(x) for x in curgroup.split(' ')]

            Bmod=Bmatch.copy()
            l=0
            for col in Bmod.columns:
                Bmod[col]=intarray[l]*Bmod[col]
                l+=1

            Bsum=pd.DataFrame(Bmod.sum(axis=1))
            ordervec=return_sorted_list(-np.array(Bsum[0]))
            ordervec=ordervec+numbering

            iarray.ix[matches,'type']=curtype
            iarray.ix[matches,'order']=ordervec
            iarray.ix[matches,'effect']=np.array(Bsum[0])

            numbering+=np.max(ordervec)+1

    return iarray

def hyper_overlap(genes1,genes2,M):

    curoverlap=genes1.intersection(genes2)
    x=len(curoverlap)
    n=len(genes1)
    N=len(genes2)
    pval=1.0-scipy.stats.hypergeom.cdf(x,M, n, N)
    return pval

def hyper_category(df_cats,genes_in):
    pvals=[]
    cat_un=np.unique(df_cats[0])
    genes2=set(genes_in).intersection(set(df_cats.index))

    for cat in cat_un:
        genes1=set(df_cats[df_cats[0]==cat].index)
        pvals.append(hyper_overlap(genes1,genes2,len(df_cats)))
    df_pvals=pd.DataFrame(-np.log10(statsmodels.sandbox.stats.multicomp.fdrcorrection0(pvals)[1]))
    df_pvals.index=cat_un
    return df_pvals

def numbins(x):
    iqr=((np.percentile(x, 75) - np.percentile(x, 25)))
    if iqr==0.0:
        return int(np.ceil(np.sqrt(len(x))))
    else:    
        bins=int(np.ceil((np.max(x)-np.min(x))/((iqr)/np.power(len(x),0.33333333))))
        return bins

def get_1sidepval(B,joint,edges,gsums,gvar,nguides):


    Bpval=B.copy()

    #create index lookup for each gene to the pairs
    genevec=np.array(range(len(gsums)))
    guidevec=np.array(range(len(nguides)))
    gsums=np.array(gsums)
    gvar=np.array(gvar)
    nguides=np.array(nguides)
    rowindex_dict={}
    colindex_dict={}

    for i in range(len(edges[0])-1):
        for j in range(len(edges[1])-1):
            logical_gsums=np.logical_and(gsums>=edges[0][i],gsums<edges[0][i+1])
            logical_gvar=np.logical_and(gvar>=edges[1][j],gvar<edges[1][j+1])

            logical_both=np.logical_and(logical_gsums,logical_gvar)

            if np.sum(logical_both)>0:
                rowindex_dict[(i,j)]=genevec[logical_both]

    for i in range(len(edges[2])-1):

        logical_nguides=np.logical_and(nguides>=edges[2][i],nguides<edges[2][i+1])
        if np.sum(logical_nguides)>0:
            colindex_dict[i]=guidevec[logical_nguides]

    maxedges=len(edges[3])-2

    for key in rowindex_dict.keys():

        for guidekey in colindex_dict.keys():

            curjoint=joint[key[0]][key[1]][guidekey]
            curjoint /= curjoint.sum()
            curjoint=pd.DataFrame(curjoint)
            curjoint.index=edges[3][:-1]
            curjoint=curjoint.cumsum()

            curmat=Bpval.ix[rowindex_dict[key],colindex_dict[guidekey]]
            
            lookup_mat=curmat.copy()

            bp=pd.DataFrame(np.searchsorted(curjoint.index,curmat))

            bpmax=bp>maxedges
            bp[bpmax]=0

            for i in range(np.shape(bp)[1]):
                lookup=1.0-np.round(np.array(curjoint)[bp.ix[:,i]],10)
                lookup_mat.ix[:,i]=lookup
                lookup_mat.ix[np.where(bpmax)]=0

            Bpval.ix[rowindex_dict[key],colindex_dict[guidekey]]=lookup_mat

    Bpval[B<=0]=1.0

    return Bpval


#create permuted coefficient matrix

def shuffle_mat(X,Xother,Y):
    flag=0
    X_shuffle=X.copy()
    X_shuffle.index=np.random.permutation(X.index)
    X_shuffle=X_shuffle.loc[Y.index]
    X3_shuffle=pd.concat([X_shuffle,Xother],axis=1)
        
    return X3_shuffle



def make_simple_shufs(X,Xother,Y,modalpha=0.005,shufnum=3):
    Be_shuffs=pd.DataFrame()
    flag=0
    for i in range(shufnum):
        print(i)
        X3_shuffle=shuffle_mat(X,Xother,Y)
        Be_shuf,X_adjust,RES=run_model(Y,X3_shuffle,modalpha=modalpha,verbose=0)
        if flag==0:
            Be_shuffs=Be_shuf
            flag=1
        else:
            Be_shuffs=pd.concat([Be_shuffs,Be_shuf])
    return Be_shuffs

def make_shufs(X,Xother,Y,shufnum=3,modalpha=0.005,verbose=1):
    Be_shuffs=pd.DataFrame()
    flag=0
    for i in range(shufnum):
        if verbose==1:
            print(i)
        X3_shuffle=shuffle_mat(X,Xother,Y)
        Be_shuf,X_adjust,RES=run_model_bycol(Y,X3_shuffle,modalpha=modalpha,EM_cols=X.columns,verbose=0)
        if flag==0:
            Be_shuffs=Be_shuf
            flag=1
        else:
            Be_shuffs=pd.concat([Be_shuffs,Be_shuf])
    return Be_shuffs

def make_shufs_linear_sub(X,Xother,Y,shufnum=3):
    Be_shuffs=pd.DataFrame()
    flag=0
    for i in range(shufnum):
        X3_shuffle=shuffle_mat(X,Xother,Y)
        from sklearn import linear_model
        lm=linear_model.Ridge(fit_intercept=True,max_iter=10000)
        lm.fit(X3_shuffle,Y)
        Be_shuf=pd.DataFrame(lm.coef_)
        Be_shuf.index=Y.columns
        Be_shuf.columns=X3_shuffle.columns
        if flag==0:
            Be_shuffs=Be_shuf
            flag=1
        else:
            Be_shuffs=pd.concat([Be_shuffs,Be_shuf])
    return Be_shuffs

def make_shufs_linear(X,Y,shufnum=3):
    Be_shuffs=pd.DataFrame()
    flag=0
    for i in range(shufnum):
        X_shuffle=X.copy()
        X_shuffle.index=np.random.permutation(X.index)
        X_shuffle=X_shuffle.loc[Y.index]
        from sklearn import linear_model
        lm=linear_model.Ridge(fit_intercept=True,max_iter=10000)
        lm.fit(X_shuffle,Y)
        Be_shuf=pd.DataFrame(lm.coef_)
        Be_shuf.index=Y.columns
        Be_shuf.columns=X_shuffle.columns
        if flag==0:
            Be_shuffs=Be_shuf
            flag=1
        else:
            Be_shuffs=pd.concat([Be_shuffs,Be_shuf])
    return Be_shuffs

#get FDR matrix
def fdr_coefs(B,B_shuf,gsums,gvar,nguides,mybins=[30,30,20,1000]):

    numshufs=(1.0*len(B_shuf))/len(B)
    if numshufs%1!=0:
        print('you screwed up permuted is not integer multiple of nonpermuted')
        return
    numshufs=int(numshufs)

    gsums_rep=np.array([list(gsums)]*numshufs).flatten()
    gvar_rep=np.array([list(gvar)]*numshufs).flatten()
    nguides=np.array(nguides)

    flag=0
    for i in range(np.shape(B_shuf)[1]):
        datas=pd.DataFrame([gsums_rep,gvar_rep,np.array([nguides[i]]*len(gsums_rep)),np.array(B_shuf.ix[:,i])]).T
        if flag==0:
            SHUFCOV=datas
            flag=1
        else:
            SHUFCOV=pd.concat([SHUFCOV,datas])

    numBins = mybins  # number of bins in each dimension
    SHUFPOS=SHUFCOV.copy()
    SHUFPOS=SHUFPOS[SHUFPOS[3]>=0]
    joint_pos, edges_pos = np.histogramdd(np.array(SHUFPOS), bins=numBins)
    joint_pos /= joint_pos.sum()

    SHUFNEG=SHUFCOV.copy()
    SHUFNEG=SHUFNEG[SHUFNEG[3]<=0]
    SHUFNEG[3]=SHUFNEG[3].abs()
    joint_neg, edges_neg = np.histogramdd(np.array(SHUFNEG), bins=numBins)
    joint_neg /= joint_neg.sum()
    
    print('Created 4D Null Distributions')

    B_sign = np.sign(B)

    Bpos=B.copy()
    Bpos[B<0]=0

    Bneg=B.copy()
    Bneg[B>0]=0
    Bneg=Bneg.abs()

    Bpval_pos=get_1sidepval(Bpos,joint_pos,edges_pos,gsums,gvar,nguides)
    print('positive pvals calculated')
    Bpval_neg=get_1sidepval(Bneg,joint_neg,edges_neg,gsums,gvar,nguides)
    print('negative pvals calculated')

    BFDR=Bpval_pos.copy()
    BFDR[Bpval_neg<1]=Bpval_neg[Bpval_neg<1]

    for col in BFDR.columns:
        curcol=BFDR[col]
        curcol_logical=curcol<1

        BFDR.ix[curcol_logical,col]=-np.log10(statsmodels.sandbox.stats.multicomp.fdrcorrection0(curcol[curcol_logical])[1])

    BFDR=np.multiply(B_sign,BFDR)

    print('FDR correction performed')

    return BFDR
  
#get FDR matrix
def fdr_colwise_coefs(B,B_shuf):

    BFDR=B.copy()

    for col in BFDR.columns:
        
        curcol=B[col]
        curfdr=BFDR[col]
        curecdf=ECDF(B_shuf[col])

        curcol_pos=curcol>0
        curcol_neg=curcol<0
        sign_col=np.sign(curcol)

        curfdr[curcol_pos]=-np.log10(statsmodels.sandbox.stats.multicomp.fdrcorrection0(1.0-curecdf(curcol[curcol_pos]))[1])
        curfdr[curcol_neg]=np.log10(statsmodels.sandbox.stats.multicomp.fdrcorrection0(curecdf(curcol[curcol_neg]))[1])

        BFDR[col]=curfdr

    print('FDR correction performed')

    return BFDR

def pointwise_p_colwisefdr(B,Bshuf):

    BFDR=B.copy()

    for col in B.columns:
        probs=[]
        sign=[]
        for ind in B.index:
            curecdf=ECDF(Bshuf[col].ix[ind])
            curval=B[col].ix[ind]
            if curval>0:
                sign.append(1)
                probs.append(1.0-curecdf(B[col].ix[ind]))
            else:
                sign.append(-1)
                probs.append(curecdf(B[col].ix[ind]))

        probs=np.array(probs)
        sign=np.array(sign)
    
        BFDR[col]=sign*(-np.log10(statsmodels.sandbox.stats.multicomp.fdrcorrection0(probs)[1]))
    return BFDR

def pointwise_p_rowwisefdr(B,Bshuf):

    BFDR=B.copy()
    SIGN=B.copy()

    for col in B.columns:
        probs=[]
        sign=[]
        for ind in B.index:
            curecdf=ECDF(Bshuf[col].ix[ind])
            curval=B[col].ix[ind]
            if curval>0:
                sign.append(1)
                probs.append(1.0-curecdf(B[col].ix[ind]))
            else:
                sign.append(-1)
                probs.append(curecdf(B[col].ix[ind]))

        probs=np.array(probs)
        sign=np.array(sign)
        SIGN[col]=sign
        BFDR[col]=probs
        
    #rowwise FDR
    for ind in B.index:
        BFDR.ix[ind,:]=SIGN.ix[ind,:]*(-np.log10(statsmodels.sandbox.stats.multicomp.fdrcorrection0(BFDR.ix[ind,:])[1]))

    return BFDR

def fregression_fdr(X,Y,B):
    FDR=pd.DataFrame()
    for i in Y.columns:
        pvals=-np.log10(sklearn.feature_selection.f_regression(X, Y[i])[1])
        FDR[i]=pvals
    FDR.index=X.columns
    FDR=FDR.T
    FDR=np.sign(B)*FDR
    return FDR
    
def compare_fdrs(BFDR,BFDR_down,thresh1=1.3,thresh2=1.3):

    COMPARE=pd.DataFrame()

    for col in BFDR.columns:
        col1=BFDR[col]
        col2=BFDR_down[col]

        a=np.sign(col1)*(col1.abs()>thresh1)
        b=np.sign(col2)*(col2.abs()>thresh2)

        CONF=sklearn.metrics.confusion_matrix(a,b,labels=[-1,0,1])
        tp=CONF[0][0]+CONF[2][2]
        fp=CONF[0][2]+CONF[2][0]+CONF[1][0]+CONF[1][2]
        fn=CONF[0][1]+CONF[2][1]
        tn=CONF[1][1]

        sensitvitiy=np.divide(1.0*tp,tp+fn)
        specificity=np.divide(1.0*tn,tn+fp)

        COMPARE[col]=[tp,tn,fp,fn,sensitvitiy,specificity]

    COMPARE.index=['TP','TN','FP','FN','Sensitivity','Specificity']

    return(COMPARE)

#############data filtering############

def fano_variable(DGEtpm,input_mean=None,meanthresh=0.5,resthresh=0.05,f=0.25,highlight_genes=None,plot=0):
    #get mean and std for each gene
    if input_mean is None:
        popmean=np.log2(np.mean(DGEtpm,axis=1)+1)
    else:
        popmean=input_mean
        
    popstd=np.std(np.log2(DGEtpm+1),axis=1)#np.divide(np.std(DGEtpm,axis=1),popmean)
    thresh=meanthresh
    
    x=popmean[np.array(popmean>thresh)]
    y=popstd[np.array(popmean>thresh)]
    
    DGE_fit=DGEtpm[np.array(popmean>thresh)]
    
    #fit line
    lowess = sm.nonparametric.lowess
    lz_pred = lowess(y, x,frac=f,return_sorted=False)

    residuals=y-lz_pred
    
    if plot==1:
        plt.scatter(x,y,c=['red' if z>resthresh else 'blue' for z in residuals])
        plt.xlabel('log2(Population Mean)')
        plt.ylabel('Standard Deviation')

    df_res=pd.DataFrame()
    df_res['residuals']=residuals
    df_res['mean']=x
    df_res['std']=y
    df_res.index=DGE_fit.index

    if highlight_genes:
        if plot==1:
            subset=df_res.loc[highlight_genes].dropna()
            for thisgene in subset.index:
                df_tmp=subset.loc[thisgene]
                plt.text(df_tmp['mean'],df_tmp['std'],thisgene,fontsize=16)
    
    return df_res
    #return variable genes

##########PCA stuff#########################

def fb_pca(DGE,k=50):
    if 'fbpca' in sys.modules:
        [Ufb,Sfb,Vfb]=fbpca.pca(DGE,k)
    else:
        pca=sklearn.decomposition.PCA(n_components=k)
        pca.fit(DGE)
        Ufb=pca.fit_transform(DGE)
        Sfb=pca.explained_variance_
        Vfb=pca.components_
    Vfb=pd.DataFrame(Vfb).T
    Vfb.index=DGE.columns
    Ufb=pd.DataFrame(Ufb)
    Ufb.index=DGE.index
    return Ufb,Sfb,Vfb

def project_ontoPC(U,DGE):

    DGE_t=DGE.copy().T
    Vfb_new=pd.DataFrame()
    #sscells=DGE_t.pow(2).sum(axis=1)
    for i in range(np.shape(U)[1]):
        Vfb_new[i]=DGE_t.dot(U[i])
    Vfb_new=Vfb_new.T
    Vfb_new=Vfb_new#/sscells
    return Vfb_new.T

def columnwise_compare_innermax(U1,U2):
    U_big=U1.copy()
    U_big=U_big.merge(U2,left_index=True,right_index=True)
    U_big.columns=range(len(U_big.columns))

    genes1=set(U1.index)
    genes2=set(U2.index)
    jac=np.divide(1.0*len(genes1.intersection(genes2)),len(genes1.union(genes2)))

    print('jaccard gene overlap =',jac)

    cols1=np.shape(U1)[1]
    cols2=np.shape(U2)[1]

    mincols=np.min([cols1,cols2])

    print(np.shape(U_big))

    comparevec=[]
    for i in range(mincols):
        corrs=[]
        for j in range(mincols):

            corrs.append(np.abs(np.corrcoef(U_big.ix[:,j],U_big.ix[:,i+cols1])[0][1]))
        comparevec.append(np.max(corrs))

    return comparevec


def PC_noise(DGEZ,noiselevels=np.linspace(-2,2,20),reps=3,sig_pcs=40):

    PC_cor=pd.DataFrame()
    [Ufb,Sfb,Vfb]= fb_pca(DGEZ,k=sig_pcs)

    for noise in noiselevels:

        df_noise=pd.DataFrame()

        for rep in range(reps):

            DGE_Z_wnoise=DGEZ+np.random.normal(0,np.power(10.0,noise),np.shape(DGEZ))

            [Ufb_noise,Sfb_noise,Vfb_noise]=fb_pca(DGE_Z_wnoise,k=sig_pcs)
            comparevec=[]
            for i in range(sig_pcs):
                corrs=[]
                for j in range(sig_pcs):
                    corrs.append(np.abs(np.corrcoef(Ufb.ix[:,j],Ufb_noise.ix[:,i])[0][1]))
                comparevec.append(np.max(corrs))

            df_noise[rep]=comparevec
        PC_cor[noise]=df_noise.mean(axis=1)

    return PC_cor

def jackstraw(DGEZ,per=0.005,sig_pcs=40,reps=100,verbose=0):
    """substitute small percentage of features with permuted versions, compare actual to permuted to obtain significance"""
    ngenes=len(DGEZ)
    [Ufb,Sfb,Vfb]= fb_pca(DGEZ,k=sig_pcs)

    Ufb_null = pd.DataFrame()
    flag=0
    #repeatedly permute and recalculate null PC distributions
    for i in range(reps):
        if (verbose==1):
            print('rep',i)
        shuf_genes=np.random.choice(range(ngenes),size=int(np.ceil(ngenes*per)),replace=False)

        DGEZ_perm=DGEZ.copy()
        DGEZ_perm.ix[shuf_genes,:]=np.array(DGEZ_perm.ix[shuf_genes,np.random.permutation(range(np.shape(DGEZ)[1]))])

        [Ufb_perm,Sfb_perm,Vfb_perm]= fb_pca(DGEZ,k=sig_pcs)
        tmp_null=Ufb.ix[shuf_genes,:]
        if flag==0:
            Ufb_null=tmp_null
            flag=1
        else:
            Ufb_null=pd.concat([Ufb_null,tmp_null])

    PVALS=Ufb.copy()
    for i in range(sig_pcs):
        curecdf=ECDF(Ufb_null.ix[:,i])
        curUfb=Ufb.ix[:,i]

        isnegative=curUfb<0.0
        ispositive=curUfb>=0.0

#statsmodels.sandbox.stats.multicomp.fdrcorrection0
        PVALS.ix[isnegative,i]=np.log10(curecdf(Ufb.ix[isnegative,i]))
        PVALS.ix[ispositive,i]=-np.log10(1-curecdf(Ufb.ix[ispositive,i]))
        PVALS[PVALS>5]=5
        PVALS[PVALS<(-5)]=-5

    return PVALS

##########significance testing##############

def ttest_DGE(DGE1,DGE2,batch=None):
    """Ttest with batchwise comparison capability"""
    FC=[]
    PVALS=[]
    
    A=DGE1.T 
    B=DGE2.T
    if batch is None:
        for gene in A:
            difmean=np.mean(A[gene])-np.mean(B[gene])
            ttest=scipy.stats.ttest_ind(A[gene],B[gene],equal_var=False)
        
            FC.append(difmean)
            PVALS.append(np.sign(difmean)*(-np.log10(ttest[1])))
    else:
        batchun=np.unique(batch[0])
        for gene in A:
            difmean=[]
            pvals=[]
            curA=A[gene]
            curB=B[gene]
            for curbatch in batchun:
                curAbatch=curA[batch[0]==curbatch]
                curBbatch=curB[batch[1]==curbatch]
                curdiff=np.mean(curAbatch)-np.mean(curBbatch)
                difmean.append(curdiff)
                ttest=scipy.stats.ttest_ind(curAbatch,curBbatch,equal_var=False)
                pvals.append(np.sign(curdiff)*(-np.log10(ttest[1])))
            FC.append(np.mean(difmean))
            PVALS.append(np.sum(pvals))

    return np.array(FC),np.array(PVALS)

def within_across(B1,B2):
    """Compare correlations between two sets of regression coefficients"""
    D1=B1.copy()
    D2=B2.copy()

    D1=D1.groupby(D1.index).mean()
    D2=D2.groupby(D2.index).mean()

    genesin=list(set(D1.index).intersection(set(D2.index)))

    D1=D1.loc[genesin]
    D2=D2.loc[genesin]

    cor_list=[]
    compare_list=[]

    for col1 in D1.columns:
        for col2 in D2.columns:

            if col1 in D2.columns:
                cur1=D1[col1]
                cur2=D2[col2]

                gene1=guide2gene(col1)
                gene2=guide2gene(col2)

                curcor=np.corrcoef(cur1,cur2)[0][1]
                cor_list.append(curcor)

                if col1==col2:
                    compare_list.append('same guide')
                elif gene1==gene2:
                    compare_list.append('within')
                else:
                    compare_list.append('across')

    df_compare=pd.DataFrame(compare_list)
    df_compare['correlations']=np.array(cor_list)
    return df_compare

def within_without(B):
    """Compare correlations between perturbations targeting the same gene"""
    COR=pd.DataFrame(np.corrcoef(B.T))
    COR.index=B.columns
    COR.columns=B.columns

    count1=0
    samesies=[]
    difsies =[]
    intsies = []
    for guide1 in COR.columns:
        count2=0
        for guide2 in COR.columns:
            if (count1>count2)&(('NTC' not in guide1)&('NTC' not in guide2))&(('INTERGENIC' not in guide1)&('INTERGENIC' not in guide2)):
                underscore_pos = []
                for x in re.finditer('_',guide1):
                    underscore_pos.append(x.span()[1])
                thisgene1=guide1[underscore_pos[0]:underscore_pos[1]-1]
                underscore_pos = []
                for x in re.finditer('_',guide2):
                    underscore_pos.append(x.span()[1])
                thisgene2=guide2[underscore_pos[0]:underscore_pos[1]-1]
                
                if thisgene1==thisgene2:
                    samesies.append(COR[guide1][guide2])
                else:
                    difsies.append(COR[guide1][guide2])
            elif (count1>count2):
                intsies.append(COR[guide1][guide2])
            count2+=1
        count1+=1
    return np.array(samesies),np.array(difsies),np.array(intsies)


#####GO analysis################
def PCA2GO(DGEZ,sigpcs=15,thresh=2,repin=100,perin=0.005,fdr_thresh=0.1,species='human'):
    
    #Jackstraw
    PVALS=jackstraw(DGEZ,sig_pcs=sigpcs,reps=repin,per=perin)
    print('done with jackstraw')

    #go analysis
    path2db='Path2obo/db/'
    obodag = GODag(path2db+"go-basic.obo")

    if (species=='human'):

        geneid2gos = read_ncbi_gene2go(path2db+"gene2go", taxids=[9606])
        print("{N:,} annotated genes".format(N=len(geneid2gos)))
        these_genes = DGEZ.index
        Xtable=pd.read_csv('Path2Xref/hg19_xref.txt',sep='\t')
        Xtable.index=Xtable['Approved Symbol']
        entrez=[int(x) for x in np.unique(Xtable.loc[these_genes].dropna()['EntrezGene ID'])]

    elif(species=='mouse'):

        geneid2gos = read_ncbi_gene2go(path2db+"gene2go", taxids=[10090])
        print("{N:,} annotated genes".format(N=len(geneid2gos)))
        these_genes = DGEZ.index
        Xtable=pd.read_csv('Path2xref/biomart_xref.mm10.txt',sep='\t')
        Xtable=Xtable[['Associated Gene Name','EntrezGene ID']].dropna()
        Xtable.index=Xtable['Associated Gene Name']
        entrez=[int(x) for x in np.unique(Xtable.loc[these_genes].dropna()['EntrezGene ID'])]

    goeaobj = GOEnrichmentStudy(
            entrez, # List of mouse protein-coding genes
            geneid2gos, # geneid/GO associations
            obodag, # Ontologies
            propagate_counts = False,
            alpha = 0.05, # default significance cut-off
            methods = ['fdr_bh']) # defult multipletest correction method

    df_bigGO=pd.DataFrame()
    count=0
    ngenes=len(PVALS)

    for pc in PVALS.columns:
        print(pc)
        df_GO=pd.DataFrame()
        

        curU=PVALS[pc]
        meanx=np.mean(curU)
        stdx=np.std(curU)
        threshlow=meanx-thresh*stdx
        threshhigh=meanx+thresh*stdx

        lookup_entrez_high=[int(x) for x in np.unique(Xtable.loc[curU[curU>threshhigh].index].dropna()['EntrezGene ID'])]
        lookup_entrez_low=[int(x) for x in np.unique(Xtable.loc[curU[curU<threshlow].index].dropna()['EntrezGene ID'])]
        # 'p_' means "pvalue". 'fdr_bh' is the multipletest method we are currently using.
        goea_results_high = goeaobj.run_study(lookup_entrez_high)
        indexlist=[]
        if len(lookup_entrez_high)>0:
            for i in range(len(goea_results_high)):
                if goea_results_high[i].p_fdr_bh<fdr_thresh:
                    df_GO[i]=[-np.log10(goea_results_high[i].p_fdr_bh)]
                    indexlist.append(goea_results_high[i].name)
        if len(lookup_entrez_low)>0:
            goea_results_low = goeaobj.run_study(lookup_entrez_low)
            for j in range(len(goea_results_high)):
                if goea_results_low[j].p_fdr_bh<fdr_thresh:
                    df_GO[j+len(goea_results_high)+1]=[np.log10(goea_results_low[j].p_fdr_bh)]
                    indexlist.append(goea_results_low[j].name)
        if(np.shape(df_GO)[0]==0):
            df_GO[0]=[0]
            df_GO.index=['NoGO']
        else:
            df_GO=df_GO.T
            df_GO.index=indexlist
        df_GO.columns=[pc]
        df_GO=df_GO.groupby(df_GO.index).first()
        if count==0:
            df_bigGO=df_GO
            count=1
        else:
            df_bigGO=df_bigGO.merge(df_GO,how='outer',left_index=True,right_index=True)
    df_bigGO=df_bigGO.fillna(0)
    if 'NoGO' in df_bigGO.index:
        df_bigGO=df_bigGO.drop('NoGO')
    return df_bigGO
    
 
def DE2GO(df_p,background,sig_thresh=3,num_genes=None,fdr_thresh=0.1,species='human'):
    #go analysis
    path2db='PATH2GOobofile/db/'
    obodag = GODag(path2db+"go-basic.obo")

    if (species=='human'):

        geneid2gos = read_ncbi_gene2go(path2db+"gene2go", taxids=[9606])
        print("{N:,} annotated genes".format(N=len(geneid2gos)))
        these_genes = background
        Xtable=pd.read_csv('PATH2hg19Xref/hg19_xref.txt',sep='\t')
        Xtable.index=Xtable['Approved Symbol']
        entrez=[int(x) for x in np.unique(Xtable.loc[these_genes].dropna()['EntrezGene ID'])]

    elif(species=='mouse'):

        geneid2gos = read_ncbi_gene2go(path2db+"gene2go", taxids=[10090])
        print("{N:,} annotated genes".format(N=len(geneid2gos)))
        these_genes = background
        Xtable=pd.read_csv('PATH2mm10xref/biomart_xref.mm10.txt',sep='\t')
        Xtable=Xtable[['Associated Gene Name','EntrezGene ID']].dropna()
        Xtable.index=Xtable['Associated Gene Name']
        entrez=[int(x) for x in np.unique(Xtable.loc[these_genes].dropna()['EntrezGene ID'])]

    goeaobj = GOEnrichmentStudy(
            entrez, # List of mouse protein-coding genes
            geneid2gos, # geneid/GO associations
            obodag, # Ontologies
            propagate_counts = False,
            alpha = 0.05, # default significance cut-off
            methods = ['fdr_bh']) # defult multipletest correction method

    df_bigGO=pd.DataFrame()
    count=0
    ngenes=len(df_p)
    xtable_genes=set(Xtable.index)

    for cluster in df_p.columns:
        print(cluster)
        df_GO=pd.DataFrame()
        

        cur_cluster=df_p[cluster]

        threshlow=(-sig_thresh)#np.percentile(cur_cluster,100.0*(1-np.divide(num_genes,ngenes)))
        threshhigh=sig_thresh#np.percentile(cur_cluster,100.0*(np.divide(num_genes,ngenes)))

        genes_high=cur_cluster[cur_cluster>threshhigh].index
        genes_high=list(xtable_genes.intersection(set(genes_high)))
        genes_low=cur_cluster[cur_cluster<threshlow].index
        genes_low=list(xtable_genes.intersection(set(genes_low)))

        lookup_entrez_high=[int(x) for x in np.unique(Xtable.loc[genes_high].dropna()['EntrezGene ID'])]
        lookup_entrez_low=[int(x) for x in np.unique(Xtable.loc[genes_low].dropna()['EntrezGene ID'])]
        # 'p_' means "pvalue". 'fdr_bh' is the multipletest method we are currently using.
        indexlist=[]
        if len(lookup_entrez_high)>2:
            goea_results_high = goeaobj.run_study(lookup_entrez_high)
            for i in range(len(goea_results_high)):
                if goea_results_high[i].p_fdr_bh<fdr_thresh:
                    df_GO[i]=[-np.log10(goea_results_high[i].p_fdr_bh)]
                    indexlist.append(goea_results_high[i].name)
            highlen=len(goea_results_high)
        else:
            highlen=0

        if len(lookup_entrez_low)>2:
            goea_results_low = goeaobj.run_study(lookup_entrez_low)
            for j in range(len(goea_results_low)):
                if goea_results_low[j].p_fdr_bh<fdr_thresh:
                    df_GO[j+highlen+1]=[np.log10(goea_results_low[j].p_fdr_bh)]
                    indexlist.append(goea_results_low[j].name)
        if(np.shape(df_GO)[0]==0):
            df_GO[0]=[0]
            df_GO.index=['NoGO']
        else:
            df_GO=df_GO.T
            df_GO.index=indexlist
        df_GO.columns=[cluster]
        df_GO=df_GO.groupby(df_GO.index).first()
        if count==0:
            df_bigGO=df_GO
            count=1
        else:
            df_bigGO=df_bigGO.merge(df_GO,how='outer',left_index=True,right_index=True)
    df_bigGO=df_bigGO.fillna(0)
    if 'NoGO' in df_bigGO.index:
        df_bigGO=df_bigGO.drop('NoGO')
    return df_bigGO


def TOPN2GO(df_p,background,num_genes=100,fdr_thresh=0.1,species='human'):
    #go analysis
    path2db='Path2obo/db/'
    obodag = GODag(path2db+"go-basic.obo")

    if (species=='human'):

        geneid2gos = read_ncbi_gene2go(path2db+"gene2go", taxids=[9606])
        print("{N:,} annotated genes".format(N=len(geneid2gos)))
        these_genes = background
        Xtable=pd.read_csv('PATH2xfref/hg19_xref.txt',sep='\t')
        Xtable.index=Xtable['Approved Symbol']
        entrez=[int(x) for x in np.unique(Xtable.loc[these_genes].dropna()['EntrezGene ID'])]

    elif(species=='mouse'):

        geneid2gos = read_ncbi_gene2go(path2db+"gene2go", taxids=[10090])
        print("{N:,} annotated genes".format(N=len(geneid2gos)))
        these_genes = background
        Xtable=pd.read_csv('PATH2xref/biomart_xref.mm10.txt',sep='\t')
        Xtable=Xtable[['Associated Gene Name','EntrezGene ID']].dropna()
        Xtable.index=Xtable['Associated Gene Name']
        entrez=[int(x) for x in np.unique(Xtable.loc[these_genes].dropna()['EntrezGene ID'])]

    goeaobj = GOEnrichmentStudy(
            entrez, # List of mouse protein-coding genes
            geneid2gos, # geneid/GO associations
            obodag, # Ontologies
            propagate_counts = False,
            alpha = 0.05, # default significance cut-off
            methods = ['fdr_bh']) # defult multipletest correction method

    df_bigGO=pd.DataFrame()
    count=0
    ngenes=len(df_p)
    xtable_genes=set(Xtable.index)

    for cluster in df_p.columns:
        print(cluster)
        df_GO=pd.DataFrame()
        

        cur_cluster=df_p[cluster]

        threshlow=np.percentile(cur_cluster,100.0*(1-np.divide(num_genes,ngenes)))
        threshhigh=np.percentile(cur_cluster,100.0*(np.divide(num_genes,ngenes)))

        genes_high=cur_cluster[cur_cluster>threshhigh].index
        genes_high=list(xtable_genes.intersection(set(genes_high)))
        genes_low=cur_cluster[cur_cluster<threshlow].index
        genes_low=list(xtable_genes.intersection(set(genes_low)))

        lookup_entrez_high=[int(x) for x in np.unique(Xtable.loc[genes_high].dropna()['EntrezGene ID'])]
        lookup_entrez_low=[int(x) for x in np.unique(Xtable.loc[genes_low].dropna()['EntrezGene ID'])]
        # 'p_' means "pvalue". 'fdr_bh' is the multipletest method we are currently using.
        indexlist=[]
        if len(lookup_entrez_high)>2:
            goea_results_high = goeaobj.run_study(lookup_entrez_high)
            for i in range(len(goea_results_high)):
                if goea_results_high[i].p_fdr_bh<fdr_thresh:
                    df_GO[i]=[-np.log10(goea_results_high[i].p_fdr_bh)]
                    indexlist.append(goea_results_high[i].name)
            highlen=len(goea_results_high)
        else:
            highlen=0

        if len(lookup_entrez_low)>2:
            goea_results_low = goeaobj.run_study(lookup_entrez_low)
            for j in range(len(goea_results_low)):
                if goea_results_low[j].p_fdr_bh<fdr_thresh:
                    df_GO[j+highlen+1]=[np.log10(goea_results_low[j].p_fdr_bh)]
                    indexlist.append(goea_results_low[j].name)
        if(np.shape(df_GO)[0]==0):
            df_GO[0]=[0]
            df_GO.index=['NoGO']
        else:
            df_GO=df_GO.T
            df_GO.index=indexlist
        df_GO.columns=[cluster]
        df_GO=df_GO.groupby(df_GO.index).first()
        if count==0:
            df_bigGO=df_GO
            count=1
        else:
            df_bigGO=df_bigGO.merge(df_GO,how='outer',left_index=True,right_index=True)
    df_bigGO=df_bigGO.fillna(0)
    if 'NoGO' in df_bigGO.index:
        df_bigGO=df_bigGO.drop('NoGO')
    return df_bigGO



def cluster_merger(DGE,cbcs,cluster_labs):
    #differential expression per cluster, each cluster against all others
    unique_labs=np.unique(cluster_labs)
    inputgenes=np.log10(len(DGE))
    numclust=len(unique_labs)

    df_numDE = pd.DataFrame(np.zeros((numclust,numclust)))

    for i in range(numclust):
        print(i)
        for j in range(numclust):
            if (i>j):
                cells1=cbcs[cluster_labs==unique_labs[i]]
                cells2=cbcs[cluster_labs==unique_labs[j]]

                fc,sl10pval=ttest_DGE(DGE[cells1],DGE[cells2])

                numDE=np.sum(np.abs(sl10pval)>inputgenes)
                df_numDE.ix[i,j]=numDE
                df_numDE.ix[j,i]=numDE

    return df_numDE
       
##Random string manipulations

# def edit_dist(barcodes1,barcodes2):
#     dist = [[Levenshtein.distance(a,b) for a in barcodes1] for b in barcodes2]
#     return np.array(dist)

# def edit_dist_closest(barcodes):
# 	DISTMAT=pd.DataFrame(np.zeros((len(barcodes),len(barcodes)))+20)
# 	for i in range(len(barcodes)):
# 		for j in range(len(barcodes)):
# 			if(i>j):
# 				DISTMAT.ix[i,j] = Levenshtein.distance(barcodes[i],barcodes[j])
# 				DISTMAT.ix[j,i] = DISTMAT.ix[i,j]
# 	return DISTMAT

#alt_map = {'ins':'0'}
#complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 

def reverse_complement(seq):    
    for k,v in alt_map.iteritems():
        seq = seq.replace(k,v)
    bases = list(seq) 
    bases = reversed([complement.get(base,base) for base in bases])
    bases = ''.join(bases)
    for k,v in alt_map.iteritems():
        bases = bases.replace(v,k)
    return bases

#######Plotting accessory function###########

#def plot_stacked_bar_chart(df_vals):

def number2color(vec,minin=None,maxin=None,cmap=None):

    if minin:
        amin=minin
    else:
        amin=np.min(vec)
    if maxin:
        amax=maxin
    else:
        amax=np.max(vec)

    import matplotlib as mpl
    import matplotlib.cm as cm
    norm=mpl.colors.Normalize(vmin=amin,vmax=amax)
    if not cmap:
        cmap=cm.viridis
    m=cm.ScalarMappable(norm=norm,cmap=cmap)
    colors_out=[m.to_rgba(x) for x in vec]

    return colors_out

def factor2color(vec,n=None):

    vec=np.array(vec)
    types=np.unique(vec)
    if n:
        n=len(types)
        max_value = 16581375 #255**3
        interval = int(max_value / n)
        colors = [hex(I)[2:].zfill(6) for I in range(0, max_value, interval)]
        cur_colors=[(int(i[:2], 16), int(i[2:4], 16), int(i[4:], 16)) for i in colors]
        cur_colors=[np.divide(1.0*x,255.0) for x in np.array(cur_colors)]
    else:
        cur_colors=sns.get_color_cycle()
    colors=pd.DataFrame([[0,0,0]]*len(vec))
    i=0
    for category in types:
        colors.ix[vec==category,:]=np.array(cur_colors[i])
        i+=1
    return np.array(colors)


#    means=df_vals.mean()
#    stds=df_vals.std()
def self_targeting(DE_COEFS):
    #get selftargeting scores
    selftarget=[]
    ourgenetargets = []

    for guide in DE_COEFS:
        curFC=DE_COEFS[guide]
        
        underscore_pos = []
        count=0
        if ('INTERGENIC' in guide)|('NTC' in guide):
            selftarget.append(0.0)
        else:
            thisgene=guide2gene(guide)
            if thisgene in DE_COEFS.index:
                selftarget.append(curFC[thisgene])
                ourgenetargets.append(thisgene)
            else:
                selftarget.append(0.0)

    #make plotting targeting colors
    import matplotlib as mpl
    import matplotlib.cm as cm

    norm = mpl.colors.Normalize(vmin=-0.3, vmax=0.3)
    cmap = cm.seismic
    m = cm.ScalarMappable(norm=norm, cmap=cmap)
    targeting_colors=[m.to_rgba(x) for x in selftarget]
    return selftarget,targeting_colors

def density_scatter(x,y,savepath,sample=None,pointsize=10,colmap='magma'):
    
    if sample:
        BIG=pd.DataFrame()
        BIG['x']=x
        BIG['y']=y
        BIG_sample=BIG.sample(sample)
        x=BIG_sample['x'].copy()
        y=BIG_sample['y'].copy()

    # Calculate the point density
    xy = np.vstack([x,y])
    z = gaussian_kde(xy)(xy)

    # Sort the points by density, so that the densest points are plotted last
    idx = z.argsort()
    x, y, z = x[idx], y[idx], z[idx]

    plt.scatter(x, y, c=z, s=pointsize, edgecolor='',cmap=colmap)
    #plt.axhline(0,c='red')
    plt.savefig(savepath,transparent=True)
    return 

def model_plot_orenplots(B,column_subset,pathout):
    DE_COEFS=B[column_subset]
    COR_COEF=pd.DataFrame(np.corrcoef(DE_COEFS.T))
    COR_COEF.index=DE_COEFS.columns
    COR_COEF.columns=DE_COEFS.columns

    selftarget,targeting_colors=self_targeting(DE_COEFS)

    cg=sns.clustermap(COR_COEF,cmap='bwr',figsize=(40,40),col_colors=targeting_colors,row_colors=targeting_colors)
    plt.setp(cg.ax_heatmap.yaxis.get_majorticklabels(),rotation=0)
    plt.savefig(pathout+'_orenplot.pdf')
    plt.clf()

    COR_on=COR_COEF[np.array(selftarget)<0.05].T[np.array(selftarget)<0.05]

    cg=sns.clustermap(COR_on,cmap='bwr',col_colors=np.array(targeting_colors)[np.array(selftarget)<0.05],row_colors=np.array(targeting_colors)[np.array(selftarget)<0.05],figsize=(20,20))
    plt.setp(cg.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
    plt.savefig(pathout+'_orenplot_nored.pdf')
    plt.clf()

def model_dump_results(B,Bshuf,Bfdr,column_subset,pathout,myspecies='human'):

    B.to_csv(pathout+'_LM_COEFS.csv')
    Bshuf.to_csv(pathout+'_LM_COEFS_SHUF.csv')
    Bfdr.to_csv(pathout+'_LM_FDR.csv')

    model_plot_orenplots(B,column_subset,pathout)
    
    df_coefGO=DE2GO(Bfdr,background=Bfdr.index,sig_thresh=1.3,species=myspecies)
    nnz_go=np.sum(df_coefGO==0,axis=1)
    max_go=np.max(np.abs(df_coefGO),axis=1)
    GO_plot=df_coefGO[np.logical_and(nnz_go<len(column_subset)-3,max_go>2)]
    print(np.shape(GO_plot))

    cg=sns.clustermap(GO_plot,figsize=(80,80))#,metric='correlation')
    plt.setp(cg.ax_heatmap.yaxis.get_majorticklabels(),rotation=0)
    plt.savefig(pathout+'_go_fdr.pdf')
    plt.clf()

    GO_plot.to_csv(pathout+'_GOenrichments.csv')
    GO_plot_gene=GO_plot.copy().T
    GO_plot_gene=GO_plot_gene.groupby([guide2gene(x) for x in GO_plot_gene.index]).mean().T
    GO_plot_gene_f=GO_plot_gene[np.abs(GO_plot_gene).max(axis=1)>1.3]
    cg=sns.clustermap(GO_plot_gene_f,figsize=(80,80))
    plt.setp(cg.ax_heatmap.yaxis.get_majorticklabels(),rotation=0)
    plt.savefig(pathout+'_go_fdr_gene.pdf')
    plt.clf()
    GO_plot_gene_f.to_csv(pathout+'_go_fdr_gene.csv')
