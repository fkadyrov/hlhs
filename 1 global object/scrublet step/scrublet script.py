import numpy as np
import pandas as pd
import scanpy as sc
import scrublet as scr
import seaborn as sns
import scipy.stats
import anndata
import os
import igraph as ig
import louvain

adata = sc.read_h5ad("./merged_postQC.h5ad")
adata

def bh(pvalues):
    '''
    Computes the Benjamini-Hochberg FDR correction.

    Input:
        * pvals - vector of p-values to correct
    '''
    n = int(pvalues.shape[0])
    new_pvalues = np.empty(n)
    values = [ (pvalue, i) for i, pvalue in enumerate(pvalues) ]
    values.sort()
    values.reverse()
    new_values = []
    for i, vals in enumerate(values):
        rank = n - i
        pvalue, index = vals
        new_values.append((n/rank) * pvalue)
    for i in range(0, int(n)-1):
        if new_values[i] < new_values[i+1]:
            new_values[i+1] = new_values[i]
    for i, vals in enumerate(values):
        pvalue, index = vals
        new_pvalues[index] = new_values[i]
    return new_pvalues


####### main
meta = adata.obs
sc.settings.verbosity = 1  # verbosity: errors (0), warnings (1), info (2), hints (3)

scorenames = ['scrublet_score','scrublet_cluster_score','bh_pval']
os.makedirs('scrublet-scores')

###

for sample in meta['sample'].unique():
    print(sample)
    #import data
    adata_sample = adata[adata.obs['sample'] == sample]
    
    #set up and run Scrublet
    scrub = scr.Scrublet(adata_sample.X)
    doublet_scores, predicted_doublets = scrub.scrub_doublets(verbose=False)
    adata_sample.obs['scrublet_score'] = doublet_scores
    #overcluster prep. run turbo basic scanpy pipeline
    sc.pp.filter_genes(adata_sample, min_cells=3)
    sc.pp.normalize_per_cell(adata_sample, counts_per_cell_after=1e4)
    sc.pp.log1p(adata_sample)
    sc.pp.highly_variable_genes(adata_sample, min_mean=0.0125, max_mean=3, min_disp=0.5)
    adata_sample = adata_sample[:, adata_sample.var['highly_variable']]
    sc.pp.scale(adata_sample, max_value=10)
    sc.tl.pca(adata_sample, svd_solver='arpack')
    sc.pp.neighbors(adata_sample)
    #eoverclustering proper - do basic clustering first, then cluster each cluster
    sc.tl.louvain(adata_sample)
    for clus in np.unique(adata_sample.obs['louvain']):
        sc.tl.louvain(adata_sample, restrict_to=('louvain',[clus]))
        adata_sample.obs['louvain'] = adata_sample.obs['louvain_R']
    #compute the cluster scores - the median of Scrublet scores per overclustered cluster
    for clus in np.unique(adata_sample.obs['louvain']):
        adata_sample.obs.loc[adata_sample.obs['louvain']==clus, 'scrublet_cluster_score'] = \
            np.median(adata_sample.obs.loc[adata_sample.obs['louvain']==clus, 'scrublet_score'])
    #now compute doublet p-values. figure out the median and mad (from above-median values) for the distribution
    med = np.median(adata_sample.obs['scrublet_cluster_score'])
    mask = adata_sample.obs['scrublet_cluster_score']>med
    mad = np.median(adata_sample.obs['scrublet_cluster_score'][mask]-med)
    #let's do a one-sided test. the Bertie write-up does not address this but it makes sense
    pvals = 1-scipy.stats.norm.cdf(adata_sample.obs['scrublet_cluster_score'], loc=med, scale=1.4826*mad)
    adata_sample.obs['bh_pval'] = bh(pvals)
    #create results data frame for single sample and copy stuff over from the adata object
    scrublet_sample = pd.DataFrame(0, index=adata_sample.obs_names, columns=scorenames)
    for meta in scorenames:
        scrublet_sample[meta] = adata_sample.obs[meta]
    #write out complete sample scores
    scrublet_sample.to_csv('scrublet-scores/'+sample+'.csv')
    
    adata.obs['sample'].unique()
    
s1 = pd.read_csv("./scrublet-scores/P0116.csv", index_col=0)
s2 = pd.read_csv("./scrublet-scores/P0156.csv", index_col=0)
s3 = pd.read_csv("./scrublet-scores/P0172.csv", index_col=0)
s4 = pd.read_csv("./scrublet-scores/P0202.csv", index_col=0)
s5 = pd.read_csv("./scrublet-scores/P0219.csv", index_col=0)
s6 = pd.read_csv("./scrublet-scores/P0254.csv", index_col=0)
s7 = pd.read_csv("./scrublet-scores/P0318.csv", index_col=0)
s8 = pd.read_csv("./scrublet-scores/P0344.csv", index_col=0)
s9 = pd.read_csv("./scrublet-scores/P0371.csv", index_col=0)
s10 = pd.read_csv("./scrublet-scores/T485.csv", index_col=0)
s11 = pd.read_csv("./scrublet-scores/T1007.csv", index_col=0)
s12 = pd.read_csv("./scrublet-scores/T1093.csv", index_col=0)
s13 = pd.read_csv("./scrublet-scores/T1126.csv", index_col=0)
s14 = pd.read_csv("./scrublet-scores/T1197.csv", index_col=0)

df_list = [s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14]
scrub_all = pd.concat(df_list)

scrub_all.to_csv('scrublet-scores/all.csv')