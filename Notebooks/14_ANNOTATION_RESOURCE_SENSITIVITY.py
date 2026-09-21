"""Non-destructive annotation/resource sensitivity analysis; see output protocol.

Run with the original Python 3.12.10 analysis environment. External source
AnnData files and classifier are read-only; caches and donor outputs stay local.
"""
from __future__ import annotations
import os
from pathlib import Path
BASE = Path(__file__).resolve().parents[1]
LOCAL = BASE / '_local_submission' / 'sensitivity_cache'
for key, folder in [('NUMBA_CACHE_DIR','numba'), ('MPLCONFIGDIR','matplotlib'), ('TEMP','tmp'), ('TMP','tmp')]:
    (LOCAL / folder).mkdir(parents=True, exist_ok=True)
    os.environ[key] = str(LOCAL / folder)
os.environ['PYTHONDONTWRITEBYTECODE'] = '1'
os.environ['NUMBA_NUM_THREADS'] = '4'
os.environ['OMP_NUM_THREADS'] = '4'
os.environ['OPENBLAS_NUM_THREADS'] = '4'
import argparse
import hashlib
import importlib.util
import json
import time
import h5py
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix, vstack
OUT = BASE / 'Results' / 'annotation_resource_sensitivity'
OUT.mkdir(parents=True, exist_ok=True)
FILES = {'child':'adata_child_processed_leiden08.h5ad', 'adult':'adata_adult_processed_selected_leiden.h5ad'}
HELPER = {'Tcm/Naive helper T cells','Tem/Effector helper T cells','Tem/Effector helper T cells PD1+',
          'Follicular helper T cells','Regulatory T cells','Treg(diff)','Type 1 helper T cells',
          'Type 17 helper T cells','Memory CD4+ cytotoxic T cells'}
CYTOTOXIC = {'Tcm/Naive cytotoxic T cells','Tem/Temra cytotoxic T cells','Tem/Trm cytotoxic T cells',
             'Trm cytotoxic T cells','CD8a/b(entry)'}
STATES = {'baseline':(False,'consensus'), 'annotation':(True,'consensus'),
          'resource':(False,'cellphonedb'), 'combined':(True,'cellphonedb')}

def decode(x):
    if isinstance(x,h5py.Group): return pd.Categorical.from_codes(x['codes'][:], decode(x['categories']))
    return x.asstr()[:] if x.dtype.kind in ('S','O') else x[:]

def metadata(handle, cohort):
    o=handle['obs']; v=handle['var']
    obs=pd.DataFrame({k:decode(o[k]) for k in ['sample','leiden_res_0.8','n_genes_by_counts','total_counts','pct_counts_mt','doublet_score']},index=decode(o[o.attrs['_index']]))
    obs.index.name='cell_id'; obs['sample']=obs['sample'].astype(str)
    obs['cluster']=obs['leiden_res_0.8'].astype(str)
    mapping=pd.read_csv(BASE/'Results'/'qc'/f'{cohort}_cluster_to_celltype_map.csv',dtype={'cluster':str}).set_index('cluster').cell_type
    obs['original_label']=obs.cluster.map(mapping)
    assert obs.index.is_unique and obs.original_label.notna().all()
    genes=np.asarray(decode(v[v.attrs['_index']])).astype(str)
    return obs,genes

def slice_csr(node, start, end, n_genes):
    ptr=node['indptr'][start:end+1]; lo,hi=int(ptr[0]),int(ptr[-1])
    return csr_matrix((node['data'][lo:hi],node['indices'][lo:hi],ptr-ptr[0]),shape=(end-start,n_genes))

def donor_matrix(handle, indices, n_genes):
    # Preserve source order; rows within a donor are usually contiguous.
    splits=np.split(indices,np.flatnonzero(np.diff(indices)!=1)+1)
    return vstack([slice_csr(handle['X'],int(i[0]),int(i[-1])+1,n_genes) for i in splits],format='csr')

def resources():
    import liana as li
    path=Path(importlib.util.find_spec('liana.resource').origin).parent/'omni_resource.csv'
    table=pd.read_csv(path)
    table.loc[((table.source_genesymbol=='GNAI2') & (table.target_genesymbol=='S1PR4')) |
              ((table.source_genesymbol=='ARF6') & (table.target_genesymbol=='SMAP1'))].to_csv(OUT/'flagged_pair_resource_membership.csv',index=False)
    for name in ['consensus','cellphonedb']:
        li.resource.select_resource(name).to_csv(OUT/f'{name}_resource.csv',index=False)
    (OUT/'resource_provenance.json').write_text(json.dumps({'liana':li.__version__,'bundled_resource_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
      'alternative':'LIANA bundled cellphonedb resource; not the current external CellPhoneDB release',
      'interpretation':'Resource membership does not identify original experimental evidence or guarantee every interaction is extracellular.'},indent=2))

def annotate(args):
    import anndata as ad
    import celltypist
    model=celltypist.models.Model.load(model=str(args.model))
    resources()
    for cohort in FILES:
        target=LOCAL/f'{cohort}_cell_review.csv.gz'
        if target.exists(): print('Using completed annotation',cohort,flush=True); continue
        source=args.archive_root/f'{cohort}_individual_h5ad'/FILES[cohort]
        with h5py.File(source,'r') as f:
            obs,genes=metadata(f,cohort); lookup={g:i for i,g in enumerate(genes)}
            ribo=np.array([g.startswith(('RPL','RPS')) for g in genes]); chunks=[]
            for start in range(0,len(obs),4000):
                end=min(start+4000,len(obs)); x=slice_csr(f['X'],start,end,len(genes)); counts=slice_csr(f['layers']['counts'],start,end,len(genes))
                part=obs.iloc[start:end].copy()
                pred=celltypist.annotate(ad.AnnData(x,obs=pd.DataFrame(index=part.index),var=pd.DataFrame(index=genes)),model=model,majority_voting=False)
                part['reference_label']=pred.predicted_labels['predicted_labels'].astype(str).to_numpy()
                part['reference_probability']=pred.probability_matrix.max(axis=1).to_numpy()
                part['ribosomal_percent_retained_genes']=100*np.asarray(counts[:,ribo].sum(axis=1)).ravel()/np.maximum(np.asarray(counts.sum(axis=1)).ravel(),1)
                for gene in ['CD3D','CD3E','CD4','CD8A','CD8B','IL7R','CCR7','TCF7','NKG7']:
                    part[f'{gene}_detected']=np.asarray(x[:,lookup[gene]].toarray()).ravel()>0 if gene in lookup else False
                part['reference_broad']=part.reference_label.map({**{k:'CD4 T Cells' for k in HELPER},**{k:'CD8 T Cells' for k in CYTOTOXIC}}).fillna('unresolved')
                part['restore']=part.original_label.eq('Ribosomal/Low-quality') & part.reference_broad.ne('unresolved') & part.reference_probability.ge(.5) & part.CD3D_detected & part.CD3E_detected
                part['sensitivity_label']=np.where(part.restore,part.reference_broad,part.original_label)
                chunks.append(part)
                print(cohort,'annotated',end,'/',len(obs),flush=True)
        review=pd.concat(chunks); review.to_csv(target,compression='gzip')
        review.groupby(['cluster','original_label','reference_label'],observed=True).size().rename('n_cells').reset_index().to_csv(OUT/f'{cohort}_reference_labels_by_cluster.csv',index=False)
        review.groupby(['sample','cluster','original_label','sensitivity_label','restore'],observed=True).size().rename('n_cells').reset_index().to_csv(OUT/f'{cohort}_restoration_counts.csv',index=False)
        qc=['n_genes_by_counts','total_counts','pct_counts_mt','doublet_score','ribosomal_percent_retained_genes','reference_probability']
        review.groupby(['cluster','original_label'],observed=True)[qc].quantile([.1,.25,.5,.75,.9]).to_csv(OUT/f'{cohort}_cell_qc_quantiles.csv')
        review.loc[review.restore,['sample','cluster','reference_label','reference_probability','sensitivity_label']].to_csv(OUT/f'{cohort}_restored_cell_ids.csv.gz',compression='gzip')
        for name,mask in [('all',review.restore),('probability_at_least_0.8',review.restore & review.reference_probability.ge(.8))]:
            print(cohort,name,int(mask.sum()),flush=True)
    (OUT/'classifier_provenance.json').write_text(json.dumps({'model_name':args.model.name,'model_sha256':hashlib.sha256(args.model.read_bytes()).hexdigest(), 'majority_voting':False,'threshold':.5,'helper_labels':sorted(HELPER),'cytotoxic_labels':sorted(CYTOTOXIC)},indent=2))

def infer(args):
    import anndata as ad
    import liana as li
    for cohort in FILES:
        if args.cohort and cohort != args.cohort: continue
        review=pd.read_csv(LOCAL/f'{cohort}_cell_review.csv.gz',index_col='cell_id',dtype={'cluster':str})
        source=args.archive_root/f'{cohort}_individual_h5ad'/FILES[cohort]
        with h5py.File(source,'r') as f:
            obs,genes=metadata(f,cohort); assert obs.index.equals(review.index)
            for sample in sorted(obs['sample'].unique()):
                if args.sample and sample!=args.sample: continue
                idx=np.flatnonzero(obs['sample'].to_numpy()==sample); donor=review.iloc[idx]; x=None
                for state in args.states:
                    dest=LOCAL/args.score_subdir/state/f'{sample}.csv.gz'; dest.parent.mkdir(parents=True,exist_ok=True)
                    if dest.exists(): continue
                    if x is None: x=donor_matrix(f,idx,len(genes))
                    restore,resource=STATES[state]; label=donor.sensitivity_label if restore else donor.original_label
                    keep=~label.isin(['Ribosomal/Low-quality','RBCs','Platelets','Mitochondrial-high/Low-quality','Mixed platelet-monocyte/Low-quality'])
                    data=ad.AnnData(x[keep.to_numpy()].copy(),obs=pd.DataFrame({'cell_type':label.loc[keep].astype(str),'sample':sample},index=donor.index[keep]),var=pd.DataFrame(index=genes))
                    start=time.time()
                    res=li.mt.rank_aggregate(data,groupby='cell_type',resource_name=resource,expr_prop=.1,min_cells=10,use_raw=False,n_perms=None if args.magnitude_only else 1000,seed=1337,inplace=False,verbose=False)
                    res['sample']=sample
                    res.to_csv(dest,index=False,compression='gzip')
                    print(state,sample,len(data),len(res),round(time.time()-start,1),'seconds',flush=True)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--archive-root',type=Path,required=True)
    p.add_argument('--model',type=Path,default=Path.home()/'.celltypist/data/models/Immune_All_Low.pkl')
    p.add_argument('--phase',choices=['annotate','infer'],required=True)
    p.add_argument('--sample')
    p.add_argument('--cohort',choices=['child','adult'])
    p.add_argument('--magnitude-only',action='store_true',help='Skip unused permutation specificity scores; verify magnitude equivalence before using.')
    p.add_argument('--score-subdir',default='scores')
    p.add_argument('--states',nargs='+',choices=list(STATES),default=list(STATES))
    args=p.parse_args()
    if args.phase=='annotate':annotate(args)
    else:infer(args)

if __name__=='__main__': main()
