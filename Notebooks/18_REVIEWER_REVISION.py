"""Reviewer-requested correction and sensitivity analyses; see analysis_protocol.txt.

Run with the recorded Python 3.12 environment. No installed package is modified.
Outputs are isolated from the submitted analysis. Stages: cached, infer, analyze.
"""
from pathlib import Path
import os, sys, json, argparse, importlib, importlib.util, hashlib, warnings
BASE=Path(__file__).resolve().parents[1]
LOCAL=BASE/'_local_submission/reviewer_revision'
OUT=BASE/'Results/reviewer_revision'
for p in [LOCAL,OUT,LOCAL/'numba',LOCAL/'mpl',LOCAL/'scores']: p.mkdir(parents=True,exist_ok=True)
os.environ.update(PYTHONDONTWRITEBYTECODE='1',NUMBA_CACHE_DIR=str(LOCAL/'numba'),MPLCONFIGDIR=str(LOCAL/'mpl'),OMP_NUM_THREADS='4',OPENBLAS_NUM_THREADS='4',NUMBA_NUM_THREADS='4')
import numpy as np
import pandas as pd
from scipy.stats import rankdata
KEYS=['source','target','ligand_complex','receptor_complex']
GROUPS=['cHD','cSLE','aHD','aSLE']
EXAMPLES=[('Classical Monocytes','NK Cells','TIMP2','ITGB1'),('CD4 T Cells','NK Cells','CD48','CD244'),('Non-classical Monocytes','CD4 T Cells','SELPLG','ITGB2'),('Classical Monocytes','NK Cells','TGFB1','ITGB1'),('Non-classical Monocytes','CD4 T Cells','CD48','CD2'),('B Cells','Non-classical Monocytes','CD52','SIGLEC10'),('CD8 T Cells','Non-classical Monocytes','CD52','SIGLEC10')]
IDS=['|'.join(k) for k in EXAMPLES]
def module(name,file):
    spec=importlib.util.spec_from_file_location(name,BASE/'Notebooks'/file)
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod);return mod
fit=module('revision_fit','15_SUMMARIZE_ANNOTATION_RESOURCE_SENSITIVITY.py').fit
def corrected_rank(frame,specs,aggregate_method):
    """Rank each distinct score once, retaining the original unique-score RRA."""
    agg=importlib.import_module('liana.method._pipe_utils._aggregate')
    unique={}
    for score,ascending in specs.values():
        if score in unique and unique[score]!=ascending: raise ValueError('Inconsistent rank directions')
        unique[score]=ascending
    mat=np.column_stack([rankdata(frame[s].to_numpy()*(1 if ascending else -1),method='average') for s,ascending in unique.items()])
    if aggregate_method=='rra': return agg._robust_rank_aggregate(mat)
    if aggregate_method=='mean': return mat.mean(axis=1)/len(mat)
    raise ValueError(aggregate_method)
def patch():
    import liana as li
    assert li.__version__=='1.8.1',li.__version__
    agg=importlib.import_module('liana.method._pipe_utils._aggregate')
    original=agg._rank_aggregate
    # Regression check: high expr_prod must get the best rank even with duplicate specs.
    toy=pd.DataFrame({'expr_prod':[1.,2.,3.]})
    specs={'Connectome':('expr_prod',False),'NATMI':('expr_prod',False)}
    assert original(toy.copy(),specs,'mean')[0]<original(toy.copy(),specs,'mean')[2]
    result=corrected_rank(toy,specs,'mean')
    assert result[2]<result[0]
    agg._rank_aggregate=corrected_rank
    return li,original
def metadata(): return pd.read_csv(BASE/'Results/severity_analysis/patient_clinical_metadata.csv')
def identify(d):
    d['interaction_id']=d[KEYS].astype(str).agg('|'.join,axis=1);return d
def load_scores(state):
    return pd.concat([pd.read_csv(p) for p in sorted((LOCAL/'scores'/state).glob('*.csv.gz'))],ignore_index=True).merge(metadata(),on='sample',validate='many_to_one')
def cached():
    li,oldrank=patch();checks=[]
    import importlib.metadata as im
    versions={p:im.version(p) for p in ['liana','scanpy','anndata','numpy','pandas','scipy','statsmodels','celltypist','harmonypy','scrublet','leidenalg','igraph','umap-learn']}
    versions['python']=sys.version
    (OUT/'software_versions.json').write_text(json.dumps(versions,indent=2))
    for state in ['baseline','annotation','resource','combined']:
        dest=LOCAL/'scores'/state;dest.mkdir(exist_ok=True)
        for p in sorted((BASE/'_local_submission/sensitivity_cache/scores'/state).glob('*.csv.gz')):
            d=pd.read_csv(p)
            old=oldrank(d.copy(),li.mt.rank_aggregate.magnitude_specs,'rra')
            assert np.allclose(old,d.magnitude_rank,atol=1e-12),p
            corrected=corrected_rank(d.copy(),li.mt.rank_aggregate.magnitude_specs,'rra')
            checks.append({'scenario':state,'sample':d['sample'].iloc[0],'rows':len(d),'pearson_r':pd.Series(old).corr(pd.Series(corrected)),'spearman_r':pd.Series(old).corr(pd.Series(corrected),method='spearman'),'max_abs_change':float(np.max(np.abs(old-corrected)))})
            d['old_magnitude_rank']=d.magnitude_rank;d['magnitude_rank']=corrected;d['score']=1-corrected
            d.to_csv(dest/p.name,index=False,compression='gzip')
        s=load_scores(state);model=fit(s);model.to_csv(OUT/f'{state}_models.csv',index=False)
        print('cached',state,len(model),int((model.fdr_q_value<.05).sum()),flush=True)
    pd.DataFrame(checks).to_csv(OUT/'donor_score_correction_comparison.csv',index=False)
    # Freeze the exact installed resource, verify against the archived snapshot.
    for resource in ['consensus','cellphonedb']:
        current=li.rs.select_resource(resource)
        saved=pd.read_csv(BASE/f'Results/annotation_resource_sensitivity/{resource}_resource.csv')
        assert set(map(tuple,current[['ligand','receptor']].values))==set(map(tuple,saved[['ligand','receptor']].values))
        saved.to_csv(OUT/f'{resource}_resource.csv',index=False)
def infer(args):
    import h5py,anndata as ad
    li,_=patch();helper=module('revision_input','14_ANNOTATION_RESOURCE_SENSITIVITY.py')
    warnings.filterwarnings('ignore',category=FutureWarning)
    warnings.filterwarnings('ignore',category=UserWarning)
    resource=pd.read_csv(OUT/'consensus_resource.csv')
    tested=pd.read_csv(BASE/'Results/revised_primary/primary_models.csv')[['interaction_id']+KEYS]
    for cohort,filename in helper.FILES.items():
        if args.cohort and args.cohort!=cohort: continue
        review=pd.read_csv(BASE/f'_local_submission/sensitivity_cache/{cohort}_cell_review.csv.gz',index_col='cell_id',dtype={'cluster':str})
        with h5py.File(args.archive_root/f'{cohort}_individual_h5ad'/filename,'r') as f:
            obs,genes=helper.metadata(f,cohort);assert obs.index.equals(review.index)
            for sample in sorted(obs['sample'].unique()):
                if args.sample and args.sample!=sample: continue
                dest=LOCAL/'scores/return_all';dest.mkdir(exist_ok=True)
                checkfile=LOCAL/f'{sample}_verification.json'
                if checkfile.exists(): continue
                idx=np.flatnonzero(obs['sample'].to_numpy()==sample)
                donor=review.iloc[idx];label=donor.sensitivity_label
                keep=~label.isin(['Ribosomal/Low-quality','RBCs','Platelets','Mitochondrial-high/Low-quality','Mixed platelet-monocyte/Low-quality'])
                x=helper.donor_matrix(f,idx,len(genes))[keep.to_numpy()]
                labels=label[keep].astype(str)
                data=ad.AnnData(x.copy(),obs=pd.DataFrame({'cell_type':labels},index=labels.index),var=pd.DataFrame(index=genes))
                counts=labels.value_counts();fractions={}
                for ct in counts.index:
                    # Count in integers before division: sparse float32 means can
                    # round an exact 10% fraction just below the threshold.
                    ct_matrix=x[labels.eq(ct).to_numpy()]
                    fractions[ct]=pd.Series(np.asarray((ct_matrix>0).sum(axis=0)).ravel()/ct_matrix.shape[0],index=genes)
                audit=tested.copy();audit['sample']=sample
                for role,gene_col in [('source','ligand_complex'),('target','receptor_complex')]:
                    audit[role+'_cells']=audit[role].map(counts).fillna(0).astype(int)
                    audit[role+'_min_fraction']=[min([fractions.get(ct,pd.Series(dtype=float)).get(g,np.nan) for g in complex_.split('_')]) for ct,complex_ in zip(audit[role],audit[gene_col])]
                    audit[role+'_under_10_cells']=audit[role+'_cells']<10
                    audit[role+'_below_10pct']=audit[role+'_min_fraction']<.1
                    audit[role+'_gene_unavailable']=[not set(complex_.split('_')).issubset(set(genes)) for complex_ in audit[gene_col]]
                audit['population_eligible']=~(audit.source_under_10_cells|audit.target_under_10_cells)
                audit['expression_eligible']=~(audit.source_below_10pct|audit.target_below_10pct|audit.source_gene_unavailable|audit.target_gene_unavailable)
                # Both calls use the same unchanged resource, X, labels and thresholds.
                kwargs=dict(groupby='cell_type',resource=resource,expr_prop=.1,min_cells=10,use_raw=False,n_perms=None,seed=1337,inplace=False,verbose=False)
                regular=identify(li.mt.rank_aggregate(data,**kwargs))
                cached=pd.read_csv(LOCAL/f'scores/annotation/{sample}.csv.gz');cached=identify(cached)
                match=regular.merge(cached,on='interaction_id',suffixes=('_rerun','_cached'),validate='one_to_one')
                assert len(match)==len(regular)==len(cached)
                delta=float((match.magnitude_rank_rerun-match.magnitude_rank_cached).abs().max())
                # CSV round trips can split floating-point ties in cached method scores.
                # Fresh inference is authoritative; record and bound reconstruction differences.
                assert delta<1e-4,delta
                fresh=LOCAL/'scores/fresh_primary';fresh.mkdir(exist_ok=True)
                regular['sample']=sample;regular['score']=1-regular.magnitude_rank
                regular.to_csv(fresh/f'{sample}.csv.gz',index=False)
                audit['observed']=audit.interaction_id.isin(regular.interaction_id)
                assert (audit.observed==(audit.population_eligible & audit.expression_eligible)).all(),sample
                audit.to_csv(LOCAL/f'{sample}_availability.csv.gz',index=False)
                allres=identify(li.mt.rank_aggregate(data,return_all_lrs=True,**kwargs))
                allres['sample']=sample;allres['score']=1-allres.magnitude_rank
                allres.to_csv(dest/f'{sample}.csv.gz',index=False)
                checkfile.write_text(json.dumps({'sample':sample,'rows':len(regular),'return_all_rows':len(allres),'max_abs_difference':delta}))
                print('inference verified',sample,len(regular),len(allres),flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['cached','infer','analyze']);p.add_argument('--archive-root',type=Path);p.add_argument('--cohort');p.add_argument('--sample');args=p.parse_args()
    if args.stage=='cached': cached()
    elif args.stage=='infer': infer(args)
    else: raise NotImplementedError('Analysis stage follows completed inference')
