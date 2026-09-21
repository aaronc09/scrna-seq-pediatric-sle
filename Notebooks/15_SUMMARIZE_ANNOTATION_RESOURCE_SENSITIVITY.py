"""Fit and compare the four completed scenarios from script 14.

Uses the closed form of the saturated four-group OLS interaction and HC3 variance.
Baseline estimates, p values, and BH q values must reproduce the saved primary
analysis before sensitivity comparisons are exported.
"""
from pathlib import Path
import json
import os
import hashlib
from datetime import datetime, timezone
BASE=Path(__file__).resolve().parents[1]
LOCAL=BASE/'_local_submission'/'sensitivity_cache'
OUT=BASE/'Results'/'annotation_resource_sensitivity'
os.environ['MPLCONFIGDIR']=str(LOCAL/'matplotlib')
import numpy as np
import pandas as pd
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests
KEYS=['source','target','ligand_complex','receptor_complex']

def fit(scores):
    summaries=scores.groupby(KEYS+['Groups'],observed=True).score.agg(['count','mean','var']).reset_index()
    count=summaries.pivot(index=KEYS,columns='Groups',values='count').reindex(columns=['cHD','cSLE','aHD','aSLE'])
    mean=summaries.pivot(index=KEYS,columns='Groups',values='mean').reindex(columns=count.columns)
    var=summaries.pivot(index=KEYS,columns='Groups',values='var').reindex(columns=count.columns)
    valid=count.ge(5).all(axis=1) & ((var.fillna(0).sum(axis=1)>0) | (mean.max(axis=1)>mean.min(axis=1)))
    count,mean,var=count.loc[valid],mean.loc[valid],var.loc[valid]
    coef=mean.cSLE-mean.cHD-mean.aSLE+mean.aHD
    # h_ii=1/n_group, so group contribution to HC3 variance is SSE/(n-1)^2.
    se=np.sqrt((var/(count-1)).sum(axis=1))
    p=2*norm.sf(np.abs(coef/se))
    res=coef.rename('interaction_coef').to_frame()
    res['interaction_ci95_low']=coef-norm.ppf(.975)*se
    res['interaction_ci95_high']=coef+norm.ppf(.975)*se
    res['n_observations']=count.sum(axis=1).astype(int)
    res['interaction_p_value']=p
    if not np.isfinite(p).all(): raise ValueError('Degenerate model requires separate inspection')
    res['fdr_q_value']=multipletests(p,method='fdr_bh')[1]
    for group in count.columns:
        res[f'n_{group}']=count[group].astype(int);res[f'mean_{group}']=mean[group]
    res=res.reset_index();res.insert(0,'interaction_id',res[KEYS].agg('|'.join,axis=1))
    return res.sort_values('fdr_q_value')

def main():
    meta=pd.read_csv(BASE/'Results/severity_analysis/patient_clinical_metadata.csv')[['sample','Groups']]
    models={}
    for state in ['baseline','annotation','resource','combined']:
        paths=list((LOCAL/'scores'/state).glob('*.csv.gz'))
        if {p.name.removesuffix('.csv.gz') for p in paths}!=set(meta['sample']):
            raise ValueError(f'{state}: incomplete donor inference ({len(paths)} files; expected 56)')
        frames=[]
        for path in paths:
            frame=pd.read_csv(path,usecols=KEYS+['sample','magnitude_rank'])
            frame['score']=1-frame.magnitude_rank;frames.append(frame)
        scores=pd.concat(frames,ignore_index=True).merge(meta,on='sample',validate='many_to_one')
        if scores.duplicated(['sample']+KEYS).any(): raise ValueError('Duplicate donor interaction')
        models[state]=fit(scores)
    original=pd.read_csv(BASE/'Results/severity_analysis/age_diseasestatus_interaction.csv')
    alternative=pd.read_csv(OUT/'cellphonedb_resource.csv')
    pairs=set(zip(alternative.ligand,alternative.receptor))
    original_sig=original[original.fdr_q_value<.05].copy()
    original_sig['in_alternative_resource']=[(a,b) in pairs for a,b in zip(original_sig.ligand_complex,original_sig.receptor_complex)]
    original_sig.to_csv(OUT/'original_significant_resource_membership.csv',index=False)
    b=models['baseline'];m=original.merge(b,on='interaction_id',suffixes=('_original','_rerun'),validate='one_to_one')
    if len(m)!=len(original) or len(b)!=len(original):raise ValueError('Baseline tested interaction set changed')
    checks={}
    for col in ['interaction_coef','interaction_ci95_low','interaction_ci95_high','interaction_p_value','fdr_q_value','n_observations']:
        delta=float((m[col+'_original']-m[col+'_rerun']).abs().max());checks[col]=delta
        if delta>1e-7:raise ValueError(f'Baseline reproduction failed for {col}: {delta}')
    (OUT/'baseline_reproduction.json').write_text(json.dumps(checks,indent=2))
    summary=[]
    for state,table in models.items():
        table.to_csv(OUT/f'{state}_models.csv',index=False)
        sig=table[table.fdr_q_value<.05]
        common=b.merge(table,on='interaction_id',suffixes=('_baseline','_scenario'),validate='one_to_one')
        common.to_csv(OUT/f'{state}_matched_comparison.csv',index=False)
        baseline_sig_common=common[common.fdr_q_value_baseline<.05]
        row={'scenario':state,'tested':len(table),'fdr05':len(sig),'fdr10':int((table.fdr_q_value<.1).sum()),
            'positive_fdr05':int((sig.interaction_coef>0).sum()),'negative_fdr05':int((sig.interaction_coef<0).sum()),
            'shared_tests':len(common),'coefficient_pearson_r':common.interaction_coef_baseline.corr(common.interaction_coef_scenario),
            'direction_agreement_shared':float((np.sign(common.interaction_coef_baseline)==np.sign(common.interaction_coef_scenario)).mean()),
            'original_significant_testable':len(baseline_sig_common),
            'original_significant_still_fdr05':int((baseline_sig_common.fdr_q_value_scenario<.05).sum()),
            'original_significant_direction_reversals':int((np.sign(baseline_sig_common.interaction_coef_baseline)!=np.sign(baseline_sig_common.interaction_coef_scenario)).sum()),
            'nonclassical_target_fdr05':int(sig.target.eq('Non-classical Monocytes').sum()),
            'min_donors':int(table.n_observations.min()),'max_donors':int(table.n_observations.max())}
        summary.append(row)
        tested=table.groupby('target').size().rename('tested')
        hit=sig.groupby('target').size().rename('fdr05')
        targets=pd.concat([tested,hit],axis=1).fillna(0);targets['fraction_fdr05']=targets.fdr05/targets.tested
        targets.to_csv(OUT/f'{state}_target_summary.csv')
    pd.DataFrame(summary).to_csv(OUT/'scenario_summary.csv',index=False)
    print(pd.DataFrame(summary).to_string(index=False))
    review_rows=[]
    group_rows=[]
    for cohort in ['child','adult']:
        review=pd.read_csv(LOCAL/f'{cohort}_cell_review.csv.gz',dtype={'cluster':str})
        grouped_review=review.merge(meta,on='sample',validate='many_to_one')
        for group,d in grouped_review.groupby('Groups'):
            original_retained=(~d.original_label.isin(['Ribosomal/Low-quality','RBCs','Platelets'])).sum()
            group_rows.append({'group':group,'original_retained':int(original_retained),'restored':int(d.restore.sum()),
                'sensitivity_retained':int(original_retained+d.restore.sum()),'donors_with_restored_cells':int(d.loc[d.restore,'sample'].nunique())})
        for cluster in (['0','1'] if cohort=='child' else ['0','4']):
            d=review[review.cluster.eq(cluster)]
            review_rows.append({'cohort':cohort,'cluster':cluster,'cells':len(d),
                'reference_CD4_or_CD8_fraction':d.reference_broad.ne('unresolved').mean(),
                'restored':int(d.restore.sum()),'restored_p08':int((d.restore & d.reference_probability.ge(.8)).sum()),
                'median_ribosomal_percent':d.ribosomal_percent_retained_genes.median(),
                'median_mito_percent':d.pct_counts_mt.median(),'median_genes':d.n_genes_by_counts.median()})
        pd.crosstab(review.original_label,review.reference_broad).to_csv(OUT/f'{cohort}_original_vs_reference_broad.csv')
        review.groupby(['original_label','restore'])[['ribosomal_percent_retained_genes','pct_counts_mt','n_genes_by_counts','doublet_score']].median().to_csv(OUT/f'{cohort}_qc_by_restoration.csv')
    pd.DataFrame(review_rows).to_csv(OUT/'excluded_cluster_review_summary.csv',index=False)
    pd.DataFrame(group_rows).to_csv(OUT/'restoration_by_group.csv',index=False)
    (OUT/'analysis_provenance.json').write_text(json.dumps({'completed_utc':datetime.now(timezone.utc).isoformat(),
        'scripts_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),BASE/'Notebooks/14_ANNOTATION_RESOURCE_SENSITIVITY.py']},
        'input_cell_count':348274,'donors':56,'statistical_model':'Four-group OLS interaction; HC3; normal Wald CI/p values; BH by scenario',
        'original_results_modified':False},indent=2))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,3,figsize=(12,4),layout='constrained')
    for ax,state in zip(axes,['annotation','resource','combined']):
        data=pd.read_csv(OUT/f'{state}_matched_comparison.csv')
        selected=data.fdr_q_value_baseline<.05
        ax.scatter(data.interaction_coef_baseline,data.interaction_coef_scenario,c=np.where(selected,'#235789','#bbbbbb'),s=7,alpha=.65,rasterized=True)
        ax.plot([-.7,.7],[-.7,.7],'k--',lw=.7);ax.axhline(0,c='grey',lw=.5);ax.axvline(0,c='grey',lw=.5)
        ax.set(xlabel='Original interaction coefficient',ylabel='Sensitivity interaction coefficient',title=state.capitalize())
    fig.savefig(OUT/'sensitivity_coefficient_comparison.png',dpi=300)
    plt.close(fig)

if __name__=='__main__':main()
