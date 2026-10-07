"""Summarize corrected LIANA inference and reviewer-requested sensitivity checks."""
import importlib.util
from pathlib import Path
spec=importlib.util.spec_from_file_location('revision',Path(__file__).with_name('18_REVIEWER_REVISION.py'))
r=importlib.util.module_from_spec(spec);spec.loader.exec_module(r)
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests
import statsmodels.formula.api as smf
import numpy as np, pandas as pd, json, warnings
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
BASE,OUT,LOCAL,KEYS,GROUPS=r.BASE,r.OUT,r.LOCAL,r.KEYS,r.GROUPS
def bh(d,p='p_value',q='fdr_q_value'):
    d[q]=np.nan
    ok=np.isfinite(d[p]);d.loc[ok,q]=multipletests(d.loc[ok,p],method='fdr_bh')[1]
    return d
def compare(a,b,label):
    m=a.merge(b,on='interaction_id',suffixes=('_primary','_sensitivity'))
    m.to_csv(OUT/f'{label}_comparison.csv',index=False)
    both=m[(m.fdr_q_value_primary<.05)&(m.fdr_q_value_sensitivity<.05)]
    return dict(tested=len(b),significant=int((b.fdr_q_value<.05).sum()),shared_tests=len(m),coefficient_r=float(m.interaction_coef_primary.corr(m.interaction_coef_sensitivity)),primary_significant_testable=int((m.fdr_q_value_primary<.05).sum()),both_significant=len(both),both_significant_same_direction=int((np.sign(both.interaction_coef_primary)==np.sign(both.interaction_coef_sensitivity)).sum()))
def robust_model(d,formula):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore');m=smf.ols(formula,d).fit(cov_type='HC3')
    rank=np.linalg.matrix_rank(m.model.exog);n,p=m.model.exog.shape
    leverage=np.max(np.sum(m.model.exog*np.linalg.pinv(m.model.exog).T,axis=1))
    if rank<p or n<=p or leverage>=1-1e-9:return None,dict(n=n,rank=rank,parameters=p,max_leverage=float(leverage),status='not_estimable')
    return m,dict(n=n,rank=rank,parameters=p,max_leverage=float(leverage),condition_number=float(np.linalg.cond(m.model.exog)),status='ok')
def main():
    print('Loading corrected donor scores and fitting primary models',flush=True)
    meta=r.metadata();old=pd.read_csv(BASE/'Results/revised_primary/primary_models.csv')
    assert {p.name.removesuffix('.csv.gz') for p in (LOCAL/'scores/fresh_primary').glob('*.csv.gz')}==set(meta['sample'])
    scores=r.identify(r.load_scores('fresh_primary'))
    primary=r.fit(scores)
    cached=r.identify(r.load_scores('annotation'))
    compared=scores[['sample','interaction_id','magnitude_rank']].merge(cached[['sample','interaction_id','old_magnitude_rank','magnitude_rank']],on=['sample','interaction_id'],suffixes=('_fresh','_cached'),validate='one_to_one')
    stats={'old_vs_corrected':compare(primary,old,'old_vs_corrected')}
    stats['old_vs_corrected']['old_significant_retained']=stats['old_vs_corrected']['both_significant']
    stats['score_pearson_r']=float(compared.old_magnitude_rank.corr(compared.magnitude_rank_fresh))
    stats['score_spearman_r']=float(compared.old_magnitude_rank.corr(compared.magnitude_rank_fresh,method='spearman'))
    stats['max_fresh_vs_cached_score_difference']=float((compared.magnitude_rank_fresh-compared.magnitude_rank_cached).abs().max())
    lost=old[~old.interaction_id.isin(primary.interaction_id)].copy()
    ranges=scores.groupby('interaction_id').score.agg(['min','max','count'])
    lost=lost.merge(ranges,on='interaction_id',how='left');lost['reason']='constant corrected scores; non-estimable model'
    lost.to_csv(OUT/'previously_tested_now_nonestimable.csv',index=False)
    primary.to_csv(OUT/'primary_models.csv',index=False)
    significant=primary[primary.fdr_q_value<.05]
    significant.to_csv(OUT/'supplementary_all_significant.csv',index=False)
    independent=[]
    for iid in list(primary.head(5).interaction_id)+r.IDS:
        d=scores[scores.interaction_id.eq(iid)]
        m=smf.ols('score ~ C(Groups)',d).fit(cov_type='HC3')
        contrast=np.zeros(len(m.params))
        for name,sign in [('C(Groups)[T.aSLE]',-1),('C(Groups)[T.cHD]',-1),('C(Groups)[T.cSLE]',1)]:contrast[list(m.params.index).index(name)]=sign
        t=m.t_test(contrast);saved=primary.set_index('interaction_id').loc[iid]
        delta=abs(float(np.asarray(t.effect).item())-saved.interaction_coef)
        pdelta=abs(float(np.asarray(t.pvalue).item())-saved.interaction_p_value)
        assert delta<1e-10 and pdelta<1e-8
        independent.append({'interaction_id':iid,'coefficient_difference':delta,'p_value_difference':pdelta})
    pd.DataFrame(independent).to_csv(OUT/'independent_model_verification.csv',index=False)
    for state in ['baseline','resource','combined']:
        stats[state]=compare(primary,pd.read_csv(OUT/f'{state}_models.csv'),state)
    female=r.fit(scores[scores.Gender.eq('F')]);female.to_csv(OUT/'female_models.csv',index=False)
    stats['female']=compare(primary,female,'female')
    original=pd.read_csv(OUT/'baseline_models.csv')
    stable=significant.merge(original,on='interaction_id',suffixes=('_primary','_original'))
    stable=stable[stable.fdr_q_value_original<.05];stable.to_csv(OUT/'annotation_shared_significant.csv',index=False)
    # The submitted overlap is historical and must not be substituted for corrected inference.
    legacy_original=pd.read_csv(BASE/'Results/annotation_resource_sensitivity/baseline_models.csv')
    legacy_ids=set(old.loc[old.fdr_q_value<.05,'interaction_id']) & set(legacy_original.loc[legacy_original.fdr_q_value<.05,'interaction_id'])
    primary[primary.interaction_id.isin(legacy_ids)].to_csv(OUT/'submitted_104_after_correction.csv',index=False)
    stats['submitted_overlap_count']=len(legacy_ids)
    stats['submitted_overlap_corrected_significant']=int(primary.loc[primary.interaction_id.isin(legacy_ids),'fdr_q_value'].lt(.05).sum())
    targets=primary.groupby('target').size().rename('tested').to_frame().join(significant.groupby('target').size().rename('significant')).fillna(0)
    targets['percent_significant']=100*targets.significant/targets.tested;targets.to_csv(OUT/'receiving_cell_denominators.csv')
    pd.crosstab(meta.Groups,meta.Batch).to_csv(OUT/'group_by_batch_counts.csv')
    # Batch-adjusted primary model, checking estimability interaction by interaction.
    print('Fitting additive batch models',flush=True)
    formula='score ~ C(Groups) + C(Batch)'
    batch=[]
    for iid,d in scores[scores.interaction_id.isin(primary.interaction_id)].groupby('interaction_id',sort=False):
        m,diag=robust_model(d,formula);item={'interaction_id':iid,**diag,'interaction_coef':np.nan,'p_value':np.nan}
        if m is not None:
            contrast=np.zeros(len(m.params))
            for name,sign in [('C(Groups)[T.aSLE]',-1),('C(Groups)[T.cHD]',-1),('C(Groups)[T.cSLE]',1)]:contrast[list(m.params.index).index(name)]=sign
            test=m.t_test(contrast);ci=np.asarray(test.conf_int()).ravel()
            item.update(interaction_coef=float(np.asarray(test.effect).item()),p_value=float(np.asarray(test.pvalue).item()),ci_low=float(ci[0]),ci_high=float(ci[1]))
        batch.append(item)
    batch=bh(pd.DataFrame(batch));batch.to_csv(OUT/'batch_adjusted_models.csv',index=False)
    stats['batch']=compare(primary,batch,'batch_adjusted')
    stats['batch']['nonestimable']=int(batch.status.ne('ok').sum())
    # Existing Table 2 examples remain fixed even when no longer significant.
    example=primary.set_index('interaction_id').loc[r.IDS].reset_index()
    example=example.merge(old[['interaction_id','interaction_coef','fdr_q_value']],on='interaction_id',suffixes=('','_submitted'))
    example['coefficient_change']=example.interaction_coef-example.interaction_coef_submitted
    example['same_coefficient_direction']=np.sign(example.interaction_coef)==np.sign(example.interaction_coef_submitted)
    example.to_csv(OUT/'table2_corrected.csv',index=False)
    selected=scores[scores.interaction_id.isin(r.IDS)].copy();selected.to_csv(OUT/'table2_individual_donor_scores.csv',index=False)
    batchchecks=[];loo=[]
    for iid,d in selected.groupby('interaction_id',sort=False):
        ped=d[d.Groups.eq('cSLE')];m,diag=robust_model(ped,'score ~ C(Batch)')
        item={'interaction_id':iid,**diag,'p_value':np.nan,'batch_sizes':json.dumps(ped.Batch.value_counts().sort_index().to_dict())}
        if m is not None:
            test=m.wald_test(np.eye(len(m.params))[1:],scalar=True);item['p_value']=float(test.pvalue)
        batchchecks.append(item)
        full=float(primary.set_index('interaction_id').loc[iid,'interaction_coef'])
        for sample in d['sample']:
            sub=d[d['sample'].ne(sample)];means=sub.groupby('Groups').score.mean()
            coef=float(means.cSLE-means.cHD-means.aSLE+means.aHD)
            loo.append({'interaction_id':iid,'omitted_sample':sample,'omitted_group':d.loc[d['sample'].eq(sample),'Groups'].iloc[0],'full_coef':full,'leave_one_out_coef':coef,'delta':coef-full,'direction_retained':np.sign(coef)==np.sign(full)})
    bh(pd.DataFrame(batchchecks)).to_csv(OUT/'table2_pediatric_SLE_batch_tests.csv',index=False)
    loo=pd.DataFrame(loo);loo.to_csv(OUT/'table2_leave_one_donor_out.csv',index=False)
    stats['table2_loo_direction_reversals']=int((~loo.direction_retained).sum())
    # Require all donor reruns before reporting missingness and return-all checks.
    print('Summarizing donor availability and return-all sensitivities',flush=True)
    verification=[json.loads(p.read_text()) for p in LOCAL.glob('*_verification.json')]
    assert {v['sample'] for v in verification}==set(meta['sample']),'Donor inference incomplete'
    pd.DataFrame(verification).to_csv(OUT/'inference_verification.csv',index=False)
    audit=pd.concat([pd.read_csv(p) for p in LOCAL.glob('*_availability.csv.gz')],ignore_index=True).merge(meta[['sample','Groups','Batch']],on='sample',validate='many_to_one')
    # This tested family required scores in all four groups, so its subunits are
    # present in both cohort matrices. No-population fractions are not missing genes.
    for role in ['source','target']:
        audit.loc[audit[role+'_cells'].eq(0),role+'_gene_unavailable']=False
    def reason(x):
        if x.observed:return 'observed'
        if not x.population_eligible:
            return '; '.join(role+' population has fewer than 10 cells' for role in ['source','target'] if x[role+'_under_10_cells'])
        return '; '.join(label for flag,label in [('source_below_10pct','ligand below 10%'),('target_below_10pct','receptor below 10%'),('source_gene_unavailable','ligand gene unavailable'),('target_gene_unavailable','receptor gene unavailable')] if x[flag])
    audit['missing_reason']=audit.apply(reason,axis=1)
    assert audit.loc[~audit.observed,'missing_reason'].str.len().gt(0).all()
    audit.to_csv(OUT/'all_primary_tests_donor_availability.csv.gz',index=False)
    relevant=set(significant.interaction_id)|set(old.loc[old.fdr_q_value<.05,'interaction_id'])
    audit[audit.interaction_id.isin(relevant)].to_csv(OUT/'significant_donor_missingness_reasons.csv.gz',index=False)
    cols=['observed','population_eligible','source_under_10_cells','target_under_10_cells','source_below_10pct','target_below_10pct','source_gene_unavailable','target_gene_unavailable']
    summary=audit.groupby(['interaction_id','Groups'])[cols].sum().reset_index()
    summary['total_donors']=summary.Groups.map(meta.Groups.value_counts());summary.to_csv(OUT/'donor_availability_by_group.csv',index=False)
    allscores=r.identify(r.load_scores('return_all'))
    print('Loaded return-all donor rows:',len(allscores),flush=True)
    # All interactions testable in the original family, retaining threshold failures where LIANA returns them.
    allfamily=allscores[allscores.interaction_id.isin(old.interaction_id)]
    allmodel=r.fit(allfamily);allmodel.to_csv(OUT/'return_all_original_family_models.csv',index=False)
    stats['return_all']=compare(primary,allmodel,'return_all')
    # Common support: population coverage and genes fixed without using disease/age labels or scores.
    counts=pd.concat([pd.read_csv(BASE/f'Results/revised_primary/qc/{c}_sample_celltype_counts.csv') for c in ['child','adult']])
    covered=counts.pivot(index='sample',columns='cell_type',values='n_cells').reindex(meta['sample']).fillna(0).ge(10).all()
    commoncells=set(covered[covered].index)
    fixed=allscores[allscores.source.isin(commoncells)&allscores.target.isin(commoncells)]
    # Require all 56 return-all observations: donor-wide absent genes are not returned by LIANA.
    complete=fixed.groupby('interaction_id')['sample'].nunique();complete=set(complete[complete.eq(56)].index)
    common=r.fit(fixed[fixed.interaction_id.isin(complete)]);common.to_csv(OUT/'common_support_models.csv',index=False)
    pd.Series(sorted(complete),name='interaction_id').to_csv(OUT/'common_support_fixed_universe.csv',index=False)
    stats['common_support']=compare(primary,common,'common_support');stats['common_support']['fixed_universe']=len(complete);stats['common_support']['cell_types']=sorted(commoncells)
    # Detection component conditional on sufficient cell population sizes; expression failure remains 0.
    det=audit[audit.population_eligible].copy();det['score']=det.observed.astype(float)
    variation=det.groupby(['interaction_id','Groups']).score.var().groupby('interaction_id').sum()
    estimable=set(variation[variation.gt(0)].index)
    detect=r.fit(det[det.interaction_id.isin(estimable)]);detect.to_csv(OUT/'two_part_detection_models.csv',index=False)
    two=old[['interaction_id']+KEYS].merge(detect[['interaction_id','interaction_coef','fdr_q_value','n_observations']],on='interaction_id',how='left').rename(columns={'interaction_coef':'detection_coef','fdr_q_value':'detection_q','n_observations':'detection_n'})
    two=two.merge(primary[['interaction_id','interaction_coef','fdr_q_value','n_observations']],on='interaction_id',how='left').rename(columns={'interaction_coef':'conditional_score_coef','fdr_q_value':'conditional_score_q','n_observations':'conditional_score_n'})
    two['detection_status']=np.where(two.detection_q.isna(),'constant_or_insufficient_support','estimated')
    two.to_csv(OUT/'two_part_models.csv',index=False)
    stats['detection']={'tested':len(detect),'significant':int((detect.fdr_q_value<.05).sum()),'constant_or_insufficient_support':int(two.detection_q.isna().sum())}
    stats.update(primary_tested=len(primary),primary_significant=len(significant),primary_fdr10=int(primary.fdr_q_value.lt(.1).sum()),positive=int(significant.interaction_coef.gt(0).sum()),negative=int(significant.interaction_coef.lt(0).sum()),min_donors=int(primary.n_observations.min()),max_donors=int(primary.n_observations.max()),incomplete_tests=int(primary.n_observations.lt(56).sum()),annotation_shared=len(stable))
    (OUT/'summary.json').write_text(json.dumps(stats,indent=2))
    import hashlib,importlib
    aggregate=Path(importlib.import_module('liana.method._pipe_utils._aggregate').__file__)
    inputs=[aggregate,BASE/'Notebooks/18_REVIEWER_REVISION.py',BASE/'Notebooks/19_SUMMARIZE_REVIEWER_REVISION.py',OUT/'consensus_resource.csv',OUT/'cellphonedb_resource.csv']+[BASE/f'Results/revised_primary/{c}_primary_cell_labels.csv.gz' for c in ['child','adult']]
    provenance={'date':'2026-10-06','patch':'Process-local unique-score ranking; no package files changed','upstream_changelog':'https://github.com/scverse/liana-py/blob/main/CHANGELOG.md','old_results_preserved':True,'primary':'Fresh inference from archived normalized log-transformed X','other_scenarios':'Corrected aggregation of unchanged cached per-method scores','sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs}}
    (OUT/'provenance.json').write_text(json.dumps(provenance,indent=2))
    figures(primary,female,selected)
    table_html(example)
    from render_table2_png import render
    render()
    print(json.dumps(stats,indent=2),flush=True)
def figures(primary,female,selected):
    figure_dir = OUT.parent / 'manuscript_figures'
    figure_dir.mkdir(parents=True, exist_ok=True)
    def save(fig,name):
        fig.savefig(figure_dir/(name+'.png'),dpi=400,bbox_inches='tight');fig.savefig(figure_dir/(name+'.svg'),bbox_inches='tight');plt.close(fig)
    fig,ax=plt.subplots(figsize=(6,4.5));hit=primary.fdr_q_value<.05
    ax.scatter(primary.interaction_coef,-np.log10(primary.fdr_q_value.clip(lower=1e-300)),c=np.where(hit,'#bd382d','#a5a5a5'),s=9)
    ax.axhline(-np.log10(.05),ls='--',c='black',lw=.7);ax.axvline(0,c='black',lw=.7)
    ax.set(xlabel='Cohort-by-SLE difference in relative rank score',ylabel='-log10(adjusted p-value)',title='Primary cohort-by-SLE analysis');save(fig,'figure2_primary_analysis')
    d=primary.merge(female,on='interaction_id',suffixes=('_p','_f'));fig,ax=plt.subplots(figsize=(6.5,4.5))
    ax.scatter(d.interaction_coef_p,d.interaction_coef_f,s=9,c=np.where(d.fdr_q_value_p<.05,'#2878a5','#aaa'))
    ax.plot([-.8,.8],[-.8,.8],'k--',lw=.7);ax.set(xlabel='All-donor coefficient',ylabel='Female-only coefficient',title=f'All-donor versus female-only: r = {d.interaction_coef_p.corr(d.interaction_coef_f):.3f}');save(fig,'figure3_female_only_comparison')
    plot_ids = pd.read_csv(OUT/'table2_corrected.csv').query('fdr_q_value < 0.05').interaction_id.tolist()
    assert len(plot_ids) == 4
    fig,axes=plt.subplots(2,2,figsize=(12,9),constrained_layout=True);rng=np.random.default_rng(0)
    for ax,iid in zip(axes.flat,plot_ids):
        d=selected[selected.interaction_id.eq(iid)];a,b,l,rec=iid.split('|')
        for i,g in enumerate(GROUPS):
            vals=d[d.Groups.eq(g)].score.to_numpy();ax.scatter(i+rng.uniform(-.12,.12,len(vals)),vals,s=48,alpha=.8)
            ax.plot([i-.18,i+.18],[vals.mean()]*2,c='black',lw=2)
        ax.set(xticks=range(4),xticklabels=['Pediatric\nhealthy','Pediatric\nSLE','Adult\nhealthy','Adult\nSLE'],ylabel='Relative rank score',title=f'{l}-{rec}\n{a} to\n{b}')
        ax.tick_params(axis='both',labelsize=16)
        ax.title.set_fontsize(18)
        ax.title.set_fontweight('bold')
        ax.yaxis.label.set_size(17)
    save(fig,'figure4_individual_donor_scores')
    fig,axes=plt.subplots(2,2,figsize=(12,9),constrained_layout=True)
    for ax,iid in zip(axes.flat,plot_ids):
        d=selected[selected.interaction_id.eq(iid)&selected.Groups.eq('cSLE')];a,b,l,rec=iid.split('|')
        for i,batch in enumerate(['B1','B2','B3','B4','B5','B6']):
            vals=d[d.Batch.eq(batch)].score.to_numpy();ax.scatter(i+rng.uniform(-.12,.12,len(vals)),vals,s=48)
            if len(vals):ax.plot([i-.18,i+.18],[vals.mean()]*2,c='black',lw=2)
        ax.set(xticks=range(6),xticklabels=['B1','B2','B3','B4','B5','B6'],ylabel='Relative rank score',title=f'{l}-{rec}\n{a} to\n{b}')
        ax.tick_params(axis='both',labelsize=16)
        ax.title.set_fontsize(18)
        ax.title.set_fontweight('bold')
        ax.yaxis.label.set_size(17)
    save(fig,'figure5_pediatric_SLE_by_batch')
def table_html(d):
    import html
    d = d.loc[d.fdr_q_value < 0.05].copy()
    assert len(d) == 4, 'Expected four significant illustrative combinations'
    a=[];b=[]
    for _,x in d.iterrows():
        cells=x.source+' to '+x.target;pair=x.ligand_complex+'-'+x.receptor_complex
        a.append([cells,pair,f'{x.interaction_coef:.3f}',f'{x.interaction_ci95_low:.3f} to {x.interaction_ci95_high:.3f}',f'{x.fdr_q_value:.5g}'])
        b.append([cells,pair]+[f'{x["mean_"+g]:.3f} (n={int(x["n_"+g])})' for g in GROUPS])
    a=pd.DataFrame(a,columns=['Source to target','Ligand-receptor','Coefficient','95% CI','FDR q'])
    b=pd.DataFrame(b,columns=['Source to target','Ligand-receptor','Pediatric healthy','Pediatric SLE','Adult healthy','Adult SLE'])
    a.to_csv(OUT/'table2_panel_A.csv',index=False);b.to_csv(OUT/'table2_panel_B.csv',index=False)
    page='<!doctype html><meta charset="utf-8"><title>Table 2</title><style>body{font:12pt Arial;margin:30px}table{border-collapse:collapse;width:100%;margin:18px 0}th,td{border:1px solid #bbb;padding:8px;text-align:left}th{background:#eee}</style><h1>Table 2. Illustrative interaction combinations</h1><h2>Panel A: interaction estimates</h2>'+a.to_html(index=False)+'<h2>Panel B: group means and donor counts</h2>'+b.to_html(index=False)+'<p>Outcome: 1 - corrected magnitude_rank, a relative rank score, not absolute signaling intensity. Four illustrative combinations meeting FDR &lt; 0.05 in the primary analysis are shown. These examples are not the complete set of significant results. Counts are contributing donors. Coefficients compare the SLE-versus-healthy change between pediatric and adult cohorts. CI: HC3 95% confidence interval. FDR q: Benjamini-Hochberg adjustment across the full eligible primary family.</p>'
    (OUT/'table2_editable.html').write_text(page,encoding='utf-8')
if __name__=='__main__':main()
