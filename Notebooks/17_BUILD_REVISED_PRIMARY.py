"""Build the revised primary analysis from reviewed script-14 annotations.

Preserves original statistical results; outputs an explicitly versioned package.
Requires script 14 annotation/consensus donor scores and script 16 review outputs.
"""
from pathlib import Path
import argparse
import importlib.util
import os
import json
import hashlib
import shutil
import sys
BASE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(BASE/'Notebooks'))
os.environ['MPLCONFIGDIR']=str(BASE/'_local_submission/sensitivity_cache/matplotlib')
os.environ['PYTHONDONTWRITEBYTECODE']='1'
import pandas as pd
import numpy as np
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
OUT=BASE/'Results/revised_primary'
FIG=OUT/'manuscript_figures'
TABLE=OUT/'manuscript_tables'
KEYS=['source','target','ligand_complex','receptor_complex']

def load_module(name,filename):
    spec=importlib.util.spec_from_file_location(name,BASE/'Notebooks'/filename)
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod);return mod

def savefig(fig,name):
    fig.savefig(FIG/f'{name}.png',dpi=600,bbox_inches='tight',facecolor='white')
    fig.savefig(FIG/f'{name}.jpg',dpi=600,bbox_inches='tight',facecolor='white',pil_kwargs={'quality':95})
    plt.close(fig)

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--archive-root',type=Path,required=True)
    parser.add_argument('--publish',action='store_true',help='Refresh visible manuscript figures/tables after preserving their previous versions locally.')
    args=parser.parse_args()
    for p in [OUT,FIG,TABLE,OUT/'qc']:p.mkdir(parents=True,exist_ok=True)
    assert (BASE/'Results/annotation_review/group_review.csv').exists()
    fit=load_module('model','15_SUMMARIZE_ANNOTATION_RESOURCE_SENSITIVITY.py').fit
    helper=load_module('annotations','14_ANNOTATION_RESOURCE_SENSITIVITY.py')
    meta=pd.read_csv(BASE/'Results/severity_analysis/patient_clinical_metadata.csv')
    frames=[]
    paths=list((helper.LOCAL/'scores/annotation').glob('*.csv.gz'))
    assert {p.name.removesuffix('.csv.gz') for p in paths}==set(meta['sample'])
    for p in paths:
        frame=pd.read_csv(p,usecols=KEYS+['sample','magnitude_rank'])
        frame['score']=1-frame.magnitude_rank;frames.append(frame)
    scores=pd.concat(frames,ignore_index=True).merge(meta,on='sample',validate='many_to_one')
    assert not scores.duplicated(['sample']+KEYS).any()
    model=fit(scores);female=fit(scores[scores.Gender.eq('F')])
    # Independent formula-based HC3 checks, including revised female-only fits.
    import statsmodels.formula.api as smf
    verification=[]
    scores['interaction_id']=scores[KEYS].agg('|'.join,axis=1)
    scores['condition']=np.where(scores.Groups.str.endswith('SLE'),'SLE','healthy')
    scores['age_group']=np.where(scores.Groups.str.startswith('c'),'pediatric','adult')
    for name,data,table in [('primary',scores,model),('female',scores[scores.Gender.eq('F')],female)]:
        selected=pd.concat([table.head(3),table.iloc[(table.fdr_q_value-.05).abs().argsort()[:3]]]).drop_duplicates('interaction_id')
        for _,row in selected.iterrows():
            result=smf.ols('score ~ C(condition, Treatment(reference="healthy")) * C(age_group, Treatment(reference="adult"))',data=data[data.interaction_id.eq(row.interaction_id)]).fit(cov_type='HC3')
            term=next(t for t in result.params.index if ':' in t)
            difference=abs(result.params[term]-row.interaction_coef)
            p_difference=abs(result.pvalues[term]-row.interaction_p_value)
            assert difference<1e-10 and p_difference<1e-8
            verification.append({'scenario':name,'interaction_id':row.interaction_id,'coefficient_difference':difference,'p_value_difference':p_difference})
    pd.DataFrame(verification).to_csv(OUT/'independent_HC3_verification.csv',index=False)
    expected=pd.read_csv(BASE/'Results/annotation_resource_sensitivity/annotation_models.csv')
    check=model.merge(expected,on='interaction_id',suffixes=('_new','_saved'),validate='one_to_one')
    assert len(check)==len(model)==len(expected)
    assert np.allclose(check.interaction_coef_new,check.interaction_coef_saved,atol=1e-12)
    assert np.allclose(check.fdr_q_value_new,check.fdr_q_value_saved,atol=1e-10)
    model.to_csv(OUT/'primary_models.csv',index=False);female.to_csv(OUT/'female_only_models.csv',index=False)
    sig=model[model.fdr_q_value<.05].copy()
    sig.to_csv(TABLE/'supplementary_all_fdr05_interactions.csv',index=False)
    means=sig.copy()
    means['pediatric_sle_minus_healthy']=means.mean_cSLE-means.mean_cHD
    means['adult_sle_minus_healthy']=means.mean_aSLE-means.mean_aHD
    means['difference_in_differences_from_means']=means.pediatric_sle_minus_healthy-means.adult_sle_minus_healthy
    means.to_csv(TABLE/'significant_interaction_group_means.csv',index=False)
    compare=model.merge(female,on='interaction_id',suffixes=('_primary','_female'),validate='one_to_one')
    compare.to_csv(OUT/'female_only_comparison.csv',index=False)
    # Freeze exact per-cell labels and exclusions, separate from archived labels.
    datasets=[];counts=[];accounting=[]
    for cohort,filename in helper.FILES.items():
        r=pd.read_csv(helper.LOCAL/f'{cohort}_cell_review.csv.gz',index_col='cell_id',dtype={'cluster':str})
        keep=~r.sensitivity_label.isin(['Ribosomal/Low-quality','RBCs','Platelets'])
        r['primary_label']=r.sensitivity_label.mask(r.sensitivity_label.eq('Ribosomal/Low-quality'),'Unresolved T-rich cluster')
        r['included_primary']=keep
        r[['sample','cluster','original_label','primary_label','included_primary','restore','reference_label','reference_probability']].to_csv(OUT/f'{cohort}_primary_cell_labels.csv.gz')
        cnt=r.loc[keep].groupby(['sample','primary_label']).size().rename('n_cells').reset_index().rename(columns={'primary_label':'cell_type'})
        cnt.to_csv(OUT/'qc'/f'{cohort}_sample_celltype_counts.csv',index=False);counts.append(cnt)
        for sample,group in r.groupby('sample'):
            accounting.append({'cohort':cohort,'sample':sample,'post_qc':len(group),'original_retained':int((group.included_primary & ~group.restore).sum()),'restored':int(group.restore.sum()),'primary_retained':int(group.included_primary.sum()),'unresolved_excluded':int(group.primary_label.eq('Unresolved T-rich cluster').sum()),'RBC_platelet_excluded':int(group.primary_label.isin(['RBCs','Platelets']).sum())})
        with h5py.File(args.archive_root/f'{cohort}_individual_h5ad'/filename,'r') as f:
            obs,genes=helper.metadata(f,cohort);assert obs.index.equals(r.index)
            coords=f['obsm']['X_umap'][:]
        datasets.append((cohort,coords[keep],r.loc[keep,'primary_label'].to_numpy()))
    pd.DataFrame(accounting).to_csv(OUT/'qc/annotation_accounting.csv',index=False)
    counts=pd.concat(counts).merge(meta[['sample','Groups']],on='sample',validate='many_to_one')
    counts.to_csv(TABLE/'supplementary_donor_cell_type_coverage.csv',index=False)
    rows=[]
    for (group,ct),part in counts.groupby(['Groups','cell_type']):
        rows.append({'group':group,'cell_type':ct,'total_donors':int(meta.Groups.eq(group).sum()),'donors_with_any_cells':part['sample'].nunique(),'donors_with_at_least_10_cells':int(part.n_cells.ge(10).sum()),'minimum_cells_among_donors_with_any':int(part.n_cells.min())})
    pd.DataFrame(rows).to_csv(TABLE/'supplementary_cell_type_coverage.csv',index=False)
    table1=pd.read_csv(BASE/'Results/manuscript_tables/table1_cohort_characteristics.csv')
    groupnames={'cHD':'Pediatric healthy','cSLE':'Pediatric SLE','aHD':'Adult healthy','aSLE':'Adult SLE'}
    totals=counts.groupby('Groups').n_cells.sum().rename(index=groupnames)
    table1['retained_cells']=table1.group.map(totals).astype(int)
    table1.to_csv(TABLE/'table1_cohort_characteristics.csv',index=False)
    shutil.copy2(BASE/'Results/manuscript_tables/supplementary_adult_sample_exclusions.csv',TABLE/'supplementary_adult_sample_exclusions.csv')
    for col in ['source','target']:
        sig[col].value_counts().rename_axis(col).reset_index(name='n_significant_interactions').to_csv(TABLE/f'significant_{col}s.csv',index=False)
    sig.groupby(['source','target']).size().sort_values(ascending=False).rename('n_significant_interactions').reset_index().to_csv(TABLE/'significant_source_target_pairs.csv',index=False)
    families=sig.assign(ligand_receptor=sig.ligand_complex+'-'+sig.receptor_complex).groupby('ligand_receptor').agg(n_significant_interactions=('interaction_id','size'),n_positive=('interaction_coef',lambda x:int((x>0).sum())),n_negative=('interaction_coef',lambda x:int((x<0).sum())))
    families.sort_values('n_significant_interactions',ascending=False).to_csv(TABLE/'significant_ligand_receptor_families.csv')
    # Preserve the original seven illustrative combinations, including any loss of significance.
    keys=[('Classical Monocytes','NK Cells','TIMP2','ITGB1'),('CD4 T Cells','NK Cells','CD48','CD244'),('Non-classical Monocytes','CD4 T Cells','SELPLG','ITGB2'),('Classical Monocytes','NK Cells','TGFB1','ITGB1'),('Non-classical Monocytes','CD4 T Cells','CD48','CD2'),('B Cells','Non-classical Monocytes','CD52','SIGLEC10'),('CD8 T Cells','Non-classical Monocytes','CD52','SIGLEC10')]
    rows=[]
    for key in keys:
        row=model.set_index('interaction_id').loc['|'.join(key)]
        ped=row.mean_cSLE-row.mean_cHD;adult=row.mean_aSLE-row.mean_aHD
        pattern='Increased in children but decreased in adults' if ped>0 and adult<0 else 'Decreased in children but increased in adults' if ped<0 and adult>0 else 'Different SLE-associated changes between cohorts'
        rows.append({'Source to target':f'{key[0]} to {key[1]}','Ligand-receptor':f'{key[2]}-{key[3]}','Coefficient':f'{row.interaction_coef:.3f}','95% CI':f'{row.interaction_ci95_low:.3f} to {row.interaction_ci95_high:.3f}','FDR':f'{row.fdr_q_value:.5f}','Pediatric healthy to SLE':f'{row.mean_cHD:.3f} to {row.mean_cSLE:.3f}','Adult healthy to SLE':f'{row.mean_aHD:.3f} to {row.mean_aSLE:.3f}','Plain-language pattern':pattern})
    pd.DataFrame(rows).to_csv(TABLE/'table2_representative_interactions.csv',index=False)
    # Primary-relative comparisons, not original-relative denominators.
    comparisons=[]
    for name,file in [('original_exclusions','baseline_models.csv'),('alternative_resource','combined_models.csv')]:
        other=pd.read_csv(BASE/'Results/annotation_resource_sensitivity'/file)
        matched=model.merge(other,on='interaction_id',suffixes=('_primary','_sensitivity'))
        matched.to_csv(OUT/f'{name}_comparison.csv',index=False)
        hit=matched[matched.fdr_q_value_primary<.05]
        comparisons.append({'comparison':name,'tested':len(other),'fdr05':int(other.fdr_q_value.lt(.05).sum()),'shared_tests':len(matched),'coefficient_r':matched.interaction_coef_primary.corr(matched.interaction_coef_sensitivity),'primary_significant_testable':len(hit),'primary_significant_still_significant':int(hit.fdr_q_value_sensitivity.lt(.05).sum()),'primary_significant_not_testable':len(sig)-len(hit)})
    pd.DataFrame(comparisons).to_csv(OUT/'sensitivity_summary.csv',index=False)
    summary={'retained_cells':int(counts.n_cells.sum()),'tested':len(model),'fdr05':len(sig),'fdr10':int(model.fdr_q_value.lt(.1).sum()),'positive':int(sig.interaction_coef.gt(0).sum()),'negative':int(sig.interaction_coef.lt(0).sum()),'female_tested':len(female),'female_fdr05':int(female.fdr_q_value.lt(.05).sum()),'female_fdr10':int(female.fdr_q_value.lt(.1).sum()),'female_r':compare.interaction_coef_primary.corr(compare.interaction_coef_female),'female_overlap':int((compare.fdr_q_value_primary.lt(.05)&compare.fdr_q_value_female.lt(.05)).sum()),'min_donors':int(model.n_observations.min()),'max_donors':int(model.n_observations.max()),'tests_fewer_than_56':int(model.n_observations.lt(56).sum())}
    (OUT/'summary.json').write_text(json.dumps(summary,indent=2));print(summary,flush=True)
    fig,ax=plt.subplots(figsize=(7,5))
    ax.scatter(model.interaction_coef,-np.log10(model.fdr_q_value.clip(lower=1e-300)),c=np.where(model.fdr_q_value<.05,'#c94242','#aaaaaa'),s=12,alpha=.75)
    ax.axhline(-np.log10(.05),c='black',ls='--',lw=.8);ax.axvline(0,c='grey',lw=.7)
    ax.set(xlabel='Age-group-by-SLE interaction coefficient',ylabel='−log10(FDR)',title='Primary analysis with reviewed T-cell restoration')
    savefig(fig,'figure2_primary_age_diseasestatus_volcano')
    fig,ax=plt.subplots(figsize=(6,5))
    ax.scatter(compare.interaction_coef_primary,compare.interaction_coef_female,c=np.where(compare.fdr_q_value_primary<.05,'#235789','#aaaaaa'),s=10,alpha=.7)
    limits=[min(compare.interaction_coef_primary.min(),compare.interaction_coef_female.min()),max(compare.interaction_coef_primary.max(),compare.interaction_coef_female.max())]
    ax.plot(limits,limits,'k--',lw=.8);ax.axhline(0,c='grey',lw=.6);ax.axvline(0,c='grey',lw=.6)
    ax.set(xlabel='All-donor interaction coefficient',ylabel='Female-only interaction coefficient',title=f"Female-only sensitivity (r = {summary['female_r']:.3f})")
    savefig(fig,'figure3_primary_vs_female_only_coefficients')
    umap=load_module('umap_style','10_GENERATE_MANUSCRIPT_UMAP.py')
    fig,axes=plt.subplots(1,2,figsize=(11.5,5.2))
    for ax,(cohort,xy,labels) in zip(axes,datasets):
        for ct in sorted(set(labels),key=lambda c:np.sum(labels==c),reverse=True):
            select=labels==ct;ax.scatter(xy[select,0],xy[select,1],s=.45,c=umap.COLORS[ct],alpha=.72,linewidths=0,rasterized=True)
        ax.set(title='A. Pediatric cohort' if cohort=='child' else 'B. Adult cohort',xlabel='UMAP 1',ylabel='UMAP 2',xticks=[],yticks=[])
        ax.spines[['top','right','bottom','left']].set_visible(False)
    handles=[Line2D([0],[0],marker='o',ls='',color=umap.COLORS[ct],label=ct) for ct in umap.CELL_TYPE_ORDER]
    fig.legend(handles=handles,loc='lower center',ncol=5,fontsize=9);fig.tight_layout(rect=[0,.14,1,1])
    savefig(fig,'figure1_pediatric_adult_celltype_umap')
    tables=load_module('table_style','09_GENERATE_MANUSCRIPT_SUMMARIES.py');tables.MANUSCRIPT=TABLE
    tables.render_main_tables()
    (OUT/'provenance.json').write_text(json.dumps({'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'annotation_source':'script 14 conservative restoration, reviewed by script 16','primary_resource':'LIANA 1.8.1 consensus','original_exclusions':'sensitivity comparison','retrospective_revision':True,'new_inference_needed':False,'explanation':'Exact reviewed label set matches completed annotation scenario; all donor scores reused, primary models refit and checked, female models newly fit.','umap':'Archived pre-exclusion cohort coordinates; no new embedding fit.'},indent=2))
    manifests={str(p.relative_to(BASE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths+[helper.LOCAL/'child_cell_review.csv.gz',helper.LOCAL/'adult_cell_review.csv.gz',BASE/'Results/annotation_resource_sensitivity/consensus_resource.csv',BASE/'Results/annotation_resource_sensitivity/classifier_provenance.json']}
    (OUT/'input_sha256.json').write_text(json.dumps(manifests,indent=2))
    if args.publish:
        backup=BASE/'_local_submission/pre_primary_revision'
        for folder in ['manuscript_figures','manuscript_tables']:
            source=OUT/folder;target=BASE/'Results'/folder;old=backup/folder
            old.mkdir(parents=True,exist_ok=True);target.mkdir(parents=True,exist_ok=True)
            for p in target.iterdir():
                if p.is_file() and not (old/p.name).exists():shutil.copy2(p,old/p.name)
            for p in source.iterdir():
                if p.is_file():shutil.copy2(p,target/p.name)
        # Old severity summaries are not part of the revised primary manuscript.
        stale=BASE/'Results/manuscript_tables/secondary_sledai_summary.csv'
        if stale.exists():
            archived=BASE/'Results/original_primary_presentation';archived.mkdir(exist_ok=True)
            if not (archived/stale.name).exists():shutil.copy2(stale,archived/stale.name)
            stale.unlink()
        journal=BASE/'Journal_Submission/Figures';journal.mkdir(parents=True,exist_ok=True)
        for i,stem in enumerate(['figure1_pediatric_adult_celltype_umap','figure2_primary_age_diseasestatus_volcano','figure3_primary_vs_female_only_coefficients'],1):
            dest=journal/f'Figure_{i}.jpg';old=backup/'Journal_Submission/Figures'/dest.name
            old.parent.mkdir(parents=True,exist_ok=True)
            if dest.exists() and not old.exists():shutil.copy2(dest,old)
            shutil.copy2(FIG/f'{stem}.jpg',dest)
        # Keep a single visible research presentation set after publication.
        # Delete only verified duplicates inside these explicitly owned folders.
        for folder in ['manuscript_figures','manuscript_tables']:
            source=(OUT/folder).resolve();target=(BASE/'Results'/folder).resolve()
            assert source.is_relative_to(OUT.resolve()) and target.is_relative_to(BASE.resolve())
            for p in source.iterdir():
                if p.is_file():
                    assert (target/p.name).is_file() and p.read_bytes()==(target/p.name).read_bytes()
                    p.unlink()
            source.rmdir()
        provenance_dir=BASE/'Results/provenance'
        bundle=provenance_dir/'run_manifests.json'
        records=json.loads(bundle.read_text(encoding='utf-8')) if bundle.exists() else {}
        individual=sorted(p for p in provenance_dir.glob('*.json') if p.name[:2].isdigit())
        for p in individual:records[p.name]=json.loads(p.read_text(encoding='utf-8'))
        bundle.write_text(json.dumps(records,indent=2,ensure_ascii=False),encoding='utf-8')
        confirmed=json.loads(bundle.read_text(encoding='utf-8'))
        for p in individual:
            assert confirmed[p.name]==json.loads(p.read_text(encoding='utf-8'))
            p.unlink()

if __name__=='__main__':main()
