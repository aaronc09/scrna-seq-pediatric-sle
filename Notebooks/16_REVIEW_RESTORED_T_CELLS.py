"""Review proposed T-cell restoration without using communication outcomes.

Diagnostic flags are not new exclusion thresholds. Reads archived matrices and
script-14 predictions; writes a donor-aware evidence review, not final labels.
"""
from pathlib import Path
import argparse
import importlib.util
import json
import hashlib
import numpy as np
import pandas as pd
import h5py

BASE = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('sensitivity', BASE/'Notebooks/14_ANNOTATION_RESOURCE_SENSITIVITY.py')
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
OUT = BASE/'Results/annotation_review'
PANELS = {
    'T': ['CD3D','CD3E','TRAC'],
    'CD4_naive_memory': ['CD4','IL7R','CCR7','TCF7','LEF1'],
    'CD8': ['CD8A','CD8B'],
    'cytotoxic': ['NKG7','CCL5','GZMK','PRF1','GNLY'],
    'B': ['MS4A1','CD79A','CD79B'],
    'myeloid': ['LYZ','LST1','FCN1','S100A8','S100A9'],
    'platelet': ['PPBP','PF4'],
    'erythroid': ['HBB','HBA1','HBA2'],
}
QC = ['n_genes_by_counts','total_counts','pct_counts_mt','doublet_score',
      'ribosomal_percent_retained_genes','reference_probability']

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--archive-root',type=Path,required=True)
    args=parser.parse_args(); OUT.mkdir(parents=True,exist_ok=True)
    results=[]; donor_results=[]; markers=[]; quantiles=[]; availability=[]; gates=[]
    for cohort, filename in module.FILES.items():
        review=pd.read_csv(module.LOCAL/f'{cohort}_cell_review.csv.gz',index_col='cell_id',dtype={'cluster':str})
        review['review_group']=np.select(
            [review.restore,review.original_label.isin(['CD4 T Cells','CD8 T Cells']),review.original_label.eq('Ribosomal/Low-quality')],
            ['restored','original_T','excluded_unrestored'],default='other')
        with h5py.File(args.archive_root/f'{cohort}_individual_h5ad'/filename,'r') as f:
            obs,genes=module.metadata(f,cohort);assert obs.index.equals(review.index)
            lookup={g:i for i,g in enumerate(genes)}
            wanted=list(dict.fromkeys(g for panel in PANELS.values() for g in panel))
            present=[g for g in wanted if g in lookup]
            for gene in wanted: availability.append({'cohort':cohort,'gene':gene,'available':gene in lookup})
            expression=np.empty((len(obs),len(present)),dtype=np.float32)
            for start in range(0,len(obs),5000):
                end=min(start+5000,len(obs))
                expression[start:end]=module.slice_csr(f['X'],start,end,len(genes))[:,[lookup[g] for g in present]].toarray()
            expr=pd.DataFrame(expression,index=review.index,columns=present)
        for panel, names in PANELS.items():
            available=[g for g in names if g in expr]
            review[f'{panel}_detected_genes']=(expr[available]>0).sum(axis=1)
        review['nonT_multimarker_flag']=review[[f'{p}_detected_genes' for p in ['B','myeloid','platelet','erythroid']]].ge(2).any(axis=1)
        review['CD8_both_detected']=review.CD8_detected_genes.eq(2)
        review['same_restoration_gate']=review.reference_broad.ne('unresolved') & review.reference_probability.ge(.5) & review.CD3D_detected & review.CD3E_detected
        review['p08']=review.reference_probability.ge(.8)
        cols=QC+['nonT_multimarker_flag','CD8_both_detected','same_restoration_gate','p08']
        for (group,label), part in review.groupby(['review_group','sensitivity_label'],observed=True):
            row={'cohort':cohort,'review_group':group,'label':label,'cells':len(part),'donors':part['sample'].nunique()}
            row.update({f'{c}_median':part[c].median() for c in QC})
            row.update({f'{c}_fraction':part[c].mean() for c in cols[len(QC):]})
            results.append(row)
            for gene in present:
                values=expr.loc[part.index,gene]
                markers.append({'cohort':cohort,'review_group':group,'label':label,'gene':gene,'cells':len(part),'fraction_detected':(values>0).mean(),'mean_log_expression':values.mean()})
            for metric in QC:
                for q,value in part[metric].quantile([.1,.25,.5,.75,.9,.95]).items():
                    quantiles.append({'cohort':cohort,'review_group':group,'label':label,'metric':metric,'quantile':q,'value':value})
        for (sample,group,label),part in review.groupby(['sample','review_group','sensitivity_label'],observed=True):
            row={'cohort':cohort,'sample':sample,'review_group':group,'label':label,'cells':len(part)}
            row.update({f'{c}_median':part[c].median() for c in QC})
            row.update({f'{c}_fraction':part[c].mean() for c in cols[len(QC):]})
            donor_results.append(row)
        for (label,ref),part in review[review.review_group.eq('original_T')].groupby(['original_label','reference_broad']):
            gates.append({'cohort':cohort,'original_label':label,'reference_broad':ref,'cells':len(part),'pass_restoration_gate':int(part.same_restoration_gate.sum())})
        print(cohort,'review complete',len(review),flush=True)
    for name,rows in [('group_review',results),('donor_review',donor_results),('marker_review',markers),('qc_quantiles',quantiles),('marker_availability',availability),('retained_T_gate_check',gates)]:
        pd.DataFrame(rows).to_csv(OUT/f'{name}.csv',index=False)
    (OUT/'provenance.json').write_text(json.dumps({'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'panels':PANELS,'purpose':'Outcome-independent diagnostic review; flags are not exclusion rules.','communication_outcomes_used':False,'source_objects_modified':False},indent=2))
    print(pd.DataFrame(results).query("review_group in ['restored','original_T']").to_string(index=False))

if __name__=='__main__':main()
