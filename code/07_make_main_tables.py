#!/usr/bin/env python3
"""Generate manuscript-ready Table 1 and Table 2 from audited tables."""
from pathlib import Path
import pandas as pd

root=Path('.').resolve(); tables=root/'04-tables'
man=pd.read_csv(tables/'GSE178411_sample_manifest_v0.2.csv')
primary=man[man.primary_analysis_role=='primary_cross_sectional'].copy()
primary['age_num']=pd.to_numeric(primary.age,errors='coerce')
rows=[]
for group,d in primary.groupby('analysis_group',sort=False):
    rows.append({
        'analysis_group':group,
        'specimens':len(d),
        'unique_source_subjects':d.subject.nunique(),
        'subjects_with_multiple_specimens_within_group':int((d.groupby('subject').size()>1).sum()),
        'median_age_years':d.age_num.median(),
        'female_specimens':int((d.sex=='female').sum()),
        'male_specimens':int((d.sex=='male').sum()),
        'paired_normal_skin_available_for_specimen':int((d.paired_normal_skin_available.astype(str).str.lower()=='true').sum())
    })
t1=pd.DataFrame(rows)
t1.to_csv(tables/'Table1_GSE178411_primary_analysis_sample_structure.csv',index=False)

s=pd.read_csv(tables/'GSE178411_stage2_sensitivity_concordance.csv')
s['shared_FDR05_genes_same_direction']=s['same_direction_overlap']
s=s[['contrast','stage_only_vs_age_sex_logFC_r','stage_only_vs_nonE_logFC_r','FDR05_overlap_with_nonE','shared_FDR05_genes_same_direction']]
s.to_csv(tables/'Table2_GSE178411_sensitivity_concordance.csv',index=False)
print('Table 1 and Table 2 generated.')
