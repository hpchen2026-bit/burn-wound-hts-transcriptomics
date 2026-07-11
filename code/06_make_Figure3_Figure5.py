#!/usr/bin/env python3
"""Generate main Figure 3–5 from the locked program and external-reference tables."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

root = Path('.').resolve(); tables=root/'04-tables'; out=root/'05-results/figures/main'; out.mkdir(parents=True,exist_ok=True)

# Figure 3: stable Reactome program signatures.
stable = pd.read_csv(tables/'GSE178411_stage3_Reactome_stable_programs.csv')
selected = [
 'REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION','REACTOME_COLLAGEN_FORMATION','REACTOME_COLLAGEN_BIOSYNTHESIS_AND_MODIFYING_ENZYMES',
 'REACTOME_NEUTROPHIL_DEGRANULATION','REACTOME_INNATE_IMMUNE_SYSTEM','REACTOME_CELL_CYCLE_MITOTIC'
]
contrasts = ['early_wound_vs_normal_skin','late_wound_vs_normal_skin','hypertrophic_scar_vs_normal_skin','hypertrophic_scar_vs_late_wound']
labels = ['Early wound\nvs normal skin','Late wound\nvs normal skin','HTS\nvs normal skin','HTS\nvs late wound']
mat=[]; sigmat=[]
for term in selected:
    r=stable[stable.gs_name==term]
    vals=[]; sigs=[]
    for c in contrasts:
        x=r[r.contrast==c].iloc[0]
        f=max(float(x.FDR_stage_only),float(x.FDR_nonE_restricted))
        v=-np.log10(max(f,1e-300))
        is_stable = bool(x.stable_FDR05)
        if x.Direction_stage_only=='Down': v=-v
        vals.append(v); sigs.append(is_stable)
    mat.append(vals); sigmat.append(sigs)
mat=np.array(mat); sigmat=np.array(sigmat)
fig,ax=plt.subplots(figsize=(9,5.8)); im=ax.imshow(mat,cmap='RdBu_r',vmin=-20,vmax=20,aspect='auto')
ax.set_xticks(range(4),labels,fontsize=9); ax.set_yticks(range(6),['ECM organization','Collagen formation','Collagen biosynthesis','Neutrophil degranulation','Innate immune system','Mitotic cell cycle'],fontsize=10)
for i in range(mat.shape[0]):
 for j in range(mat.shape[1]):
  text=f'{mat[i,j]:.1f}' + ('*' if sigmat[i,j] else '')
  ax.text(j,i,text,ha='center',va='center',fontsize=8,color='white' if abs(mat[i,j])>11 else 'black')
cb=fig.colorbar(im,ax=ax,pad=.02); cb.set_label('Signed −log10(maximum FDR across two models); red = increased in named contrast')
ax.set_title('Figure 3. Reactome program direction and sensitivity robustness',loc='left',fontweight='bold')
fig.text(.01,.01,'* FDR < 0.05 with the same direction in both stage-only and non-E restricted models. HTS, hypertrophic scar.',fontsize=8)
fig.tight_layout();fig.savefig(out/'Figure3_stable_stage_programs.png',dpi=300,bbox_inches='tight');plt.close(fig)

# Figure 4: cell-source localization in human scar reference.
sc = pd.read_csv(tables/'GSE156326_stage4C_module_score_by_sample_compartment.csv')
sc=sc[(sc.cells>=30) & sc.reviewed_cell_compartment.str.contains('fibroblast')].copy()
order=['activated_fibroblast_POSTN_ADAM12','fibroblast_PCOLCE2_IGFBP6','fibroblast_PI16_ELN','fibroblast_CXCL12_APOE','fibroblast_CXCL14_CFD','fibroblast_CILP','fibroblast_TWIST1']
sc=sc[sc.reviewed_cell_compartment.isin(order)]; sc['compartment']=pd.Categorical(sc.reviewed_cell_compartment,order,ordered=True);sc=sc.sort_values('compartment')
short={'activated_fibroblast_POSTN_ADAM12':'POSTN/ADAM12','fibroblast_PCOLCE2_IGFBP6':'PCOLCE2/IGFBP6','fibroblast_PI16_ELN':'PI16/ELN','fibroblast_CXCL12_APOE':'CXCL12/APOE','fibroblast_CXCL14_CFD':'CXCL14/CFD','fibroblast_CILP':'CILP','fibroblast_TWIST1':'TWIST1'}
fig,ax=plt.subplots(figsize=(9,5.5)); pos={c:i for i,c in enumerate(order)}; rng=np.random.default_rng(20260711)
for cond,color in [('normal_skin','#457B9D'),('hypertrophic_scar','#D1495B')]:
 d=sc[sc.condition==cond]
 x=np.array([pos[v] for v in d.reviewed_cell_compartment])+rng.uniform(-.13,.13,len(d))
 ax.scatter(x,d.ECM_organization_score1,s=48,alpha=.85,color=color,label='Normal skin' if cond=='normal_skin' else 'HTS')
ax.set_xticks(range(len(order)),[short[x] for x in order],rotation=25,ha='right');ax.set_ylabel('ECM program score (descriptive)');ax.set_title('Figure 4. ECM program localizes to multiple human scar fibroblast compartments',loc='left',fontweight='bold');ax.legend(frameon=False);ax.axhline(0,color='#BDBDBD',lw=.7);fig.tight_layout();fig.savefig(out/'Figure4_fibroblast_ECM_localization.png',dpi=300,bbox_inches='tight');plt.close(fig)

# Figure 5: donor-level spatial trajectory.
sp=pd.read_csv(tables/'GSE241124_stage4E_module_scores_by_sample.csv'); order=['Skin','Wound1','Wound7','Wound30']; x=np.arange(4)
fig,axes=plt.subplots(1,2,figsize=(11,4.7),sharex=True)
for ax,module,title in zip(axes,['ECM_organization','collagen_formation'],['ECM program','Collagen-formation program']):
 for donor,d in sp.groupby('Donor'):
  d=d.set_index('Condition').reindex(order)
  ax.plot(x,d[module],marker='o',lw=1.8,label=donor)
 ax.set_xticks(x,order);ax.set_ylabel('Mean spot score within sample');ax.set_title(title,fontweight='bold');ax.axhline(0,color='#BDBDBD',lw=.7)
axes[0].legend(frameon=False,title='Donor');fig.suptitle('Figure 5. Human wound spatial reference: donor-level remodeling context',x=.1,ha='left',fontweight='bold');fig.tight_layout();fig.savefig(out/'Figure5_spatial_wound_remodeling_context.png',dpi=300,bbox_inches='tight');plt.close(fig)
print('Figure 3–5 generated.')
