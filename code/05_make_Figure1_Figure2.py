#!/usr/bin/env python3
"""Generate Figure 1 and Figure 2 from locked, audited project tables."""
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

root = Path('.').resolve()
tables = root / '04-tables'
outdir = root / '05-results/figures/main'
outdir.mkdir(parents=True, exist_ok=True)

# ---------- Figure 1: study design ----------
fig, ax = plt.subplots(figsize=(12, 7))
ax.set_xlim(0, 12); ax.set_ylim(0, 8); ax.axis('off')

def box(x, y, w, h, text, fc, ec='#253746', fs=10):
    patch = FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0.03,rounding_size=0.1', facecolor=fc, edgecolor=ec, linewidth=1.3)
    ax.add_patch(patch)
    ax.text(x+w/2, y+h/2, text, ha='center', va='center', fontsize=fs, wrap=True)

def arrow(x1,y1,x2,y2):
    ax.annotate('', xy=(x2,y2), xytext=(x1,y1), arrowprops=dict(arrowstyle='->', lw=1.4, color='#253746'))

box(0.4, 5.5, 3.0, 1.55, 'Primary bulk cohort\nGSE178411 · 108 human specimens\nnormal skin 24 | early wound 22\nlate wound 29 | HTS 28', '#DDEBF7', fs=10)
box(4.45, 5.55, 3.0, 1.45, 'Subject-aware bulk models\nTMM + limma-voom\nage/sex and non-E sensitivity checks', '#FFF2CC', fs=10)
box(8.55, 5.55, 3.0, 1.45, 'Frozen stage-associated programs\nECM/collagen · inflammation\ncell cycle', '#E2F0D9', fs=10)
arrow(3.4,6.25,4.45,6.25); arrow(7.45,6.25,8.55,6.25)
box(1.1, 2.5, 3.2, 1.5, 'Human HTS scRNA reference\nGSE156326\n3 normal skin + 3 scar samples\ncell-source localization only', '#FCE4D6', fs=10)
box(7.7, 2.5, 3.2, 1.5, 'Human wound spatial reference\nGSE241124\n4 donors × Skin/Wound1/Wound7/Wound30\ncontext only', '#EDEDED', fs=10)
arrow(10.05,5.55,9.3,4.0); arrow(10.05,5.55,2.7,4.0)
box(3.75, 0.55, 4.5, 1.1, 'Interpretation boundary\nCross-sectional stage association; no individual trajectory, prediction, causal mechanism or therapeutic claim', '#F4CCCC', fs=10)
arrow(2.7,2.5,5.2,1.65); arrow(9.3,2.5,6.8,1.65)
ax.text(0.4,7.55,'Figure 1. Study design and prespecified evidence chain', fontsize=14, fontweight='bold')
fig.tight_layout(); fig.savefig(outdir/'Figure1_study_design.png', dpi=300, bbox_inches='tight'); plt.close(fig)

# ---------- Figure 2: PCA and sensitivity concordance ----------
qc = pd.read_csv(tables / 'GSE178411_stage1_QC_metrics.csv')
concord = pd.read_csv(tables / 'GSE178411_stage2_sensitivity_concordance.csv')
label = {'normal_skin':'Normal skin','early_wound':'Early wound','late_wound':'Late wound','hypertrophic_scar':'HTS'}
colors = {'normal_skin':'#E76F51','early_wound':'#8AB17D','late_wound':'#2A9D8F','hypertrophic_scar':'#7B2CBF'}
fig, axes = plt.subplots(1,2,figsize=(13,5.2),gridspec_kw={'width_ratios':[1.15,1]})
for group in ['normal_skin','early_wound','late_wound','hypertrophic_scar']:
    d = qc.loc[qc.analysis_group == group]
    axes[0].scatter(d.PC1, d.PC2, s=36, color=colors[group], alpha=.85, label=f'{label[group]} (n={len(d)})', edgecolor='none')
axes[0].axhline(0,color='#BBBBBB',lw=.7); axes[0].axvline(0,color='#BBBBBB',lw=.7)
axes[0].set_xlabel('PC1 (45.3% variance)'); axes[0].set_ylabel('PC2 (15.5% variance)')
axes[0].set_title('A  Bulk RNA-seq stage separation',loc='left',fontweight='bold')
axes[0].legend(frameon=False,fontsize=8,loc='best')

names = [x.replace('_vs_',' vs ').replace('hypertrophic_scar','HTS').replace('normal_skin','normal skin').replace('early_wound','early wound').replace('late_wound','late wound') for x in concord.contrast]
y = range(len(concord))
axes[1].scatter(concord['stage_only_vs_age_sex_logFC_r'], y, color='#457B9D', label='Age/sex-adjusted', s=55)
axes[1].scatter(concord['stage_only_vs_nonE_logFC_r'], y, color='#E76F51', label='Non-E restricted', s=55)
for i, row in concord.iterrows():
    axes[1].plot([row['stage_only_vs_age_sex_logFC_r'],row['stage_only_vs_nonE_logFC_r']],[i,i], color='#C7C7C7',lw=1)
axes[1].set_yticks(list(y)); axes[1].set_yticklabels(names,fontsize=8)
axes[1].set_xlim(.90,1.002); axes[1].set_xlabel('Correlation of gene-level logFC with stage-only model')
axes[1].set_title('B  Prespecified sensitivity concordance',loc='left',fontweight='bold')
axes[1].legend(frameon=False,fontsize=8,loc='lower left')
fig.tight_layout(); fig.savefig(outdir/'Figure2_bulk_stage_separation_and_robustness.png', dpi=300, bbox_inches='tight'); plt.close(fig)
print('Figure 1 and Figure 2 generated.')
