#!/usr/bin/env python3
"""Assemble manuscript source and DOCX from reviewed modular draft files."""
from pathlib import Path
import csv
import subprocess

root = Path('.').resolve()
man = root / '06-manuscript'
out_md = man / 'manuscript_full_draft_v0.1.md'
out_docx = man / 'manuscript_full_draft_v0.1.docx'

front = (man / 'front_matter_introduction_discussion_draft_v0.1.md').read_text()
results_methods = (man / 'results_methods_draft_v0.1.md').read_text()
legends = (man / 'figure_legends_v0.1.md').read_text()
references = (man / 'references_v0.1.md').read_text()

def between(text, start, end):
    return text.split(start, 1)[1].split(end, 1)[0].strip()

def csv_markdown(path):
    with path.open(newline='') as handle:
        rows = list(csv.reader(handle))
    header, data = rows[0], rows[1:]
    lines = ['| ' + ' | '.join(header) + ' |', '| ' + ' | '.join(['---'] * len(header)) + ' |']
    for row in data:
        lines.append('| ' + ' | '.join(row) + ' |')
    return '\n'.join(lines)

title = between(front, '## Title', '## Running title')
running = between(front, '## Running title', '## Abstract')
abstract = between(front, '## Abstract', '## Keywords')
keywords = between(front, '## Keywords', '## Introduction')
introduction = between(front, '## Introduction', '## Discussion')
discussion = front.split('## Discussion', 1)[1].strip()
methods = between(results_methods, '## Methods', '## Results')
results = results_methods.split('## Results', 1)[1].strip()
refs = references.split('\n', 1)[1].split('## Draft citation mapping', 1)[0].strip()
table1 = '''| Analysis group | Specimens, n | Unique source participants, n | Participants with >1 specimen within group, n | Median age, years | Female specimens, n | Male specimens, n | Specimens with paired normal skin, n |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Early wound | 22 | 19 | 3 | 27.0 | 6 | 16 | 11 |
| Late wound | 29 | 23 | 5 | 35.5 | 11 | 18 | 14 |
| Normal skin | 24 | 24 | 0 | 30.5 | 8 | 16 | 0 |
| Hypertrophic scar | 28 | 26 | 1 | 19.0 | 12 | 16 | 2 |'''
table2 = '''| Contrast | Stage-only vs age/sex-adjusted logFC correlation | Stage-only vs non-E-restricted logFC correlation | Shared FDR<0.05 genes in non-E analysis, n | Shared FDR<0.05 genes with same direction, n |
| --- | ---: | ---: | ---: | ---: |
| Early wound vs normal skin | 0.9999 | 0.9996 | 15,031 | 15,031 |
| Late wound vs normal skin | 0.9997 | 0.9990 | 14,930 | 14,930 |
| HTS vs normal skin | 0.9994 | 0.9719 | 10,260 | 10,260 |
| HTS vs late wound | 0.9981 | 0.9460 | 7,179 | 7,179 |'''
legend_body = legends.split('\n', 1)[1].strip()

body = f'''# {title.strip('*')}

**Running title:** {running.strip('*')}

**Authors:** Lining Chan, M.S.¹; Huiping Chen, M.D.¹*  
¹ Department of Plastic Surgery, Xinhua Hospital, Shanghai Jiao Tong University School of Medicine, No. 1665 Kangjiang Road, Yangpu District, Shanghai 200092, People’s Republic of China.  
*Corresponding author: Huiping Chen, M.D.; hpchen2005@163.com

# Abstract

{abstract}

**Keywords:** {keywords}

# Introduction

{introduction}

# Methods

{methods}

# Results

{results}

# Discussion

{discussion}

# Declarations

**Ethics approval:** This secondary analysis used publicly available, deidentified GEO datasets; no new participants were recruited and no identifiable individual-level data were accessed.

**Consent to participate/publication:** Not applicable.

**Competing interests:** The proposed declaration is that the authors have no known competing financial interests or personal relationships that could have influenced this work; confirm before submission.

**Funding:** The proposed declaration is no specific external grant funding; confirm before submission.

**Author contributions:** Lining Chan: Data curation, Investigation, Formal analysis, Visualization, Writing – original draft. Huiping Chen: Conceptualization, Methodology, Supervision, Project administration, Writing – review & editing. Confirm roles and final approval before submission.

**Data availability:** GSE178411, GSE156326 and GSE241124 are publicly available from GEO. Analysis scripts, manifests, checksums, derived tables and figures are available at https://github.com/hpchen2026-bit/burn-wound-hts-transcriptomics. A permanent DOI archive will be created before publication.

# References

{refs}

# Table 1. Primary analysis sample structure

{table1}

# Table 2. Sensitivity concordance of stage-associated gene-level effects

{table2}

# Figure legends

{legend_body}
'''
out_md.write_text(body)
subprocess.run(['pandoc', str(out_md), '-o', str(out_docx)], check=True)
print(f'Wrote {out_md.name} and {out_docx.name}')
