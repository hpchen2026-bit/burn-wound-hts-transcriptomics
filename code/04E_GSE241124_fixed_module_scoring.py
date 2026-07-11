#!/usr/bin/env python3
"""Stage 4E: descriptive fixed-module scoring in annotated GSE241124 spatial spots."""
from pathlib import Path
import csv, gzip
import h5py
import numpy as np
import pandas as pd
from scipy.sparse import csc_matrix

root = Path(".").resolve()
h5_dir = root / "05-results/intermediate_references/GSE241124_h5"
meta_path = root / "02-data-raw/GEO/GSE241124/metadata/GSE241124_spatialseq_metadata_acutewound.txt.gz"
soft_path = root / "02-data-raw/GEO/GSE241124/metadata/GSE241124_family.soft.gz"
modules_path = root / "04-tables/stage4_fixed_reactome_module_genes.csv"
out_spot = root / "04-tables/GSE241124_stage4E_annotated_spot_module_scores.csv"
out_summary = root / "04-tables/GSE241124_stage4E_module_scores_by_condition_annotation.csv"

soft = gzip.open(soft_path, "rt", encoding="utf-8", errors="replace").read()
gsm_to_title = {}
for block in soft.split("^SAMPLE = ")[1:]:
    lines = block.splitlines()
    gsm = lines[0].strip()
    title = next((line.split(" = ", 1)[1] for line in lines if line.startswith("!Sample_title = ")), "")
    gsm_to_title[gsm] = title.replace("-", "_")
metadata = pd.read_csv(meta_path, sep="\t")
module_genes = pd.read_csv(modules_path).groupby("module")["gene_symbol"].apply(set).to_dict()
rows = []
for path in sorted(h5_dir.glob("*.h5")):
    gsm = path.name.split("_", 1)[0]
    sample_name = gsm_to_title[gsm]
    with h5py.File(path, "r") as handle:
        mat = handle["matrix"]
        genes = np.array([x.decode() for x in mat["features"]["name"][:]])
        barcodes = np.array([x.decode() for x in mat["barcodes"][:]])
        counts = csc_matrix((mat["data"][:], mat["indices"][:], mat["indptr"][:]), shape=tuple(mat["shape"][:]))
    library_size = np.asarray(counts.sum(axis=0)).ravel()
    scores = {}
    for module, geneset in module_genes.items():
        idx = np.flatnonzero(np.isin(genes, list(geneset)))
        if len(idx) < 5:
            raise RuntimeError(f"Insufficient feature match for {module} in {path.name}")
        values = counts[idx, :].toarray()
        scores[module] = np.log1p(values * 10000.0 / np.maximum(library_size, 1)).mean(axis=0)
    sample_meta = metadata.loc[metadata["Sample_name"] == sample_name].copy()
    index = {f"{sample_name}_{barcode}": i for i, barcode in enumerate(barcodes)}
    sample_meta["matrix_index"] = sample_meta["barcode"].map(index)
    if sample_meta["matrix_index"].isna().any():
        raise RuntimeError(f"Annotated spot mapping failed for {sample_name}")
    for module, values in scores.items():
        sample_meta[module] = [float(values[int(i)]) for i in sample_meta["matrix_index"]]
    rows.append(sample_meta)
spot = pd.concat(rows, ignore_index=True)
spot.to_csv(out_spot, index=False)
summary = spot.groupby(["Sample_name", "Condition", "Donor", "Seq_Batch", "AnnoType"], as_index=False).agg(
    annotated_spots=("barcode", "size"), **{m: (m, "mean") for m in module_genes}
)
summary.to_csv(out_summary, index=False)
print(f"GSE241124 Stage 4E complete: {len(spot)} annotated spots scored across {spot.Sample_name.nunique()} samples.")
