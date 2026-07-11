#!/usr/bin/env python3
"""Stage 4D: audit public GSE241124 H5 matrices against public spot metadata."""
from pathlib import Path
import csv
import gzip
import h5py
import pandas as pd

root = Path(".").resolve()
h5_dir = root / "05-results/intermediate_references/GSE241124_h5"
metadata_path = root / "02-data-raw/GEO/GSE241124/metadata/GSE241124_spatialseq_metadata_acutewound.txt.gz"
soft_path = root / "02-data-raw/GEO/GSE241124/metadata/GSE241124_family.soft.gz"
out = root / "04-tables/GSE241124_stage4D_spatial_matrix_audit.csv"

soft = gzip.open(soft_path, "rt", encoding="utf-8", errors="replace").read()
gsm_to_title = {}
for block in soft.split("^SAMPLE = ")[1:]:
    lines = block.splitlines()
    gsm = lines[0].strip()
    title = next((line.split(" = ", 1)[1] for line in lines if line.startswith("!Sample_title = ")), "")
    gsm_to_title[gsm] = title.replace("-", "_")

meta = pd.read_csv(metadata_path, sep="\t")
rows = []
for path in sorted(h5_dir.glob("*.h5")):
    gsm = path.name.split("_", 1)[0]
    sample_name = gsm_to_title[gsm]
    with h5py.File(path, "r") as handle:
        matrix = handle["matrix"]
        barcodes = [item.decode() for item in matrix["barcodes"][:]]
        n_genes, n_spots = map(int, matrix["shape"][:])
        nnz = int(len(matrix["data"]))
    expected = {f"{sample_name}_{barcode}" for barcode in barcodes}
    sample_meta = meta.loc[meta["Sample_name"] == sample_name]
    rows.append({
        "GSM": gsm,
        "sample_name": sample_name,
        "condition": sample_meta["Condition"].iloc[0] if len(sample_meta) else "",
        "donor": sample_meta["Donor"].iloc[0] if len(sample_meta) else "",
        "seq_batch": sample_meta["Seq_Batch"].iloc[0] if len(sample_meta) else "",
        "genes": n_genes,
        "h5_spots": n_spots,
        "metadata_spots": len(sample_meta),
        "barcode_metadata_matches": len(expected.intersection(set(sample_meta["barcode"]))),
        "matrix_nonzero_entries": nnz,
    })

with out.open("w", newline="", encoding="utf-8") as fh:
    writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)

if not all(row["metadata_spots"] == row["barcode_metadata_matches"] for row in rows):
    raise RuntimeError("Some public annotated spots do not map back to H5 barcodes")
print(f"GSE241124 Stage 4D passed: {len(rows)} matrices; the public metadata annotates a subset of each H5 matrix's spots.")
