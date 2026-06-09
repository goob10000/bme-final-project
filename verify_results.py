import numpy as np
from pathlib import Path

viruses = ["HSV1", "HSV2", "ChickenPox_VZV", "Measles"]
results = []

for virus in viruses:
    npz_path = Path(f"Results2/{virus}/stuff.npz")
    if not npz_path.exists():
        results.append(f"{virus}: MISSING")
        continue
    
    file_size_mb = npz_path.stat().st_size / (1024 * 1024)
    npz = np.load(npz_path)
    
    keys = list(npz.files)
    matrix = npz['matrix']
    gene_ids = npz['gene_ids']
    
    results.append(
        f"{virus}: ✓ file={file_size_mb:.2f}MB, matrix_shape={matrix.shape}, "
        f"n_gene_ids={len(gene_ids)}, keys={keys}"
    )

for r in results:
    print(r)
