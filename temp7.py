import numpy as np
import json

counts = json.load(open("test/human/count.npy"))
matrix = np.load("test/human/matrix.npz")

matrix["matrix"]
matrix["gene_ids"]
matrix["codon_order"]