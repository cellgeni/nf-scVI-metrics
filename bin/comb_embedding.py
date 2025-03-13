#!/usr/bin/env python3

import sys
import numpy as np
import anndata as ad

adata = ad.read_h5ad(sys.argv[1])

for i in sorted(sys.argv[2].split()):
    adata.obsm[i] = np.load(i)

adata.write_h5ad("adata_embedding.h5ad")