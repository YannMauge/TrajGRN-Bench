import sys; sys.path += ['../']
import numpy as np
import anndata as ad
import getopt
import glob
import os
import re

import warnings
for module in ["anndata", "scipy"]:
    warnings.filterwarnings("ignore", module=module)


def _extract_dataset_id(file_path: str) -> str:
    filename = os.path.basename(file_path)
    match = re.match(r"^data_\d+_ko_([A-Za-z0-9_-]+)\.h5ad$", filename)
    if match:
        return match.group(1)
    return filename.replace(".h5ad", "")

def main(argv):
    try:
        opts, args = getopt.getopt(argv, "hi:g:d:o:n:", ["ifile=", "igene=", "degradation=", "output", "name"])
    except getopt.GetoptError:
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-o", "--output"):
            output_path = '{}'.format(arg)

    # Only process KO files (pattern: data_X_ko_*)
    file_list = glob.glob(f'{output_path}/data_*_ko_*.h5ad')
    if not file_list:
        raise FileNotFoundError(f"No KO .h5ad files found in {output_path}")
    
    # Group files by replicate number
    replicates = {}
    for file_path in file_list:
        filename = os.path.basename(file_path)
        # Extract replicate number (e.g., "1" from "data_1_ko_Gene0.h5ad")
        match = re.match(r"^data_(\d+)_ko_", filename)
        if match:
            replicate_num = match.group(1)
            if replicate_num not in replicates:
                replicates[replicate_num] = []
            replicates[replicate_num].append(file_path)
    
    # Merge files for each replicate
    for replicate_num in sorted(replicates.keys()):
        files = replicates[replicate_num]
        adata_list = []
        kos = set()
        all_uns = {}
        
        for file_path in files:
            adata = ad.read_h5ad(file_path)
            dataset_id = _extract_dataset_id(file_path)
            adata.obs['dataset_id'] = dataset_id
            if dataset_id != "WT":
                kos.add(dataset_id)
            # Preserve individual uns data
            all_uns.update(adata.uns)
            adata_list.append(adata)
        
        # Concatenate with join='outer' to preserve all obs/var columns
        ad_concat = ad.concat(adata_list, join='outer')
        print(ad_concat)
        
        # Merge all uns data and add kos list
        ad_concat.uns = all_uns
        ad_concat.uns['kos'] = sorted(kos)
        
        output_file = f'{output_path}/data_{replicate_num}.h5ad'
        ad_concat.write(output_file)
        print(f'Merged {len(files)} datasets into {output_file}')
    
    # Remove merged KO files
    for file_path in file_list:
        os.remove(file_path)
if __name__ == "__main__":
   main(sys.argv[1:])
