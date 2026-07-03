import sys; sys.path += ['../']
import numpy as np
import pandas as pd
import anndata as ad
import getopt

verb = 1

def main(argv):
    try:
        opts, args = getopt.getopt(argv, "hi:g:d:o:n:s:p:", ["ifile=", "igene=", "degradation=", "output", "name", "dataset_id", "simulation_name="])
    except getopt.GetoptError:
        sys.exit(2)
    dataset_id = "WT"
    simulation_name = "Harissa"
    for opt, arg in opts:
        if opt in ("-i", "--ifile"):
            fpath_panel = '{}'.format(arg)
        if opt in ("-g", "--igene"):
            fname_genes = '{}'.format(arg)
        if opt in ("-d", "--degradation"):
            fname_deg = '{}'.format(arg)
        if opt in ("-o", "--output"):
            output_path = '{}'.format(arg)
        if opt in ("-n", "--name"):
            fname_panel = '/{}'.format(arg)
        if opt in ("-s", "--dataset_id"):
            dataset_id = '{}'.format(arg)
        if opt in ("-p", "--simulation_name"):
            simulation_name = '{}'.format(arg)

    raw_matrix = np.loadtxt(fpath_panel, delimiter='\t').astype(np.int64)

    time = raw_matrix[0, 1:]
    data_rna = raw_matrix[1:, 1:].T 

    genes_df = pd.read_csv(fname_genes, sep='\t', header=None, names=['gene_id', 'gene_name'])
    assert data_rna.shape[1] == genes_df.shape[0], "Number of genes do not match between expression data and gene list."

    adata = ad.AnnData(X=data_rna)
    adata.var_names = genes_df['gene_name'].values.astype(str)

    adata.obs['time'] = time.astype(int)
    adata.obs['timepoint'] = time.astype(int)
    adata.obs['dataset_id'] = dataset_id
    degradations = np.loadtxt(fname_deg, delimiter='\t')

    adata.var['d0'] = degradations[:, 0]
    adata.var['d1'] = degradations[:, 1]
    adata.uns['simulation'] = simulation_name
    print(adata)
    adata.write(output_path+fname_panel+'.h5ad')


if __name__ == "__main__":
   main(sys.argv[1:])
