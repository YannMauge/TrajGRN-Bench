# Methods Input/Output Format

The method I/O contract is documented in
[`docs/methods_io.md`](docs/methods_io.md) (also on the
[docs site](https://yannmauge.github.io/TrajGRN-Bench/methods_io/)), and the
step-by-step integration guide in
[`docs/development_add_method.md`](docs/development_add_method.md).

In short:

- **Call**: `python <script> -i <inputs>/data_<i>.h5ad -o <out>/ [-u output_mode] [-k ko_genes] [-p perturbation_training]`
- **Inputs**: `data_<i>.h5ad` (WT + training perturbations), `train_tps.npy`,
  `test_tps.npy`, optional `subsample_train_ids.npy`, next to the input file.
- **Outputs**: `data_<i>_GRN.npy` — 2-D, `(n_genes+1, n_genes+1)`, Stimulus at
  index 0, `grn[regulator, target]`; `data_<i>_adata.h5ad` — predicted cells
  with `obs['time']`; `data_<i>_ko_<gene>_adata.h5ad` for KO predictions.
- **Helpers**: `methods/_common.py` (loading/filtering, normalisation, GRN saving, seeding).
- **Environment**: a pixi environment in `pixi.toml` → Docker image `benchmark/<pixi_env>`.
