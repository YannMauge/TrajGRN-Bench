
# Methods included in this benchmark

- **CardamomOT**
  - GitHub: https://github.com/eliasventre/CardamomOT
  - Publications: bioRxiv (2026) DOI: [10.64898/2026.03.31.715390](https://doi.org/10.64898/2026.03.31.715390); earlier works In PLOS Comput Biol (2023) DOI: [10.1371/journal.pcbi.1010962](https://doi.org/10.1371/journal.pcbi.1010962) and Silico Biology (2021)
  - Description: Mechanistic GRN inference that fits stochastic models to time-course scRNA-seq and then simulates the calibrated model. Designed to produce interpretable kinetic parameters and synthetic trajectories for benchmarking.

- **FLeCS**
  - GitHub: https://github.com/Bertinus/FLeCS
  - Preprint: [arXiv:2503.20027](https://arxiv.org/abs/2503.20027) (2025) — preprint available
  - Description: Functional and Learnable Cell dynamicS — a scalable gene-network-based differential-equation model that incorporates network structure to infer regulatory dynamics and simulate single-cell trajectories.

- **Reference fitting**
  - GitHub: https://github.com/zsteve/referencefitting
  - Preprint: [arXiv:2409.06879](https://arxiv.org/abs/2409.06879) (2024)
  - Description: Joint trajectory and network inference via reference fitting; leverages dynamical and perturbational single-cell data to jointly infer trajectories and directed/signed networks using min-entropy estimation.

- **RENGE**
  - GitHub: https://github.com/masastat/RENGE
  - Paper: Communications Biology (2023) DOI: [10.1038/s42003-023-05594-4](https://doi.org/10.1038/s42003-023-05594-4)
  - Description: REgulatory Network inference using GEne perturbation data; infers GRNs from time-series single-cell data using a kinetic model with hyperparameter optimization. Optionally leverages CRISPR perturbation (KO) data when available. Predicts expression changes under gene perturbations.

- **scNODE**
  - GitHub: https://github.com/rsinghlab/scNODE
  - Preprint: bioRxiv 2023 DOI: [10.1101/2023.11.22.568346](https://doi.org/10.1101/2023.11.22.568346)
  - Description: A VAE + Neural ODE generative model that learns continuous-time dynamics on a latent manifold to predict and simulate single-cell transcriptomic trajectories.

- **TrajectoryNet**
  - GitHub: https://github.com/krishnaswamylab/TrajectoryNet
  - Paper / preprint: [arXiv:2002.04461](https://arxiv.org/abs/2002.04461) (TrajectoryNet: A Dynamic Optimal Transport Network for Modeling Cellular Dynamics; ICML 2020)
  - Description: Continuous normalizing-flow / optimal-transport-based model that learns continuous paths between distributions for modeling cellular dynamics and sampling predictions at arbitrary timepoints.

- **GENIE3**
  - GitHub (R package): https://github.com/aertslab/GENIE3
  - Paper: PLoS ONE (2010) DOI: [10.1371/journal.pone.0012776](https://doi.org/10.1371/journal.pone.0012776)
  - Description: Tree-ensemble regression approach (Random Forests) for gene regulatory network inference from expression data; widely used baseline for GRN reconstruction.

- **Pearson (coexpression baseline)**
  - Local script: `methods/pearson/inference_pearson.py`
  - Citation: /
  - Description: Computes pairwise Pearson correlation between genes as a naive GRN scoring baseline.

- **Waddington-OT (WOT)**
  - GitHub: https://github.com/broadinstitute/wot
  - Resources / documentation: https://broadinstitute.github.io/wot/
  - Citation: Cell (2019) DOI: [10.1016/j.cell.2019.01.006](https://doi.org/10.1016/j.cell.2019.01.006)
  - Description: Optimal-transport-based approach (Waddington-OT) that computes transport maps between timepoints and interpolates gene-expression clouds to reconstruct trajectories.

