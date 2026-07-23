## Joint GRN + Trajectory Methods

These methods jointly infer both the gene regulatory network and cellular trajectories from time-course scRNA-seq data.

### Reference Fitting

[:fontawesome-brands-github: Repository](https://github.com/zsteve/referencefitting){ .md-button .md-button--primary }

Joint trajectory and network inference via reference fitting; leverages dynamical and perturbational single-cell data to jointly infer trajectories and directed/signed networks using min-entropy estimation.

- **Publication:** arXiv:2409.06879 (2024)
- **GRN:** :material-check:{ .green } **Trajectory:** :material-check:{ .green } **Perturbation:** :material-check:{ .green }

---

### CardamomOT

[:fontawesome-brands-github: Repository](https://github.com/eliasventre/CardamomOT){ .md-button .md-button--primary }

Mechanistic GRN inference that fits stochastic models to time-course scRNA-seq and then simulates the calibrated model. Designed to produce interpretable kinetic parameters and synthetic trajectories for benchmarking.

- **Publication:** bioRxiv (2026) DOI:10.64898/2026.03.31.715390
- **GRN:** :material-check:{ .green } **Trajectory:** :material-check:{ .green } **Perturbation:** :material-check:{ .green }

---

### RENGE

[:fontawesome-brands-github: Repository](https://github.com/masastat/RENGE){ .md-button .md-button--primary }

REgulatory Network inference using GEne perturbation data; infers GRNs from time-series single-cell data using a kinetic model with hyperparameter optimization. Optionally leverages CRISPR perturbation (KO) data when available.

- **Publication:** Communications Biology (2023) DOI:10.1038/s42003-023-05594-4
- **GRN:** :material-check:{ .green } **Trajectory:** :material-check:{ .green } **Perturbation:** :material-check:{ .green }

---

### FLeCS

[:fontawesome-brands-github: Repository](https://github.com/Bertinus/FLeCS){ .md-button .md-button--primary }

Functional and Learnable Cell dynamicS — a scalable gene-network-based differential-equation model that incorporates network structure to infer regulatory dynamics and simulate single-cell trajectories.

- **Publication:** arXiv:2503.20027 (2025)
- **GRN:** :material-check:{ .green } **Trajectory:** :material-check:{ .green } **Perturbation:** :material-close:{ .red }
- **Variant:** `FLeCS-TPs`

---

### FLeCS-TPs

[:fontawesome-brands-github: Repository](https://github.com/Bertinus/FLeCS){ .md-button .md-button--primary }

FLeCS variant that forces timepoint alignment during training.

- **Publication:** arXiv:2503.20027 (2025)
- **GRN:** :material-check:{ .green } **Trajectory:** :material-check:{ .green } **Perturbation:** :material-close:{ .red }

---

## GRN-Only Methods



### GENIE3

[:fontawesome-brands-github: Repository](https://github.com/aertslab/GENIE3){ .md-button .md-button--primary }

Tree-ensemble regression approach (Random Forests) for gene regulatory network inference from expression data; widely used baseline for GRN reconstruction.

- **Publication:** PLoS ONE (2010) DOI:10.1371/journal.pone.0012776
- **GRN:** :material-check:{ .green } **Trajectory:** :material-close:{ .red } **Perturbation:** :material-close:{ .red }

---

### Pearson (coexpression)

Local script: `methods/pearson/inference_pearson.py`

Computes pairwise Pearson correlation between genes as a naive GRN scoring baseline.

- **Publication:** —
- **GRN:** :material-check:{ .green } **Trajectory:** :material-close:{ .red } **Perturbation:** :material-close:{ .red }

---

## Trajectory-Only Methods



### scNODE

[:fontawesome-brands-github: Repository](https://github.com/rsinghlab/scNODE){ .md-button .md-button--primary }

A VAE + Neural ODE generative model that learns continuous-time dynamics on a latent manifold to predict and simulate single-cell transcriptomic trajectories.

- **Publication:** bioRxiv 2023 DOI:10.1101/2023.11.22.568346
- **GRN:** :material-close:{ .red } **Trajectory:** :material-check:{ .green } **Perturbation:** :material-close:{ .red }

---

### TrajectoryNet

[:fontawesome-brands-github: Repository](https://github.com/krishnaswamylab/TrajectoryNet){ .md-button .md-button--primary }

Continuous normalizing-flow / optimal-transport-based model that learns continuous paths between distributions for modeling cellular dynamics and sampling predictions at arbitrary timepoints.

- **Publication:** ICML 2020, arXiv:2002.04461
- **GRN:** :material-close:{ .red } **Trajectory:** :material-check:{ .green } **Perturbation:** :material-close:{ .red }

---

### Waddington-OT

[:fontawesome-brands-github: Repository](https://github.com/broadinstitute/wot){ .md-button .md-button--primary }

Optimal-transport-based approach (Waddington-OT) that computes transport maps between timepoints and interpolates gene-expression clouds to reconstruct trajectories.

- **Publication:** Cell (2019) DOI:10.1016/j.cell.2019.01.006
- **GRN:** :material-close:{ .red } **Trajectory:** :material-check:{ .green } **Perturbation:** :material-close:{ .red }