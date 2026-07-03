import getopt
import os
import re
import sys
from typing import List, Tuple

import anndata as ad
import warnings
for module in ["anndata", "scipy", "torchdiffeq"]:
    warnings.filterwarnings("ignore", module=module)


def _parse_gene_tokens(raw: str) -> List[str]:
    if not raw:
        return []
    token = str(raw).strip()
    if token.lower() in {"", "0", "none", "nan"}:
        return []
    return [g.strip() for g in token.split(",") if g.strip()]


def _extract_ko_ov_from_dataset_id(dataset_id: str) -> Tuple[List[str], List[str]]:
    did = str(dataset_id).strip()
    low = did.lower()
    if low in {"", "wt", "wildtype", "wild_type", "control", "ctrl", "none"}:
        return [], []

    kos: List[str] = []
    ovs: List[str] = []

    for match in re.findall(r"(?i)ko_([A-Za-z0-9_,.-]+)", did):
        kos.extend(_parse_gene_tokens(match))
    for match in re.findall(r"(?i)ov_([A-Za-z0-9_,.-]+)", did):
        ovs.extend(_parse_gene_tokens(match))

    if not kos and not ovs:
        kos = [did]

    return kos, ovs


def _normalize_genes(genes: List[str], gene_map_upper):
    normalized = []
    for gene in genes:
        g = str(gene).strip()
        if not g:
            continue
        matched = gene_map_upper.get(g.upper())
        if matched is not None and matched not in normalized:
            normalized.append(matched)
    return normalized


def _write_basal_ref_kov(
    data_dir: str,
    adata: ad.AnnData,
) -> None:
    out_path = os.path.join(data_dir, "basal_ref_KOV.txt")
    if "dataset_id" not in adata.obs:
        raise ValueError("perturbation_training=true requires adata.obs['dataset_id']")

    dataset_ids = sorted(
        [str(v) for v in adata.obs["dataset_id"].astype(str).unique()],
        key=lambda x: (x != "WT", x),
    )
    gene_map_upper = {str(g).upper(): str(g) for g in adata.var_names}

    print(f"[prepare_benchmark_perturbations] {adata.obs['dataset_id'].nunique()} unique dataset_id values found")
    with open(out_path, "w", encoding="utf-8") as f:
        f.write("sample_id\tKO\tOV\n")
        for dataset_id in dataset_ids:
            kos_raw, ovs_raw = _extract_ko_ov_from_dataset_id(dataset_id)
            kos = _normalize_genes(kos_raw, gene_map_upper)
            ovs = _normalize_genes(ovs_raw, gene_map_upper)
            f.write(f"{dataset_id}\t{','.join(kos)}\t{','.join(ovs)}\n")

    print(f"[prepare_benchmark_perturbations] Wrote {out_path}")


def _parse_ko_output_genes(raw: str, adata: ad.AnnData) -> List[str]:
    value = str(raw).strip()
    if value.lower() in {"", "none"}:
        return []

    allowed_genes = [str(g) for g in adata.var_names if str(g) != "Stimulus"]
    gene_map_upper = {g.upper(): g for g in allowed_genes}

    if value.lower() == "all":
        return allowed_genes

    requested = [g.strip() for g in value.split(",") if g.strip()]
    normalized = _normalize_genes(requested, gene_map_upper)
    missing = [g for g in requested if g.upper() not in gene_map_upper]
    if missing:
        raise ValueError(f"Unknown KO output gene(s): {missing}")
    return normalized


def _write_ko_ov_list(data_dir: str, ko_targets: List[str]) -> None:
    out_path = os.path.join(data_dir, "KO_OV_list.txt")
    with open(out_path, "w", encoding="utf-8") as f:
        f.write("KO\tOV\n")
        for gene in ko_targets:
            f.write(f"{gene}\t\n")
    print(f"[prepare_benchmark_perturbations] Wrote {out_path}")


def _safe_remove(path: str) -> None:
    if os.path.exists(path):
        os.remove(path)
        print(f"[prepare_benchmark_perturbations] Removed {path}")


def main(argv):
    inputfile = ""
    perturbation_training = False
    ko_output_genes = "none"
    try:
        opts, _ = getopt.getopt(
            argv,
            "hi:p:k:",
            ["input=", "perturbation_training=", "ko_genes="],
        )
    except getopt.GetoptError:
        print(
            "Usage: prepare_benchmark_perturbations.py "
            "-i <project_path> -p <true|false> -k <none|all|g1,g2>"
        )
        sys.exit(2)

    for opt, arg in opts:
        if opt in ("-i", "--input"):
            inputfile = str(arg)
        elif opt in ("-p", "--perturbation_training"):
            perturbation_training = arg.lower() == 'true'
        elif opt in ("-k", "--ko_genes"):
            ko_output_genes = str(arg)

    if not inputfile:
        raise ValueError("Missing required --input argument")

    data_dir = os.path.join(inputfile, "Data")
    data_path = os.path.join(data_dir, "data_full.h5ad")
    if not os.path.exists(data_path):
        raise FileNotFoundError(f"Data file not found: {data_path}")
    adata = ad.read_h5ad(data_path)

    basal_ref_path = os.path.join(data_dir, "basal_ref_KOV.txt")
    ko_ov_path = os.path.join(data_dir, "KO_OV_list.txt")

    if perturbation_training:
        _write_basal_ref_kov(data_dir, adata)
    else:
        _safe_remove(basal_ref_path)

    ko_targets = _parse_ko_output_genes(ko_output_genes, adata)
    if ko_targets:
        _write_ko_ov_list(data_dir, ko_targets)
    else:
        _safe_remove(ko_ov_path)


if __name__ == "__main__":
    main(sys.argv[1:])
