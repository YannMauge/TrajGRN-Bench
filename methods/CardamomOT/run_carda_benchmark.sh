#!/bin/bash

set -euo pipefail

source_adata="$1"
file="$2"
output_mode="$3"
perturbation_training="${4:-false}"
ko_output_genes="${5:-none}"
runner_mode="${BENCHMARK_RUNNER_MODE:-conda}"
cardamom_env="${BENCHMARK_CARDAMOMOT_CONDA_ENV:-cardamom_env}"

split="full"

run_python() {
  if [ "$runner_mode" = "conda" ]; then
    conda run -n "$cardamom_env" --live-stream python "$@"
  else
    python "$@"
  fi
}

if [ ! -d "$file/project_$(basename "${source_adata%.h5ad}")" ]; then
  mkdir -p "$file/project_$(basename "${source_adata%.h5ad}")"
  mkdir -p "$file/project_$(basename "${source_adata%.h5ad}")/Data"
  mkdir -p "$file/project_$(basename "${source_adata%.h5ad}")/CardamomOT"
  mkdir -p "$file/project_$(basename "${source_adata%.h5ad}")/Check"
  cp "$source_adata" "$file/project_$(basename "${source_adata%.h5ad}")/Data"
  mv "$file/project_$(basename "${source_adata%.h5ad}")/Data/$(basename "${source_adata%.h5ad}").h5ad" "$file/project_$(basename "${source_adata%.h5ad}")/Data/data_orig.h5ad"
fi

folder="$file/project_$(basename "${source_adata%.h5ad}")"

split="full"

run_python methods/CardamomOT/split_tps_extrapol.py -i "${folder}" -a "$1" -u "$output_mode" -p "$perturbation_training"
run_python methods/CardamomOT/prepare_benchmark_perturbations.py -i "${folder}" -p "$perturbation_training" -k "$ko_output_genes"

echo "Read depth inference"
run_python methods/CardamomOT/infer_rd.py -i "${folder}" 

echo "Inference mixture"
run_python methods/CardamomOT/infer_mixture.py -i "${folder}" -s "${split}"

echo "Infer network structure"
run_python methods/CardamomOT/infer_network_structure.py -i "${folder}" -s "${split}"

echo "Adapt network to simulate and degradation rates"
run_python methods/CardamomOT/infer_network_simul.py -i "${folder}" -s "${split}"

if [ "$output_mode" != "no_traj" ]; then

    if [ -f "$(dirname "$1")/subsample_train_ids.npy" ]; then
    split="test"
    fi

    echo "Simulate network"
    run_python methods/CardamomOT/simulate_network.py -i "${folder}" -s "${split}"
    
    echo "Check simulation"
    run_python methods/CardamomOT/check_sim_to_data.py -i "${folder}" -s "${split}"

    run_python methods/CardamomOT/concat_pred_train.py -i "${folder}" -a "$1" -k "$ko_output_genes" -p "$perturbation_training"
    
    echo "Saving adata_sim_final"
    cp "${folder}/CardamomOT/adata_sim_final.h5ad" "${file}/$(basename "${source_adata%.h5ad}")_adata.h5ad"

    if [[ "${ko_output_genes,,}" != "none" && -f "${folder}/Data/KO_OV_list.txt" ]]; then
        echo "Simulate KO/OV perturbations from Data/KO_OV_list.txt"
        run_python methods/CardamomOT/simulate_network_KOV.py -i "${folder}" -s "${split}"
        run_python methods/CardamomOT/check_KOV_to_sim.py -i "${folder}" -s "${split}"

        while IFS=$'\t' read -r ko_gene ov_genes; do
            if [[ "${ko_gene}" == "KO" ]]; then
                continue
            fi
            ko_gene="$(echo "${ko_gene}" | xargs)"
            ov_genes="$(echo "${ov_genes}" | xargs)"
            if [[ -z "${ko_gene}" || "${ko_gene}" == *","* || -n "${ov_genes}" ]]; then
                continue
            fi

            shopt -s nullglob
            ko_adata_matches=("${folder}/CardamomOT/adata_sim_KO_${ko_gene}_OV_none_stim"*"_prior"*.h5ad)
            shopt -u nullglob
            if [ ${#ko_adata_matches[@]} -gt 0 ]; then
                cp "${ko_adata_matches[0]}" "${file}/$(basename "${source_adata%.h5ad}")_ko_${ko_gene}_adata.h5ad"
            fi
        done < "${folder}/Data/KO_OV_list.txt"
    fi
else
    echo "Skipping simulation as output_mode is set to no_traj"
fi

echo "Saving inter_benchmark.npy as GRN"
if [ -f "${folder}/CardamomOT/inter_benchmark.npy" ]; then
    cp "${folder}/CardamomOT/inter_benchmark.npy" "${file}/$(basename "${source_adata%.h5ad}")_GRN.npy"
elif [ -f "${folder}/CardamomOT/inter.npy" ]; then
    cp "${folder}/CardamomOT/inter.npy" "${file}/$(basename "${source_adata%.h5ad}")_GRN.npy"
else
    echo "Warning: no GRN file found in ${folder}/CardamomOT (expected inter_benchmark.npy or inter.npy)"
fi

echo "All scripts executed !"
