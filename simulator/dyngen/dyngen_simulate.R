#!/usr/bin/env Rscript
# =============================================================================
# dyngen_simulate.R — Chain-based dyngen simulation (gold standard pipeline).
#
# Each benchmark gene = chain of CHAIN_LEN dyngen modules (intra-chain
# strength=1 for propagation delay). Cross-chain edges from benchmark GRN
# (regulator_last → target_first). num_tfs=50 per module for dimred stability.
# Expression aggregated per gene, binned by cumulative global time.
# =============================================================================

suppressPackageStartupMessages({
  library(dyngen); library(dynutils); library(tibble); library(Matrix); library(dplyr)
})

CHAIN_LEN <- 5

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  get_arg <- function(name, default = NULL) {
    idx <- which(args == name)
    if (length(idx) == 0) return(default)
    if (idx == length(args)) stop(sprintf("Missing value for %s", name))
    args[idx + 1]
  }
  list(output_folder=get_arg("--output_folder"), n_runs=as.integer(get_arg("--n_runs","2")),
       n_genes=as.integer(get_arg("--n_genes","8")), n_cells=as.integer(get_arg("--n_cells","1000")),
       n_time_bins=as.integer(get_arg("--n_time_bins","20")), ko_genes=get_arg("--ko_genes","none"),
       seed=as.integer(get_arg("--seed","42")), benchmark_grn=get_arg("--benchmark_grn",""))
}

build_backbone <- function(grn_str, n_genes) {
  parts <- strsplit(grn_str, ";")[[1]]
  stimulus_targets <- integer(0)
  reg_from <- reg_to <- integer(0); reg_eff <- integer(0); reg_str <- numeric(0)
  for (p in parts) {
    if (nchar(p)==0) next
    m <- regmatches(p, regexec("^([0-9]+)>([0-9]+):(-?[0-9]+):([0-9]+)$", p))[[1]]
    if (length(m)!=5) next
    r <- as.integer(m[2]); t <- as.integer(m[3])
    if (r==0) { stimulus_targets <- c(stimulus_targets, t) }
    else { reg_from<-c(reg_from,r); reg_to<-c(reg_to,t); reg_eff<-c(reg_eff,as.integer(m[4])); reg_str<-c(reg_str,as.numeric(m[5])) }
  }

  module_ids <- character(0); basal_vec <- numeric(0); burn_vec <- logical(0)
  gene_first <- character(n_genes); gene_last <- character(n_genes)
  for (g in seq_len(n_genes)) {
    chain_ids <- paste0("G", g-1, "_", seq_len(CHAIN_LEN))
    module_ids <- c(module_ids, chain_ids)
    is_stim <- g %in% stimulus_targets
    basal_vec <- c(basal_vec, if(is_stim) c(1, rep(0, CHAIN_LEN-1)) else rep(0, CHAIN_LEN))
    burn_vec  <- c(burn_vec,  rep(is_stim, CHAIN_LEN))
    gene_first[g] <- chain_ids[1]; gene_last[g] <- chain_ids[CHAIN_LEN]
  }
  module_info <- tibble(module_id=module_ids, basal=basal_vec, burn=burn_vec, independence=rep(1,length(module_ids)))

  intra_from <- character(0); intra_to <- character(0)
  for (g in seq_len(n_genes))
    for (j in seq_len(CHAIN_LEN-1)) {
      intra_from <- c(intra_from, paste0("G",g-1,"_",j))
      intra_to   <- c(intra_to,   paste0("G",g-1,"_",j+1))
    }
  intra_net <- tibble(from=intra_from, to=intra_to, effect=rep(1L,length(intra_from)),
                       strength=rep(1,length(intra_from)), hill=rep(2,length(intra_from)))
  cross_net <- if (length(reg_from) > 0)
    tibble(from=gene_last[reg_from], to=gene_first[reg_to], effect=reg_eff, strength=reg_str, hill=rep(2,length(reg_from)))
  else tibble(from=character(0),to=character(0),effect=integer(0),strength=numeric(0),hill=numeric(0))
  module_network <- bind_rows(intra_net, cross_net)

  non_stim <- setdiff(seq_len(n_genes), stimulus_targets)
  if (length(non_stim) > 0) {
    epf<-"sBurn"; ept<-"sG1"; epp<-""; eps<-TRUE; epb<-TRUE; eptime<-CHAIN_LEN*20
    for (i in seq_along(non_stim)) {
      g <- non_stim[i]; chain_mods <- paste0("G",g-1,"_",seq_len(CHAIN_LEN))
      epf<-c(epf,paste0("sG",i)); ept<-c(ept,paste0("sG",i+1))
      epp<-c(epp,paste0("+",paste(chain_mods,collapse=",+")))
      eps<-c(eps,FALSE); epb<-c(epb,FALSE); eptime<-c(eptime, CHAIN_LEN*20)
    }
    epf<-c(epf,paste0("sG",length(non_stim)+1)); ept<-c(ept,"sEnd")
    epp<-c(epp,""); eps<-c(eps,FALSE); epb<-c(epb,FALSE); eptime<-c(eptime,10)
  } else {
    epf<-c("sBurn","sBurn"); ept<-c("sG1","sEnd"); epp<-c("","")
    eps<-c(TRUE,FALSE); epb<-c(TRUE,FALSE); eptime<-c(CHAIN_LEN*20, CHAIN_LEN*20*8)
  }
  ep <- tibble(from=epf,to=ept,module_progression=epp,start=eps,burn=epb,time=eptime)
  backbone(module_info=module_info, module_network=module_network, expression_patterns=ep)
}

# ── Run full SSA pipeline ────────────────────────────────────────────────────
run_model <- function(n_genes, n_cells, seed, benchmark_grn) {
  bb <- build_backbone(benchmark_grn, n_genes)
  backbone_time <- simtime_from_backbone(bb)
  set.seed(seed)
  cfg <- initialise_model(
    backbone = bb, num_cells = n_cells,
    num_tfs = 10, num_targets = 5, num_hks = 0,
    gold_standard_params = gold_standard_default(tau=100/3600, census_interval=1),
    verbose = FALSE, download_cache_dir = tools::R_user_dir("dyngen", "data"),
    simulation_params = simulation_default(
      total_time = backbone_time, census_interval = 10,
      ssa_algorithm = ssa_etl(tau = 100 / 3600),
      experiment_params = simulation_type_wild_type(num_simulations = 100),
      compute_dimred = FALSE   # skip dimred to avoid null-space error
    )
  )
  model <- generate_tf_network(cfg)
  model <- generate_feature_network(model)
  model <- generate_kinetics(model)
  model <- generate_gold_standard(model)   # dimred: warning only
  model <- generate_cells(model)            # dimred skipped
  model <- generate_experiment(model)
  model
}

# ── Extract expression from SSA experiment, aggregate chains, bin by step_ix ─
extract_expression <- function(model, n_genes, n_cells, n_time_bins) {
  expr_mat <- model$experiment$counts_mrna
  if (is.null(expr_mat)) stop("No expression matrix")
  if (inherits(expr_mat, "sparseMatrix")) expr_mat <- as.matrix(expr_mat)

  fi <- model$feature_info; tfs <- fi[fi$is_tf,]
  id_col <- if ("feature_id" %in% colnames(tfs)) "feature_id" else "gene_id"

  gene_expr <- matrix(0, nrow=n_genes, ncol=nrow(expr_mat))
  gene_names <- paste0("Gene", seq(0, n_genes-1))
  for (g in seq_len(n_genes)) {
    tf_ids <- tfs[[id_col]][grepl(paste0("^G", g-1, "_"), tfs$module_id)]
    cols <- intersect(tf_ids, colnames(expr_mat))
    if (length(cols) > 0) gene_expr[g, ] <- rowMeans(expr_mat[, cols, drop=FALSE])
  }

  # Bin by trajectory progression: compute a global progression score
  # that increases from 0 to 1 as cells advance through the trajectory.
  # We use the gold standard edge order: each unique (from,to) edge gets
  # a sequential index, and within each edge we use the normalized time.
  cell_info <- model$experiment$cell_info
  edges <- unique(cell_info[, c("from","to")])
  edge_order <- setNames(seq_len(nrow(edges)) - 1, paste(edges$from, edges$to, sep="->"))
  edge_key <- paste(cell_info$from, cell_info$to, sep="->")
  edge_idx <- edge_order[edge_key]
  # Global progression: edge_index contributes its offset, time contributes within-edge position
  global_prog <- (edge_idx + cell_info$time) / max(edge_idx + 1, 1)
  prog_order <- order(global_prog)
  n_avail <- length(prog_order)
  cells_per_bin <- max(1, floor(n_avail / n_time_bins))
  bin_centers <- seq(2.5, 97.5, length.out=n_time_bins)

  total_cells <- min(n_cells, cells_per_bin * n_time_bins)
  counts <- matrix(0, nrow=n_genes, ncol=total_cells)
  time_vec <- numeric(total_cells)
  cell_idx <- 1
  for (b in seq_len(n_time_bins)) {
    si <- (b-1)*cells_per_bin + 1; ei <- min(b*cells_per_bin, n_avail)
    for (ci in prog_order[si:ei]) {
      if (cell_idx > total_cells) break
      counts[, cell_idx] <- pmax(0, round(gene_expr[, ci]))
      time_vec[cell_idx] <- bin_centers[b]
      cell_idx <- cell_idx + 1
    }
    if (cell_idx > total_cells) break
  }
  list(counts=counts[,1:(cell_idx-1),drop=FALSE], time=time_vec[1:(cell_idx-1)],
       gene_names=gene_names, model=model)
}

run_ko <- function(model_base, ko_gene_idx, n_genes, n_cells, n_time_bins, seed) {
  fi <- model_base$feature_info; tfs <- fi[fi$is_tf,]
  id_col <- if("feature_id" %in% colnames(tfs)) "feature_id" else "gene_id"
  ko_pattern <- paste0("^G", ko_gene_idx-1, "_")
  ko_tf_ids <- tfs[[id_col]][grepl(ko_pattern, tfs$module_id)]
  set.seed(seed)
  model_ko <- model_base
  model_ko$simulation_params$experiment_params <- simulation_type_knockdown(
    num_simulations = 100, timepoint = 0, genes = ko_tf_ids,
    num_genes = length(ko_tf_ids), multiplier = 0
  )
  model_ko <- generate_cells(model_ko)
  model_ko <- generate_experiment(model_ko)
  extract_expression(model_ko, n_genes, n_cells, n_time_bins)
}

write_outputs <- function(counts, time_vec, gene_names, output_dir, run_idx, ko_label=NULL) {
  ng <- length(gene_names); nc <- ncol(counts)
  data <- matrix(0, nrow=nc+1, ncol=ng+2)
  data[1,1]<-0; data[1,2]<-0; data[1,3:(ng+2)]<-seq_len(ng)
  data[2:(nc+1),1]<-time_vec; data[2:(nc+1),2]<-ifelse(time_vec>0,100,0)
  data[2:(nc+1),3:(ng+2)]<-pmax(0,round(counts)); data<-t(data)
  suffix<-if(is.null(ko_label))""else paste0("_ko_",ko_label)
  write.table(data,file=file.path(output_dir,sprintf("data_%d%s.txt",run_idx,suffix)),
              sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
}
write_gene_panel <- function(output_dir, gene_names) {
  panel<-cbind(c(0,seq_along(gene_names)),c("Stimulus",gene_names))
  write.table(panel,file=file.path(output_dir,"panel_genes.txt"),sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
}
write_degradation_rates <- function(output_dir, n_genes) {
  dir.create(file.path(output_dir,"Rates"),showWarnings=FALSE,recursive=TRUE)
  write.table(cbind(rep(0.25,n_genes+1),rep(0.05,n_genes+1)),
              file=file.path(output_dir,"Rates","degradation_rates.txt"),sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
}

main <- function() {
  params <- parse_args()
  if(is.null(params$output_folder)) stop("--output_folder is required")
  dir.create(params$output_folder, showWarnings=FALSE, recursive=TRUE)
  ko_spec <- trimws(params$ko_genes); simulate_ko <- !(ko_spec %in% c("","none"))
  for (run_idx in seq_len(params$n_runs)) {
    run_seed <- params$seed + run_idx - 1
    message(sprintf("\n=== Run %d/%d (seed=%d) ===", run_idx, params$n_runs, run_seed))
    model <- run_model(params$n_genes, params$n_cells, run_seed, params$benchmark_grn)
    result <- extract_expression(model, params$n_genes, params$n_cells, params$n_time_bins)
    if (run_idx==1) {
      write_gene_panel(params$output_folder, result$gene_names)
      write_degradation_rates(params$output_folder, params$n_genes)
    }
    ko_suffix <- if(simulate_ko) "WT" else NULL
    write_outputs(result$counts, result$time, result$gene_names, params$output_folder, run_idx, ko_label=ko_suffix)
    if (simulate_ko) {
      ko_gene_list <- if(ko_spec=="all") result$gene_names else intersect(trimws(strsplit(ko_spec,",")[[1]]),result$gene_names)
      for (ko_gene in ko_gene_list) {
        message(sprintf("  KO: %s",ko_gene))
        gene_idx <- match(ko_gene, result$gene_names)
        if(is.na(gene_idx)){warning(sprintf("Cannot map '%s'.",ko_gene));next}
        ko_expr <- run_ko(model, gene_idx, params$n_genes, params$n_cells, params$n_time_bins)
        if(!is.null(ko_expr)) write_outputs(ko_expr$counts, ko_expr$time, result$gene_names, params$output_folder, run_idx, ko_label=ko_gene)
      }
    }
  }
  message("\ndyngen simulation completed successfully.")
}
main()
