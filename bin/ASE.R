#!/usr/bin/env Rscript

# Fit the copy-number-aware BayesASE model to one shard of genes from a single contrast.
# The model is compiled once by COMPILE_ASE_MODEL; a .stan path also works when
# the script is run on its own.

suppressPackageStartupMessages({
  library(data.table)
  library(rstan)
})

argv <- commandArgs(trailingOnly = TRUE)
opt <- function(name, default = NULL, type = c("character", "integer")) {
  type <- match.arg(type)
  index <- match(name, argv)
  value <- if (is.na(index) || index == length(argv)) default else argv[index + 1]
  if (type == "integer") return(as.integer(value))
  value
}

input_file <- opt("--input")
model_file <- opt("--model")
output_file <- opt("--output")
nsim <- opt("--iterations", 100000, "integer")
nburnin <- opt("--warmup", 10000, "integer")
seed <- opt("--seed", 1, "integer")

if (is.null(input_file) || is.null(model_file) || is.null(output_file)) {
  stop("--input, --model and --output are required")
}
if (nburnin >= nsim) stop("warmup must be smaller than iterations")

rstan_options(auto_write = FALSE)
ASE_df <- fread(input_file)
required <- c("sample", "group", "ID", "H1_counts", "H2_counts", "NonHS_counts",
  "CN_H1", "CN_H2", "H1_prior", "H2_prior", "test", "groupA", "groupB")
missing <- setdiff(required, names(ASE_df))
if (length(missing)) stop(sprintf("ASE shard missing columns: %s", paste(missing, collapse = ", ")))
if (uniqueN(ASE_df$test) != 1 || uniqueN(ASE_df$groupA) != 1 || uniqueN(ASE_df$groupB) != 1) {
  stop("each shard must contain exactly one test")
}

group_levels <- c(unique(ASE_df$groupA), unique(ASE_df$groupB))
ASE_df <- ASE_df[group %in% group_levels]
ASE_df[, group := factor(group, levels = group_levels)]
setorder(ASE_df, ID, group)

gam.mles.data = function(x){
  n = length(x)
  xb = mean(x)
  xd = mean(log(x))
  s = log(xb)-xd
  a0 = ifelse(s>0,(3.0-s+sqrt((s-3.0)^2+24.0*s))/12.0/s,50)
  l = 1
  repeat{
    ans = (log(a0)-digamma(a0)-s)
    a1 = a0-ans/(1.0/a0-trigamma(a0))
    if(abs(ans) <= 1.0e-7 | l >= 30){break}
    a0 = a1
    l = l+1
  }
  ah = a1; bh = xb/a1
  return(c(ah,bh))
}

prior_empBayes_forbeta <- function(xs, ys, zs, cnv1, cnv2){
  bbeta_est = (xs + ys + zs) / (cnv1 + cnv2)
  bbeta_est[which(bbeta_est==0)]=0.1
  tem=gam.mles.data(bbeta_est)
  tem[1]=min(max(tem[1],10^(-3)),10^5)
  tem[2]=min(max(tem[2],10^(-3)),10^5)
  a_beta=tem[1]
  a_b_beta=2*tem[2]^(-1)
  b_b_beta=2
  return(list(a_beta=a_beta,a_b_beta=a_b_beta,b_b_beta=b_b_beta))
}

process_gene <- function(gene_data, group_levels, compiled_model, nsim, nburnin, seed) {
  tryCatch({
    test_name <- paste(group_levels, collapse = "_VS_")
    K <- nrow(gene_data)
    n_groups <- length(group_levels)
    xenv <- as.numeric(gene_data$group)
    smallestGroupSize <- min(table(xenv))
    xs <- as.integer(gene_data$H1_counts)
    ys <- as.integer(gene_data$H2_counts)
    zs <- as.integer(gene_data$NonHS_counts)
    cnv1 <- as.numeric(gene_data$CN_H1)
    cnv2 <- as.numeric(gene_data$CN_H2)
    r <- gene_data[, .(H1_prior=mean(H1_prior), H2_prior=mean(H2_prior)), by=group][, .(H1_prior, H2_prior)] |> as.matrix()
    for (grp_idx in 1:n_groups) {
      grp_mask <- which(xenv == grp_idx)
      current_xs <- xs[grp_mask]
      current_ys <- ys[grp_mask]
      current_zs <- zs[grp_mask]
      if (sum(current_zs) == 0) {zs[grp_mask[which.max(current_xs + current_ys)]] <- 1}
      if (sum(current_xs) == 0) {xs[grp_mask[which.max(zs[grp_mask])]] <- 1}
      if (sum(current_ys) == 0) {ys[grp_mask[which.max(zs[grp_mask])]] <- 1}
    }
    hyper_beta <- prior_empBayes_forbeta(xs, ys, zs, cnv1, cnv2)
    datastan <- list(K = K, n_environment = n_groups, xenv = xenv, xs = xs, ys = ys, zs = zs, r = r, cnv1 = cnv1, cnv2 = cnv2,
      a_beta = hyper_beta$a_beta, a_b_beta = hyper_beta$a_b_beta, b_b_beta = hyper_beta$b_b_beta,
      a_overdispersion = 2.01, b_overdispersion = 0.05)
    starting_values <- function() {list(overdispersion = 0.01, bbeta = (datastan$xs + datastan$ys + datastan$zs) / (cnv1 + cnv2), alpha = rep(1.0, datastan$n_environment))}
    total_xs_vec <- sapply(1:n_groups, function(i) sum(datastan$xs[datastan$xenv == i]))
    total_ys_vec <- sapply(1:n_groups, function(i) sum(datastan$ys[datastan$xenv == i]))
    total_zs_vec <- sapply(1:n_groups, function(i) sum(datastan$zs[datastan$xenv == i]))
    mean_xs_vec <- sapply(1:n_groups, function(i) mean(datastan$xs[datastan$xenv == i]))
    mean_ys_vec <- sapply(1:n_groups, function(i) mean(datastan$ys[datastan$xenv == i]))
    mean_zs_vec <- sapply(1:n_groups, function(i) mean(datastan$zs[datastan$xenv == i]))
    total_counts_sample <- datastan$xs + datastan$ys + datastan$zs
    cond1 <- sum(total_counts_sample >= 10, na.rm = TRUE) >= smallestGroupSize
    cond2 <- sum((r - r^2) != 0) > 0
    cond3 <- smallestGroupSize >= 3
    if (cond1 && cond2 && cond3) {
      max_try <- 15
      n_try <- 0
      repeat{
        n_try <- n_try + 1
        fit1 <- sampling(object=compiled_model, data=datastan, chains=1, warmup=nburnin, iter=nsim, refresh=0, init=starting_values, pars=c("alpha", "theta"), cores=1, open_progress=FALSE, seed=seed + n_try)
        theta <- rstan::extract(fit1, pars = "theta")$theta
        alpha <- rstan::extract(fit1, pars = "alpha")$alpha
        pvalue_test <- 2 * min(mean(alpha[, 1] > alpha[, 2]), mean(alpha[, 1] < alpha[, 2]))
        analysis_flag <- "Success"
        if(pvalue_test > 0){break}
        if(n_try >= max_try){analysis_flag <- "pvalue0"; break}
      }
    }else{theta <- matrix(NA_real_,n_groups,n_groups); alpha <- matrix(NA_real_,n_groups,n_groups); analysis_flag <- "Skipped_Low_Counts"}
    delta_threshold <- 0.15
    diff_distribution <- theta[, 1] - theta[, 2]
    gene_res <- list(test = test_name,
      groupA = group_levels[1],
      groupB = group_levels[2],
      groupA_priorH1 = datastan$r[1,1],
      groupA_priorH2 = datastan$r[1,2],
      groupA_totalH1 = as.numeric(total_xs_vec[1]),
      groupA_totalH2 = as.numeric(total_ys_vec[1]),
      groupA_totalNonHS = as.numeric(total_zs_vec[1]),
      groupA_meanH1 = as.numeric(mean_xs_vec[1]),
      groupA_meanH2 = as.numeric(mean_ys_vec[1]),
      groupA_meanNonHS = as.numeric(mean_zs_vec[1]),
      groupB_priorH1 = datastan$r[2,1],
      groupB_priorH2 = datastan$r[2,2],
      groupB_totalH1 = as.numeric(total_xs_vec[2]),
      groupB_totalH2 = as.numeric(total_ys_vec[2]),
      groupB_totalNonHS = as.numeric(total_zs_vec[2]),
      groupB_meanH1 = as.numeric(mean_xs_vec[2]),
      groupB_meanH2 = as.numeric(mean_ys_vec[2]),
      groupB_meanNonHS = as.numeric(mean_zs_vec[2]),
      groupA_alpha_mean = mean(alpha[, 1]),
      groupA_theta_mean = mean(theta[, 1]),
      groupA_theta_q025 = as.numeric(quantile(theta[, 1], 0.025, na.rm=TRUE)),
      groupA_theta_q975 = as.numeric(quantile(theta[, 1], 0.975, na.rm=TRUE)),
      groupB_alpha_mean = mean(alpha[, 2]),
      groupB_theta_mean = mean(theta[, 2]),
      groupB_theta_q025 = as.numeric(quantile(theta[, 2], 0.025, na.rm=TRUE)),
      groupB_theta_q975 = as.numeric(quantile(theta[, 2], 0.975, na.rm=TRUE)),
      groupA_alphaAI_pvalue = 2 * min(mean(alpha[, 1] > 1), mean(alpha[, 1] < 1)),
      groupA_thetaAI_pvalue = 2 * min(mean(theta[, 1] > 0.5), mean(theta[, 1] < 0.5)),
      groupB_alphaAI_pvalue = 2 * min(mean(alpha[, 2] > 1), mean(alpha[, 2] < 1)),
      groupB_thetaAI_pvalue = 2 * min(mean(theta[, 2] > 0.5), mean(theta[, 2] < 0.5)),
      delta_theta_mean = mean(theta[, 1]) - mean(theta[, 2]),
      delta_alpha_mean = mean(alpha[, 1]) - mean(alpha[, 2]),
      diffAI_pvalue = 2 * min(mean(alpha[,1] > alpha[,2]), mean(alpha[,1] < alpha[,2])),
      rope_value = mean(abs(diff_distribution) < delta_threshold),
      analysis_flag = analysis_flag
    )
    return(gene_res)
  }, error = function(e) {
    gene_res <- data.table(test=NA_character_, groupA=NA_character_, groupB=NA_character_,
      groupA_priorH1=NA_real_, groupA_priorH2=NA_real_, groupA_totalH1=NA_real_, groupA_totalH2=NA_real_, groupA_totalNonHS=NA_real_, groupA_meanH1=NA_real_, groupA_meanH2=NA_real_, groupA_meanNonHS=NA_real_,
      groupB_priorH1=NA_real_, groupB_priorH2=NA_real_, groupB_totalH1=NA_real_, groupB_totalH2=NA_real_, groupB_totalNonHS=NA_real_, groupB_meanH1=NA_real_, groupB_meanH2=NA_real_, groupB_meanNonHS=NA_real_,
      groupA_alpha_mean=NA_real_, groupA_theta_mean=NA_real_, groupA_theta_q025=NA_real_, groupA_theta_q975=NA_real_,
      groupB_alpha_mean=NA_real_, groupB_theta_mean=NA_real_, groupB_theta_q025=NA_real_, groupB_theta_q975=NA_real_,
      groupA_alphaAI_pvalue=NA_real_, groupA_thetaAI_pvalue=NA_real_, groupB_alphaAI_pvalue=NA_real_, groupB_thetaAI_pvalue=NA_real_,
      delta_theta_mean=NA_real_, delta_alpha_mean=NA_real_, diffAI_pvalue=NA_real_, rope_value=NA_real_, analysis_flag=paste("Error:", e$message))
    return(gene_res)
  })
}

compiled_model <- if (grepl("\\.rds$", model_file)) readRDS(model_file) else stan_model(file = model_file)
results_df <- ASE_df[, process_gene(.SD, group_levels, compiled_model, nsim, nburnin, seed), by = ID]
fwrite(results_df, output_file, sep="\t", quote=FALSE, na="NA")
