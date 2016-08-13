#' Filter l2fc table to 'best' sparsity cutoffs & density windows
#'
#' \code{filter_l2fc} filters a l2fc table to 'best' sparsity cutoffs & density windows.
#'
#' @param df_l2fc  data.frame of log2 fold change values
#' @param padj_cutoff  Adjusted p-value cutoff for rejecting the null hypothesis
#'   that l2fc values were not greater than the l2fc_threshold.
#' @param padj_method  Method for global p-value adjustment (See \code{p.adjust})
#'
filter_l2fc = function(df_l2fc, padj_cutoff=0.1, padj_method='BH'){
  padj_cutoff = as.numeric(padj_cutoff)
  if(is.null(df_l2fc$xxGUIDxx)){
    stop('l2fc table already filtered: no unique identifier found')
  }

  df_l2fc_f = df_l2fc %>%
    mutate(padj = p.adjust(p, method=padj_method)) %>%
    group_by(xxGUIDxx, density_min, density_max) %>%
    mutate(rej_hypo = padj < padj_cutoff) %>%
    group_by(xxGUIDxx, density_min, density_max, sparsity_threshold) %>%
    mutate(n_rej_hypo = sum(rej_hypo)) %>%
    group_by(xxGUIDxx, density_min, density_max) %>%
    filter(n_rej_hypo == max(n_rej_hypo)) %>%
    filter(sparsity_threshold == min(sparsity_threshold)) %>%
    ungroup() %>%
    dplyr::select(-rej_hypo, -n_rej_hypo)

  # filtering OTUs to just density window with the highest l2fc
  df_l2fc_f = df_l2fc_f %>%
    group_by(xxGUIDxx, OTU) %>%
    filter(log2FoldChange == max(log2FoldChange)) %>%
    ungroup() %>%
    dplyr::select(-xxGUIDxx)

  return(df_l2fc_f)
}



#' (MW)-HR-SIP analysis
#'
#' \code{HRSIP} conducts high resolution stable isotope probing (HR-SIP)
#'   analysis on phyloseq object.
#'
#' @inheritParams DESeq2_l2fc_multi
#' @param padj_cutoff  Adjusted p-value cutoff for rejecting the null hypothesis
#'   that l2fc values were not greater than the l2fc_threshold.
#' @param padj_method  Method for global p-value adjustment (See \code{p.adjust})
#'
#' @export
#'
#' @examples
#' data(physeq)
#' # This will compare all 13C treatments to the 12C controls (all time points combined)
#' density_windows = data.frame(start=c(1.70), end=c(1.75))
#' df_l2fc = HRSIP(physeq, density_windows=density_windows, design=~Substrate)
#' head(df_l2fc)
#'
#' # This will conduct individual comparisons between 12C controls & 13C treatments (by time point)
#' ## Get a listing of the parameters used for subsetting the phyloseq object for each comparison
#' params = get_treatment_params(physeq, c('Substrate', 'Day'), "Substrate != '12C-Con'")
#' ## Make an expression for subsetting the phyloseq object.
#' ## The values in braces will be replaced by values in the param.
#' ## The first "()" is selecting the control; the second "()" is selecting a treatment.
#' ex = "(Substrate=='12C-Con' & Day=='${Day}') | (Substrate=='${Substrate}' & Day == '${Day}')"
#' ## HRSIP run on each comparison
#' density_windows = data.frame(start=c(1.70), end=c(1.75))
#' df_l2fc = HRSIP(physeq, pairwise_expr=ex, pairwise_params=params, density_windows=density_windows, design=~Substrate)
#' head(df_l2fc)
#'
#' # Same as above, using multiple buoyant density windows (MW-HR-SIP)
#' density_windows = data.frame(start=c(1.70, 1.72), end=c(1.73, 1.75))
#' df_l2fc = HRSIP(physeq, pairwise_expr=ex, pairwise_params=params, density_windows=density_windows, design=~Substrate)
#' head(df_l2fc)
#'
#' # Same as above, but running in parallel
#' library(doParallel)
#' ncores=4
#' registerDoParallel(ncores)
#' df_l2fc = HRSIP(physeq, pairwise_expr=ex, pairwise_params=params, density_windows=density_windows, design=~Substrate, .parallel=TRUE)
#' head(df_l2fc)
#'
HRSIP = function(physeq,
                  density_windows,
                  design,
                  sparsity_threshold=seq(0, 0.3, 0.1),
                  sparsity_apply='all',
                  l2fc_threshold=0.25,
                  pairwise_expr=NULL,
                  pairwise_params=NULL,
                  padj_cutoff=0.1,
                  padj_method='BH',
                  parallel=FALSE,
                  ...){

  # assertions
  if(is.factor(density_windows)){
    density_windows = as.vector(density_windows)
  }
  if(is.vector(density_windows)){
    stopifnot(length(density_windows)>=2)
    density_windows = data.frame(start=c(density_windows[1]),
                                 end=c(density_windows[2]))
  }
  stopifnot(is.numeric(l2fc_threshold))
  stopifnot(is.character(sparsity_apply))
  stopifnot(all(sapply(sparsity_threshold, function(x) x>=0 & x<=1))==TRUE)
  stopifnot(is.data.frame(density_windows))
  stopifnot(ncol(density_windows) >= 2)

  # pairwise or not
  if(is.null(pairwise_expr)){
      df_l2fc = DESeq2_l2fc_multi(physeq=physeq,
                        sparsity_threshold=sparsity_threshold,
                        sparsity_apply=sparsity_apply,
                        density_windows=density_windows,
                        design=design,
                        l2fc_threshold=l2fc_threshold,
                        parallel=parallel,
                        ...)
  } else {
    # pairwise subsetting to just control & 1 treatment
    pairwise_expr = as.character(pairwise_expr)
    pairwise_params = as.data.frame(pairwise_params)
    pairwise_params_l = apply(pairwise_params, 1, as.list)

    df_l2fc = plyr::ldply(pairwise_params_l, DESeq2_l2fc_multi_sub,
                          ex=pairwise_expr,
                          physeq=physeq,
                          sparsity_threshold=sparsity_threshold,
                          sparsity_apply=sparsity_apply,
                          density_windows=density_windows,
                          design=design,
                          l2fc_threshold=l2fc_threshold,
                          parallel=parallel,
                          ...)
  }

  # filtering l2fc table
  if(!is.null(padj_cutoff)){
    df_l2fc = filter_l2fc(df_l2fc, padj_cutoff=padj_cutoff, padj_method=padj_method)
  } else {
    df_l2fc = dplyr::select(df_l2fc, -xxGUIDxx)
  }

  return(df_l2fc)
}
