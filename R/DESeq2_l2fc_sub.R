#' Get parameters for subsetting the phyloseq dataset
#'
#' This function is needed if you want to make multiple
#' subsets of the phyloseq object in order to make specific
#' comparisons between isotopically labeled-treatments and
#'
#'  their corresponding controls (eg., from the same time point).
#'
#' Makes a data.frame of all of the parameter values that differ
#' among the treatment-control comparisons.
#'
#' For example, if you want to compare the gradient fractions from
#' each labeled-treatment to its corresponding unlabeled-Control (both from the
#' same time point).
#'
#' @param physeq  Phyloseq object
#' @param exp_params  a vector listing the columns in the phyloseq sample_data
#'   table that can subset the phyloseq dataset in order to make the specific
#'   labeled-treatment vs labeled-control comparisons that you would like to make.
#' @param treatment  This is an expression used to filter out the
#'   control-specific parameters (if needed).
#'
#' @export
#'
#' @examples
#' data(physeq)
#' # Here, the treatment/controls (12C & 13C) are listed in substrate,
#' # and should be matched by 'Day'. The 13C-treatments can be identified by
#' # the expression: "Substrate != '12C-Con'"
#' get_treatment_params(physeq, c('Substrate', 'Day'), "Substrate != '12C-Con'")
#'
get_treatment_params = function(physeq, exp_params, treatment=NULL){
  physeq.m = phyloseq2df(physeq, phyloseq::sample_data)
  # filter out control (if needed; depending on subsetting expression)
  if(!is.null(treatment)){
    physeq.m = filter_(physeq.m, treatment)
  }
  # all pairwise params
  params = dplyr::distinct(physeq.m[,exp_params])
  return(params)
}


#' \code{DESeq2_l2fc_multi}: subset comparison
#'
#' Calls \code{DESeq2_l2fc_multi} with subsets of the dataset:
#' (eg., just gradient fractions from 1 treatment & its matching
#' control gradient samples)
#'
#' The phyloseq object subsetting is performed with an expression
#' that generates boolean values (eg., 1:10 < 3). Since each subset
#' requires a specific expression, the user provides a 'generalized'
#' expression in which specific values are inserted into the specific
#' expression for each subset.
#' For example: expression="(Substrate=='12C-Con') |
#' (Substrate=='${Substrate}')"
#' The "${Substrate}" value is changed for each subset, and the values
#' used for the parameter are set by the \code{params} argument.
#'
#'
#' @param params  data.frame of parameters supplies to \code{ex}
#' @param ex  Expression for subsetting the phyloseq object
#' @param physeq  Phyloseq object
#' @param density_min  Minimum buoyant density of the 'heavy' gradient fractions
#' @param density_max  Maximum buoyant density of the 'heavy' gradient fractions
#' @param design  \code{design} parameter used for DESeq2 analysis.
#'   See \code{DESeq2::DESeq} for more details.
#' @param l2fc_threshold  log2 fold change (l2fc) values must be significantly above this
#'   threshold in order to reject the hypothesis of equal counts.
#' @param sparsity_threshold  All OTUs observed in less than this portion (fraction: 0-1)
#'   of gradient fraction samples are pruned. A a form of indepedent filtering,
#'   The sparsity cutoff with the most rejected hypotheses is used.
#' @param sparsity_apply  Apply sparsity threshold to all gradient fraction samples ('all')
#'   or just heavy fraction samples ('heavy')
#' @return dataframe of HRSIP results
#'
#' @examples
#' data(physeq)
#' params = get_treatment_params(physeq, c('Substrate', 'Day'), "Substrate != '12C-Con'")
#' ex = "(Substrate=='12C-Con' & Day=='${Day}') | (Substrate=='${Substrate}' & Day == '${Day}')"
#' #
#' df_l2fc = plyr::ldply(params_l, DESeq2_l2fc_multi_sub, ex=ex, physeq=physeq)
#' head(df_l2fc)
#'
DESeq2_l2fc_multi_sub = function(params, ex, physeq, density_windows,
                                 design=~Substrate,
                                 l2fc_threshold=0.25,
                                 sparsity_threshold=seq(0,0.3,0.1),
                                 sparsity_apply='all',
                                 parallel=parallel){
  # Pruning phyloseq to pairwise comparison
  exx = stringterpolate(ex, params)
  cat('Subset:', exx, '\n')
  physeq.m = phyloseq2df(physeq, phyloseq::sample_data)
  bool = mutate_(physeq.m, exx)[,ncol(physeq.m)+1]
  ## assertions
  if(sum(bool) <= 1){
    stop('<2 samples selected based on selection expression')
  }
  stopifnot(length(bool)==nrow(physeq.m))
  ## pruning
  physeq.p = phyloseq::prune_samples(bool, physeq)
  ## HRSIP
  df_l2fc = DESeq2_l2fc_multi(physeq.p,
                              density_windows=density_windows,
                              design=design,
                              l2fc_threshold=l2fc_threshold,
                              sparsity_threshold=sparsity_threshold,
                              sparsity_apply=sparsity_apply,
                              parallel=parallel)

  # adding comparisons
  n = names(params)
  for(nn in n){
    df_l2fc[,nn] = params[nn]
  }
  return(df_l2fc)
}
