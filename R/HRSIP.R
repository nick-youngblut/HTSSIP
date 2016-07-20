#' phyloseq data object conversion to data.frame
#'
#' Conducts conversion of 1 of the data objects
#' in a phyloseq object (eg., tax_table) to a dataframe
#'
#' @param physeq  Phyloseq object
#' @param table_func  See \code{Phyloseq::phyloseq-class} for options
#' @return dataframe
#'
#' @export
#'
#' @examples
#' data(physeq)
#' df_otu = phyloseq2df(physeq, table_func=otu_table)
#' head(df_otu)
phyloseq2df = function(physeq, table_func){
  physeq.md = table_func(physeq)
  physeq.md = suppressWarnings(as.data.frame(as.matrix(physeq.md)))
  physeq.md = as.matrix(data.frame(lapply(physeq.md, as.character)))
  physeq.md = as.data.frame(apply(physeq.md, 2, trimws))
  return(physeq.md)
}

#' Format phyloseq & call DESeq2
#'
#' Conducts the actual trimming of the phyloseq dataset & calling of DESeq2
#'
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
#' df_l2fc = DESeq2_l2fc(physeq, density_min=1.71, density_max=1.75, design=~Substrate)
#' head(df_l2fc)
DESeq2_l2fc = function(physeq,
                  density_min,
                  density_max,
                  design,
                  l2fc_threshold=0.25,
                  sparsity_threshold=0.25,
                  sparsity_apply='all', ...){
  # assertions
  l2fc_threshold = as.numeric(l2fc_threshold)
  stopifnot(l2fc_threshold >= 0 & l2fc_threshold <= 1)
  sparsity_apply = tolower(sparsity_apply)
  stopifnot(sparsity_apply %in% c('all', 'heavy'))
  physeq.md = phyloseq::sample_data(physeq)
  stopifnot(!is.null(physeq.md$Buoyant_density))

  # status
  cat('Sparsity threshold:', sparsity_threshold, '\n')

  # sparcity cutoff applied to all gradient fractions
  prn = function(x) sum(x > 0) > sparsity_threshold * length(x)
  if(sparsity_apply=='all'){
    physeq = phyloseq::filter_taxa(physeq, prn, TRUE)
  }

  # window selection
  ## applying 'heavy' window pruning
  physeq = phyloseq::prune_samples((physeq.md$Buoyant_density >= density_min) &
                                   (physeq.md$Buoyant_density <= density_max),  physeq)

  # removing 0-abundance taxa
  physeq = phyloseq::filter_taxa(physeq, function(x) sum(x > 0) > 0 * length(x), TRUE)

  # sparsity cutoff applied to just heavy fractions
  if(sparsity_apply=='heavy'){
    physeq = phyloseq::filter_taxa(physeq, prn, TRUE)
  }

  # deseq
  dds = phyloseq::phyloseq_to_deseq2(physeq, design)   # design=~Substrate
  dds = DESeq2::DESeq(dds, quiet = TRUE, fitType = "local")
  theta = l2fc_threshold

  # results
  res = DESeq2::results(dds, independentFiltering=FALSE)
  res$OTU = rownames(res)

  # p-value
  beta = res$log2FoldChange
  betaSE = res$lfcSE
  p = pnorm(beta, theta, betaSE, lower.tail=FALSE)
  res$p = p
  d = data.frame(res[, c('OTU','log2FoldChange', 'p')])

  # p-value adjust
  d$padj = p.adjust(p, method = 'BH')

  # taxonomy data
  TT = phyloseq::tax_table(physeq)
  if(!is.null(TT)){
    TT = as.data.frame(as.matrix(TT))
    TT$OTU = rownames(TT)
    d = dplyr::left_join(d, TT, c('OTU'))
  }

  # setting pruning info
  d$density_min = density_min
  d$density_max = density_max
  d$sparsity_threshold = sparsity_threshold
  d$sparsity_apply = sparsity_apply
  return(d)
}


#' \code{DESeq2_l2fc}: multiple sparsity thresholds
#'
#' Calls \code{DESeq2_l2fc} with different sparsity thresholds
#' and uses the l2fc results of the sparsity threshold that produced
#' the most rejected hypotheses
#'
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
#' @param padj_cutoff  Adjusted P-value cutoff for identifying rejected hypotheses,
#'   which is in turn used to select the sparsity cutoff with the most rejected hypotheses.
#' @param .parallel  Run in parallel. See plyr::adply for more information.
#' @return dataframe of HRSIP results
#'
#' @examples
#' data(physeq)
#' df_l2fc = DESeq2_l2fc_MS(physeq, density_min=1.71, density_max=1.75, design=~Substrate)
#' head(df_l2fc)
DESeq2_l2fc_MS = function(physeq,
                          density_min,
                          density_max,
                          design,
                          l2fc_threshold=0.25,
                          sparsity_threshold=seq(0, 0.5, 0.1),
                          sparsity_apply='all',
                          padj_cutoff=0.1,
                          .parallel=FALSE,
                          ...){
  df_l2fc = plyr::adply(sparsity_threshold, 1, function(x){
    DESeq2_l2fc(physeq=physeq,
           density_min=density_min,
           density_max=density_max,
           design=design,
           l2fc_threshold=l2fc_threshold,
           sparsity_threshold=x,
           sparsity_apply=sparsity_apply,
           ...)
    }, .parallel=.parallel)

  # summing number of rejected hypotheses
  ## selecting lowest sparsity with the most rej_hypo
  df_l2fc = df_l2fc %>%
      mutate(rej_hypo = padj < padj_cutoff) %>%
      group_by(sparsity_threshold) %>%
      mutate(n_rej_hypo = sum(rej_hypo)) %>%
      ungroup() %>%
      filter(n_rej_hypo == max(n_rej_hypo)) %>%
      filter(sparsity_threshold == min(sparsity_threshold)) %>%
      dplyr::select(-n_rej_hypo)
  return(df_l2fc)
}


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


#' \code{DESeq2_l2fc_MS}: subset comparison
#'
#' Calls \code{DESeq2_l2fc_MS} with subsets of the dataset:
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
#' df_l2fc = plyr::ldply(params_l, DESeq2_l2fc_MS_sub, ex=ex, physeq=physeq)
#' head(df_l2fc)
DESeq2_l2fc_MS_sub = function(params, ex, physeq,
                           density_min=1.71,
                           density_max=1.75,
                           design=~Substrate,
                           l2fc_threshold=0.25,
                           padj_cutoff=0.1,
                           sparsity_threshold=seq(0,0.3,0.1),
                           sparsity_apply='all'){
  # Pruning phyloseq to pairwise comparison
  exx = stringterpolate(ex, params)
  cat(exx, '\n')
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
  df_l2fc = DESeq2_l2fc_MS(physeq.p,
                    density_min=density_min,
                    density_max=density_max,
                    design=design,
                    l2fc_threshold=l2fc_threshold,
                    padj_cutoff=padj_cutoff,
                    sparsity_threshold=sparsity_threshold,
                    sparsity_apply=sparsity_apply)

  # adding comparisons
  n = names(params)
  for(nn in n){
    df_l2fc[,nn] = params[nn]
  }
  return(df_l2fc)
}



#' HR-SIP analysis
#'
#' \code{HRSIP} conducts high resolution stable isotope probing (HR-SIP)
#'   analysis on phyloseq object.
#'
#' @inheritParams DESeq2_l2fc_MS
#'
#' @export
#'
#' @examples
#' data(physeq)
#' # This will compare all 13C treatments to the 12C controls (all time points combined)
#' df_l2fc = HRSIP(physeq, density_min=1.71, density_max=1.75, design=~Substrate)
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
#' df_l2fc = HRSIP(physeq, pairwise_expr=ex, pairwise_params=params, density_min=1.71, density_max=1.75, design=~Substrate)
#' head(df_l2fc)
#'
#' # Same as above, but running in parallel
#' library(doParallel)
#' ncores=4
#' registerDoParallel(ncores)
#' df_l2fc = HRSIP(physeq, pairwise_expr=ex, pairwise_params=params, density_min=1.71, density_max=1.75, design=~Substrate, .parallel=TRUE)
#' head(df_l2fc)
HRSIP = function(physeq,
                  density_min,
                  density_max,
                  design,
                  sparsity_threshold=seq(0, 0.3, 0.1),
                  sparsity_apply='all',
                  l2fc_threshold=0.25,
                  padj_cutoff=0.1,
                  pairwise_expr=NULL,
                  pairwise_params=NULL,
                  .parallel=FALSE,
                  ...){

  # assertions
  stopifnot(is.numeric(l2fc_threshold))
  stopifnot(is.numeric(density_min))
  stopifnot(is.numeric(density_max))
  stopifnot(is.character(sparsity_apply))
  stopifnot(all(sapply(sparsity_threshold, function(x) x>=0 & x<=1))==TRUE)

  # pairwise or not
  if(is.null(pairwise_expr)){
      df_l2fc = DESeq2_l2fc_MS(physeq=physeq,
                        sparsity_threshold=sparsity_threshold,
                        sparsity_apply=sparsity_apply,
                        density_min=density_min,
                        density_max=density_max,
                        design=design,
                        l2fc_threshold=l2fc_threshold,
                        padj_cutoff=padj_cutoff,
                        .parallel=.parallel,
                        ...)
  } else {
    # pairwise subsetting to just control & 1 treatment
    pairwise_expr = as.character(pairwise_expr)
    pairwise_params = as.data.frame(pairwise_params)
    pairwise_params_l = apply(pairwise_params, 1, as.list)

    plyr::ldply(pairwise_params_l, DESeq2_l2fc_MS_sub,
                      ex=pairwise_expr,
                      physeq=physeq,
                      sparsity_threshold=sparsity_threshold,
                      sparsity_apply=sparsity_apply,
                      density_min=density_min,
                      density_max=density_max,
                      design=design,
                      l2fc_threshold=l2fc_threshold,
                      padj_cutoff=padj_cutoff,
                      .parallel=.parallel,
                      ...)
  }

  return(df_l2fc)
}

