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
#' phyloseq2df(physeq, table_func=otu_table)
#'
# conversion of phyloseq data table to dataframe
phyloseq2df = function(physeq, table_func){
  physeq.md = table_func(physeq)
  physeq.md = suppressWarnings(as.data.frame(as.matrix(physeq.md)))
  physeq.md = as.matrix(data.frame(lapply(physeq.md, as.character)))
  physeq.md = as.data.frame(apply(physeq.md, 2, trimws))
  return(physeq.md)
}


DESeq2_l2fc = function(physeq,
                  sparsity_threshold,
                  density_min,
                  density_max,
                  design,
                  l2fc_threshold=0.25,
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


#--- HR-SIP: testing multiple sparsity thresholds ---#
#' HR-
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
#' phyloseq2df(physeq, table_func=otu_table)
#'
# conversion of phyloseq data table to dataframe
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


#--- HR-SIP: all pairwise ---#
# Making a df of all of the HRSIP tests to run (boolean vectors)
get_treatment_params = function(physeq, exp_params, control){
  physeq.m = phyloseq2df(physeq, phyloseq::sample_data)
  # filter out control
  physeq.m = filter_(physeq.m, control)
  # all pairwise params
  params = dplyr::distinct(physeq.m[,exp_params])
  return(params)
}


# HRSIP
DESeq2_l2fc_MS_PW = function(params, ex, physeq,
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



#-- full HR-SIP
#' HR-SIP analysis
#'
#' \code{HRSIP} conducts high resolution stable isotope probing (HR-SIP)
#'   analysis on phyloseq object.
#'
#' @param physeq  Phyloseq object
#' @param density_min  Minimum buoyant density of the 'heavy' gradient fractions
#' @param density_max  Maximum buoyant density of the 'heavy' gradient fractions
#' @param design  \code{design} parameter used for DESeq2 analysis.
#'   See \code{DESeq2::DESeq} for more details.
#' @param l2fc_threshold  DESeq2 log2 fold change threshold
#' @param sparsity_threshold  A vector of values (range=0-1).
#'   All OTUs observed in less than this portion (fraction: 0-1)
#'   of gradient fraction samples are pruned. A a form of indepedent filtering,
#'   The sparsity cutoff with the most rejected hypotheses is used.
#' @param sparsity_apply  Apply sparsity threshold to all gradient fraction samples ('all')
#'   or just heavy fraction samples ('heavy')
#' @param l2fc_threshold  log2 fold change (l2fc) values must be significantly above this
#'   threshold in order to reject the hypothesis of equal counts.
#' @return dataframe of HRSIP results
#'
#' @export
#'
#' @examples
#' HRSIP()
#'
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

    plyr::ldply(pairwise_params_l, DESeq2_l2fc_MS_PW,
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

