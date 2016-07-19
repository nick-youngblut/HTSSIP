#' phyloseq data object conversion to data.frame
#'
#' Conducts conversion of 1 of the data objects
#' in a phyloseq object (eg., tax_table) to a dataframe
#'
#' @param physeq Phyloseq object
#' @param table_func See Phyloseq::phyloseq-class for options
#'   (example: table_func=otu_table)
#' @return dataframe
#'
#' @examples
#' phyloseq2df(physeq)
#'
# conversion of phyloseq data table to dataframe
phyloseq2df = function(physeq, table_func){
  physeq.md = table_func(physeq)
  physeq.md = suppressWarnings(as.data.frame(as.matrix(physeq.m)))
  physeq.md = as.matrix(data.frame(lapply(physeq.md,as.character)))
  physeq.md = as.data.frame(apply(physeq.md, 2, trimws))
  return(physeq.md)
}


#' HR-SIP analysis
#'
#' \code{HRSIP} conducts HR-SIP analysis on phyloseq object
#'
#' blah blah blah
#'
#' @param physeq Phyloseq object
#' @param l2fc_threshold DESeq2 log2 fold change threshold
#' @param sparsity_threshold All OTUs observed in less than this portion (fraction: 0-1)
#'   of gradient fraction samples are pruned
#' @param sparsity_apply Apply sparsity threshold to all gradient fraction samples ('all')
#'   or just heavy fraction samples ('heavy')
#' @return dataframe of HRSIP results
#'
#' @examples
#' HRSIP()
#'
.HRSIP = function(physeq,
                   sparsity_threshold,
                   density_min, density_max,
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


#-- HR-SIP: pairwise of comparisions--#
# Making a df of all of the HRSIP tests to run (boolean vectors)
get_treatment_params = function(physeq, exp_params, control){
  physeq.m = phyloseq2df(physeq, sample_data)
  # filter out control
  physeq.m = filter_(physeq.m, control)
  # all pairwise params
  params = dplyr::distinct(physeq.m[,exp_params])
  return(params)
}


#---- HR-SIP: all pairwise ---#
.HRSIP_pairwise = function(params, ex, physeq,
                           density_min=1.71,
                           density_max=1.75,
                           design=~Substrate,
                           l2fc_threshold=0.25,
                           sparsity_threshold=0.25,
                           sparsity_apply='all'){
  # Pruning phyloseq to pairwise comparison
  exx = stringterpolate(ex, params)
  cat(exx, '\n')
  physeq.m = phyloseq2df(physeq, phyloseq::sample_data)
  bool = mutate_(physeq.m, exx)[,ncol(physeq.m)+1]
  physeq.p = phyloseq::prune_samples(bool, physeq)
  ## HRSIP
  df_l2fc = .HRSIP(physeq.p,
                    density_min=density_min,
                    density_max=density_max,
                    design=design,
                    l2fc_threshold=l2fc_threshold,
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
HRSIP = function(physeq,
                  density_min,
                  density_max,
                  design,
                  sparsity_threshold=seq(0, 0.5, 0.1),
                  sparsity_apply='all',
                  l2fc_threshold=0.25,
                  pairwise_expr=NULL,
                  pairwise_params=NULL,
                  .parallel=FALSE,
                  ...){
  sparsity_threshold = as.array(sparsity_threshold)

  # pairwise or not
  if(is.null(pairwise_expr)){
    df_l2fc = plyr::adply(sparsity_threshold, 1, function(x){
                  .HRSIP(physeq=physeq,
                        sparsity_threshold=x,
                        sparsity_apply=sparsity_apply,
                        density_min=density_min,
                        density_max=density_max,
                        design=design,
                        l2fc_threshold=l2fc_threshold,
                        ...)
    }, .parallel=.parallel)
  } else {
    # pairwise subsetting to just control & 1 treatment
    pairwise_expr = as.character(pairwise_expr)
    pairwise_params = as.data.frame(pairwise_params)
    pairwise_params_l = apply(pairwise_params, 1, as.list)

    df_l2fc = plyr::adply(sparsity_threshold, 1, function(x){
          plyr::ldply(pairwise_params_l, .HRSIP_pairwise,
                      ex=pairwise_expr,
                      physeq=physeq,
                      sparsity_threshold=x,
                      sparsity_apply=sparsity_apply,
                      density_min=density_min,
                      density_max=density_max,
                      design=design,
                      l2fc_threshold=l2fc_threshold,
                      ...)
    }, .parallel=.parallel)
  }

  # TODO: selecting cutoff with most hypotheses rejected
  return(df_l2fc)
}



