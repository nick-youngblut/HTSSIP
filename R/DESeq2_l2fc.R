#' Calculating log2 fold change for HTS-SIP data.
#'
#' The phyloseq object will be filtered to 1) just OTUs
#' that pass the sparsity cutoff 2) just samples in the user-defined
#' 'heavy' fractions. The log2 fold change (l2fc) is calculated
#' between labeled treatment and control gradients.
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
#' @export
#'
#' @examples
#' data(physeq_S2D2)
#' \dontrun{
#' df_l2fc = DESeq2_l2fc(physeq_S2D2, density_min=1.71, density_max=1.75, design=~Substrate)
#' head(df_l2fc)
#' }
#'
DESeq2_l2fc = function(physeq, density_min, density_max, design,
                       l2fc_threshold=0.25, sparsity_threshold=0.25,
                       sparsity_apply='all'){
  # assertions
  l2fc_threshold = as.numeric(l2fc_threshold)
  stopifnot(l2fc_threshold >= 0 & l2fc_threshold <= 1)
  sparsity_apply = tolower(sparsity_apply)
  stopifnot(sparsity_apply %in% c('all', 'heavy'))
  physeq.md = phyloseq::sample_data(physeq)
  stopifnot(!is.null(physeq.md$Buoyant_density))

  # status
  cat('Sparsity threshold:', sparsity_threshold, '\n')
  cat('Density window:', paste(c(density_min, density_max), collapse='-'), '\n')

  # sparsity cutoff applied to all gradient fractions
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
  p = stats::pnorm(beta, theta, betaSE, lower.tail=FALSE)
  res$p = p
  d = data.frame(res[, c('OTU','log2FoldChange', 'p')])

  # p-value adjust
  d$padj = stats::p.adjust(p, method = 'BH')

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
  d$l2fc_threshold = l2fc_threshold
  return(d)
}

