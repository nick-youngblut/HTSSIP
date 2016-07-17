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

  # window selection
  prn = function(x) sum(x > 0) > sparsity_threshold * length(x)

  # sparcity cutoff applied to all gradient fractions
  if(sparsity_apply=='all'){
    physeq = phyloseq::filter_taxa(physeq, prn, TRUE)
  }

  # applying 'heavy' window pruning
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
  if(!is.null(phyloseq::tax_table(physeq))){
    TT = phyloseq::tax_table(physeq.m)
    TT = as.data.frame(tax_table)
    TT$OTU = rownames(TT)
    d = dplyr::left_join(d, TT, c('OTU'))
  }

  # setting pruning info
  d$density_min = density_min
  d$density_max = density_max
  d$sparsity_threshold = sparsity_threshold
  return(d)
}

# example
physeq = readRDS('/home/nick/dev/HTSSIP/data-raw/fullCyc_con-cel-xyl.rds')
res = .HRSIP(physeq, sparsity_threshold=0.25,
         density_min=1.71, density_max=1.75,
        design=~Substrate,
        l2fc_threshold=0.25,
        sparsity_apply='all')
head(res)
