# conversion of phyloseq data table to dataframe
phyloseq2df = function(physeq, table_func){
  physeq.md = table_func(physeq)
  physeq.md = suppressWarnings(as.data.frame(as.matrix(physeq.m)))
  physeq.md = as.matrix(data.frame(lapply(physeq.md,as.character)))
  physeq.md = as.data.frame(apply(physeq.md, 2, trimws))
  return(physeq.md)
}
# test
df = phyloseq2df(physeq, sample_data)
head(df)


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

# example
## dataset
physeq = readRDS('/home/nick/dev/HTSSIP/data-raw/fullCyc_con-cel-xyl_r1000.rds')
physeq.m = sample_data(physeq)
physeq.p = prune_samples(physeq.m$Substrate %in% c('12C-Con', '13C-Cel') &
                         physeq.m$Day == 6,
                         physeq)
physeq.p = filter_taxa(physeq.p, function(x) sum(x) > 0, TRUE)
## HR-SIP on 1 sample & substrate
res = .HRSIP(physeq, sparsity_threshold=0.25,
         density_min=1.71, density_max=1.75,
        design=~Substrate,
        l2fc_threshold=0.25,
        sparsity_apply='all')
head(res)


#' All comparisons between controls and treatments
#'
#' \code{comparisons} conducts HR-SIP analysis on phyloseq object
#'
#' blah blah blah

.pairwise = function(physeq, control_ids, sample_id_column, sel_columns=c('Substrate', 'Day'), join_columns=c('Day')){
  physeq.m = suppressWarnings(as.data.frame(as.matrix(sample_data(physeq))))

  # control / treatment
  subject = as.vector(physeq.m[,sample_id_column])
  physeq.m.con = physeq.m[subject %in% control_ids,sel_columns]
  physeq.m.trt = physeq.m[!subject %in% control_ids,sel_columns]
  physeq.j = dplyr::inner_join(physeq.m.con, physeq.m.trt, join_columns)
  return(physeq.j)
}

# test
#physeq = readRDS('/home/nick/dev/HTSSIP/data-raw/fullCyc_con-cel-xyl_r1000.rds')
physeq.m = sample_data(physeq)
ctr_ids = physeq.m$X.Sample[grepl('12C-Con', physeq.m$X.Sample)]
ctr_ids = unique(as.vector(ctr_ids))
res = .pairwise (physeq, control_ids=ctr_ids, sample_id_column='X.Sample')
res %>% tail


#-- all pairwise of comparisions of HR-SIP --#
## TODO: figure out how to make a matrix of logical vectors, which correspond to each HR-SIP comparison

# NS evalulation; formatting query
.filter_ = function(df, x, invert=FALSE){
  #x = c('Substrate'='12C-Con', 'Day'=6)
  ## query
  y = as.character(names(x))
  z = as.character(x)
  z = gsub('(.+)', '"\\1"', z)
  if(invert==TRUE){
    xx = apply(rbind(y,z), 2, function(x) paste(x, collapse=' != '))
  } else {
    xx = apply(rbind(y,z), 2, function(x) paste(x, collapse=' == '))
  }
  yy = paste(xx, collapse=' & ')

  ##
  print(yy)
  dplyr::filter_(df, yy)
}

## test
physeq.md = phyloseq2df(physeq, sample_data)
x = .filter_(physeq.md, c('Substrate'='12C-Con', 'Day'=6))
head(x)


# Making a matrix/df of all of the HRSIP tests to run (boolean vectors)
comparisons = function(physeq, exp_params, control){
  physeq.m = phyloseq2df(physeq, sample_data)

  # all pairwise params
  params = dplyr::distinct(physeq.m[,exp_params])
  # filter out control
  params = .filter_(params, control, invert=TRUE)

  return(params)
  #physeq.m %>% head %>% print
}

# test
#physeq = readRDS('/home/nick/dev/HTSSIP/data-raw/fullCyc_con-cel-xyl_r1000.rds')
res = comparisons(physeq, c('Substrate', 'Day'), c('Substrate'='12C-Con'))
res %>% tail






#-- full HR-SIP
HRSIP = function(physeq.control,
                 physeq.treatment,
                 sparsity_threshold,
                 density_min,
                 density_max,
                  design,
                  pairwise=FALSE,
                  l2fc_threshold=0.25,
                  sparsity_apply='all', ...){

  # merging

  if(pairwise=TRUE){
    # pairwise subsetting to just control & 1 treatment

    #physeq.m = merge_phyloseq(physeq.control, physeq.treatment)
  } else {

  }

}

