#' Total sum scaling
#'
#' @param x  data.frame of numeric values
#' @param MARGIN  table margin (1=rows, 2=columns)
#' @param na.rm  remove NAs?
#'
#' @return data.frame of qPCR values
#'
#' @export
#'
#' @examples
#' # making functions for simulating values
#' df = data.frame(1:5, 5:9)
#' df_t = tss(df)
#' apply(df_t, 2, sum)
#'
tss = function(x, MARGIN=2, na.rm=FALSE){
  k = min(x, na.rm = na.rm)
  tmp = pmax(k, apply(x, MARGIN, sum, na.rm = na.rm))
  x = sweep(x, MARGIN, tmp, "/")
  return(x)
}

#' Transform OTU counts based on qPCR data
#'
#' OTU counts in the phyloseq otu_table object will be normalized
#' to sample totals (total sum scaling), then multiplied by the
#' qPCR value associated with each sample. Thus, the qPCR table
#' should have ONE value matching the OTU count table. Value
#' matching between the OTU table & qPCR value table to set by
#' \code{sample_idx()}.
#'
#'
#' @param physeq  A phyloseq object
#' @param qPCR  A list of qPCR data from \code{qPCR_sim()}
#' @param sample_idx  The qPCR table column index for
#' matching to otu table samples.
#' @param value_idx  The qPCR table column index for qPCR values.
#'
#' @return A phyloseq object with transformed OTU counts
#'
#' @export
#'
#' @examples
#' # qPCR data simulation
#' data(physeq_rep3)
#' data(phsyeq_rep3_qPCR)
#' physeq_rep3_t = OTU_qPCR_trans(physeq_rep3, physeq_rep3_qPCR)
#'
OTU_qPCR_trans = function(physeq, qPCR,
                          sample_idx='Sample',
                          value_idx='qPCR_tech_rep_mean'){
  # means of qPCR (if needed)
  if(!is.null(qPCR$summary)){
    qPCR = qPCR$summary
  }
  stopifnot(class(qPCR)=='data.frame' | class(qPCR)=='matrix')
  stopifnot(!is.null(qPCR$Sample))

  # OTU table
  df_OTU_col = physeq %>% otu_table %>% colnames
  df_OTU = phyloseq2df(physeq, otu_table)
  df_OTU_rn = rownames(df_OTU)
  df_OTU = as.data.frame(apply(df_OTU, 2, as.Num))
  rownames(df_OTU) = df_OTU_rn

  # sum scale transformation
  df_OTU = tss(df_OTU)

  # qPCR multiplication
  ## ordering of qPCR table
  rownames(qPCR) = make.names(qPCR[,sample_idx])
  qPCR = qPCR[colnames(df_OTU),]
  qPCR_vals = qPCR[,value_idx]
  if(length(qPCR_vals) != ncol(df_OTU)){
    stop('length qPCR_vals (', length(qPCR_vals),
         ') != ncol df_OTU (', ncol(df_OTU), ')')
  }
  ## transform
  df_OTU = sweep(df_OTU %>% as.data.frame, 2, qPCR_vals, "*")
  df_OTU = apply(df_OTU, 2, function(x) round(x, 0))
  colnames(df_OTU) = df_OTU_col

  # nan to zero
  df_OTU[is.na(df_OTU)] = 0

  # making new physeq object
  tree = phy_tree(physeq, errorIfNULL=FALSE)
  tax  = tax_table(physeq, errorIfNULL=FALSE)
  sam  = sample_data(physeq, errorIfNULL=FALSE)

  physeq2 = phyloseq(otu_table(df_OTU, taxa_are_rows=TRUE),
                    phy_tree(tree, errorIfNULL=FALSE),
                    tax_table(tax, errorIfNULL=FALSE),
                    sample_data(sam, errorIfNULL=FALSE))

  return(physeq2)
}

