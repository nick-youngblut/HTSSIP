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
#' @param physeq  A phyloseq object
#' @param qPCR  A list of qPCR data from \code{qPCR_sim()}
#'
#' @return A phyloseq object with transformed OTU counts
#'
#' @export
#'
#' @examples
#' # qPCR data simulation
#' ## making functions for simulating values
#' control_mean_fun = function(x) dnorm(x, mean=1.70, sd=0.01) * 1e8
#' control_sd_fun = function(x) control_mean_fun(x) / 3
#' treat_mean_fun = function(x) dnorm(x, mean=1.75, sd=0.01) * 1e8
#' treat_sd_fun = function(x) treat_mean_fun(x) / 3
#' ## simulating qPCR values
#' qPCR = qPCR_sim(physeq_S2D2,
#'                 control_expr='Substrate=="12C-Con"',
#'                 control_mean_fun=control_mean_fun,
#'                 control_sd_fun=control_sd_fun,
#'                 treat_mean_fun=treat_mean_fun,
#'                 treat_sd_fun=treat_sd_fun)
#' ## transforming OTU values
#' OTU_qPCR_trans(physeq, qPCR)
#'
OTU_qPCR_trans = function(physeq, qPCR){
  # means of qPCR (if needed)
  stopifnot(class(qPCR)=='list')
  stopifnot(!is.null(qPCR$summary))

  # OTU table
  df_OTU = phyloseq2df(physeq, otu_table)
  df_OTU_rn = rownames(df_OTU)
  df_OTU = as.data.frame(apply(df_OTU, 2, as.Num))
  rownames(df_OTU) = df_OTU_rn

  # sum scale transformation
  df_OTU = tss(df_OTU)

  # qPCR multiplication
  ## ordering of qPCR table
  qPCR_s = qPCR$summary
  stopifnot(!is.null(qPCR_s$Sample))
  rownames(qPCR_s) = make.names(qPCR_s$Sample)
  qPCR_s = qPCR_s[colnames(df_OTU),]
  qPCR_vals = qPCR_s$qPCR_tech_rep_mean
  if(length(qPCR_vals) != ncol(df_OTU)){
    stop('length qPCR_vals (', length(qPCR_vals),
         ') != ncol df_OTU (', ncol(df_OTU), ')')
  }
  ## transform
  df_OTU = sweep(df_OTU %>% as.data.frame, 2, qPCR_s$qPCR_tech_rep_mean, "*")

  # making new physeq object
  tree = phy_tree(physeq, errorIfNULL=FALSE)
  tax  = tax_table(physeq, errorIfNULL=FALSE)
  sam  = sample_data(physeq, errorIfNULL=FALSE)
  rownames(sam) = make.names(rownames(sam))

  physeq2 = phyloseq(otu_table(df_OTU, taxa_are_rows=TRUE),
                    phy_tree(tree, errorIfNULL=FALSE),
                    tax_table(tax, errorIfNULL=FALSE),
                    sample_data(sam, errorIfNULL=FALSE))

  return(physeq2)
}

