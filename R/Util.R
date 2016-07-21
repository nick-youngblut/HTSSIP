#' phyloseq data object conversion to data.frame
#'
#' Conducts conversion of 1 of the data objects
#' in a phyloseq object (eg., tax_table) to a dataframe
#'
#' @param physeq  Phyloseq object
#' @param table_func  See \code{Phyloseq::phyloseq-class} for options
#' @return data.frame
#'
#' @export
#'
#' @examples
#' data(physeq)
#' df_otu = phyloseq2df(physeq, table_func=otu_table)
#' head(df_otu)
#'
phyloseq2df = function(physeq, table_func){
  physeq.md = table_func(physeq)
  physeq.md = suppressWarnings(as.data.frame(as.matrix(physeq.md)))
  physeq.md = as.matrix(data.frame(lapply(physeq.md, as.character)))
  physeq.md = as.data.frame(apply(physeq.md, 2, trimws))
  return(physeq.md)
}


#' Phyloseq OTU table conversion to a qSIP table
#'
#'
#' @param physeq  Phyloseq object
#' @param table_func  See \code{Phyloseq::phyloseq-class} for options
#' @return data.frame
#'
#' @export
#'
#' @examples
#' data(physeq)
#'
phyloseq2qSIP = function(physeq, col.keep=NULL){
  # OTU table
  df_OTU = otu_table(physeq)
  df_OTU = suppressWarnings(as.data.frame(as.matrix(df_OTU)))
  df_OTU$OTU = rownames(df_OTU)
  df_OTU = gather(df_OTU, SAMPLE_JOIN, Count, 1:(ncol(df_OTU)-1))

  # sample metdata
  df_meta = sample_data(physeq)
  df_meta = suppressWarnings(as.data.frame(as.matrix(df_meta)))
  df_meta$SAMPLE_JOIN = rownames(df_meta)
  ## trimming
  if(!is.null(col.keep)){
    col.keep = c('SAMPLE_JOIN', col.keep)
    print(col.keep)
    df_meta = dplyr::select_(df_meta, col.keep)
    head(df_meta) %>% print
  }

  # join
  df_OTU = inner_join(df_OTU, df_meta, c('SAMPLE_JOIN'))

  return(df_OTU)
}

df_OTU = phyloseq2qSIP(physeq, c('OTU', 'Count', 'Buoyant_density', 'Substrate', 'Day'))
head(df_OTU)
