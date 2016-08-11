as.Num = function(x) x %>% as.character %>% as.numeric

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
#' df_sample = phyloseq2df(physeq, table_func=sample_data)
#' head(df_sample)
#'
phyloseq2df = function(physeq, table_func){
  physeq.md = table_func(physeq)
  physeq.md = suppressWarnings(as.data.frame(as.matrix(physeq.md)))
  physeq.md = as.matrix(data.frame(lapply(physeq.md, as.character)))
  physeq.md = as.data.frame(apply(physeq.md, 2, trimws))
  rownames(physeq.md) = rownames(table_func(physeq))
  return(physeq.md)
}


#' Phyloseq conversion to a ggplot-formatted table
#'
#' Convert the OTU table (+ metadata) to a format that can be
#' easily plotted with phyloseq
#'
#' @param physeq  Phyloseq object
#' @param include_sample_data  Include \code{sample_table} information?
#' @param sample_col_keep  Which columns in the \code{sample_data} table to keep?
#'   Use \code{NULL} to keep all columns.
#' @param include_tax_table  Include \code{tax_table} information?
#' @param tax_col_keep  Which columns in the \code{tax_table} table to keep?
#'   Use \code{NULL} to keep all columns.
#' @return data.frame
#'
#' @export
#'
#' @examples
#' data(physeq)
#' # Including some columns from sample metadata
#' df_OTU = phyloseq2table(physeq, include_sample_data=TRUE, sample_col_keep=c('Buoyant_density', 'Substrate', 'Day'))
#' head(df_OTU)
#'
#' # Including some columns from sample metadata & taxonomy
#' df_OTU = phyloseq2table(physeq, include_sample_data=TRUE, sample_col_keep=c('Buoyant_density', 'Substrate', 'Day'),include_tax_table=TRUE)
#' head(df_OTU)
#'
phyloseq2table = function(physeq, include_sample_data=FALSE,
                          sample_col_keep=NULL,
                          include_tax_table=FALSE, tax_col_keep=NULL){
  # OTU table
  df_OTU = otu_table(physeq)
  df_OTU = suppressWarnings(as.data.frame(as.matrix(df_OTU)))
  df_OTU$OTU = rownames(df_OTU)
  df_OTU = gather(df_OTU, SAMPLE_JOIN, Count, 1:(ncol(df_OTU)-1))

  # sample metdata
  if(include_sample_data==TRUE){
    df_meta = sample_data(physeq)
    df_meta = suppressWarnings(as.data.frame(as.matrix(df_meta)))
    df_meta$SAMPLE_JOIN = rownames(df_meta)
    ## trimming
    if(!is.null(sample_col_keep)){
      sample_col_keep = c('SAMPLE_JOIN', sample_col_keep)
      df_meta = dplyr::select_(df_meta, .dots=as.list(sample_col_keep))
    }
    # join
    df_OTU = inner_join(df_OTU, df_meta, c('SAMPLE_JOIN'))
  }

  # sample metdata
  if(include_tax_table==TRUE){
    df_tax = tax_table(physeq)
    df_tax = suppressWarnings(as.data.frame(as.matrix(df_tax)))
    df_tax$OTU = rownames(df_tax)
    ## trimming
    if(!is.null(tax_col_keep)){
      tax_col_keep = c('OTU', tax_col_keep)
      df_tax = dplyr::select_(df_tax, .dots=as.list(tax_col_keep))
    }
    # join
    df_OTU = inner_join(df_OTU, df_tax, c('OTU'))
  }

  return(df_OTU)
}

#' Make a list of phyloseq object subsets
#'
#' Create a list of phyloseq object subsets based on phyloseq
#' sample data parameters (e.g., a phyloseq subset for each treatment)
#'
#' @param physeq  Phyloseq object
#' @param params  data.frame of parameters supplies to \code{ex}
#' @param ex  Expression for subsetting the phyloseq object
#'
#' @return A list of Phyloseq objects
#'
#' @export
#'
#' @examples
#' data(physeq)
#' # making subsets by substrate and time point
#' params = get_treatment_params(physeq, c('Substrate', 'Day'))
#' params = dplyr::filter(params, Substrate!='12C-Con')
#' ex = "(Substrate=='12C-Con' & Day=='${Day}') | (Substrate=='${Substrate}' & Day == '${Day}')"
#' physeq_l = phyloseq_subset(physeq, params, ex)
#'
phyloseq_subset = function(physeq, params, ex){
  if(is.data.frame(params)){
    params = apply(params, 1, as.list)
  }
  # subsetting phyloseq object
  physeq_l = lapply(params, function(x){
    # x should be a list of parameters
    exx = stringterpolate(ex, x)
    physeq.m = phyloseq2df(physeq, phyloseq::sample_data)
    bool = mutate_(physeq.m, exx)[,ncol(physeq.m)+1]
    phyloseq::prune_samples(bool, physeq)
  })
  # adding names to list
  ## names based on subsetting expression
  n = lapply(params, function(x){
    # x should be a list of parameters
    exx = stringterpolate(ex, x)
  })
  names(physeq_l) = do.call(rbind, n)

  return(physeq_l)
}
