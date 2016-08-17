# makes 1 gradient
gradient_sim = function(locs, params,
                        responseModel='gaussian',
                        countModel='poisson',
                        #sample_prefix ='sample',
                        ...){
  df_OTU = coenocliner::coenocline(locs,
                                   params=params,
                                   responseModel=responseModel,
                                   countModel=countModel,
                                   ...)
  df_OTU = as.data.frame(df_OTU)
  colnames(df_OTU) = gsub('^', 'OTU.', 1:M)
  df_OTU$Buoyant_density = locs #as.character(round(locs, digits=4))
  return(df_OTU)
}



# making whole dataset
phyloseq_sim = function(locs, params,
                   responseModel='gaussian',
                   countModel='poisson',
                   parallel=FALSE,
                   ...){

  #samples = data.frame(sample_prefix=samples)
  df_OTU = plyr::ldply(params, gradient_sim,
                       locs=locs,
                       responseModel=responseModel,
                       countModel=countModel,
                       .parallel=parallel,
                       .id='Gradient',
                       ...)

  # vary the BDs a bit
  x = rnorm(nrow(df_OTU), mean=0, sd=0.002)
  df_OTU$Buoyant_density = as.Num(df_OTU$Buoyant_density) + x
  df_OTU$Buoyant_density = round(df_OTU$Buoyant_density, digits=4)

  # metadata
  df_meta = df_OTU[,c('Gradient', 'Buoyant_density')]

  # formatting OTU table
  rownames(df_OTU) = apply(df_meta, 1, paste, collapse="_")
  df_OTU$Gradient = NULL
  df_OTU$Buoyant_density = NULL
  df_OTU = t(df_OTU)

  # sample matching between metdata &  OTU table
  rownames(df_meta) = colnames(df_OTU)

  # making phyloseq object
  physeq = phyloseq::phyloseq(
    phyloseq::otu_table(df_OTU, taxa_are_rows=TRUE),
    phyloseq::sample_data(df_meta)
  )

  return(physeq)
}





