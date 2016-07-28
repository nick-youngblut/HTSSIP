# calculating pairwise beta-div between fraction communities

fraction_beta_div = function(physeq, params, ex){
  #if(is.null(subset)){
  #  df_OTU = phyloseq2table(physeq)
 # } else {
  #  df_OTU = phyloseq2table(physeq, include_sample_data=TRUE, sample_col_keep=physeq_subset)
    #df_OTU = group_by_(df_OTU, .dots=as.list(physeq_subset))
  #}

  # making a list of phyloseq subsets
  bool_l = lapply(params, function(x){
    # x should be a list of parameters
    exx = stringterpolate(ex, x)
    physeq.m = phyloseq2df(physeq, phyloseq::sample_data)
    bool = mutate_(physeq.m, exx)[,ncol(physeq.m)+1]
    phyloseq::prune_samples(bool, physeq)
  })

  # subsetting phyloseq


  # wunif.dist = distance(physeq,
  #                       method = "unifrac",
  #                       weighted = TRUE,
  #                       fast = TRUE,
  #                       parallel = TRUE,
  #                       normalized = FALSE)
  #plyr::dlply(df_OTU, physeq_subset, )
}

# data(physeq)
# params = get_treatment_params(physeq, c('Substrate', 'Day'))
# params_l = apply(params, 1, as.list)
# ex = "Substrate=='${Substrate}' & Day == '${Day}'"
#
# res = fraction_beta_div(physeq, params_l, ex)
# head(res)

# calc.wunif.dist = function(physeq, cores=cores){
#   message('Processing sample...')
#   registerDoParallel(cores=cores)
#   wunif.dist = distance(physeq,
#                         method = "unifrac",
#                         weighted = TRUE,
#                         fast = TRUE,
#                         parallel = TRUE,
#                         normalized = FALSE)
#   return(wunif.dist)
# }

# for (d in as.character(u.day)){
#   # filter by Day and removing h2o controls & mock communities
#   tmp = subset_samples(physeq.SIP, Day == d)
#   physeq.SIP.l[[d]] = subset_samples(tmp, ! is.na(Buoyant_density))
#   print(d)
#   physeq.SIP.l[[d]] %>% print
# }
