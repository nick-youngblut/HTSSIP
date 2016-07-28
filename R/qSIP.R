# simulating qPCR data for a phyloseq object
OTU_SIP_rank_abunds = function(df_OTU, exp_design){
  exp_design_plus = c('OTU', exp_design)
  df_OTU %>%
    group_by_(.dots=as.list(exp_design_plus)) %>%
    mutate(Max_count = max(Count)) %>%
    group_by_(.dots=as.list(exp_design)) %>%
    mutate(Rank_count = -Max_count %>% as.factor %>% as.numeric) %>%
    ungroup()
}

# df_OTU = phyloseq2table(physeq, c('Buoyant_density', 'Substrate', 'Day'))
# df_OTU = OTU_SIP_rank_abunds(df_OTU, c('Substrate', 'Day'))
# df_OTU %>% head

qPCR_sim = function(physeq, mean_fun, sd_fun, ...){
  # using OTU max count of any gradient fraction as mean

  sd_fun(...)
}

#mean_fun = function(x) mu + 10
#sd_fun = function(mu) 5889 + mu + 0.714 * mu**2
#qPCR_sim(physeq, mean_fun, sd_fun, mu=10)
