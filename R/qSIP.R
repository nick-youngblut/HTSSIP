# simulating qPCR data for a phyloseq object

OTU_SIP_rank_abunds = function(OTU, exp_design){
  OTU %>%
    group_by_(exp_design) %>%
    summarize(n=n())
}

OTU = phyloseq2df(physeq, otu_table)
OTU_SIP_rank_abunds(OTU, c('Substrate', 'Day'))


qPCR_sim = function(physeq, mean_fun, sd_fun, ...){
  # using OTU max count of any gradient fraction as
  sd_fun(...)
}

mean_fun = function(x) mu + 10
sd_fun = function(mu) 5889 + mu + 0.714 * mu**2
qPCR_sim(physeq, mean_fun, sd_fun, mu=10)
