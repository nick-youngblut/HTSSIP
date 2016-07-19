# dataset
physeq = readRDS('/home/nick/dev/HTSSIP/data-raw/fullCyc_con-cel-xyl_r1000.rds')
physeq.m = sample_data(physeq)
physeq.p = prune_samples(physeq.m$Substrate %in% c('12C-Con', '13C-Cel') &
                           physeq.m$Day == 6,
                         physeq)
physeq.p = filter_taxa(physeq.p, function(x) sum(x) > 0, TRUE)


# phyloseq2df
test_that('phyloseq data object can be converted to dataframe',{
  cat('\n.phyloseq2df=')
  df1 = phyloseq2df(physeq, sample_data)
  df2 = suppressWarnings(as.data.frame(as.matrix(sample_data(physeq))))
  expect_is(df1, 'data.frame')
  expect_equal(class(df1), class(df2))
  expect_equal(nrow(df1), nrow(df2))
  expect_equal(ncol(df1), ncol(df2))
  })


# .HRSIP
test_that('.HRSIP runs properly',{
  cat('\n.HRSIP=')
  res = .HRSIP(physeq,
               sparsity_threshold=0.25,
               density_min=1.71,
               density_max=1.75,
               design=~Substrate,
               l2fc_threshold=0.25,
               sparsity_apply='all')

  expect_equal(nrow(res), 461)
  expect_equal(ncol(res), 16)
  expect_equal(unique(res$density_min)[1], 1.71)
  expect_equal(unique(res$density_max)[1], 1.75)

  cat('\n.HRSIP (sparsity=heavy)=')
  res = .HRSIP(physeq,
               sparsity_threshold=0.25,
               density_min=1.71,
               density_max=1.75,
               design=~Substrate,
               l2fc_threshold=0.25,
               sparsity_apply='heavy')
  expect_equal(nrow(res), 472)
  expect_equal(ncol(res), 16)
  expect_equal(unique(res$density_min)[1], 1.71)
  expect_equal(unique(res$density_max)[1], 1.75)
})

#
#
# # .HRSIP_pairwise
# ## params
# params = get_treatment_params(physeq, c('Substrate', 'Day'), "Substrate != '12C-Con'")
# ## HRSIP call
# ex = "(Substrate=='12C-Con' & Day=='${Day}') | (Substrate=='${Substrate}' & Day == '${Day}')"
# params_l = apply(params, 1, as.list)
# df_l2fc = plyr::ldply(params_l, .HRSIP_pairwise,
#                       ex=ex,
#                       physeq=physeq,
#                       density_min=1.71,
#                       density_max=1.75,
#                       design=~Substrate,
#                       l2fc_threshold=0.25,
#                       sparsity_threshold=0.25,
#                       sparsity_apply='all')
#
#
#
# # HR-SIP
# ## basic call
# df_l2fc = HRSIP(physeq.p,
#                 density_min=1.71,
#                 density_max=1.75,
#                 design=~Substrate)
# print(nrow(df_l2fc))
# head(df_l2fc)
#
#
# ## HR-SIP in parallel
# library(doParallel)
# registerDoParallel(4)
# df_l2fc = HRSIP(physeq.p,
#                 density_min=1.71,
#                 density_max=1.75,
#                 design=~Substrate,
#                 .parallel=TRUE)
# print(nrow(df_l2fc))
# head(df_l2fc)
#
# ## pairwise HRSIP
# params = get_treatment_params(physeq, c('Substrate', 'Day'), "Substrate != '12C-Con'")
# ex = "(Substrate=='12C-Con' & Day=='${Day}') | (Substrate=='${Substrate}' & Day == '${Day}')"
# df_l2fc = HRSIP(physeq,
#                 pairwise_expr=ex,
#                 pairwise_params=params,
#                 sparsity_threshold=seq(0.1, 0.3, 0.1),
#                 density_min=1.71,
#                 density_max=1.75,
#                 design=~Substrate)
# head(df_l2fc)
#
# ## pairwise HRSIP w/ multiple cores
# library(doParallel)
# registerDoParallel(4)
# df_l2fc = HRSIP(physeq,
#                 pairwise_expr=ex,
#                 pairwise_params=params,
#                 sparsity_threshold=seq(0.1, 0.3, 0.1),
#                 density_min=1.71,
#                 density_max=1.75,
#                 design=~Substrate,
#                 .parallel=TRUE)
# head(df_l2fc)
