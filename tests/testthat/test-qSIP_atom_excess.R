test_that('qSIP_BD_shift working', {
  # control_mean_fun = function(x) dnorm(x, mean=1.70, sd=0.01) * 1e8
  # control_sd_fun = function(x) control_mean_fun(x) / 3
  # treat_mean_fun = function(x) dnorm(x, mean=1.75, sd=0.01) * 1e8
  # treat_sd_fun = function(x) treat_mean_fun(x) / 3
  #
  # physeq = physeq_S2D2_l[[1]]
  # qPCR = qPCR_sim(physeq,
  #                 control_expr='Substrate=="12C-Con"',
  #                 control_mean_fun=control_mean_fun,
  #                 control_sd_fun=control_sd_fun,
  #                 treat_mean_fun=treat_mean_fun,
  #                 treat_sd_fun=treat_sd_fun)
  # OTU table transformation
  physeq_rep3_t = OTU_qPCR_trans(physeq_rep3, physeq_rep3_qPCR)

  # BD shift (Z) & atom excess (A)
  atomX = qSIP_atom_excess(physeq_rep3_t,
                           control_expr='Treatment=="12C-Con"',
                           treatment_rep='Replicate')

  #df_atomX %>% head
  expect_is(atomX, 'list')
  expect_false(is.null(atomX$A$Z))
  expect_false(is.null(atomX$A$A))
  expect_equal(atomX$A$OTU %>% length,
               atomX$A$OTU %>% unique %>% length)
})

test_that('bootstrap iteration is working', {
  atomX_boot = .qSIP_bootstrap(atomX)
  atomX_boot %>% head
})

test_that('bootstrap in parallel working', {
  doParallel::registerDoParallel(10)
  df_atomX_boot = qSIP_bootstrap(atomX, parallel=TRUE)
  df_atomX_boot %>% head(n=3) %>% as.data.frame
  df_atomX_boot %>% dplyr::filter(OTU=='OTU.1')
})

# calculating CIs
# a=0.1
# df_CI = df_atomX_boot %>%
#   dplyr::group_by(OTU) %>%
#   summarize(A_CI_low = quantile(A, a / 2, na.rm=TRUE),
#             A_CI_high = quantile(A, 1 - a/2, na.rm=TRUE))
#
# df_atomX_j = dplyr::inner_join(atomX$A, df_CI, c('OTU'='OTU'))
# df_atomX_j %>% head(n=3)
# df_atomX_j %>% nrow

#atomX_shuf %>% head(n=10)

#df = atomX$W
#n_sample = c(3,3)



# making shuffled dataset (shuffling between control & treatment)
# physeq = physeq_S2D2_l[[1]]
# physeq_t = OTU_qPCR_trans(physeq, qPCR)
# cols = c('IS_CONTROL', 'Buoyant_density', gradient_rep)
# control_expr='Substrate=="12C-Con"'
# df_OTU = phyloseq2table(physeq_t,
#                         include_sample_data=TRUE,
#                         sample_col_keep=cols,
#                         control_expr=control_expr)
# #head(df_OTU)
# # how many to subsample?
# df_OTU_s = df_OTU %>%
#   distinct(IS_CONTROL, Microcosm_replicate) %>%
#   group_by(IS_CONTROL) %>%
#   summarize(n_tech_reps = n()) %>%
#   as.data.frame
# df_OTU_j = dplyr::inner_join(df_OTU, df_OTU_s, c('IS_CONTROL'='IS_CONTROL'))
#
# df_OTU_j %>% head

# sampling with replacement among replicates
## sampling weighted mean densities (W)
# df_OTU_f = df_OTU %>%
#   mutate(IS_CONTROL = sample(IS_CONTROL, size=nrow(df_OTU)))


