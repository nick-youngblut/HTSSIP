test_that('DESeq2_l2fc_multi_sub working with basic parameters',{
  windows = data.frame(start=c(1.70, 1.72), end=c(1.73, 1.75))
  params = get_treatment_params(physeq, c('Substrate', 'Day'), "Substrate != '12C-Con'")
  params_l = apply(params, 1, as.list)
  ## HRSIP call
  ex = "(Substrate=='12C-Con' & Day=='${Day}') | (Substrate=='${Substrate}' & Day == '${Day}')"
  df_l2fc = plyr::ldply(params_l, DESeq2_l2fc_multi_sub,
                        ex=ex,
                        physeq=physeq,
                        sparsity_threshold = c(0, 0.1),
                        density_windows=windows,
                        design=~Substrate,
                        l2fc_threshold=0.25,
                        sparsity_apply='all')

  expect_is(df_l2fc, 'data.frame')
  expect_false(is.null(df_l2fc$Substrate))
  expect_false(is.null(df_l2fc$Day))
  expect_equal(df_l2fc$Substrate %>% unique %>% length, 2)
  expect_equal(df_l2fc$Day %>% unique %>% length, 1)
})
