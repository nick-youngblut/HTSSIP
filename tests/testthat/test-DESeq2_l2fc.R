test_that('DESeq2_l2fc runs with default params',{
  df_l2fc = DESeq2_l2fc(physeq,
                        sparsity_threshold=0.25,
                        density_min=1.71,
                        density_max=1.75,
                        design=~Substrate)

  expect_is(df_l2fc, 'data.frame')
  expect_equal(ncol(df_l2fc), 16)
  expect_equal(unique(df_l2fc$density_min)[1], 1.71)
  expect_equal(unique(df_l2fc$density_max)[1], 1.75)
})

test_that('DESeq2_l2fc runs with sparsity_apply=heavy',{
  df_l2fc = DESeq2_l2fc(physeq,
                        sparsity_threshold=0.25,
                        density_min=1.71,
                        density_max=1.75,
                        design=~Substrate,
                        sparsity_apply='heavy')

  expect_is(df_l2fc, 'data.frame')
  expect_equal(ncol(df_l2fc), 16)
  expect_equal(unique(df_l2fc$density_min)[1], 1.71)
  expect_equal(unique(df_l2fc$density_max)[1], 1.75)
})

test_that('DESeq2_l2fc_multi working with basic parameters',{
  windows = data.frame(start=c(1.70, 1.72), end=c(1.73, 1.75))
  df_l2fc = DESeq2_l2fc_multi(physeq,
                              sparsity_threshold = c(0, 0.1),
                              density_windows=windows,
                              padj_cutoff=0.1,
                              design=~Substrate)

  expect_is(df_l2fc, 'data.frame')
  expect_true(is.null(df_l2fc$rej_hypo))
  expect_true(is.null(df_l2fc$n_rej_hypo))
  expect_equal(length(unique(df_l2fc$sparsity_threshold)), 2)
})


test_that('DESeq2_l2fc_multi_sub working with basic parameters',{
  windows = data.frame(start=c(1.70, 1.72), end=c(1.73, 1.75))
  params = get_treatment_params(physeq, c('Substrate', 'Day'), "Substrate != '12C-Con'")
  params_l = apply(params, 1, as.list)
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
  expect_equal(df_l2fc$sparsity_threshold %>% unique %>% length, 2)
  expect_equal(df_l2fc$Substrate %>% unique %>% length, 2)
  expect_equal(df_l2fc$Day %>% unique %>% length, 1)
})
