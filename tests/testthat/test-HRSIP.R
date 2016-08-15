test_that('HRSIP runs with default',{
  ## basic call
  df_l2fc = HRSIP(physeq_S2D2_l[[1]],
                  design=~Substrate)
  expect_is(df_l2fc, 'data.frame')
  expect_gt(nrow(df_l2fc), 0)

  df_l2fc_s = df_l2fc %>%
    dplyr::group_by(sparsity_threshold,
                    density_min, density_max) %>%
    dplyr::summarize(n = n())
  expect_gt(nrow(df_l2f_s), 1)
})

test_that('MW-HRSIP runs with default & parallel ',{
  ncores = 8
  doParallel::registerDoParallel(ncores)
  windows = data.frame(density_min=c(1.70, 1.72), density_max=c(1.73, 1.75))
  df_l2fc = HRSIP(physeq_S2D2_l[[1]],
                  design=~Substrate,
                  density_windows=windows,
                  parallel=TRUE)
  expect_is(df_l2fc, 'data.frame')
  expect_gt(nrow(df_l2fc), 0)

  df_l2fc_s = df_l2fc %>%
    dplyr::group_by(sparsity_threshold,
                    density_min, density_max) %>%
    dplyr::summarize(n = n())
  expect_gt(nrow(df_l2f_s), 1)
})

