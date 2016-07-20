test_that('HRSIP runs with default',{
  ## basic call
  windows = data.frame(start=c(1.70, 1.72), end=c(1.73, 1.75))
  df_l2fc = HRSIP(physeq,
                  sparsity_threshold = c(0, 0.1),
                  density_windows=windows,
                  design=~Substrate)
  #print(nrow(df_l2fc))
  #head(df_l2fc)
  #df_l2fc %>% .$rej_hypo %>% table %>% print
  #df_l2fc %>% .$n_rej_hypo %>% unique
  #df_l2fc %>% .$sparsity_threshold %>% table
})

test_that('HRSIP runs on multiple dataset subsets with default params',{
  windows = data.frame(start=c(1.70, 1.72), end=c(1.73, 1.75))
  params = get_treatment_params(physeq, c('Substrate', 'Day'), "Substrate != '12C-Con'")
  ex = "(Substrate=='12C-Con' & Day=='${Day}') | (Substrate=='${Substrate}' & Day == '${Day}')"
  df_l2fc = HRSIP(physeq,
                  pairwise_expr=ex,
                  pairwise_params=params,
                  sparsity_threshold=c(0, 0.1),
                  density_windows=windows,
                  design=~Substrate)
  expect_is(df_l2fc, 'data.frame')
  # df_l2fc %>% .$xxGUIDxx %>% table
  # df_l2fc %>% .$rej_hypo %>% table
  # df_l2fc %>% .$n_rej_hypo %>% table
  # df_l2fc %>% group_by(Substrate) %>% summarize(sum(rej_hypo))
})


test_that('HRSIP runs pairwise with multiple cores',{
  ## pairwise HRSIP w/ multiple cores
  ncores = 10
  doParallel::registerDoParallel(ncores)
  params = get_treatment_params(physeq, c('Substrate', 'Day'), "Substrate != '12C-Con'")
  ex = "(Substrate=='12C-Con' & Day=='${Day}') | (Substrate=='${Substrate}' & Day == '${Day}')"
  windows = data.frame(start=c(1.70, 1.72), end=c(1.73, 1.75))

  df_l2fc = HRSIP(physeq,
                pairwise_expr=ex,
                pairwise_params=params,
                sparsity_threshold=seq(0, 0.5, 0.05),
                density_windows=windows,
                design=~Substrate,
                parallel=TRUE)

  expect_is(df_l2fc, 'data.frame')
  expect_true(is.null(df_l2fc$xxGUIDxx))
  # df_l2fc %>% .$rej_hypo %>% table
  # df_l2fc %>% .$n_rej_hypo %>% table
  # df_l2fc %>% group_by(Substrate) %>% summarize(sum(rej_hypo))
})



