# dataset
physeq = readRDS('/home/nick/dev/HTSSIP/data-raw/fullCyc_con-cel-glu.RData')
physeq

# phyloseq2df
test_that('phyloseq data object can be converted to dataframe',{
  cat('\n.phyloseq2df()\n')
  df1 = phyloseq2df(physeq, sample_data)
  df2 = suppressWarnings(as.data.frame(as.matrix(sample_data(physeq))))
  expect_is(df1, 'data.frame')
  expect_equal(class(df1), class(df2))
  expect_equal(nrow(df1), nrow(df2))
  expect_equal(ncol(df1), ncol(df2))
  })

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

test_that('DESeq2_l2fc_MS working with basic parameters',{
  df_l2fc = DESeq2_l2fc_MS(physeq,
                        density_min=1.71,
                        density_max=1.75,
                        padj_cutoff=0.1,
                        design=~Substrate)

  expect_is(df_l2fc, 'data.frame')
  expect_false(is.null(df_l2fc$rej_hypo))
  expect_false(is.null(df_l2fc$n_rej_hypo))
  expect_equal(length(unique(df_l2fc$sparsity_threshold)), 1)
})


test_that('DESeq2_l2fc_MS_PW working with basic parameters',{
  params = get_treatment_params(physeq, c('Substrate', 'Day'), "Substrate != '12C-Con'")
  params_l = apply(params, 1, as.list)
  ## HRSIP call
  ex = "(Substrate=='12C-Con' & Day=='${Day}') | (Substrate=='${Substrate}' & Day == '${Day}')"
  df_l2fc = plyr::ldply(params_l, DESeq2_l2fc_MS_PW,
                        ex=ex,
                        physeq=physeq,
                        density_min=1.71,
                        density_max=1.75,
                        design=~Substrate,
                        l2fc_threshold=0.25,
                        sparsity_apply='all')

  expect_is(df_l2fc, 'data.frame')
  expect_false(is.null(df_l2fc$Substrate))
  expect_false(is.null(df_l2fc$Day))
  expect_equal(df_l2fc$Substrate %>% unique %>% length, 2)
  expect_equal(df_l2fc$Day %>% unique %>% length, 1)
})


test_that('HRSIP runs with default',{
  ## basic call
  df_l2fc = HRSIP(physeq,
                  density_min=1.71,
                  density_max=1.75,
                  padj_cutoff=0.2,
                  design=~Substrate)
  print(nrow(df_l2fc))
  head(df_l2fc)
  df_l2fc %>% .$rej_hypo %>% table %>% print
  df_l2fc %>% .$n_rej_hypo %>% unique
  df_l2fc %>% .$sparsity_threshold %>% table
})

test_that('HRSIP runs pairwise with default',{
  cat('\nHRSIP(ex, params)\n')
  params = get_treatment_params(physeq, c('Substrate', 'Day'), "Substrate != '12C-Con'")
  ex = "(Substrate=='12C-Con' & Day=='${Day}') | (Substrate=='${Substrate}' & Day == '${Day}')"
  df_l2fc = HRSIP(physeq,
                  pairwise_expr=ex,
                  pairwise_params=params,
                  sparsity_threshold=seq(0, 0.5, 0.1),
                  density_min=1.71,
                  density_max=1.75,
                  padj_cutoff=0.1,
                  design=~Substrate)
  expect_is(df_l2fc, 'data.frame')
  # df_l2fc %>% .$rej_hypo %>% table
  # df_l2fc %>% .$n_rej_hypo %>% table
  # df_l2fc %>% group_by(Substrate) %>% summarize(sum(rej_hypo))
})


test_that('HRSIP runs pairwise with multiple cores',{
  ## pairwise HRSIP w/ multiple cores
  ncores = 20
  doParallel::registerDoParallel(ncores)
  params = get_treatment_params(physeq, c('Substrate', 'Day'), "Substrate != '12C-Con'")
  ex = "(Substrate=='12C-Con' & Day=='${Day}') | (Substrate=='${Substrate}' & Day == '${Day}')"
  df_l2fc = HRSIP(physeq,
                pairwise_expr=ex,
                pairwise_params=params,
                sparsity_threshold=seq(0, 0.5, 0.05),
                density_min=1.71,
                density_max=1.75,
                design=~Substrate,
                .parallel=TRUE)
  expect_is(df_l2fc, 'data.frame')
  # df_l2fc %>% .$rej_hypo %>% table
  # df_l2fc %>% .$n_rej_hypo %>% table
  # df_l2fc %>% group_by(Substrate) %>% summarize(sum(rej_hypo))
})



