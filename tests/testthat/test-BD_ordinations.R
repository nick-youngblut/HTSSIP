test_that('Beta diversity from a list of phyloseq objects',{
  physeq_l_d = physeq_list_betaDiv(physeq_S2D2_l)
  expect_is(physeq_l_d, 'list')
  expect_equal(length(physeq_l_d), 4)
})

test_that('Beta diversity from a list of phyloseq objects (parallel)',{
  doParallel::registerDoParallel(4)
  physeq_l_d = physeq_list_betaDiv(physeq_S2D2_l, parallel=TRUE)
  expect_is(physeq_l_d, 'list')
  expect_equal(length(physeq_l_d), 4)
})


test_that('Can make a data.frame to plot',{
  # params for subseting
  doParallel::registerDoParallel(4)
  physeq_l_d = physeq_list_betaDiv(physeq_S2D2_l, parallel=TRUE)
  physeq_l_d_ord = physeq_list_ord(physeq_S2D2_l, physeq_l_d)
  physeq_l_d_ord_df = phyloseq_list_ord_dfs(physeq_S2D2_l, physeq_l_d_ord)
  expect_true('phyloseq_subsets' %in% colnames(physeq_l_d_ord_df))
})

test_that('Plots created from phyloseq object',{
  # params for subseting
  params = get_treatment_params(physeq_S2D2, c('Substrate', 'Day'))
  expect_is(params, 'data.frame')
  expect_equal(nrow(params), 6)
  ## filtering params
  params = dplyr::filter(params, Substrate!='12C-Con')
  expect_equal(nrow(params), 4)

  # subsetting phyloseq
  ex = "(Substrate=='12C-Con' & Day=='${Day}') | (Substrate=='${Substrate}' & Day == '${Day}')"
  physeq_S2D2_l = phyloseq_subset(physeq_S2D2, params, ex)
  expect_is(physeq_S2D2_l, 'list')
  expect_equal(length(physeq_S2D2_l), 4)

  # calculating beta diversity
  physeq_S2D2_p = SIP_betaDiv_ord(physeq_S2D2_l)
  expect_is(physeq_S2D2_p, 'data.frame')

  ## plot output
  physeq_S2D2_p = SIP_betaDiv_ord(physeq_S2D2_l, plot=TRUE)
  expect_is(physeq_S2D2_p, 'ggplot')

  #expect_equal(length(physeq_S2D2_l_p), 4)
  #lapply(physeq_S2D2_l_p, expect_is, class='ggplot')
  #expect_is(physeq_S2D2_l_p[[1]], 'ggplot')
  #expect_is(physeq_S2D2_l_p[[2]], 'ggplot')
})
