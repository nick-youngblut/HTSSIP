test_that('Beta diversity from a list of phyloseq objects',{
  physeq_l_d = physeq_list_betaDiv(physeq_l)
  expect_is(physeq_l_d, 'list')
  expect_equal(length(physeq_l_d), 2)
})

test_that('Plots created from phyloseq object',{
  # params for subseting
  params = get_treatment_params(physeq, c('Substrate', 'Day'))
  params = dplyr::filter(params, Substrate!='12C-Con')
  expect_equal(nrow(params), 2)

  # subsetting phyloseq
  ex = "(Substrate=='12C-Con' & Day=='${Day}') | (Substrate=='${Substrate}' & Day == '${Day}')"
  physeq_l = phyloseq_subset(physeq, params, ex)
  expect_is(physeq_l, 'list')
  expect_equal(length(physeq_l), 2)

  # calculating beta diversity
  physeq_l_p = SIP_betaDiv_ord(physeq_l)
  expect_is(physeq_l_p, 'list')
  expect_equal(length(physeq_l_p), 2)
  expect_is(physeq_l_p[[1]], 'ggplot')
  expect_is(physeq_l_p[[2]], 'ggplot')
})
