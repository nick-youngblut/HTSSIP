test_that('Plots created from phyloseq list',{
  physeq_l_d = physeq_list_betaDiv(physeq_l)
  params = get_treatment_params(physeq, c('Substrate', 'Day'))
  params = dplyr::filter(params, Substrate!='12C-Con')
  ex = "(Substrate=='12C-Con' & Day=='${Day}') | (Substrate=='${Substrate}' & Day == '${Day}')"
  physeq_l = phyloseq_subset(physeq, params, ex)
  physeq_l_p = SIP_betaDiv(physeq_l)

  expect_is(physeq_l, 'list')
  expect_equal(length(physeq_l), 2)
  expect_is(physeq_l_p, 'list')
  expect_equal(length(physeq_l_p), 2)
  expect_is(physeq_l_p[[1]], 'ggplot')
})
