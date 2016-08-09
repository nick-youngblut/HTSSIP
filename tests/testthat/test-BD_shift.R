test_that('Percent overlap is working',{
  x = perc_overlap(0, 1, 0, 0.5)
  expect_equal(x, 50)
  x = perc_overlap(0, 0.5, 0, 1)
  expect_equal(x, 100)
  x = perc_overlap(0, 0.5, 0.5, 1)
  expect_equal(x, 0)
})


expect_wmean = function(wmean){
  expect_is(wmean, 'data.frame')
  expect_gt(wmean %>% nrow, 0)

  wmean_min = wmean$distance %>% min
  expect_gte(wmean_min, 0)
  expect_lt(wmean_min, 1)
  wmean_max = wmean$distance %>% max
  expect_gt(wmean_max, 0)
  expect_lte(wmean_max, 1)

  cat('\n\n--Weighted mean distance summary--\n')
  wmean$distance %>% summary %>% print
  cat('\n\n')
}

test_that('BD_shift runs w/ default',{
  ## basic call
  data(physeq_l)

  # dataset 1
  wmean = BD_shift(physeq_l[[1]])
  expect_wmean(wmean)

  # dataset 2
  wmean = BD_shift(physeq_l[[2]])
  expect_wmean(wmean)

  # ggplot
  p = ggplot(wmean, aes(BD_min.x, wmean_dist)) +
        geom_point()
  expect_is(p, 'ggplot')
})

test_that('BD_shift runs w/ Bray-Curtis',{
  ## basic call
  data(physeq_l)

  # dataset 1
  wmean = BD_shift(physeq_l[[1]], method='bray')
  expect_wmean(wmean)

  # ggplot
  p = ggplot(wmean, aes(BD_min.x, wmean_dist)) +
    geom_point()
  expect_is(p, 'ggplot')
})
