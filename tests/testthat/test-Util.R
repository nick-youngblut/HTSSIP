test_that('Expression parameter extraction regardless of quotes',{
  x = "(Substrate=='12C-Con' & Day=='14') | (Substrate=='13C-Cel' & Day == '14')"
  y = '(Substrate=="12C-Con" & Day=="14") | (Substrate=="13C-Cel" & Day == "14")'
  xx = expr_param_extract(x)
  expect_true(all(c( "12C-Con", "14", "13C-Cel","14") %in% xx))
  yy = expr_param_extract(y)
  expect_true(all(c( "12C-Con", "14", "13C-Cel","14") %in% yy))
})

test_that('Expression parameter extraction regardless of expr vector length',{
# returns a matrix
  x = c('(Substrate=="12C-Con" & Day=="14")',
        '(Substrate=="13C-Cel" & Day == "14")')
  xx = expr_param_extract(x)
  expect_is(xx, 'matrix')

  # returns a list
  y = c('(Substrate=="12C-Con" & Day=="14")',
        '(Substrate=="13C-Cel" & Day == "14")',
        '(Substrate=="13C-Cel")')
  yy = expr_param_extract(y)
  expect_is(yy, 'list')
})

expect_fun = function(df1, df2){
  expect_is(df1, 'data.frame')
  expect_equal(class(df1), class(df2))
  expect_equal(nrow(df1), nrow(df2))
  expect_equal(ncol(df1), ncol(df2))
}

test_that('phyloseq sample_data can be converted to dataframe',{
  df1 = phyloseq2df(physeq_S2D2, phyloseq::sample_data)
  df2 = suppressWarnings(as.data.frame(as.matrix(sample_data(physeq_S2D2))))

  expect_fun(df1, df2)
})


test_that('phyloseq tax_table can be converted to dataframe',{
  df1 = phyloseq2df(physeq_S2D2, phyloseq::tax_table)
  df2 = suppressWarnings(as.data.frame(as.matrix(tax_table(physeq_S2D2))))

  expect_fun(df1, df2)
})


test_that('phyloseq otu_table can be converted to dataframe',{
  df1 = phyloseq2df(physeq_S2D2, phyloseq::otu_table)
  df2 = suppressWarnings(as.data.frame(as.matrix(otu_table(physeq_S2D2))))

  expect_fun(df1, df2)
})


test_that('Subsetting phyloseq object',{
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
})


test_that('Naming phyloseq object',{
  # params for subseting
  params = get_treatment_params(physeq_S2D2, c('Substrate', 'Day'))
  ## filtering params
  params = dplyr::filter(params, Substrate!='12C-Con')

  # subsetting phyloseq
  ex = "(Substrate=='12C-Con' & Day=='${Day}') | (Substrate=='${Substrate}' & Day == '${Day}')"
  physeq_S2D2_l = phyloseq_subset(physeq_S2D2, params, ex)
})

