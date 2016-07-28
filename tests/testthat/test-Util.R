#load('/home/nick/dev/HTSSIP/data/physeq.RData')

test_that('phyloseq sample_data can be converted to dataframe',{
  df1 = phyloseq2df(physeq, phyloseq::sample_data)
  df2 = suppressWarnings(as.data.frame(as.matrix(sample_data(physeq))))

  expect_is(df1, 'data.frame')
  expect_equal(class(df1), class(df2))
  expect_equal(nrow(df1), nrow(df2))
  expect_equal(ncol(df1), ncol(df2))
})


test_that('phyloseq tax_table can be converted to dataframe',{
  df1 = phyloseq2df(physeq, phyloseq::tax_table)
  df2 = suppressWarnings(as.data.frame(as.matrix(tax_table(physeq))))

  expect_is(df1, 'data.frame')
  expect_equal(class(df1), class(df2))
  expect_equal(nrow(df1), nrow(df2))
  expect_equal(ncol(df1), ncol(df2))
})


test_that('phyloseq otu_table can be converted to dataframe',{
  df1 = phyloseq2df(physeq, phyloseq::otu_table)
  df2 = suppressWarnings(as.data.frame(as.matrix(otu_table(physeq))))

  expect_is(df1, 'data.frame')
  expect_equal(class(df1), class(df2))
  expect_equal(nrow(df1), nrow(df2))
  expect_equal(ncol(df1), ncol(df2))
})
