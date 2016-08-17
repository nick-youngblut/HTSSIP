# setting parameters for tests
set.seed(2)
M = 10                                  # number of species
ming = 1.67                             # gradient minimum...
maxg = 1.78                                # ...and maximum
nfrac = 24                                 # number of gradient fractions
locs = seq(ming, maxg, length=nfrac)       # gradient locations
tol  = rep(0.005, M)                       # species tolerances
h    = ceiling(rlnorm(M, meanlog=11))    # max abundances

opt1 = rnorm(M, mean=1.7, sd=0.005)      # species optima
params1 = cbind(opt=opt1, tol=tol, h=h)  # put in a matrix
opt2 = rnorm(M, mean=1.7, sd=0.005)      # species optima
params2 = cbind(opt=opt2, tol=tol, h=h)  # put in a matrix
opt3 = rnorm(M, mean=1.7, sd=0.005)      # species optima
params3 = cbind(opt=opt3, tol=tol, h=h)  # put in a matrix
opt4 = rnorm(M, mean=1.72, sd=0.008)      # species optima
params4 = cbind(opt=opt4, tol=tol, h=h)  # put in a matrix
opt5 = rnorm(M, mean=1.72, sd=0.008)      # species optima
params5 = cbind(opt=opt5, tol=tol, h=h)  # put in a matrix
opt6 = rnorm(M, mean=1.72, sd=0.008)      # species optima
params6 = cbind(opt=opt6, tol=tol, h=h)  # put in a matrix


param_l = list(
  '12C-Con_rep1' = params1,
  '12C-Con_rep2' = params2,
  '12C-Con_rep3' = params3,
  '13C-Glu_rep1' = params4,
  '13C-Glu_rep2' = params5,
  '13C-Glu_rep3' = params6
)


test_that('Gradient sim', {
  df_OTU = gradient_sim(locs, params)
  #df_OTU %>% head
  expect_is(df_OTU, 'data.frame')
})

test_that('phyloseq sim',{
  physeq = phyloseq_sim(locs, param_l)
  #physeq
  expect_is(physeq, 'phyloseq')
})


