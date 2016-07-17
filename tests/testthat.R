library(testthat)
library(HTSSIP)

test_check('HTSSIP')

physeq = readRDS('/home/nick/dev/HTSSIP/data-raw/fullCyc_con-cel-xyl.rds')


.HRSIP(physeq, sparsity_threshold=0.25,
  density_min=1.71, density_max=1.75,
  design=~Substrate,
  l2fc_threshold=0.25,
  sparsity_apply='all')
