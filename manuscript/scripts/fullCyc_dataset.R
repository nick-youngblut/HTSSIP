# ---Description---
# Generating a dataset for HTS-SIP manuscript example
# Using data from Full C-Cycle SIP experiment
# Also simulating a q-SIP compatible dataset

# Variables
baseDir = '/var/seq_data/fullCyc/MiSeq_16SrRNA/515f-806r/lib1-7rs/'
biomFile = file.path(baseDir, 'OTU_binning/otu_table_wtax.biom')
metadataFile = file.path(baseDir, 'metadata_SIP.txt')
treeFile = file.path(baseDir, 'fasttree/otusn_pick.tree')
outDir = '~/dev/HTSSIP/manuscript/data/'
rarefy_depth = 500

# Init
library(ggplot2)
library(dplyr)
library(tidyr)
library(phyloseq)
#library(HTSSIP)
devtools::reload()

# Loading data
## biom file
physeq_all = import_biom(biomFile, treeFile)

## loading sample metadata
sample.data = import_qiime_sample_data(metadataFile)
physeq_all = merge_phyloseq(physeq_all,sample.data)
physeq_all_m = physeq_all %>% sample_data

#--3 time points--#
## subset
exp_types = c('SIP')
sub_types = c('12C-Con', '13C-Cel', '13C-Xyl')
days = c(3, 6, 14)
physeq_S2D3 = prune_samples(physeq_all_m$Exp_type %in% exp_types &
                              physeq_all_m$core_dataset == TRUE &
                              physeq_all_m$Substrate %in% sub_types &
                              physeq_all_m$Day %in% days &
                              physeq_all_m$Sample_type == 'unknown',
                            physeq_all) %>%
  filter_taxa(function(x) sum(x) > 0, TRUE)
## rarefaction
rarefy_even_depth(physeq_S2D3, sample.size=rarefy_depth)
## saving dataset
outF = file.path(outDir, 'physeq_S2D3')
saveRDS(physeq_S2D3, outF)

# list of subsetted phyloseq objects
params = get_treatment_params(physeq_S2D3, c('Substrate', 'Day'))
params = dplyr::filter(params, Substrate!='12C-Con')
ex = "(Substrate=='12C-Con' & Day=='${Day}') | (Substrate=='${Substrate}' & Day == '${Day}')"
physeq_S2D3_l = phyloseq_subset(physeq_S2D3, params, ex)
outF = file.path(outDir, 'physeq_S2D3_l')
saveRDS(physeq_S2D3_l, outF)


#--2 time points--#
## subset
days = c(3, 14)
physeq_S2D2 = prune_samples(physeq_all_m$Exp_type %in% exp_types &
                              physeq_all_m$core_dataset == TRUE &
                              physeq_all_m$Substrate %in% sub_types &
                              physeq_all_m$Day %in% days &
                              physeq_all_m$Sample_type == 'unknown',
                            physeq_all) %>%
  filter_taxa(function(x) sum(x) > 0, TRUE)
## rarefaction
rarefy_even_depth(physeq_S2D2, sample.size=rarefy_depth)
## saving dataset
outF = file.path(outDir, 'physeq_S2D2')
saveRDS(physeq_S2D2, outF)

# list of subsetted phyloseq objects
params = get_treatment_params(physeq_S2D2, c('Substrate', 'Day'))
params = dplyr::filter(params, Substrate!='12C-Con')
ex = "(Substrate=='12C-Con' & Day=='${Day}') | (Substrate=='${Substrate}' & Day == '${Day}')"
physeq_S2D2_l = phyloseq_subset(physeq_S2D2, params, ex)
outF = file.path(outDir, 'physeq_S2D2_l')
saveRDS(physeq_S2D2_l, outF)



#-----------------------------#
# simulating a q-SIP dataset
#library(HTSSIP)
## simulation params
set.seed(2)
M = 100                                  # number of species
ming = 1.67                             # gradient minimum...
maxg = 1.78                                # ...and maximum
nfrac = 24                                 # number of gradient fractions
locs = seq(ming, maxg, length=nfrac)       # gradient locations
tol  = rep(0.005, M)                       # species tolerances
h    = ceiling(rlnorm(M, meanlog=11))    # max abundances

opt = rnorm(M, mean=1.7, sd=0.005)      # species optima
params = cbind(opt=opt, tol=tol, h=h)  # put in a matrix

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
  '13C-Xyl_rep1' = params4,
  '13C-Xyl_rep2' = params5,
  '13C-Xyl_rep3' = params6
)

meta = data.frame(
  'Gradient' = c('12C-Con_rep1', '12C-Con_rep2', '12C-Con_rep3',
                 '13C-Xyl_rep1', '13C-Xyl_rep2', '13C-Xyl_rep3'),
  'Treatment' = c(rep('12C-Con', 3), rep('13C-Xyl', 3)),
  'Replicate' = c(1:3, 1:3)
)

## physeq object simulation
physeq_rep3 = HTSSIP_sim(locs, param_l, meta=meta)
outF = file.path(outDir, 'physeq_rep3')
saveRDS(physeq_rep3, outF)
## qPCR data
control_mean_fun = function(x) dnorm(x, mean=1.70, sd=0.01) * 1e8
control_sd_fun = function(x) control_mean_fun(x) / 3
treat_mean_fun = function(x) dnorm(x, mean=1.75, sd=0.01) * 1e8
treat_sd_fun = function(x) treat_mean_fun(x) / 3
physeq_rep3_qPCR = qPCR_sim(physeq_rep3,
                            control_expr='Treatment=="12C-Con"',
                            control_mean_fun=control_mean_fun,
                            control_sd_fun=control_sd_fun,
                            treat_mean_fun=treat_mean_fun,
                            treat_sd_fun=treat_sd_fun)
outF = file.path(outDir, 'physeq_rep3_qPCR')
saveRDS(physeq_rep3_qPCR, outF)
