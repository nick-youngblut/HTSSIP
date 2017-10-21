#--- Notes ---#
# Simulating data with HTSSIP functions and conducting analyses
# using the simulated data


#--- library loading ---#
library(dplyr)
library(ggplot2)
devtools::load_all('.')

#--- setting variables ---#
# output
outDir = '~/dev/HTSSIP/manuscript/figures/'
## adjusted P-value cutoff
padj_cutoff = 0.1
## number of cores for parallel processing (increase depending on your computational hardware)
ncores = 10
## number of bootstrap replicates
nboot = 100

#--- simulations ---#
set.seed(2)
M = 10                                     # number of species
ming = 1.67                                # gradient minimum...
maxg = 1.78                                # ...and maximum
nfrac = 24                                 # number of gradient fractions
locs = seq(ming, maxg, length=nfrac)       # gradient locations
tol  = rep(0.013, M)                       # species tolerances
h    = ceiling(rlnorm(M, meanlog=19, sdlog=1.1))      # max abundances

opt = rnorm(M, mean=1.7, sd=0.005)        # species optima
params = cbind(opt=opt, tol=tol, h=h)     # put in a matrix

## control
opt1 = rnorm(M, mean=1.7, sd=0.006)       # species optima
params1 = cbind(opt=opt1, tol=tol, h=h)   # put in a matrix
opt2 = rnorm(M, mean=1.7, sd=0.006)       # species optima
params2 = cbind(opt=opt2, tol=tol, h=h)   # put in a matrix
opt3 = rnorm(M, mean=1.7, sd=0.006)       # species optima
params3 = cbind(opt=opt3, tol=tol, h=h)   # put in a matrix
## treatment
opt4 = rnorm(M, mean=1.72, sd=0.014)      # species optima
params4 = cbind(opt=opt4, tol=tol, h=h)   # put in a matrix
opt5 = rnorm(M, mean=1.72, sd=0.014)      # species optima
params5 = cbind(opt=opt5, tol=tol, h=h)   # put in a matrix
opt6 = rnorm(M, mean=1.72, sd=0.014)      # species optima
params6 = cbind(opt=opt6, tol=tol, h=h)   # put in a matrix

param_l = list(
  '12C-Control_rep1' = params1,
  '12C-Control_rep2' = params2,
  '12C-Control_rep3' = params3,
  '13C-Treatment_rep1' = params4,
  '13C-Treatment_rep2' = params5,
  '13C-Treatment_rep3' = params6
)

meta = data.frame(
  'Gradient' = c('12C-Control_rep1', '12C-Control_rep2', '12C-Control_rep3',
                 '13C-Treatment_rep1', '13C-Treatment_rep2', '13C-Treatment_rep3'),
  'Treatment' = c(rep('12C-Control', 3), rep('13C-Treatment', 3)),
  'Replicate' = c(1:3, 1:3)
)


## physeq object simulation
physeq_rep3 = HTSSIP_sim(locs, param_l, meta=meta, sim_tree=TRUE)
## qPCR data
control_mean_fun = function(x) dnorm(x, mean=1.70, sd=0.01) * 1e8
control_sd_fun = function(x) control_mean_fun(x) / 3
treat_mean_fun = function(x) dnorm(x, mean=1.75, sd=0.01) * 1e8
treat_sd_fun = function(x) treat_mean_fun(x) / 3
physeq_rep3_qPCR = qPCR_sim(physeq_rep3,
                            control_expr='Treatment=="12C-Control"',
                            control_mean_fun=control_mean_fun,
                            control_sd_fun=control_sd_fun,
                            treat_mean_fun=treat_mean_fun,
                            treat_sd_fun=treat_sd_fun)

#physeq to list of comparisons
ex = "(Treatment=='12C-Control') | (Treatment=='${Treatment}')"
params = get_treatment_params(physeq_rep3, c('Treatment'), "Treatment != '12C-Control'")
physeq_rep3_l = phyloseq_subset(physeq_rep3, params, ex)
physeq_rep3_l %>% names


#--- plotting distributions ---#
OTUs = sapply(1:10, function(x) paste(c('OTU', x), collapse='.'))
# OTU data
physeq_rep3_otu = phyloseq2table(physeq_rep3, include_sample_data=TRUE)
physeq_rep3_otu = physeq_rep3_otu %>%
  dplyr::mutate(Buoyant_density = Buoyant_density %>% as.character %>% as.numeric)

# plotting
p_count = ggplot(physeq_rep3_otu %>% dplyr::filter(Replicate==1, OTU %in% OTUs),
                 aes(Buoyant_density, Count, fill=Treatment)) +
  geom_area(position='dodge', alpha=0.3) +
  facet_wrap(~ OTU) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle=45, hjust=1)
  )
p_count

## normalized abundances
physeq_rep3_otu = physeq_rep3_otu %>%
  dplyr::group_by(SAMPLE_JOIN, Gradient, Buoyant_density, Treatment) %>%
  dplyr::mutate(Total_count = sum(as.numeric(Count), na.rm=TRUE)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(Rel_abund = Count / Total_count)

## transforming by qPCR
physeq_rep3_otu_j = physeq_rep3_otu %>%
  dplyr::inner_join(physeq_rep3_qPCR[['summary']],
                    c('SAMPLE_JOIN'='Sample',
                      'Buoyant_density'='Buoyant_density',
                      'Treatment'='Treatment',
                      'Replicate'='Replicate')) %>%
  dplyr::mutate(Rel_abund_qPCR = Rel_abund * qPCR_tech_rep_mean)

# all together
physeq_rep3_otu_g = physeq_rep3_otu_j %>%
  tidyr::gather(Abund_type, Abundance, Count, Rel_abund, Rel_abund_qPCR) %>%
  dplyr::mutate(OTU_num = gsub('OTU.', '', OTU) %>% as.character %>% as.numeric,
                OTU = OTU %>% reorder(OTU_num),
                Treatment = Treatment %>% as.character,
                Treatment = ifelse(Treatment == '12C-Control', '12C-Con', Treatment),
                Treatment = ifelse(Treatment == '13C-Treatment', '13C-Treat', Treatment),
                Abund_type = ifelse(Abund_type == 'Rel_abund', 'Rel. Abund.', Abund_type),
                Abund_type = ifelse(Abund_type == 'Rel_abund_qPCR', 'Rel. Abund.\nqPCR-trans.', Abund_type))

## plotting abundances
x_lab = bquote('Buoyant density (g '* ml^-1*')')
p_abund = ggplot(physeq_rep3_otu_g %>% dplyr::filter(Replicate==1, OTU %in% OTUs),
                 aes(Buoyant_density, Abundance, fill=Treatment)) +
  geom_area(position='dodge', alpha=0.3) +
  scale_x_continuous(limits=c(1.67, 1.77), breaks=c(1.68, 1.72, 1.76)) +
  facet_grid(Abund_type ~ OTU, scales='free_y') +
  labs(x=x_lab) +
  theme_bw() +
  theme(
    axis.title.y = element_blank(),
    axis.text.x = element_text(angle=45, hjust=1)
  )
p_abund


#--- HR-SIP ---#
doParallel::registerDoParallel(ncores)
df_l2fc_HRSIP = HRSIP(physeq_rep3,
                design = ~Treatment,
                padj_cutoff = padj_cutoff,
                data.frame(density_min = c(1.72),
                           density_max = c(1.75)),
                sparsity_threshold = seq(0,0.5,0.05),
                parallel = TRUE)
df_l2fc_HRSIP %>% head(n=3)

# number of incorporators
df_l2fc_HRSIP %>%
  filter(padj < padj_cutoff) %>%
  group_by() %>%
  summarize(n_incorp_OTUs = OTU %>% unique %>% length)


#--- MW-HR-SIP ---#
# BD windows
windows = data.frame(density_min=c(1.70, 1.72, 1.74),
                     density_max=c(1.73, 1.75, 1.77))

# MW-HR-SIP call
doParallel::registerDoParallel(ncores)
df_l2fc_MWHRSIP = plyr::ldply(physeq_rep3_l,
                      HRSIP,
                      density_windows = windows,
                      design = ~Treatment,
                      padj_cutoff = padj_cutoff,
                      sparsity_threshold = seq(0,0.5,0.05),
                      parallel = TRUE)
df_l2fc_MWHRSIP %>% head(n=3)

# number of incorporators
df_l2fc_MWHRSIP %>%
  filter(padj < padj_cutoff) %>%
  group_by(.id, sparsity_threshold) %>%
  summarize(n_incorp_OTUs = OTU %>% unique %>% length) %>%
  as.data.frame


#--- q-SIP ---#
# summing qPCR
physeq_rep3_qPCR_sum = physeq_rep3_qPCR$summary
physeq_rep3_qPCR_sum %>% head(n=4)
# transforming OTU counts
physeq_rep3_t = OTU_qPCR_trans(physeq_rep3, physeq_rep3_qPCR_sum)
# calculating atom fraction excess
atomX = qSIP_atom_excess(physeq_rep3_t,
                         control_expr='Treatment=="12C-Control"',
                         treatment_rep='Replicate')
# calc CI
doParallel::registerDoParallel(ncores)
df_atomX_boot = qSIP_bootstrap(atomX, n_boot=100, parallel=TRUE)
df_atomX_boot %>% head(n=4)


#--- delta_BD ---#
df_dBD = delta_BD(physeq_rep3, control_expr='Treatment=="12C-Control"')
df_dBD %>% head(n=4)


#--- summarizing incorporators ---#
HRSIP_res = df_l2fc_HRSIP %>%
  dplyr::mutate(HRSIP_log2FoldChange = log2FoldChange,
                HRSIP_padj = padj,
                HRSIP_incorp = padj < padj_cutoff) %>%
  dplyr::select(OTU, HRSIP_log2FoldChange, HRSIP_padj, HRSIP_incorp)

MWHRSIP_res = df_l2fc_MWHRSIP %>%
  dplyr::mutate(MWHRSIP_log2FoldChange = log2FoldChange,
                MWHRSIP_padj = padj,
                MWHRSIP_incorp = padj < padj_cutoff) %>%
  dplyr::select(OTU, MWHRSIP_log2FoldChange, MWHRSIP_padj, MWHRSIP_incorp)

qSIP_res = df_atomX_boot %>%
  dplyr::mutate(qSIP_A = A,
                qSIP_incorp = A_CI_low > 0) %>%
  dplyr::select(OTU, qSIP_A, qSIP_incorp)

incorp_res = HRSIP_res %>%
  dplyr::inner_join(MWHRSIP_res, c('OTU')) %>%
  dplyr::inner_join(qSIP_res, c('OTU')) %>%
  dplyr::select(OTU,HRSIP_incorp, MWHRSIP_incorp, qSIP_incorp) %>%
  tidyr::gather(method, incorp_status, HRSIP_incorp, MWHRSIP_incorp, qSIP_incorp) %>%
  dplyr::mutate(OTU_num = gsub('OTU.', '', OTU) %>% as.character %>% as.numeric,
                OTU = OTU %>% reorder(OTU_num),
                method = ifelse(method=='qSIP_incorp', 'q-SIP', method),
                method = ifelse(method=='HRSIP_incorp', 'HR-SIP', method),
                method = ifelse(method=='MWHRSIP_incorp', 'MW-HR-SIP', method))

p_incorp = ggplot(incorp_res, aes(OTU, method, fill=incorp_status)) +
  geom_tile() +
  scale_fill_brewer('Incorporator?', palette="Dark2") +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )
p_incorp

#--- table/plots ---#
# plot of BD
## deltaBD + qSIP
df_j = dplyr::inner_join(df_atomX_boot, df_dBD, c('OTU'='OTU'))
df_j = df_j %>%
  dplyr::mutate(OTU = reorder(OTU, -delta_BD)) %>%
  dplyr::select(OTU, Z, delta_BD) %>%
  tidyr::gather(Method, value, -OTU) %>%
  dplyr::mutate(Method = ifelse(Method=='Z', 'q-SIP', 'delta BD'),
                OTU_num = gsub('OTU.', '', OTU) %>% as.character %>% as.numeric,
                OTU = OTU %>% reorder(OTU_num))

# adding 'true' BD shift


## plotting
p_shift = ggplot(df_j, aes(OTU, value, color=Method)) +
  geom_hline(yintercept=0.00, linetype='dashed', alpha=0.5) +
  geom_point(alpha=0.3, size=3) +
  geom_point(alpha=0.8, size=3, shape='O') +
  scale_color_manual(values=c('Blue', 'Orange')) +
  labs(x='OTU', y='BD shift (Z)') +
  theme_bw() +
  theme(
    axis.title.x = element_blank()
  )
p_shift

## saving
#outF = file.path(outDir, 'deltaBD_qSIP_Z.pdf')
#ggsave(outF, p_shift, width=6, height=3)
#outF = file.path(outDir, 'deltaBD_qSIP_Z.png')
#ggsave(outF, p_shift, width=6, height=3)


#--- combining plots ---#
p_comb = cowplot::ggdraw() +
  cowplot::draw_plot(p_abund,  0.025, 0.35, 0.975, 0.65) +
  cowplot::draw_plot(p_incorp, 0, 0.20, 0.96, 0.15) +
  cowplot::draw_plot(p_shift,  0.02, 0.00, 0.92, 0.20) +
  cowplot::draw_plot_label(c('A)', 'B)', 'C)'),
                           c(0, 0, 0),
                           c(1, 0.36, 0.20), size = 11)
p_comb

## saving
outF = file.path(outDir, 'sim-data_abund_incorp_shift.pdf')
ggsave(outF, p_comb, width=8, height=7)
outF = file.path(outDir, 'sim-data_abund_incorp_shift.png')
ggsave(outF, p_comb, width=8, height=7)

