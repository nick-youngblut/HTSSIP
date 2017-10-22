#--- Notes ---#
#  Used SIPSim to simulate phyloseq object for 10 genomes
#


#--- init ---#
library(dplyr)
library(tidyr)
library(ggplot2)
devtools::load_all('.')

#--- setting variables ---#
# input files
## phyloseq
#physeq_file = '/home/nick/notebook/SIPSim/dev/bac_genome10/HTSSIP/OTU_n3_abs1e9_PCR_subNorm.physeq'
physeq_file = '/home/nick/notebook/SIPSim/dev/bac_genome10/HTSSIP/OTU_n3_abs1e9.physeq'
## qPCR
qPCR_file = '/home/nick/notebook/SIPSim/dev/bac_genome10/HTSSIP/OTU_n3_abs1e9_PCR_subNorm_qSIP.txt'
## true BD shifts from isotope incorp
BD_shift_file = '/home/nick/notebook/SIPSim/dev/bac_genome10/HTSSIP/genome10_ampFrag_BD-shift.txt'
# output
outDir = '~/dev/HTSSIP/manuscript/figures/'
## adjusted P-value cutoff
padj_cutoff = 0.1
## number of cores for parallel processing (increase depending on your computational hardware)
ncores = 10
## number of bootstrap replicates
nboot = 100

#--- loading ---#
physeq = readRDS(physeq_file)
x = sample_data(physeq)
colnames(x)[c(1,2)] = c('Library', 'Fraction')
physeq = phyloseq(otu_table(physeq), x)
qPCR = read.delim(qPCR_file, sep='\t')
colnames(qPCR)[c(1,2)] = c('Library', 'Fraction')
qPCR %>% head(n=3)
BD_shift = read.delim(BD_shift_file, sep='\t')
colnames(BD_shift)[1] = 'Library'
BD_shift %>% head(n=3)
## formatting physeq
physeq_m = physeq %>% sample_data() %>% as.data.frame
physeq_m$Treatment = ifelse(physeq_m$Library %% 2 == 0, '13C-Treat', '12C-Con')
physeq_m$Buoyant_density = physeq_m$BD_mid
physeq = phyloseq(physeq %>% otu_table,
                  physeq_m %>% sample_data)

#--- plotting distributions ---#
# 'raw' counts
physeq_otu = phyloseq2table(physeq, include_sample_data=TRUE) %>%
  dplyr::filter(Library %in% c(1,2)) %>%
  dplyr::mutate(Buoyant_density = Buoyant_density %>% as.character %>% as.numeric,
                Count = Count %>% as.numeric) %>%
  dplyr::filter(!is.infinite(abs(Buoyant_density)))
## plotting
p_count = ggplot(physeq_otu, aes(Buoyant_density, Count, fill=Treatment)) +
  geom_area(alpha=0.3, position=position_dodge(width=0.001)) +
  facet_wrap(~ OTU) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)
  )
p_count

# relative abundances
physeq_otu = physeq_otu %>%
   dplyr::group_by(SAMPLE_JOIN, Buoyant_density, Treatment) %>%
   dplyr::mutate(Total_count = sum(as.numeric(Count), na.rm=TRUE)) %>%
   dplyr::ungroup() %>%
   dplyr::mutate(Rel_abund = Count / Total_count)

## plotting
p_norm = ggplot(physeq_otu, aes(Buoyant_density, Rel_abund, fill=Treatment)) +
  geom_area(alpha=0.3, position=position_dodge(width=0.001)) +
  facet_wrap(~ OTU) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)
  )
p_norm


# qPCR-transformed
qPCR_tmp = qPCR %>%
  dplyr::mutate(Library = Library %>% as.character %>% as.factor) %>%
  dplyr::select(Library, BD_mid, Fraction, taxon, prop_abs_abund)
physeq_otu_j = physeq_otu %>%
  dplyr::inner_join(qPCR_tmp,
                    c('Library'='Library',
                      'Buoyant_density'='BD_mid',
                      'Fraction'='Fraction',
                      'OTU' = 'taxon'))

stopifnot(nrow(physeq_otu) >= nrow(physeq_otu_j))

## plotting
p_qPCR = ggplot(physeq_otu_j, aes(Buoyant_density, prop_abs_abund, fill=Treatment)) +
  geom_area(alpha=0.3, position=position_dodge(width=0.001)) +
  facet_wrap(~ OTU) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)
  )
p_qPCR


# combined distribution
## renaming
taxon_names = c(
  'Bacillus_cereus_F837_76',
  'Cyanobium_gracile_PCC_6307',
  'Escherichia_coli_APEC_O78',
  'Geobacter_lovleyi_SZ',
  'Idiomarina_loihiensis_GSL_199',
  'Leuconostoc_citreum_KM20',
  'Micrococcus_luteus_NCTC_2665',
  'Nitrosomonas_europaea_ATCC_19718',
  'Roseobacter_litoralis_Och_149',
  'Xanthomonas_axonopodis_Xac29-1'
)
renames = data.frame(taxon_names = taxon_names,
                     otu_names = sapply(1:length(taxon_names),
                                        function(x) paste0('OTU.', x)))

## adding qPCR values
physeq_otu_j_g = physeq_otu_j %>%
  tidyr::gather(Abund_type, Abundance, Count, Rel_abund, prop_abs_abund) %>%
  dplyr::mutate(Treatment = Treatment %>% as.character,
                Treatment = ifelse(Treatment == '12C', '12C-Con', Treatment),
                Treatment = ifelse(Treatment == '13C', '13C-Treat', Treatment),
                Abund_type = ifelse(Abund_type == 'Rel_abund', 'Rel. Abund.', Abund_type),
                Abund_type = ifelse(Abund_type == 'prop_abs_abund', 'Rel. Abund.\nqPCR-trans.', Abund_type)) %>%
  inner_join(renames, c('OTU'='taxon_names')) %>%
  dplyr::select(-OTU) %>%
  rename('OTU' = otu_names) %>%
  mutate(OTU = OTU %>% reorder(gsub('OTU.', '', OTU) %>% as.numeric))

## plotting
x_lab = bquote('Buoyant density (g '* ml^-1*')')
p_abund = ggplot(physeq_otu_j_g, aes(Buoyant_density, Abundance, fill=Treatment)) +
  geom_area(alpha=0.3, position=position_dodge(width=0.001)) +
  scale_x_continuous(limits=c(1.67, 1.77), breaks=c(1.68, 1.72, 1.76)) +
  scale_fill_manual('Treatment', values=c('red', 'blue')) +
  facet_grid(Abund_type ~ OTU, scales='free_y') +
  labs(x=x_lab) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle=45, hjust=1)
  )
p_abund



#------- ID incorporators & quant. BD shift ------#
# loading physeq
physeq_file = '/home/nick/notebook/SIPSim/dev/bac_genome10/HTSSIP/OTU_n3_abs1e9_PCR_subNorm.physeq'
physeq = readRDS(physeq_file)
x = sample_data(physeq)
colnames(x)[c(1,2)] = c('Library', 'Fraction')
physeq = phyloseq(otu_table(physeq), x)

## formatting physeq
physeq_m = physeq %>% sample_data() %>% as.data.frame
physeq_m$Treatment = ifelse(physeq_m$Library %% 2 == 0, '13C-Treat', '12C-Con')
physeq_m$Buoyant_density = physeq_m$BD_mid
physeq = phyloseq(physeq %>% otu_table,
                  physeq_m %>% sample_data)


#--- heavy-SIP ---#
df_heavy = heavy_SIP(physeq,
                     ex="Treatment=='12C-Con'",
                     rep='Library',
                     light_window=c(1.68, 1.70),
                     heavy_window=c(1.73, 1.75),
                     comparison='H-v-H',
                     hypo_test='wilcox')
df_heavy$OTU = rownames(df_heavy)

# number of incorporators
df_heavy %>%
  filter(padj < padj_cutoff) %>%
  group_by() %>%
  summarize(n_incorp_OTUs = OTU %>% unique %>% length)

#--- HR-SIP ---#
doParallel::registerDoParallel(ncores)
df_l2fc_HRSIP = HRSIP(physeq,
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
  nrow

#--- MW-HR-SIP ---#
# BD windows
windows = data.frame(density_min=c(1.72, 1.74),
                     density_max=c(1.75, 1.77))

# MW-HR-SIP call
doParallel::registerDoParallel(ncores)
df_l2fc_MWHRSIP = plyr::ldply(list('control_vs_treat' = physeq),
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
# transforming OTU counts
qPCR_sum = qPCR %>%
  dplyr::distinct(Library, Fraction, BD_mid, total_qPCR_copies) %>%
  tidyr::unite('Sample', Library, Fraction, sep='__', remove=FALSE) %>%
  dplyr::rename('qPCR_tech_rep_mean' = total_qPCR_copies)
physeq_t = OTU_qPCR_trans(physeq, qPCR_sum)
# calculating atom Fraction excess
atomX = qSIP_atom_excess(physeq_t,
                         control_expr='Treatment=="12C-Con"',
                         treatment_rep='Library')
# calc CI
doParallel::registerDoParallel(ncores)
df_atomX_boot = qSIP_bootstrap(atomX, n_boot=100, parallel=TRUE)
df_atomX_boot %>% head(n=4) %>% as.data.frame

#--- delta_BD ---#
df_dBD = delta_BD(physeq, control_expr='Treatment=="12C-Con"')
df_dBD %>% head(n=4)

#--- summarizing incorporators ---#
hSIP_res = df_heavy %>%
  dplyr::mutate(OTU = rownames(.),
                heavySIP_statistic = statistic,
                heavySIP_padj = padj,
                heavySIP_incorp = padj < padj_cutoff) %>%
  dplyr::select(OTU, heavySIP_statistic, heavySIP_padj, heavySIP_incorp)

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
                qSIP_incorp = A_CI_low > 0.05) %>%
  dplyr::select(OTU, qSIP_A, qSIP_incorp)

incorp_res = HRSIP_res %>%
  dplyr::inner_join(MWHRSIP_res, c('OTU')) %>%
  dplyr::inner_join(qSIP_res, c('OTU')) %>%
  dplyr::inner_join(hSIP_res, c('OTU')) %>%
  dplyr::select(OTU, heavySIP_incorp, HRSIP_incorp, MWHRSIP_incorp, qSIP_incorp) %>%
  inner_join(renames, c('OTU'='taxon_names')) %>%
  dplyr::select(-OTU) %>%
  rename('OTU' = otu_names) %>%
  tidyr::gather(method, incorp_status, heavySIP_incorp,
                HRSIP_incorp, MWHRSIP_incorp, qSIP_incorp) %>%
  dplyr::mutate(OTU_num = gsub('OTU.', '', OTU) %>% as.character %>% as.numeric,
                OTU = OTU %>% reorder(OTU_num),
                method = ifelse(method=='heavySIP_incorp', 'Mann-Whitney', method),
                method = ifelse(method=='qSIP_incorp', 'q-SIP', method),
                method = ifelse(method=='HRSIP_incorp', 'HR-SIP', method),
                method = ifelse(method=='MWHRSIP_incorp', 'MW-HR-SIP', method),
                method = factor(method, levels=c('q-SIP', 'MW-HR-SIP',
                                                 'HR-SIP', 'Mann-Whitney')))

p_incorp = ggplot(incorp_res, aes(OTU, method, fill=incorp_status)) +
  geom_tile() +
  scale_fill_manual('Incorporator?', values=c('grey', 'orange')) +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )
p_incorp


#--- plot of BD shift ---#
## deltaBD + qSIP
df_j = dplyr::inner_join(df_atomX_boot, df_dBD, c('OTU'='OTU'))
df_j = df_j %>%
  dplyr::mutate(OTU = reorder(OTU, -delta_BD)) %>%
  dplyr::select(OTU, Z, delta_BD) %>%
  tidyr::gather(Method, value, -OTU) %>%
  inner_join(renames, c('OTU'='taxon_names')) %>%
  dplyr::select(-OTU) %>%
  rename('OTU' = otu_names) %>%
  dplyr::mutate(Method = ifelse(Method=='Z', 'q-SIP', 'delta BD'),
                OTU_num = gsub('OTU.', '', OTU) %>% as.character %>% as.numeric,
                OTU = OTU %>% reorder(OTU_num))

## adding true BD shift
BD_shift_s = BD_shift %>%
  dplyr::mutate(treatment = ifelse(Library %% 2 == 0, 'X13C', 'X12C')) %>%
  group_by(taxon, treatment) %>%
  summarize(mean_BD_shift = mean(mean)) %>%
  ungroup() %>%
  spread(treatment, mean_BD_shift) %>%
  dplyr::mutate(mean_BD_shift = X13C - X12C,
                mean_BD_shift_txt = round(mean_BD_shift, 3) %>% as.character,
                incorporator = X13C > 0.001) %>%
  dplyr::rename('OTU' = taxon) %>%
  inner_join(renames, c('OTU'='taxon_names')) %>%
  dplyr::select(-OTU) %>%
  rename('OTU' = otu_names)
BD_shift_s

df_j_j = df_j %>%
  inner_join(BD_shift_s, c('OTU')) %>%
  mutate(OTU = OTU %>% reorder(gsub('OTU.', '', OTU) %>% as.numeric))

## plotting
p_shift = ggplot(df_j_j, aes(OTU, value, color=Method)) +
  geom_hline(yintercept=0.00, linetype='dashed', alpha=0.5) +
  geom_point(alpha=0.3, size=3) +
  geom_point(alpha=0.8, size=3, shape='O') +
  geom_point(aes(y=mean_BD_shift), size=10, shape='-', color='red') +
  scale_color_manual(values=c('darkgreen', 'purple')) +
  labs(x='OTU', y='BD shift (Z)') +
  theme_bw() +
  theme(
    axis.title.x = element_blank()
  )
p_shift

## saving
# #outF = file.path(outDir, 'deltaBD_qSIP_Z.pdf')
# #ggsave(outF, p_shift, width=6, height=3)
# #outF = file.path(outDir, 'deltaBD_qSIP_Z.png')
# #ggsave(outF, p_shift, width=6, height=3)

#--- combining plots ---#
p_comb = cowplot::ggdraw() +
  cowplot::draw_plot(p_abund,  0.015, 0.35, 0.98, 0.65) +
  cowplot::draw_plot(p_incorp, 0, 0.20, 0.965, 0.15) +
  cowplot::draw_plot(p_shift,  0.035, 0.00, 0.91, 0.20) +
  cowplot::draw_plot_label(c('A)', 'B)', 'C)'),
                           c(0, 0, 0),
                           c(1, 0.37, 0.20), size = 11)
p_comb

## saving
outF = file.path(outDir, 'sim-data_abund_incorp_shift.pdf')
ggsave(outF, p_comb, width=9, height=7)
cat('File written:', outF, '\n')
outF = file.path(outDir, 'sim-data_abund_incorp_shift.png')
ggsave(outF, p_comb, width=9, height=7)
cat('File written:', outF, '\n')

