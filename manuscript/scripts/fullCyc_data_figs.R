#--- Notes ---#
# using fullCyc dataset subset for generating figures


#--- init ---#
library(dplyr)
library(ggplot2)
devtools::load_all('.')

#--- setting variables ---#
# input
inDir = '~/dev/HTSSIP/manuscript/data/'
# output
outDir = '~/dev/HTSSIP/manuscript/figures/'
## adjusted P-value cutoff
padj_cutoff = 0.1
## number of cores for parallel processing (increase depending on your computational hardware)
ncores = 20

#--- loading data ---#
inF = file.path(inDir, 'physeq_S2D2_l')
physeq_S2D2_l = readRDS(inF)
inF = file.path(inDir, 'physeq_CR2_l')
physeq_CR2_l = readRDS(inF)


#--- BD ordination ---#
# treat-vs-control
## calculating ordination
doParallel::registerDoParallel(ncores)
physeq_S2D2_l_df = SIP_betaDiv_ord(physeq_S2D2_l, parallel=TRUE)
## formatting
physeq_S2D2_l_df_e = physeq_S2D2_l_df %>%
  dplyr::mutate(phyloseq_subset = gsub(' \\| ', '\n', phyloseq_subset),
                phyloseq_subset = gsub('\\(Substrate==\'', '', phyloseq_subset),
                phyloseq_subset = gsub('\' & Day *== *\'', ', D', phyloseq_subset),
                phyloseq_subset = gsub('\\)|\'', '', phyloseq_subset),
                phyloseq_subset = gsub('D3', 'D03', phyloseq_subset)) %>%
  dplyr::filter(Buoyant_density >= 1.69)
## plotting
p_ord = phyloseq_ord_plot(physeq_S2D2_l_df_e)
p_ord = p_ord +
  scale_size_continuous('Buoyant\ndensity',
                        range=c(0.5,5),
                        breaks=round(seq(1.68, 1.772, (1.772-1.69)/5), 3))
p_ord

# replicate controls
## calculating ordination
doParallel::registerDoParallel(ncores)
physeq_CR2_l_df = SIP_betaDiv_ord(physeq_CR2_l, parallel=TRUE)
## formatting
physeq_CR2_l_df_e = physeq_CR2_l_df %>%
  dplyr::mutate(phyloseq_subset = gsub(' \\| ', '\n', phyloseq_subset),
                phyloseq_subset = gsub('\\(Substrate==\'', '', phyloseq_subset),
                phyloseq_subset = gsub('Replicate==', 'Replicate=', phyloseq_subset),
                phyloseq_subset = gsub('\\)|\'', '', phyloseq_subset),
                Microcosm_replicate = Microcosm_replicate - 1) %>%
  dplyr::filter(Buoyant_density >= 1.69)
## plotting
p_ord_cont = phyloseq_ord_plot(physeq_CR2_l_df_e, point_shape='Microcosm_replicate')
p_ord_cont = p_ord_cont +
  scale_size_continuous('Buoyant\ndensity',
                        range=c(0.5,5),
                        breaks=round(seq(1.68, 1.772, (1.772-1.69)/5), 3)) +
  scale_shape_discrete('Microcosm\nreplicate')
p_ord_cont


#--- BD shift ---#
# treatment vs control
## calculating BD shift
doParallel::registerDoParallel(ncores)
wmean = plyr::ldply(physeq_S2D2_l, BD_shift,
                    perm_method='control',
                    parallel_perm=TRUE)
## formatting
wmean_m = wmean %>%
  mutate(Substrate = gsub('.+(13C-[A-z]+).+', '\\1', .id),
         Day = gsub('.+Day ==[ \']*([0-9]+).+', 'Day \\1', .id),
         Day = Day %>% reorder(gsub('Day ', '', Day) %>% as.numeric),
         BD_shift = wmean_dist > wmean_dist_CI_high) %>%
  arrange(Day, Substrate, BD_min.x) %>%
  group_by(Day, Substrate) %>%
  mutate(window = (BD_shift == TRUE & lag(BD_shift) == TRUE & lag(BD_shift, 2) == TRUE) |
                  (BD_shift == TRUE & lag(BD_shift) == TRUE & lead(BD_shift) == TRUE) |
                  (BD_shift == TRUE & lead(BD_shift) == TRUE & lead(BD_shift, 2) == TRUE),
         BD_shift = BD_shift == TRUE & window == TRUE,
         BD_shift = ifelse(is.na(BD_shift), FALSE, BD_shift))

## plotting, with facetting by 13C-treatment
x_lab = bquote('Buoyant density (g '* ml^-1*')')
y_lims = c(0, 0.38)
y_lims = c(0, round(max(wmean_m$wmean_dist_CI_high)+0.005, 2))
p_shift = ggplot(wmean_m, aes(BD_min.x, wmean_dist)) +
  geom_line(alpha=0.3) +
  geom_linerange(aes(ymin=wmean_dist_CI_low,
                     ymax=wmean_dist_CI_high),
                 size=0.7, alpha=0.3) +
  geom_point(aes(color=BD_shift)) +
  scale_y_continuous(limits=y_lims) +
  scale_color_discrete('Gradient\nfraction\nin BD shift\nwindow?') +
  labs(x=x_lab, y='Weighted Unifrac') +
  facet_grid(Day ~ Substrate) +
  theme_bw() +
  theme(title = element_text(size=14),
        axis.text.x = element_text(angle=45, hjust=1))
p_shift


# control vs control
## calculating BD shift
doParallel::registerDoParallel(ncores)
wmean = plyr::ldply(physeq_CR2_l, BD_shift,
                    perm_method='control',
                    parallel_perm=TRUE,
                    ex="Microcosm_replicate=='2'")
## formatting
wmean_m = wmean %>%
  mutate(Substrate = '12C-Con, Rep2 vs\n12C-Con, Rep1',
         BD_shift = wmean_dist > wmean_dist_CI_high) %>%
  arrange(Substrate, BD_min.x) %>%
  group_by(Substrate) %>%
  mutate(window = (BD_shift == TRUE & lag(BD_shift) == TRUE & lag(BD_shift, 2) == TRUE) |
                  (BD_shift == TRUE & lag(BD_shift) == TRUE & lead(BD_shift) == TRUE) |
                  (BD_shift == TRUE & lead(BD_shift) == TRUE & lead(BD_shift, 2) == TRUE),
         BD_shift = BD_shift == TRUE & window == TRUE,
         BD_shift = ifelse(is.na(BD_shift), FALSE, BD_shift))

## plotting
x_lab = bquote('Buoyant density (g '* ml^-1*')')
p_shift_cont = ggplot(wmean_m, aes(BD_min.x, wmean_dist)) +
  geom_line(alpha=0.3) +
  geom_linerange(aes(ymin=wmean_dist_CI_low,
                     ymax=wmean_dist_CI_high),
                 alpha=0.2) +
  geom_point(aes(color=BD_shift)) +
  scale_y_continuous(limits=y_lims) +
  scale_color_discrete('Gradient\nfraction\nin BD shift\nwindow?') +
  labs(x=x_lab, y='Weighted Unifrac') +
  facet_grid(. ~ Substrate) +
  theme_bw() +
  theme(title = element_text(size=14),
        axis.text.x = element_text(angle=45, hjust=1))
p_shift_cont


#--- combined figure ---#
p_explr = cowplot::ggdraw() +
  cowplot::draw_plot(p_ord_cont,   0.00, 0.55, 0.38, 0.33) +
  cowplot::draw_plot(p_ord,        0.38, 0.40, 0.62, 0.60) +
  cowplot::draw_plot(p_shift_cont, 0.00, 0.05, 0.38, 0.30) +
  cowplot::draw_plot(p_shift,      0.38, 0.00, 0.62, 0.40) +
  cowplot::draw_plot_label(c('A)', 'B)', 'C)', 'D)'),
                           c(0.00, 0.33, 0.00, 0.33),
                           c(1.00, 1.00, 0.40, 0.40), size = 13)
p_explr

outF = file.path(outDir, 'BD-ord_BD-shift.pdf')
ggsave(outF, p_explr, width=10, height=8.75)
cat('File written:', outF, '\n')
outF = file.path(outDir, 'BD-ord_BD-shift.png')
ggsave(outF, p_explr, width=10, height=8.75)
cat('File written:', outF, '\n')


# p_explr = cowplot::ggdraw() +
#   cowplot::draw_plot(p_ord, 0, .4, 1, .6) +
#   cowplot::draw_plot(p_shift, 0, 0, 1, .4) +
#   cowplot::draw_plot_label(c("A", "B"), c(0, 0), c(1, 0.4), size = 15)
# p_explr
# outF = file.path(outDir, 'BD-ord_BD-shift.pdf')
# ggsave(outF, p_explr, width=6, height=7)
# outF = file.path(outDir, 'BD-ord_BD-shift.png')
# ggsave(outF, p_explr, width=6, height=7)




#--- HR-SIP ---#
## all treatment-control comparisons
doParallel::registerDoParallel(ncores)
df_HRSIP = plyr::ldply(physeq_S2D2_l,
                      HRSIP,
                      design = ~Substrate,
                      padj_cutoff = padj_cutoff,
                      .parallel=TRUE)

## formatting
df_HRSIP = df_HRSIP %>%
  mutate(.id = gsub(' \\| ', '\n', .id))
df_HRSIP %>% .$.id %>% unique %>% print


## How many incorporators (rejected hypotheses) & which sparsity cutoff was used for each comparison?
HRSIP_nincorp = df_HRSIP %>%
  filter(padj < padj_cutoff) %>%
  group_by(.id, sparsity_threshold) %>%
  summarize(n_incorp_OTUs = OTU %>% unique %>% length) %>%
  as.data.frame

## A breakdown of incorporators for each phylum in each comparision.
### summarizing
df_HRSIP_s = df_HRSIP %>%
  filter(padj < padj_cutoff) %>%
  mutate(Rank2 = gsub('^__', '', Rank2)) %>%
  group_by(.id, Rank2) %>%
  summarize(n_incorp_OTUs = OTU %>% unique %>% length) %>%
  ungroup()

### plotting
df_HRSIP_s_p = ggplot(df_HRSIP_s, aes(Rank2, n_incorp_OTUs)) +
    geom_bar(stat='identity') +
    labs(x='Phylum', y='Number of incorporators') +
    facet_wrap(~ .id, scales='free') +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle=55, hjust=1)
    )
df_HRSIP_s_p


#--- MW-HR-SIP ---#
## windows
windows = data.frame(density_min=c(1.70, 1.72, 1.74),
                     density_max=c(1.73, 1.75, 1.77))

## run (parallel)
doParallel::registerDoParallel(ncores)
df_MWHRSIP = plyr::ldply(physeq_S2D2_l,
                      HRSIP,
                      density_windows = windows,
                      design = ~Substrate,
                      padj_cutoff = padj_cutoff,
                      .parallel = TRUE)
#df_MWHRSIP %>% head(n=3)


## how many incorps?
MWHRSIP_nincorp = df_MWHRSIP %>%
  filter(padj < padj_cutoff) %>%
  group_by(.id, sparsity_threshold) %>%
  summarize(n_incorp_OTUs = OTU %>% unique %>% length) %>%
  as.data.frame
MWHRSIP_nincorp

## summarizing & plotting
df_MWHRSIP_s = df_MWHRSIP %>%
  mutate(.id = gsub(' \\| ', '\n', .id)) %>%
  filter(padj < padj_cutoff) %>%
  mutate(density_range = paste(density_min, density_max, sep='-')) %>%
  group_by(.id, sparsity_threshold, density_range) %>%
  summarize(n_incorp_OTUs = OTU %>% unique %>% length)

df_MRHRggplot(df_MWHRSIP_s, aes(.id, n_incorp_OTUs, fill=density_range)) +
    geom_bar(stat='identity', position='fill') +
    labs(x='Control-treatment comparision', y='Fraction of incorporators') +
    scale_y_continuous(expand=c(0,0)) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle=55, hjust=1)
    )



## A breakdown of incorporators for each phylum in each comparision.
## summarizing & plotting
df_MWHRSIP_s = df_MWHRSIP %>%
  mutate(.id = gsub(' \\| ', '\n', .id)) %>%
  filter(padj < padj_cutoff) %>%
  mutate(Rank2 = gsub('^__', '', Rank2)) %>%
  group_by(.id, Rank2) %>%
  summarize(n_incorp_OTUs = OTU %>% unique %>% length) %>%
  ungroup()

ggplot(df_MWHRSIP_s, aes(Rank2, n_incorp_OTUs)) +
    geom_bar(stat='identity') +
    labs(x='Phylum', y='Number of incorporators') +
    facet_wrap(~ .id, scales='free') +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle=55, hjust=1)
    )


#--- qSIP ---#
## dataset
physeq_rep3_qPCR_sum = physeq_rep3_qPCR$summary
physeq_rep3_t = OTU_qPCR_trans(physeq_rep3, physeq_rep3_qPCR)
## calculating atom fraction excess
atomX = qSIP_atom_excess(physeq_rep3_t,
                         control_expr='Treatment=="12C-Con"',
                         treatment_rep='Replicate')
## bootstrapping
doParallel::registerDoParallel(ncores)
df_atomX_boot = qSIP_bootstrap(atomX, n_boot=nboot, parallel=TRUE)
## applying CI threshold
CI_threshold = 0
df_atomX_boot = df_atomX_boot %>%
  mutate(Incorporator = A_CI_low > CI_threshold,
         OTU = reorder(OTU, -A))
## number of incorporators
qSIP_nincorp = df_atomX_boot %>%
  filter(Incorporator == TRUE) %>%
  nrow
## plotting
ggplot(df_atomX_boot, aes(OTU, A, ymin=A_CI_low, ymax=A_CI_high, color=Incorporator)) +
  geom_pointrange(size=0.25) +
  geom_linerange() +
  geom_hline(yintercept=0, linetype='dashed', alpha=0.5) +
  labs(x='OTU', y='Atom fraction excess') +
  theme_bw() +
  theme(
    axis.text.x = element_blank()
  )


