# Variables
baseDir = '/var/seq_data/fullCyc/MiSeq_16SrRNA/515f-806r/lib1-7rs/'
biomFile = file.path(baseDir, 'OTU_binning/otu_table_wtax.biom')
metadataFile = file.path(baseDir, 'metadata_SIP.txt')
treeFile = file.path(baseDir, 'fasttree/otusn_pick.tree')

# Init
library(ggplot2)
library(dplyr)
library(tidyr)
library(phyloseq)

# Loading data
## biom file
physeq_all = import_biom(biomFile) #treeFile)

## loading sample metadata
sample.data = import_qiime_sample_data(metadataFile)
physeq_all = merge_phyloseq(physeq_all,sample.data)
physeq_all_m = physeq_all %>% sample_data


# Partitioning dataset
## con + cel + xyl
exp_types = c('SIP')
sub_types = c('12C-Con', '13C-Cel', '13C-Xyl')
physeq_all_p = prune_samples(physeq_all_m$Exp_type %in% exp_types &
                         physeq_all_m$core_dataset == TRUE &
                         physeq_all_m$Substrate %in% sub_types,
                         physeq_all) %>%
  filter_taxa(function(x) sum(x) > 0, TRUE)
physeq_all_p
### saving
setwd(file.path(getwd(), 'data-raw'))
save(physeq_all_p, file='fullCyc_con-cel-xyl.RData')

# making smaller dataset
physeq_all_p_r = rarefy_even_depth(physeq_all_p, 1000)
save(physeq_all_p_r, file='fullCyc_con-cel-xyl_r1000.RData')

# small sample size dataset
exp_types = c('SIP')
sub_types = c('12C-Con', '13C-Cel', '13C-Glu')
days = c(14)
physeq = prune_samples(physeq_all_m$Exp_type %in% exp_types &
                       physeq_all_m$core_dataset == TRUE &
                       physeq_all_m$Substrate %in% sub_types &
                       physeq_all_m$Day %in% days,
                       physeq_all) %>%
  filter_taxa(function(x) sum(x) > 0, TRUE)
save(physeq, file='fullCyc_con-cel-glu.RData')
unlink('../data/physeq.RData')
file.copy('fullCyc_con-cel-glu.RData', '../data/physeq.RData')
