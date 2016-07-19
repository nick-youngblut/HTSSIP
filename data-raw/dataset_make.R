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
physeq = import_biom(biomFile) #treeFile)

## loading sample metadata
sample.data = import_qiime_sample_data(metadataFile)
physeq = merge_phyloseq(physeq,sample.data)
physeq.m = physeq %>% sample_data


# Partitioning dataset
## con + cel + xyl
exp_types = c('SIP')
sub_types = c('12C-Con', '13C-Cel', '13C-Xyl')
physeq.p = prune_samples(physeq.m$Exp_type %in% exp_types &
                         physeq.m$core_dataset == TRUE &
                         physeq.m$Substrate %in% sub_types,
                         physeq) %>%
  filter_taxa(function(x) sum(x) > 0, TRUE)
physeq.p
### saving
setwd(file.path(getwd(), 'data-raw'))
saveRDS(physeq.p, 'fullCyc_con-cel-xyl.rds')

# making smaller dataset
physeq.p.r = rarefy_even_depth(physeq.p, 1000)
saveRDS(physeq.p.r, 'fullCyc_con-cel-xyl_r1000.rds')

