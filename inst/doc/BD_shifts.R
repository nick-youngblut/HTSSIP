## ---- message=FALSE, warning=FALSE---------------------------------------
library(dplyr)
library(tidyr)
library(ggplot2)
library(HTSSIP)

## ---- message=FALSE, warning=FALSE---------------------------------------
physeq_S2D2

## ------------------------------------------------------------------------
params = get_treatment_params(physeq_S2D2, c('Substrate', 'Day'))
params = dplyr::filter(params, Substrate!='12C-Con')
ex = "(Substrate=='12C-Con' & Day=='${Day}') | (Substrate=='${Substrate}' & Day == '${Day}')"
physeq_S2D2_l = phyloseq_subset(physeq_S2D2, params, ex)
physeq_S2D2_l

## ------------------------------------------------------------------------
wmean1 = BD_shift(physeq_S2D2_l[[2]])
cat('Subset:', names(physeq_S2D2_l)[2], '\n')
wmean1 %>% head(n=3)

## ---- fig.height=3.5, fig.width=7----------------------------------------
x_lab = 'Buoyant density (g/ml^-1)'
y_lab = 'Weighted mean of\nweighted-Unifrac distances'
ggplot(wmean1, aes(BD_min.x, wmean_dist)) +
  geom_point() +
  labs(x=x_lab, y=y_lab, title='Beta diversity of 13C-treatment relative to 12C-Con') +
  theme_bw() 

## ------------------------------------------------------------------------
wmean = plyr::ldply(physeq_S2D2_l, BD_shift)
wmean %>% head(n=3)

## ---- fig.height=5, fig.width=7------------------------------------------
# formatting the treatment names to look a bit better as facet labels
wmean = wmean %>%
  mutate(Substrate = gsub('.+(13C-[A-z]+).+', '\\1', .id),
         Day = gsub('.+Day ==[ \']*([0-9]+).+', '\\1', .id),
         Day = Day %>% reorder(Day %>% as.numeric))

# plotting, with facetting by 13C-treatment
ggplot(wmean, aes(BD_min.x, wmean_dist)) +
  geom_point() +
  labs(x=x_lab, y=y_lab, 
       title='Beta diversity of 13C-treatments relative to 12C-Con') +
  facet_grid(Day ~ Substrate) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1))

## ------------------------------------------------------------------------
sessionInfo()

