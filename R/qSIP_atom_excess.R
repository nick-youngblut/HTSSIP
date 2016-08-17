# input: phyloseq (OTU counts transformed)
## must be able to distinguish experimental design
### experimental design table? (sample = control/treatment); or expression

# calculate weighted average (W)
# calculating mean W for control & treatment
# calculate BD shift (Z)

# calculating density shifts
# calculating atom fraction excess
# bootstrapping & calculating CIs



# W: weighted average (BD) of taxon, with weights by abundance
## calculated on a per-gradient basis
# mean W (labeled vs control)
# Z = Wm_lab - Wm_light
# Z = BD shift


calc_Gi = function(Wlight){
  if(length(Wlight) > 1){
    Gi = sapply(Wlight, calc_Gi)
  } else{
    Gi = 1 / 0.083506 * (Wlight - 1.646057)
  }
  return(Gi)
}


calc_Mheavymax = function(Mlight, isotope='13C', Gi=NA){
  isotope = toupper(isotope)
  if(isotope=='13C'){
    Mhm = -0.4987282 * Gi + 9.974564 + Mlight
  } else
    if(isotope=='18O'){
      Mhm = 12.07747 + Mlight
    } else {
      stop('isotope not recognized')
    }
  return(Mhm)
}


calc_atom_excess = function(Mlab, Mlight, Mheavymax, isotope='13C'){
  isotope=toupper(isotope)
  if(isotope=='13C'){
    x = 0.01111233
  } else
    if(isotope=='18O'){
      x = 0.002000429
    } else {
      stop('isotope not recognized')
    }
  A = (Mlab - Mlight) / (Mheavymax - Mlight) * (1 - x)
  return(A)
}

qSIP_atom_excess_format = function(physeq, control_expr, gradient_rep){
  # formatting input
  cols = c('IS_CONTROL', 'Buoyant_density', gradient_rep)
  df_OTU = phyloseq2table(physeq,
                          include_sample_data=TRUE,
                          sample_col_keep=cols,
                          control_expr=control_expr)
  return(df_OTU)
}


#' Calculate mean buoyant density shift
#'
#' @param physeq  A phyloseq object
#' @param qPCR  A list of qPCR data from \code(qPCR_sim())
#'
#' @return A data.frame object
#'
#' @export
#'
#' @examples
#'
#'
qSIP_atom_excess = function(physeq,
                            control_expr,
                            gradient_rep=NULL,
                            isotope='13C',
                            df_OTU_W=NULL){
  # formatting input
  if(is.null(df_OTU_W)){
    no_boot = TRUE
  } else {
    no_boot = FALSE
  }

  if(no_boot){
    df_OTU = qSIP_atom_excess_format(physeq, control_expr, gradient_rep)

    # BD shift (Z)
    df_OTU_W = df_OTU %>%
      # weighted mean buoyant density (W)
      dplyr::group_by_('IS_CONTROL', 'OTU', gradient_rep) %>%
      dplyr::mutate(Buoyant_density = Buoyant_density %>% as.Num,
                    Count = Count %>% as.Num) %>%
      dplyr::summarize(W = weighted.mean(Buoyant_density, Count, na.rm=TRUE))
  }


  df_OTU_s = df_OTU_W %>%
    # mean W of replicate gradients
    dplyr::group_by_('IS_CONTROL', 'OTU') %>%
    dplyr::summarize(Wm = mean(W, na.rm=TRUE)) %>%
    # BD shift (Z)
    dplyr::group_by_('OTU') %>%
    dplyr::mutate(IS_CONTROL = ifelse(IS_CONTROL==TRUE, 'Wlight', 'Wlab')) %>%
    tidyr::spread(IS_CONTROL, Wm) %>%
    dplyr::mutate(Z = Wlab - Wlight) %>%
    dplyr::ungroup()

  # atom excess (A)
  MoreArgs = list(isotope=isotope)
  df_OTU_s = df_OTU_s %>%
    dplyr::mutate(Gi = calc_Gi(Wlight),
                  Mlight = 0.496 * Gi + 307.691,
                  Mheavymax = mapply(calc_Mheavymax,
                                     Mlight=Mlight,
                                     Gi=Gi,
                                     MoreArgs=MoreArgs),
                  Mlab = (Z / Wlight + 1) * Mlight,
                  A = mapply(calc_atom_excess,
                             Mlab=Mlab,
                             Mlight=Mlight,
                             Mheavymax=Mheavymax,
                             MoreArgs=MoreArgs))

  if(no_boot){
    return(list(W=df_OTU_W, A=df_OTU_s))
  } else {
    return(df_OTU_s)
  }
}


# sampling with replacement from control & treatment for each OTU
sample_W = function(df, n_sample){
  n_light = n_sample[1]
  n_lab = n_sample[2]
  # parsing df
  df_light = df[df$IS_CONTROL==TRUE,]
  df_lab = df[df$IS_CONTROL==FALSE,]
  # sampling
  if(length(df_light$W) > 1){
    W_light = base::sample(df_light$W, n_light, replace=TRUE)
  } else {
    W_light = rep(df_light$W, n_light)
  }
  if(length(df_lab$W) > 1){
    W_lab = base::sample(df_lab$W, n_lab, replace=TRUE)
  } else {
    W_lab = rep(df_lab$W, n_lab)
  }
  # creating new data.frames
  df_light = data.frame(IS_CONTROL=TRUE,
                        #OTU=rep(otu, n_light),
                        W = W_light)
  df_lab = data.frame(IS_CONTROL=FALSE,
                      #OTU=rep(otu, n_lab),
                      W = W_lab)
  return(rbind(df_light, df_lab))
}


# shuffling weighted mean densities (W)
.qSIP_bootstrap = function(atomX,
                           isotope='13C',
                           n_sample=c(3,3),
                           bootstrap_id = 1){
  # making a new (subsampled with replacement) dataset
  n_sample = c(3,3)  # control, treatment
  df_OTU_W = atomX$W %>%
    dplyr::group_by(OTU) %>%
    tidyr::nest() %>%
    dplyr::mutate(ret = lapply(data, sample_W, n_sample=n_sample)) %>%
    dplyr::select(-data) %>%
    tidyr::unnest()

  # calculating atom excess
  atomX = qSIP_atom_excess(physeq=NULL,
                   df_OTU_W=df_OTU_W,
                   control_expr=NULL,
                   gradient_rep=NULL,
                   isotope=isotope)
  atomX$bootstrap_id = bootstrap_id
  return(atomX)
}


qSIP_bootstrap = function(atomX, isotope='13C', n_sample=c(3,3),
                          n_boot=10, parallel=FALSE){
  # atom excess for each bootstrap replicate
  df_boot_id = data.frame(bootstrap_id = 1:n_boot)
  df_boot = plyr::mdply(df_boot_id, .qSIP_bootstrap,
                        atomX = atomX,
                        isotope=isotope,
                        n_sample=n_sample,
                        .parallel=parallel)

  # calculating atomX CIs for each OTU
  df_boot = df_boot %>%
    dplyr::group_by(OTU) %>%
    summarize(A_CI_low = quantile(A, a / 2, na.rm=TRUE),
              A_CI_high = quantile(A, 1 - a/2, na.rm=TRUE))

  # combining with atomX summary data
  df_boot = dplyr::inner_join(atomX$A, df_boot, c('OTU'='OTU'))
  return(df_boot)
}


