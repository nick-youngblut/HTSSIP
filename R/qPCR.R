# simulating qPCR data for a phyloseq object
# OTU_SIP_rank_abunds = function(df_OTU, exp_design){
#   exp_design_plus = c('OTU', exp_design)
#   df_OTU %>%
#     group_by_(.dots=as.list(exp_design_plus)) %>%
#     mutate(Max_count = max(Count)) %>%
#     group_by_(.dots=as.list(exp_design)) %>%
#     mutate(Rank_count = -Max_count %>% as.factor %>% as.numeric) %>%
#     ungroup()
# }


.qPCR_sim = function(Buoyant_density,
                     IS_CONTROL,
                     control_mean_fun,
                     control_sd_fun,
                     treat_mean_fun,
                     treat_sd_fun,
                     n_tech_rep=3,
                     ...){
  #Buoyant_density = Buoyant_density %>% as.Num
  if(IS_CONTROL==TRUE){
    M = control_mean_fun(Buoyant_density, ...)
    S = control_sd_fun(Buoyant_density, ...)
    X = rnorm(n=n_tech_rep, mean=M, sd=S)
  } else {
    M = treat_mean_fun(Buoyant_density, ...)
    S = treat_sd_fun(Buoyant_density, ...)
    X = rnorm(n=n_tech_rep, mean=M, sd=S)
  }
  X = ifelse(X < 0, 0, X)
  return(X)
}

#' Simulate qPCR values
#'
#' @param physeq  Object of class "phyloseq"
#' @param control_expr  Expression used to identify control samples based on sample_data.
#' @param control_mean_fun  Function used for simulating the qPCR normal distribution mean
#' for control samples.
#' @param control_sd_fun  Function used for simulating the qPCR normal distribution
#' standard deviation for control samples.
#' @param treat_mean_fun  Function used for simulating the qPCR normal distribution mean
#' for treatment samples.
#' @param treat_sd_fun  Function used for simulating the qPCR normal distribution
#' standard deviation for treatment samples.
#' @param n_tech_rep  Number of technical replicates.
#'
#' @return data.frame of qPCR values
#'
qPCR_sim = function(physeq,
                    control_expr,
                    control_mean_fun,
                    control_sd_fun,
                    treat_mean_fun,
                    treat_sd_fun,
                    n_tech_rep=3,
                    ...){
  # sample_metadata
  m = phyloseq2df(physeq_S2D2, sample_data) %>%
    dplyr::mutate_(IS_CONTROL = control_expr)

  if(is.null(m$Buoyant_density)){
    stop('Buoyant_density column not found in phyloseq object sample_data')
  }
  m_f = m %>%
    dplyr::mutate(Buoyant_density = as.Num(Buoyant_density)) %>%
    dplyr::select(Buoyant_density, IS_CONTROL)


  df_qPCR = plyr::mdply(m_f, .qPCR_sim,
                        control_mean_fun=control_mean_fun,
                        control_sd_fun=control_sd_fun,
                        treat_mean_fun=treat_mean_fun,
                        treat_sd_fun=treat_sd_fun,
                        n_tech_rep=n_tech_rep)

  colnames(df_qPCR)[3:(2+n_tech_rep)] = gsub('^', 'qPCR_tech_rep', 1:n_tech_rep)

  return(df_qPCR)
}



