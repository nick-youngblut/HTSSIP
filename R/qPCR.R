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

# Simulate qPCR values
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
#' @examples
#' # making functions for simulating values
#' control_mean_fun = function(x) dnorm(x, mean=1.70, sd=0.01) * 1e8
#' control_sd_fun = function(x) control_mean_fun(x) / 3
#' treat_mean_fun = function(x) dnorm(x, mean=1.75, sd=0.01) * 1e8
#' treat_sd_fun = function(x) treat_mean_fun(x) / 3
#' # simulating qPCR values
#' df_qPCR = qPCR_sim(physeq_S2D2,
#'                 control_expr='Substrate=="12C-Con"',
#'                 control_mean_fun=control_mean_fun,
#'                 control_sd_fun=control_sd_fun,
#'                 treat_mean_fun=treat_mean_fun,
#'                 treat_sd_fun=treat_sd_fun)
#'
qPCR_sim = function(physeq,
                    control_mean_fun,
                    control_sd_fun,
                    treat_mean_fun,
                    treat_sd_fun,
                    n_tech_rep=3,
                    control_expr=NULL,
                    ...){
  # sample_metadata
  m = phyloseq2df(physeq, sample_data)
  if(is.null(m$IS_CONTROL)){
    if(is.null(control_expr)){
      stop('sample_data in phyloseq object must have IS_CONTROL column (logical), or you must provide the control_expr parameter')
    }
    m = m %>%
      dplyr::mutate_(IS_CONTROL = control_expr)
  }

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
  df_qPCR$Sample = physeq %>% sample_data %>% rownames

  # gather & summarize
  df_qPCR_s = df_qPCR %>%
    tidyr::gather(qPCR_tech_rep_id, qPCR_tech_rep_value, starts_with('qPCR_tech_rep')) %>%
    dplyr::group_by(IS_CONTROL, Sample, Buoyant_density) %>%
    dplyr::summarize(qPCR_tech_rep_mean = mean(qPCR_tech_rep_value),
           qPCR_tech_rep_sd = sd(qPCR_tech_rep_value)) %>%
    dplyr::ungroup() %>%
    as.data.frame

  return(list(raw=df_qPCR, summary=df_qPCR_s))
}



