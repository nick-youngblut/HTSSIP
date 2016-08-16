test_that('qPCR simulation (control/treatment)', {
  control_mean_fun = function(x) dnorm(x, mean=1.70, sd=0.01) * 1e8
  control_sd_fun = function(x) control_mean_fun(x) / 3
  treat_mean_fun = function(x) dnorm(x, mean=1.75, sd=0.01) * 1e8
  treat_sd_fun = function(x) treat_mean_fun(x) / 3


  df_qPCR = qPCR_sim(physeq_S2D2,
         control_expr='Substrate=="12C-Con"',
         control_mean_fun=control_mean_fun,
         control_sd_fun=control_sd_fun,
         treat_mean_fun=treat_mean_fun,
         treat_sd_fun=treat_sd_fun)

  df_qPCR_s = df_qPCR %>%
    gather(qPCR_tech_rep_id, qPCR_tech_rep_value, starts_with('qPCR_tech_rep')) %>%
    group_by(IS_CONTROL, Buoyant_density) %>%
    summarize(mean_value = mean(qPCR_tech_rep_value),
            sd_value = sd(qPCR_tech_rep_value)) %>%
    ungroup()


  ggplot(df_qPCR_s, aes(Buoyant_density, mean_value,
                      ymin=mean_value-sd_value,
                      ymax=mean_value+sd_value,
                      color=IS_CONTROL)) +
    geom_pointrange() +
    theme_bw()
})
