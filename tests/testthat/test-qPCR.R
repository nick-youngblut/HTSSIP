test_that('qPCR simulation (control/treatment)', {
  control_mean_fun = function(x) dnorm(x, mean=1.70, sd=0.01) * 1e8
  control_sd_fun = function(x) control_mean_fun(x) / 3
  treat_mean_fun = function(x) dnorm(x, mean=1.75, sd=0.01) * 1e8
  treat_sd_fun = function(x) treat_mean_fun(x) / 3

  qPCR = qPCR_sim(physeq_S2D2,
         control_expr='Substrate=="12C-Con"',
         control_mean_fun=control_mean_fun,
         control_sd_fun=control_sd_fun,
         treat_mean_fun=treat_mean_fun,
         treat_sd_fun=treat_sd_fun)

  p = ggplot(qPCR$summary, aes(Buoyant_density, qPCR_tech_rep_mean,
                      ymin=qPCR_tech_rep_mean-qPCR_tech_rep_sd,
                      ymax=qPCR_tech_rep_mean+qPCR_tech_rep_sd,
                      color=IS_CONTROL)) +
    geom_pointrange() +
    theme_bw()

  expect_is(p, 'ggplot')
})
