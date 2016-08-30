data(physeq_rep3)

df = data.frame(
  Buoyant_density = seq(1.68, 1.76, 0.01),
  Count = rnorm(9, mean=1000, sd=100)
)

BD_interp = lin_interp(df)

dBD = delta_BD(physeq_rep3,
               control_expr='Treatment=="12C-Con"',
               treatment_rep='Replicate')

dBD = delta_BD(physeq_rep3,
               n=20,
               control_expr='Treatment=="12C-Con"',
               treatment_rep='Replicate')

#dBD %>% .$Buoyant_density %>% table %>% length
dBD %>% head
dBD %>% filter(OTU=='OTU.1')
