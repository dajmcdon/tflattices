library(tflattices)
library(tidyverse)
data("world_temperatures")
data("canada_map")

date_to_use <- "2010-01-01"
canada <- world_temperatures %>%
  filter(longitude > -142, longitude < -51, latitude > 42) %>%
  select(longitude, latitude, !!date_to_use) %>%
  rename(temperature = !!date_to_use) %>%
  mutate(temp_c = temperature - mean(temperature)) %>% 
  arrange(longitude, desc(latitude))





ca_v0 <- etfgrid(canada$temp_c, exp_fam = "gauss_var", dims = c(20,17), 
                 lambda_max = 1000, maxiter = 15000, tol = 1e-6) 
ca_v1 <- etfgrid(canada$temp_c, exp_fam = "gauss_var", dims = c(20,17),
                 k = c(1,1), lambda_max = 1000, maxiter = 15000, tol = 1e-6, 
                 estimate_df = TRUE)
ca_v2 <- etfgrid(canada$temp_c, exp_fam = "gauss_var", dims = c(20,17),
                 k = c(2,2), lambda_max = 1000, maxiter = 15000, tol = 1e-6)
# ca_m0 <- etfgrid(-canada$temp_c^2, dims = c(20,17),
#                  k = c(0,0), lambda_max = 1000, maxiter = 15000, tol = 1e-6) 

# ggplot(data.frame(x = canada$longitude, y = canada$latitude, 
#                   z = abs(canada$temp_c))) + 
#                   #z = ca_v1$theta[,35])) +
#   geom_tile(aes(x,y,fill=z)) +
#   scale_fill_viridis_c() 


tf_plots <- tibble(long = canada$longitude, lat = canada$latitude,
                   temp = canada$temp_c^2, 
                   #tfm = ca_m0$theta[,35], #looks good
                   tfv0a = ca_v0$theta[,20], # should be sqrt(n)/log(n)* lam[k=0]?
                   tfv1_wig = ca_v1$theta[,6],
                   tfv1 = ca_v1$theta[,26], # best by kl
                   tfv1_smooth = ca_v1$theta[,45],
                   tfv0b = ca_v2$theta[,20] # should be same as mean
                   ) %>%
  #mutate(tfm = -pmin(tfm,0)) %>%
  mutate(across(temp:tfv0b, function(x) x)) %>%
  pivot_longer(temp:tfv0b)

labs = c(temp = "|y - mean(y)|", 
         #tfm = "Mean TF (k = 0)", 
         tfv0a = "k = 0", 
         tfv1_wig = "k = 1, df = 330",
         tfv1 = "k = 1, kl optimal", 
         tfv1_smooth = "k = 1, df = 125",
         tfv0b = "k = 2")

ggplot(tf_plots, aes(long, lat)) + 
  geom_tile(aes(fill = sqrt(value))) +
  geom_path(data=canada_map,
            aes(x=long,y=lat,group=group),
            color="grey80", size = .2) +
  facet_wrap(~name, labeller = labeller(name=labs), nrow = 2, strip.position = "bottom") +
  theme_void(base_family = "Times", base_size = 16) +
  scale_fill_viridis_c(name = "ºC") +
  coord_map(projection = 'lambert', parameters=c(49,77)) +
  guides(fill = guide_colourbar(barwidth = 20, barheight = 0.5)) +
  theme(legend.position = "bottom",
        strip.text = element_text(margin = margin(b = 5)))

ggsave("../figs/temp-ex.pdf", width = 8, height = 5)
