library(tidyverse)
library(tsibble)
library(tflattices)
set.seed(12345)
hosp <- read_csv("data/hosp_age_counts.csv", col_types = "_Dii")
hosp <- hosp %>%
  arrange(desc(age)) %>%
  mutate(ddd = yearweek(adm_date),
         age_bin = cut_interval(age, length = 5)) %>%
  filter(adm_date > "2020-01-01")

hosp_grid <- hosp %>%
  group_by(age_bin,ddd) %>%
  summarise(count = sum(count)) %>% 
  pivot_wider(id_cols = ddd, names_from = age_bin, values_from = count,
              values_fill = 0L) %>%
  arrange(ddd)

hosp_mat <- t(as.matrix(hosp_grid[,-1]))
sc <- function(x) (x - min(x)) / (max(x) - min(x))

t1 <- etfgrid(hosp_mat, exp_fam = "poisson", korder = c(1L, 1L), 
              lambda_max = 1000, maxiter = 10000L, estimate_df = TRUE)
kl <- kl_estimate(t1)
  
data.frame(lambda = kl$lambda, kl = kl$kl_estimate) %>%
  ggplot(aes(lambda, kl)) + 
  geom_point() +
  scale_x_log10() +
  theme_bw()

best <- 12 # minimizes kl subject to underregularization as lambda -> 0

name_vec = c(lambda = "estimated rate", observed = "observed")
peaks <- tibble(date = yearweek(c("2020 W33", "2020 W30", "2020 W52", "2020 W53")), 
                age = c(60, 55, 60, 85),
                name = "lambda")

tibble(date = rep(hosp_grid$ddd, each = 15),
       age = rep(3:17*5, times = 48),
       lambda = t1$theta[,best],
       observed = c(hosp_mat)) %>%
  #mutate(date = lubridate::ymd(date)) %>%
  pivot_longer(lambda:observed) %>%
  ggplot(aes(date, age)) +
  geom_raster(aes(fill = value)) +
  facet_wrap(~ name, labeller = labeller(name = name_vec)) +
  scale_x_yearweek(date_breaks = "3 months", date_labels = "%b",
                   expand = expansion()) +
  scale_y_continuous(expand = expansion()) +
  theme_bw(base_family = "Times", base_size = 16) + xlab("") +
  guides(fill = guide_colourbar(barheight = 12)) +
  theme(legend.title = element_blank()) +
  scale_fill_viridis_c(trans = "pseudo_log", breaks = c(0,2,5,10,20,35,60)) +
  geom_point(data = peaks)

ggsave("gfx/hospitals.pdf", width = 8, height = 4)
