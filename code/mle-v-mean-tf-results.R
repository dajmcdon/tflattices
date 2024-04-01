library(batchtools)
library(tidyverse)
library(cowplot)
library(scales)
theme_set(theme_bw(base_size = 12, base_family = "serif"))
loadRegistry("mle-v-mean-tf")
cols <- c("#0b62a4", "#ff9200")



# collect results ---------------------------------------------------------


res <- reduceResultsDataTable()
wide_res <- unnest_wider(res, result)
pars <- getJobPars()
res <- inner_join(unwrap(pars), wide_res, by = "job.id")



# Some summary plots ------------------------------------------------------

GeoMean <- function(x, na.rm = TRUE) {
  exp(mean(log(x), na.rm = na.rm))
}



mtv <- res %>% 
  filter(control_tv == "mean")
# sum_mtv <- mtv %>%
#   group_by(n, tf_version, distribution) %>%
#   summarise(mse = GeoMean(mse))
mtv_plot <- ggplot(mtv, aes(n, mse, color = tf_version)) +
  geom_jitter(alpha = .25, size = .5) +
  #geom_line(data = sum_mtv) +
  geom_smooth(na.rm = TRUE, se = FALSE) +
  facet_wrap(~ distribution, ncol = 1, scales = "free_y") +
  ylab("mean squared error") +
  scale_y_log10(breaks = c(.0001,.001,.01,.1, 1),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10() +
  scale_color_manual(values = cols, name = "trend filter") +
  theme(legend.position = "bottom")

kltv <- res %>% 
  filter(control_tv == "nat_par")
# sum_kltv <- kltv %>%
#   group_by(n, tf_version, distribution) %>%
#   summarise(kl = GeoMean(kl))
kltv_plot <- ggplot(kltv, aes(n, kl, color = tf_version)) +
  geom_jitter(alpha = .25, size = .5) +
  #geom_line(data = sum_mtv) +
  geom_smooth(na.rm = TRUE, se = FALSE) +
  facet_wrap(~ distribution, ncol = 1) +
  ylab("KL divergence") +
  coord_cartesian(ylim = c(10e-4, 1)) +
  scale_y_log10(breaks = c(.0001, .001, .01, .1, 1),#) +
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10() +
  scale_color_manual(values = cols, name = "trend filter") +
  theme(legend.position = "bottom")

leg <- get_legend(mtv_plot)

err <- plot_grid(
  mtv_plot + theme(legend.position = "none"),
  kltv_plot + theme(legend.position = "none"),
  nrow = 1,
  labels = LETTERS[1:2])
err_plot <- plot_grid(err, leg, ncol = 1, rel_heights = c(1,.1))
ggsave("gfx/sim-errors.pdf", err_plot, width = 6, height = 5)


# Show some estimates -----------------------------------------------------

nn <- 104

make_signal_v <- function(n) {
  # V shaped (tf has trouble at the boundaries, so put the hard part 
  # in the middle)
  n <- as.integer(round(n))
  n2 <- n %/% 2
  r <- seq(1, n2 + 1, length.out = n2 + 1)
  l <- rev(r)
  r <- r[-1]
  if (n %% 2 == 0L) l <- l[-(n2 + 1)]
  return(c(l,r) / n)
}

est <- res %>% 
  filter(n == nn) %>%
  select(job.id, distribution, control_tv, tf_version, thhat_mse, thhat_kl) %>%
  pivot_longer(starts_with("thhat")) %>%
  mutate(name = ifelse(str_detect(name, "mse"), "mean", "nat_par")) %>%
  filter(control_tv == name) %>%
  unnest(cols = "value") %>%
  rename(control = control_tv) %>%
  mutate(control = ifelse(control == "nat_par", "natural parameter", "mean")) %>%
  mutate(tf_version = paste0(tf_version, " trend filter"))
est$x <- rep(1:nn / (nn + 1), length.out = nrow(est))

signal <- tibble(
  x  = rep(1:nn / (nn + 1), 4),
  sig = rep(make_signal_v(nn), 4),
  distribution = rep(c("exponential", "poisson"), each = 2*nn),
  control = rep(rep(c("mean", "natural parameter"), each = nn), 2)) %>%
  mutate(value = case_when(
    distribution == "exponential" & control == "natural parameter" ~ 1 / sig,
    distribution == "poisson" & control == "mean" ~ 0.5 - sig + log(nn),
    distribution == "poisson" & control == "natural parameter" ~ exp(0.5 - sig + log(nn)),
    TRUE ~ sig
  ),
  job.id = -1,
  tf_version = "true signal")

ggplot(bind_rows(est, signal), aes(x, value)) + 
  geom_line(aes(group = job.id, color = tf_version, size = tf_version, alpha = tf_version)) +
  scale_size_manual(values = c(.4,.4,1.2), guide = "none") +
  scale_alpha_manual(values = c(.4,.4,1), guide = "none") +
  facet_wrap(~ distribution + control, scales = "free_y", nrow = 1,
             labeller = label_both) +
  # scale_color_brewer(palette = "Dark2", name = "trend filter", 
  #                    guide = guide_legend(override.aes = list(alpha = 1))) +
  scale_color_manual(values = c(cols, "black"), name = "", 
                     guide = guide_legend(override.aes = list(alpha = 1, size = 1))) +
  xlab("") +
  ylab(expression(hat(theta))) +
  theme(legend.position = "bottom") 
ggsave("gfx/show-estimate.pdf", width = 8, height = 4)
