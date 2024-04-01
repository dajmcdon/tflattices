library(tidyverse)
library(cowplot)

orange <- '#FF9200'
darkblue <- "#053B64"
green = '#00884E'


temp <- readRDS("data/northern_hemisphere_temperatures.rds")
ests <- readRDS("data/estimated_northern_hemisphere_variance.rds")
toronto_dist <- readRDS("data/toronto_distributions.rds")
world_df <- readRDS("data/world_map_lines.rds")
toronto_loc <- tibble(lats = 43.67, lons = -79.22)


# plot the CHANGE in means for 1960 vs 2000
g1 <- ggplot() +
  geom_tile(data = temp, aes(x = lons, y = lats, fill = temp)) +
  facet_wrap(~ period) +
  scale_fill_gradient2(low = darkblue, high=orange, limits=c(-13, 13),
                       breaks = c(-12,-4,-1,0,1,4,12),
                       guide = guide_colorbar(barheight = 10, barwidth = 0.5),
                       name = "ºC", trans = scales::modulus_trans(0)) + 
  geom_path(data = world_df, aes(x = long, y = lat, group = group)) +
  geom_point(data = toronto_loc, aes(lons, lats), stroke = 1, color = green,
             size = 3, shape = 21, fill=NA) +
  coord_map("ortho", orientation=c(90, 0, 0)) + 
  theme_void(base_family = "Times", base_size = 16) +
  theme(strip.text = element_text(margin = margin(b = 5)))

# ggsave("gfx/average_chg_globes.pdf", g1, width = 8, height = 5)



# plot the CHANGE in std for 1960 vs 2000
g2 <- ggplot() +
  geom_tile(data = ests, aes(x = lons, y = lats, fill = temp)) +
  facet_wrap(~ period) +
  scale_fill_gradient2(low=darkblue, high=orange, limits=c(-2.81,2.81),
                       breaks = c(-2,-1,0,1,2),
                       guide = guide_colorbar(barheight = 10, barwidth = 0.5),
                       name = "ºC", trans = scales::modulus_trans(-1)) + 
  geom_path(data = world_df, aes(x = long, y = lat, group = group)) +
  geom_point(data = toronto_loc, aes(lons, lats), stroke = 1, color = orange,
             size = 3, shape = 21, fill=NA) +
  coord_map("ortho", orientation=c(90, 0, 0)) + 
  theme_void(base_family = "Times", base_size = 16) +
  theme(strip.text = element_text(margin = margin(b = 5)))

# ggsave("gfx/std_chg_globes.pdf", g2, width = 8, height = 5)
  
g3 <- cowplot::plot_grid(
  g1, g2, nrow=2, labels="AUTO",
  label_fontfamily = "Times", label_size = 16
)

ggsave("gfx/both_globes.pdf", g3, width = 8, height = 6)

toronto_dist %>% 
  mutate(summer = fct_relevel(summer, "summer")) %>%
  ggplot(aes(x = temp-273.15, color = period1, fill = period1)) + 
  geom_density(alpha=.5) + xlab("") +
  facet_wrap(~summer, scales="free_x") +
  scale_fill_manual(values=c(orange, darkblue), labels=c("2000s","1960s")) +
  scale_color_manual(values=c(orange, darkblue), labels=c("2000s","1960s")) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_x_continuous(expand = expansion(mult = c(0, 0))) +
  theme_bw(base_family = "Times", base_size = 16) +
  theme(legend.position = "right",
        legend.title = element_blank())#,
        #legend.justification = c(0,1))

ggsave("gfx/toronto_densities.pdf", width = 8, height = 3)

