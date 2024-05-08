library(tidyverse)
library(purrr)
library(stringr)
library(ggthemes)
library(patchwork)

pTSG <- function(q, eps){
  # Function to compute the CDF of a Two-sided Geometric Distribution 
  # with parameter exp(-eps)
  alpha <- exp(-eps)
  (q < 0)*alpha^(-q)/(1+alpha) + (q >= 0)*(1 - alpha^(q + 1)/(1+alpha))
}

find_eps <- function(r_star, p_grid = seq(1e-3, 1, 1e-3), 
                     q_grid = seq(1e-3, 1, 1e-3)){
  # Function to search over a grid to find the minimal epsilon satisfying
  # the risk profile r_star
  min_eps <- Inf
  for(p in p_grid){
    for(q in q_grid){
      if(q == 1){ 
        if(p > 1/r_star(p,1)){ next }
        eps <- log((1-p)/(1/r_star(p,1) - p)) 
      }
      else{ 
        if(4*p*(1-q)*(1/r_star(p,q) - p*q) < 0){ next }
        eps <- log(2*p*(1-q)/(sqrt((1-p)^2 + 4*p*(1-q)*(1/r_star(p,q) - p*q)) - (1-p)))
      }
      if(min_eps > eps){ min_eps <- eps; min_p <- p; min_q <- q }
    }
  }
  return(list(eps = min_eps, p = min_p, q = min_q))
}


##################################### 
###   CODE TO PRODUCE FIGURE 1    ###
##################################### 


ta <- 0.1; tr <- 3; tq <- 0.999
df <- data.frame(p = seq(0.001, 0.999, 0.001)) %>%
  mutate(r_star = case_when(p <= ta/(tr*tq) ~ ta/(p*tq),
                            p <= 1 ~ tr)) %>%
  mutate(eps = log((2*p*(1-tq))/(sqrt((1-p)^2 + 4*p*(1-tq)*(1/r_star - p*tq)) - (1-p))))

p1 <- df %>%
  ggplot(aes(x = p, y = r_star)) + 
  geom_line(linewidth = 1) +
  scale_y_log10(limits = c(1,30)) +
  xlim(0, 1/3) +
  labs(x = expression("p"[i]), y = expression("r*(p"[i]*", 1)"), subtitle = "Agency 1") +
  theme_tufte() +
  theme(plot.subtitle = element_text(hjust = 0.5),
        panel.grid.major.y = element_line(color = "lightgray", linewidth = 0.5))

p2 <- df %>%
  ggplot(aes(x = p, y = eps)) +
  geom_line(linewidth = 1) +
  ylim(0,3) +
  xlim(0, 1/3) +
  labs(x = expression("p"[i]), y = expression(epsilon[i]*"(p"[i]*", 1)")) +
  theme_tufte() +
  theme(panel.grid.major.y = element_line(color = "lightgray", linewidth = 0.5))

tp <- 0.05
df <- data.frame(q = seq(0.001, 0.999, 0.001)) %>%
  mutate(r_star = case_when(q <= ta/(tr*tp) ~ ta/(tp*q),
                            q <= 1 ~ tr)) %>%
  mutate(eps = log((2*tp*(1-q))/(sqrt((1-tp)^2 + 4*tp*(1-q)*(1/r_star - tp*q)) - (1-tp))))

p3 <- df %>%
  ggplot(aes(x = q, y = r_star)) + 
  geom_line(linewidth = 1) +
  scale_y_log10(limits = c(1,30)) +
  labs(x = expression("q"[i]), y = expression("r*(0.05, q"[i]*")"), subtitle = "Agency 2") +
  theme_tufte() +
  theme(panel.grid.major.y = element_line(color = "lightgray", linewidth = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

p4 <- df %>%
  ggplot(aes(x = q, y = eps)) +
  geom_line(linewidth = 1) +
  ylim(0,3) +
  labs(x = expression("q"[i]), y = expression(epsilon[i]*"(0.05, q"[i]*")")) +
  theme_tufte() +
  theme(panel.grid.major.y = element_line(color = "lightgray", linewidth = 0.5))

(p1 + p3)/(p2 + p4) + plot_layout(axes = "collect")

ggsave("Figures/choosing_eps_risks.png", width = 6.5, height = 3, dpi=600, units = "in")



##################################### 
###   CODE TO PRODUCE FIGURE 2    ###
##################################### 


df <- expand.grid(ta = c(0, 0.1, 0.25, 0.5), tr = c(1.2, 2, 5)) %>%
  mutate(eps = log(tr - ta) - log(1 - ta),
         RMSE = sqrt(2*exp(-eps))/(1-exp(-eps)),
         prob_0.25 =  pTSG(24-25, eps),
         prob_0.5 =  pTSG(24-26, eps),
         prob_1 =  pTSG(24-28, eps),
         prob_2 =  pTSG(24-32, eps)) %>%
  mutate(ta = factor(ta, levels = c(0, 0.1, 0.25, 0.5), 
                     labels = c("0", ".1", ".25", ".5")))

labels = paste0("\u0072\u0303 = ", df$tr, "\n\u0061\u0303 = ", df$ta)
df = mutate(df, label = factor(labels, levels = labels))

p1 <- df %>%
  ggplot(aes(x = label, y = RMSE)) + 
  geom_col() +
  labs(x = "") +
  theme_tufte() +
  theme(panel.grid.major.y = element_line(color = "lightgray", linewidth = 0.5),
        axis.title.y = element_text(angle = 0, vjust = 0.5))

p2 <- df %>%
  dplyr::select(-c(ta, tr, eps, RMSE)) %>%
  pivot_longer(cols = -ncol(.), names_to = "deviation", values_to = "probability") %>%
  mutate(deviation = str_sub(deviation, start = 6)) %>%
  ggplot(aes(x = label, y = probability, fill = deviation)) +
  geom_col(position = "dodge") +
  labs(x = "", y = "Probability\nDecision\nChanges", fill = "Deviation") +
  theme_tufte() +
  theme(panel.grid.major.y = element_line(color = "lightgray", linewidth = 0.5),
        axis.title.y = element_text(angle = 0, vjust = 0.5))

p3 <- df %>%
  ggplot(aes(x = label, y = eps)) + 
  geom_col() +
  labs(x = "", y = "Implied \u03B5") +
  theme_tufte() +
  theme(panel.grid.major.y = element_line(color = "lightgray", linewidth = 0.5),
        axis.title.y = element_text(angle = 0, vjust = 0.5))

p2 / p1 / p3 + plot_layout(axes = "collect")

ggsave("Figures/choosing_eps_real_dat.png", width = 6.5, height = 4, dpi=600, units = "in")



##################################### 
###   CODE TO PRODUCE FIGURE 3    ###
##################################### 


r_a <- function(p, tr, ta){
  return(case_when(p <= ta/tr ~ ta/p, 
                   p <= 1 ~ tr))
}

p1 <- data.frame(p = 10^seq(-2, 0, length.out = 1e3)) %>%
  mutate(`Baseline 1` = r_a(p, 1.5, 0),
         `Baseline 2` = r_a(p, 3, 0),
         `Baseline 3` = r_a(p, 6, 0)) %>%
  pivot_longer(-1, names_to = "Agency", values_to = "r_star") %>%
  mutate(eps = log((1-p)/(1/r_star - p)),
         eps = if_else(eps <= 4, eps, 4)) %>%
  ggplot(aes(x = p, y = r_star, color = eps)) +
  geom_line(linewidth = 2) +
  scale_x_log10(limits = c(0.01, 1/1.5),
                breaks = c(0.01, 0.03, 0.1, 0.3), 
                labels = c(".01", ".03", ".1", ".3")) +
  scale_y_log10(limits = c(1,30)) +
  scale_color_viridis_c(direction = -1) +
  labs(x = expression("p"[i]), y = expression("r*(p"[i]*", 1)"),
       color = expression("  "*epsilon)) +
  facet_grid(~Agency) +
  theme_tufte()  +
  theme(panel.grid.major.y = element_line(color = "lightgray", linewidth = 0.5),
        legend.position = "none", 
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

p2 <- data.frame(p = 10^seq(-2, 0, length.out = 1e3)) %>%
  mutate(`Agency 1` = r_a(p, 1.5, 0.25),
         `Agency 2` = r_a(p, 3, 0.25),
         `Agency 3` = r_a(p, 6, 0.25)) %>%
  pivot_longer(-1, names_to = "Agency", values_to = "r_star") %>%
  mutate(eps = log((1-p)/(1/r_star - p)),
         eps = if_else(eps <= 4, eps, 4)) %>%
  ggplot(aes(x = p, y = r_star, color = eps)) +
  geom_line(linewidth = 2) +
  scale_x_log10(limits = c(0.01, 1/1.5),
                breaks = c(0.01, 0.03, 0.1, 0.3), 
                labels = c(".01", ".03", ".1", ".3")) +
  scale_y_log10(limits = c(1,30)) +
  scale_color_viridis_c(direction = -1) +
  labs(x = expression("p"[i]), y = expression("r*(p"[i]*", 1)"),
       color = expression(epsilon[i]*"(p"[i]*", 1)")) +
  facet_grid(~Agency) +
  theme_tufte()  +
  theme(panel.grid.major.y = element_line(color = "lightgray", linewidth = 0.5))

p1/p2 + plot_layout(guides = "collect") 

ggsave("Figures/choosing_eps_ex1.png", width = 6.5, height = 4, dpi=600, units = "in")



##################################### 
###   CODE TO PRODUCE FIGURE 4    ###
##################################### 


r_b <- function(q, tr, ta, tp){
  return(case_when(q <= ta/(tr*tp) ~ ta/(tp*q),
                   q <= 1 ~ tr))
}

tp <- 0.05
data.frame(q = 10^seq(-2, 0, length.out = 1e3)) %>%
  #data.frame(p = seq(0.01, 0.999, 0.001)) %>%
  mutate(`Baseline` = r_b(q, 3, 0, 0.05),
         `Agency 1` = r_b(q, 3, 0.025, 0.05),
         `Agency 2` = r_b(q, 3, 0.15, 0.05),
         `Agency 3` = r_b(q, 3, 0.3, 0.05)) %>%
  pivot_longer(-1, names_to = "Agency", values_to = "r_star") %>%
  mutate(eps = log((2*tp*(1-q))/(sqrt((1-tp)^2 + 4*tp*(1-q)*(1/r_star - tp*q)) - (1-tp))),
         eps = if_else(eps <= 4, eps, 4),
         Agency = relevel(factor(Agency), ref = "Baseline")) %>%
  ggplot(aes(x = q, y = r_star, color = eps)) +
  geom_line(linewidth = 2) +
  scale_x_log10(limits = c(0.01, 1),
                breaks = c(0.01, 0.03, 0.1, 0.3, 1), 
                labels = c(".01", ".03", ".1", ".3", "1")) +
  scale_y_log10(limits = c(1,30)) +
  scale_color_viridis_c(direction = -1, limits = c(0.41, 4), oob = scales::squish) +
  labs(x = expression("q"[i]), y = expression("r*(0.05, q"[i]*")"),
       color = expression(epsilon[i]*"(0.05, q"[i]*")")) +
  facet_grid(~Agency) +
  theme_tufte()  +
  theme(panel.grid.major.y = element_line(color = "lightgray", linewidth = 0.5))

ggsave("Figures/choosing_eps_ex2.png", width = 6.5, height = 2.25, dpi=600, units = "in")



##################################### 
###   CODE TO PRODUCE FIGURE 5    ###
#####################################

r_star <- function(p, q){
  max(0.25/(p*q), 3)
}

# Find the optimal epsilon
optimum <- find_eps(r_star)
optimum

p1 <- expand.grid(p = 10^seq(-2, 0, length.out = 2e2),
                  q = 10^seq(-2, 0, length.out = 2e2)) %>%
  mutate(r_star = map2_dbl(.x = p, .y = q, .f = r_star)) %>%
  mutate(eps = log((2*tp*(1-q))/(sqrt((1-tp)^2 + 4*tp*(1-q)*
                                        (1/r_star - tp*q)) - (1-tp))),
         r_star = if_else(r_star <= 100, r_star, 100),
         eps = if_else(eps <= 4, eps, 4)) %>%
  ggplot(aes(x = p, y = q, fill = r_star)) + 
  geom_tile() + 
  scale_fill_viridis_c(option = "B", direction=-1, trans = "log",
                       breaks = c(3, 10, 30, 100)) +
  scale_x_log10() + 
  scale_y_log10() +
  geom_point(aes(x = optimum$p, y = optimum$q), color = "red", size = 2.5) +
  labs(x = expression("p"[i]), y = expression("q"[i]), 
       fill = expression("r*(p"[i]*", q"[i]*")")) +
  theme_tufte()

p2 <- expand.grid(p = 10^seq(-2, 0, length.out = 2e2),
                  q = 10^seq(-2, 0, length.out = 2e2)) %>%
  mutate(r_star = map2_dbl(.x = p, .y = q, .f = r_star)) %>%
  mutate(eps = if_else(q == 1, log((1-p)/(1/r_star - p)),
                       log((2*p*(1-q))/(sqrt((1-p)^2 + 4*p*(1-q)*
                                               (1/r_star - p*q)) - (1-p)))),
         eps = if_else(eps <= 4, eps, 4)) %>%
  ggplot(aes(x = p, y = q, fill = eps)) + 
  geom_tile() + 
  scale_x_log10() + 
  scale_y_log10() +
  geom_point(aes(x = optimum$p, y = optimum$q), color = "red", size = 2.5) +
  scale_fill_viridis_c(direction=-1) +
  labs(x = expression("p"[i]), y = expression("q"[i]), 
       fill = expression(epsilon[i]*"(p"[i]*", q"[i]*")")) +
  theme_tufte()

p1 / p2 + plot_layout(axis_titles = "collect_x")

ggsave("Figures/choosing_eps_ex_2D.png", width = 6.5, height = 6, dpi=600, units = "in")