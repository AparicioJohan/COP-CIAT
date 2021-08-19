

library(plotly)
library(tidyverse)
library(desplot)
library(broom)
library(emmeans)
library(multcomp)
library(multcompView) # mean comparisons
library(ggrepel)      # labels ggplot
library(ggpubr)       # stats with ggplot

datos <- read.csv("data/Example_rice.csv") %>% type.convert()
head(datos)
str(datos)
summary(datos)


# Grafico usando ggplot2
g1 <- datos %>% 
  ggplot(aes(x = col, y = row, fill = Bloq )) +
  geom_tile(color = "black")+
  geom_text(aes(label = GEN))+
  theme_bw()+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  facet_wrap(~ TRAT, scales = "free_x")
g1

ggplotly(g1)

# -------------------------------------------------------------------------
# Chequeo de la variable de respuesta
# -------------------------------------------------------------------------

# summary by GEN
raw_m <- datos %>% 
  group_by(GEN, TRAT) %>% 
  summarize(mean    = mean(REND, na.rm = T),
            std.dev = sd(REND, na.rm = T),
            cv      = std.dev/mean ) %>% 
  print(n=Inf)
raw_m

# summary by Bloq
datos %>% 
  group_by(TRAT, Bloq) %>% 
  summarize(mean    = mean(REND),
            std.dev = sd(REND),
            cv      = std.dev/mean )

# summary by condition
datos %>% 
  group_by(TRAT) %>% 
  summarize(mean    = mean(REND),
            std.dev = sd(REND),
            cv      = std.dev/mean )

# points 
g2 <- datos %>% 
  ggplot(aes(x = GEN,  y = REND))+
  geom_point(size = 4)+
  theme_bw()+
  theme(axis.text.x = element_text(angle=75, vjust=0.5))
g2

g3 <- datos %>% 
  ggplot(aes(x = GEN,  y = REND, color = TRAT))+
  geom_point(size = 4, alpha = 0.5)+
  theme_bw()+
  theme(axis.text.x = element_text(angle=75, vjust=0.5), legend.position = "top")
g3

# means
raw_m %>% 
  ggplot(aes(x = GEN, y = mean, group = TRAT, color = TRAT))+
  geom_point()+
  geom_line()+
  theme_bw()+
  theme(axis.text.x = element_text(angle=90, vjust=0.5), legend.position = "top")+
  labs( title =  "Interaction Plot")

# boxplots by gen
g4 <- datos %>% 
  ggplot(aes(x = GEN,  y = REND, fill = TRAT))+
  geom_boxplot()+
  theme_bw()+
  labs(title = "Yield distribution by gen")+
  theme(axis.text.x = element_text(angle=90, vjust=0.5))
g4

ggplotly(g4)

# -------------------------------------------------------------------------
# modelling
# -------------------------------------------------------------------------


mod.fb <- lm(REND ~ GEN + TRAT + GEN:TRAT + TRAT:Bloq, data = datos)


# summary
summary(mod.fb)

# ANOVA
ao_1 <- anova(mod.fb)
ao_1

# goodness of fit
gf_1 <- glance(mod.fb)
gf_1

# withing TRAT
mean_comparisons_1 <- mod.fb %>% 
  emmeans(pairwise ~ GEN|TRAT, adjust="tukey") %>% 
  pluck("emmeans") %>% 
  cld(details=TRUE, Letters=letters) # add letter display

mean_comparisons_1$emmeans # adjusted cultivar means


# withing TRAT
mean_comparisons_2 <- mod.fb %>% 
  emmeans(pairwise ~ TRAT|GEN, adjust="tukey") %>% 
  pluck("emmeans") %>% 
  cld(details=TRUE, Letters=letters) # add letter display

mean_comparisons_2$emmeans # adjusted cultivar means


# Plot1 -------------------------------------------------------------------


gf <- ggplot() +
  facet_grid(~TRAT) +
  # black dots representing the raw data
  geom_point(
    data = datos,
    aes(y = REND, x = reorder(GEN, -REND))
  ) +
  # red dots representing the adjusted means
  geom_point(
    data = mean_comparisons_1$emmeans,
    aes(y = emmean, x = GEN),
    color = "red",
    position = position_nudge(x = 0.1)
  ) +
  # red error bars representing the confidence limits of the adjusted means
  geom_errorbar(
    data = mean_comparisons_1$emmeans,
    aes(ymin = lower.CL, ymax = upper.CL, x = GEN),
    color = "red",
    width = 0.1,
    position = position_nudge(x = 0.1)
  ) +
  # red letters 
  geom_text(
    data = mean_comparisons_1$emmeans,
    aes(y = lower.CL, x = GEN, label = .group),
    color = "red",
    angle = 90,
    hjust = 1,
    position = position_nudge(y = -0.5)
  ) + 
  ylim(0, NA) + # force y-axis to start at 0
  ylab("Yield in t/ha") + # label y-axis
  xlab("Genotype") +      # label x-axis
  labs(caption = "
       Black dots represent raw data
       Red dots and error bars represent adjusted mean with 95% confidence limits per genotype-treat level combination
       Per genotype, means followed by a common letter are not significantly different according to the Tukey-test and within") +
  theme_bw() + # clearer plot format 
  theme(axis.text.x = element_text(angle=75, vjust=0.5)) # rotate x-axis label

gf

ggsave(gf , filename = "RCBD_TREATMENT.png", units = "in", dpi = 300,  width = 10,height = 5)


# plot 2 ------------------------------------------------------------------



ggplot() +
  facet_grid(~GEN) +
  # black dots representing the raw data
  geom_point(
    data = datos,
    aes(y = REND, x = reorder(TRAT, -REND))
  ) +
  # red dots representing the adjusted means
  geom_point(
    data = mean_comparisons_2$emmeans,
    aes(y = emmean, x = TRAT),
    color = "red",
    position = position_nudge(x = 0.1)
  ) +
  # red error bars representing the confidence limits of the adjusted means
  geom_errorbar(
    data = mean_comparisons_2$emmeans,
    aes(ymin = lower.CL, ymax = upper.CL, x = TRAT),
    color = "red",
    width = 0.1,
    position = position_nudge(x = 0.1)
  ) +
  # red letters 
  geom_text(
    data = mean_comparisons_2$emmeans,
    aes(y = lower.CL, x = TRAT, label = .group),
    color = "red",
    angle = 90,
    hjust = 1,
    position = position_nudge(y = -0.5)
  ) + 
  ylim(0, NA) + # force y-axis to start at 0
  ylab("Yield in t/ha") + # label y-axis
  xlab("TREATMENT") +      # label x-axis
  labs(caption = "
       Black dots represent raw data
       Red dots and error bars represent adjusted mean with 95% confidence limits per genotype-treat level combination
       Per genotype, means followed by a common letter are not significantly different according to the Tukey-test and within") +
  theme_bw() + # clearer plot format 
  theme(axis.text.x = element_text(angle=90, vjust=0.5)) # rotate x-axis label


# Correlation  ------------------------------------------------------------


gff <- mean_comparisons_1$emmeans %>%
  as.data.frame() %>% 
  pivot_wider(names_from = TRAT, values_from = c(emmean:.group)) %>% 
  ggplot(aes(x = emmean_Irrigación, y = emmean_Sequía))+
  geom_point(size = 4)+
  theme_bw()+
  geom_smooth(method = "lm", se =  F, color = "red", linetype = 2)+
  geom_label_repel(aes(label = GEN))+
  labs(title = "Correlation")+
  stat_cor(size= 5, color = "grey8")+
  stat_regline_equation(
    label.y = 5.3, size = 5
  )
  
gff

ggsave(gff , filename = "Correlation.png", units = "in", dpi = 300,  width = 6,height = 6)



