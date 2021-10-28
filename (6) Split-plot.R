
library(lmerTest)
library(plotly)
library(tidyverse)
library(desplot)
library(broom)
library(emmeans)
library(multcomp)
library(multcompView) # mean comparisons
library(ggrepel)      # labels ggplot
library(ggpubr)       # stats with ggplot


datos <- read.csv("data/splitplot.csv") %>% type.convert()
head(datos)
str(datos)


desplot(data = datos,
        form = rep ~ col + row | rep, # fill color per rep, headers per rep
        text = G, cex = 1, shorten = "no", # show genotype names per plot
        col = N, # color of genotype names for each N-level
        out1 = mainplot, out1.gpar = list(col = "black"), # lines between mainplots
        out2 = row, out2.gpar = list(col = "darkgrey"), # lines between rows
        main = "Field layout", show.key = TRUE, key.cex = 0.7) # formatting


datos %>% 
  ggplot(aes(x = col, y = row, fill = mainplot ))+
  geom_tile(color = "black")+
  geom_text(aes(label = G))+
  theme_bw()+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  facet_wrap(~ rep, scales = "free_x")

# -------------------------------------------------------------------------
# Chequeo de la variable de respuesta
# -------------------------------------------------------------------------

# summary by GEN
raw_m <- datos %>% 
  group_by(G, N) %>% 
  summarize(mean    = mean(yield, na.rm = T),
            std.dev = sd(yield, na.rm = T),
            cv      = std.dev/mean ) %>% 
  print(n=Inf)
raw_m

# means
raw_m %>% 
  ggplot(aes(x = N, y = mean, group = G, color = G))+
  geom_point()+
  geom_line()+
  theme_bw()+
  theme(axis.text.x = element_text(angle=90, vjust=0.5), legend.position = "top")+
  labs( title =  "Interaction Plot")


datos %>% 
  group_by(G) %>% 
  summarize(mean    = mean(yield),
            std.dev = sd(yield)) %>% 
  arrange(desc(mean)) %>% # sort
  print(n=Inf) # print full table


datos %>% 
  group_by(N) %>% 
  summarize(mean    = mean(yield),
            std.dev = sd(yield))  %>% 
  arrange(desc(mean)) %>% # sort
  print(n=Inf) # print full table


ggplot(data = datos, 
       aes(y = yield, 
           x = N,
           color = N)) +
  facet_grid(~G) + # facette per N level
  geom_point() +  # scatter plot observed
  theme_bw() + # clearer plot format 
  theme(legend.position = "top") # legend on top


# Modelling ---------------------------------------------------------------


mod <- lmer(yield ~ G + N + G:N + 
              rep + (1|rep:mainplot), 
            data=datos)

mod %>% 
  VarCorr() %>% 
  as.data.frame() %>% 
  dplyr::select(grp, vcov)


mod %>% anova(ddf="Kenward-Roger")


all_mean_comparisons <- mod %>%
  emmeans(pairwise ~ N:G,
          adjust = "tukey",
          lmer.df = "kenward-roger") %>%
  pluck("emmeans") %>%
  cld(details = TRUE, Letters = letters) # add letter display

all_mean_comparisons$emmeans # adjusted means



withinG_mean_comparisons <- mod %>%
  emmeans(pairwise ~ N | G,
          adjust = "tukey",
          lmer.df = "kenward-roger") %>%
  pluck("emmeans") %>%
  cld(details = TRUE, Letters = letters) # add letter display

withinG_mean_comparisons$emmeans # adjusted means


formatted_emmeans <- withinG_mean_comparisons$emmeans %>% 
  as_tibble()

write.csv(formatted_emmeans, file = "emmean_NG.csv", quote = F , row.names = F)

g1 <- ggplot() +
  facet_grid(~G) +
  # black dots representing the raw data
  geom_point(
    data = datos,
    aes(y = yield, x = N)
  ) +
  # red dots representing the adjusted means
  geom_point(
    data = formatted_emmeans,
    aes(y = emmean, x = N),
    color = "red",
    position = position_nudge(x = 0.1)
  ) +
  # red error bars representing the confidence limits of the adjusted means
  geom_errorbar(
    data = formatted_emmeans,
    aes(ymin = lower.CL, ymax = upper.CL, x = N),
    color = "red",
    width = 0.1,
    position = position_nudge(x = 0.1)
  ) +
  # red letters 
  geom_text(
    data = formatted_emmeans,
    aes(y = lower.CL, x = N, label = .group),
    color = "red",
    angle = 90,
    hjust = 1,
    position = position_nudge(y = - 2)
  ) + 
  ylim(0, NA) + # force y-axis to start at 0
  ylab("Yield in t/ha") + # label y-axis
  xlab("Nitrogen Level") +      # label x-axis
  labs(caption = "The four facettes represent genotypes A, B, C and D
       Black dots represent raw data
       Red dots and error bars represent adjusted mean with 95% confidence limits per genotype-nitrogen level combination
       Per genotype, means followed by a common letter are not significantly different according to the Tukey-test and within") +
  theme_bw() + # clearer plot format 
  theme(axis.text.x = element_text(angle=90, vjust=0.5)) # rotate x-axis label

