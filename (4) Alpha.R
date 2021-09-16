

####################################
# Alpha-lattice
####################################

library(lmerTest)
library(plotly)
library(tidyverse)
library(desplot)
library(broom)
library(emmeans)
library(multcomp)
library(multcompView) # mean comparisons

datos <- read.csv("data/Alpha.csv")
head(datos)
str(datos)
summary(datos)

# La variedad la rep y el bloque incompleto deben ser definidos como factors
datos$gen   <- as.factor(datos$gen)
datos$rep <- as.factor(datos$rep)
datos$inc.block <- as.factor(datos$inc.block)

# Grafico usando desplot
desplot(data = datos, flip = FALSE,
        form = gen ~ col + row | rep,          # fill color per genotype, headers per replicate
        text = gen, cex = 0.7, shorten = "no", # show genotype names per plot
        out1 = rep,                            # lines between complete blocks/replicates
        out2 = inc.block,                      # lines between incomplete blocks
        main = "Field layout", show.key = F)   # formatting

# Grafico usando ggplot2
g1 <- datos %>% 
  ggplot(aes(x = col, y = row, fill = inc.block )) +
  geom_tile(color = "black")+
  geom_text(aes(label = gen))+
  theme_bw()+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  facet_wrap(~ rep, scales = "free_x")
g1

# Grafico interactivo usando plotly
ggplotly(g1)

# -------------------------------------------------------------------------
# Chequeo de la variable de respuesta
# -------------------------------------------------------------------------

# summary by gen
datos %>% 
  group_by(gen) %>% 
  summarize(mean    = mean(yield),
            std.dev = sd(yield),
            cv      = std.dev/mean, 
            n = n() ,
            n_mis = sum(is.na(yield))) %>% 
  print(n = Inf)

# summary by rep
datos %>% 
  group_by(rep) %>% 
  summarize(mean    = mean(yield),
            std.dev = sd(yield),
            cv      = std.dev/mean )

# summary by inblock
datos %>% 
  group_by(rep, inc.block) %>% 
  summarize(mean    = mean(yield),
            std.dev = sd(yield),
            cv      = std.dev/mean )

# points 
g2 <- datos %>% 
  ggplot(aes(x = gen,  y = yield))+
  geom_point()+
  theme_bw()
g2

# points + mean 
plotdata <- datos %>% 
  group_by(gen) %>% 
  mutate(mean_yield = mean(yield)) %>% # add column with mean yield per gen
  ungroup() %>% 
  mutate(gen = fct_reorder(.f = gen, .x = mean_yield)) # sort factor variable by mean yield

g3 <- ggplot(data = plotdata, 
       aes(x = gen)) +
  geom_point(aes(y = yield, shape = rep)) +  # scatter plot observed
  geom_point(aes(y = mean_yield), color = "cornflowerblue") + # scatter plot mean
  ylim(0, NA) +   # force y-axis to start at 0
  labs(caption = "Blue dots represent arithmetic mean per genotype") +
  theme_bw() + # clearer plot format 
  theme(axis.text.x = element_text(angle=90, vjust=0.5)) # rotate x-axis labels
g3

# boxplots by gen
g4 <- datos %>% 
  ggplot(aes(x = gen,  y = yield, fill = gen))+
  geom_boxplot()+
  theme_bw()+
  labs(title = "Yield distribution by gen")
g4

# boxplots by gen
g5 <- datos %>% 
  ggplot(aes(x = rep,  y = yield, fill = rep))+
  geom_boxplot()+
  theme_bw()+
  labs(title = "Yield distribution by gen")
g5


# interactive plots
ggplotly(g3)


# -------------------------------------------------------------------------
# modelling
# -------------------------------------------------------------------------


# blocks as fixed (linear model)
mod.fb <- lm(yield ~ gen + rep + rep:inc.block,
             data = datos)

mod.fb %>%
  emmeans(pairwise ~ "gen",
          adjust = "tukey") %>%
  pluck("contrasts") %>% # extract diffs
  as_tibble %>% # format to table
  pull("SE") %>% # extract s.e.d. column
  mean() # get arithmetic mean


# blocks as random (linear mixed model)
mod.rb <- lmer(yield ~ gen + rep + (1 | rep:inc.block),
               data = datos)

mod.rb %>%
  emmeans(pairwise ~ "gen",
          adjust = "tukey",
          lmer.df = "kenward-roger") %>%
  pluck("contrasts") %>% # extract diffs
  as_tibble %>% # format to table
  pull("SE") %>% # extract s.e.d. column
  mean() # get arithmetic mean


mod.rb %>% anova(ddf="Kenward-Roger")



mean_comparisons <- mod.rb %>%
  emmeans(pairwise ~ "gen",
          adjust = "tukey",
          lmer.df = "kenward-roger") %>%
  pluck("emmeans") %>%
  cld(details = TRUE, Letters = letters) # add letter display

mean_comparisons$emmeans # adjusted genotype means



# -------------------------------------------------------------------------



gf <- ggplot() +
  # black dots representing the raw data
  geom_point(
    data = plotdata,
    aes(y = yield, x = gen)
  ) +
  # red dots representing the adjusted means
  geom_point(
    data = mean_comparisons$emmeans,
    aes(y = emmean, x = gen),
    color = "red",
    position = position_nudge(x = 0.1)
  ) +
  # red error bars representing the confidence limits of the adjusted means
  geom_errorbar(
    data = mean_comparisons$emmeans,
    aes(ymin = lower.CL, ymax = upper.CL, x = gen),
    color = "red",
    width = 0.1,
    position = position_nudge(x = 0.1)
  ) +
  # red letters 
  geom_text(
    data = mean_comparisons$emmeans,
    aes(y = lower.CL, x = gen, label = .group),
    color = "red",
    angle = 90,
    hjust = 1,
    position = position_nudge(y = - 0.1)
  ) + 
  ylim(0, NA) + # force y-axis to start at 0
  ylab("Yield in t/ha") + # label y-axis
  xlab("Genotype") +      # label x-axis
  labs(caption = "Black dots represent raw data
       Red dots and error bars represent adjusted mean with 95% confidence limits per genotype
       Means followed by a common letter are not significantly different according to the Tukey-test") +
  theme_bw() + # clearer plot format 
  theme(axis.text.x = element_text(angle=90, vjust=0.5)) # rotate x-axis

gf

ggsave(paste0("alpha_lattice.png"), plot = gf, units = "in", dpi = 300, width = 10, height = 6)


# LSD mixed model ---------------------------------------------------------

t1<-(table(mod.rb@frame["gen"]))
ns<-cbind(median(t1),median(t1))

MS<-sigma(mod.rb)**2
DF<-(anova(mod.rb))["gen","DenDF"]
NREP1<-ns[,1]
NREP2<-ns[,2]
LSD<-qt(1-((1-0.95)/2),DF)*sqrt(MS*(1/NREP1+1/NREP2))
LSD

LSD2 <- agricolae::LSD.test(mod.fb, "gen", alpha = 0.05)
LSD2$statistics




