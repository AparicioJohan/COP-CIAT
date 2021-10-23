

####################################
# Row-Col
####################################

library(lmerTest)
library(plotly)
library(tidyverse)
library(desplot)
library(broom)
library(emmeans)
library(multcomp)
library(multcompView) # mean comparisons

datos <- read.csv("data/row_col.csv")
head(datos)
str(datos)
summary(datos)

# La variedad la rep y el bloque incompleto deben ser definidos como factors
datos$gen  <- as.factor(datos$gen)
datos$rep <- as.factor(datos$rep)

# Grafico usando desplot
desplot(data = datos, flip = FALSE,
        form = gen ~ col + row | rep,            # fill color per genotype, headers per replicate
        text = gen, cex = 0.7, shorten = "no",   # show genotype names per plot
        out1 = row, out1.gpar=list(col="black"), # lines between rows
        out2 = col, out2.gpar=list(col="black"), # lines between columns
        main = "Field layout", show.key = F)     # formatting

g0 <- ggdesplot(data = datos, flip = FALSE,
        form = gen ~ col + row | rep,            # fill color per genotype, headers per replicate
        text = gen, cex = 0.7, shorten = "no",   # show genotype names per plot
        out1 = row, out1.gpar=list(col="black"), # lines between rows
        out2 = col, out2.gpar=list(col="black"), # lines between columns
        main = "Field layout", show.key = F)     # formatting
g0

ggplotly(g0)

# Grafico usando ggplot2
g1 <- datos %>% 
  ggplot(aes(x = col, y = row, fill = gen )) +
  geom_tile(color = "black")+
  geom_text(aes(label = gen))+
  theme_bw()+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  facet_wrap(~ rep, scales = "free_x")
g1

# Grafico interactivo usando plotly
ggplotly(g0)
ggplotly(g1)

# -------------------------------------------------------------------------
# Chequeo de la variable de respuesta
# -------------------------------------------------------------------------

# summary by gen
datos %>% 
  group_by(gen) %>% 
  summarize(mean    = mean(yield, na.rm = T),
            std.dev = sd(yield, na.rm = T),
            cv      = std.dev/mean,
            n_missing  = sum(is.na(yield))
  ) %>% 
  arrange(desc(mean)) %>% 
  print(n=Inf) # print full table

# summary by rep
datos %>% 
  group_by(rep) %>% 
  summarize(mean    = mean(yield, na.rm = T),
            std.dev = sd(yield, na.rm = T),
            cv      = std.dev/mean )

# summary by col
datos %>% 
  group_by(col) %>% 
  summarize(mean    = mean(yield, na.rm = T),
            std.dev = sd(yield, na.rm = T),
            cv      = std.dev/mean )

# summary by row
datos %>% 
  group_by(row) %>% 
  summarize(mean    = mean(yield, na.rm = T),
            std.dev = sd(yield, na.rm = T),
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
  mutate(mean_yield = mean(yield, na.rm = T)) %>% # add column with mean yield per gen
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

# boxplots by rep
g5 <- datos %>% 
  ggplot(aes(x = rep,  y = yield, fill = rep))+
  geom_boxplot()+
  theme_bw()+
  labs(title = "Yield distribution by rep")
g5


# interactive plots
ggplotly(g3)

ggplot(data = datos,
       aes(y = yield, x = gen, color=rep)) +
  geom_point() +  # scatter plot
  ylim(0, NA) +   # force y-axis to start at 0
  theme_classic() + # clearer plot format 
  theme(axis.text.x = element_text(angle=90, vjust=0.5), # rotate x-axis label
        panel.grid.major.x = element_line(), # add vertikal grid lines
        legend.position = "top") # legend on top 


# -------------------------------------------------------------------------
# modelling
# -------------------------------------------------------------------------

datos <- datos %>% 
  mutate(row_fct = as.factor(row),
         col_fct = as.factor(col))


# rows and cols fixed (linear model)
mod.fb <- lm(yield ~ gen + rep + rep:row_fct + rep:col_fct,
             data = datos)

mod.fb %>%
  emmeans(pairwise ~ "gen",
          adjust = "tukey") %>%
  pluck("contrasts") %>% # extract diffs
  as_tibble %>% # format to table
  pull("SE") %>% # extract s.e.d. column
  mean() # get arithmetic mean


# rows and cols random (linear mixed model)
mod.rb <- lmer(yield ~ gen + rep + (1|rep:row_fct) + (1|rep:col_fct),
               data = datos)

mod.rb %>%
  emmeans(pairwise ~ "gen",
          adjust = "tukey",
          lmer.df = "kenward-roger") %>%
  pluck("contrasts") %>% # extract diffs
  as_tibble %>% # format to table
  pull("SE") %>% # extract s.e.d. column
  mean() # get arithmetic mean



mod.rb %>% 
  VarCorr() %>% 
  as.data.frame() %>% 
  dplyr::select(grp, vcov)


# ANOVA
mod.rb %>% anova(ddf="Kenward-Roger")

# RANOVA
ranova(mod.rb)


# Mean comparison
mean_comparisons <- mod.rb %>%
  emmeans(pairwise ~ "gen",
          adjust = "tukey",
          lmer.df = "kenward-roger") %>%
  pluck("emmeans") %>%
  cld(details = TRUE, Letters = letters) # add letter display

mean_comparisons$emmeans # adjusted genotype means


# -------------------------------------------------------------------------

# sort gen factors according to emmean
# for adjusted means
mean_comparisons$emmeans <- mean_comparisons$emmeans %>% 
  mutate(gen = fct_reorder(.f = gen, .x = emmean))

# for raw data
datos <- datos %>% 
  mutate(gen = fct_relevel(.f = gen, 
                           as.character(mean_comparisons$emmeans$gen)))

ggplot() +
  # black dots representing the raw data
  geom_point(
    data = datos,
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
  theme(axis.text.x = element_text(angle=90, vjust=0.5)) # rotate x-axis label


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




