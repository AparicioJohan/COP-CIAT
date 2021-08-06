

####################################
# Randomized complete block design
####################################

library(tidyverse)
library(desplot)
library(plotly)
library(broom)
library(emmeans)
library(multcomp)
library(multcompView) # mean comparisons

datos <- read.csv("data/RCBD.csv")
head(datos)
str(datos)
summary(datos)

# La variedad y el bloque deben ser definidos como factor
datos$cultivar <- as.factor(datos$cultivar)
datos$block    <- as.factor(datos$block)

table(datos$block, datos$cultivar)

# Grafico usando desplot
desplot(data = datos, flip = FALSE,
        form = cultivar ~ col + row,              # fill color per genotype
        out1 = block,
        text = cultivar, cex = 1, shorten = "no", # show genotype names per plot
        main = "Field layout", show.key = TRUE)     

# Grafico usando ggplot2
g1 <- datos %>% 
  ggplot(aes(x = col, y = row, fill = cultivar )) +
  geom_tile(color = "black")+
  geom_text(aes(label = cultivar))+
  theme_bw()+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  facet_wrap(~ block, scales = "free_x")
g1

# Grafico interactivo usando plotly
ggplotly(g1)


# -------------------------------------------------------------------------
# Chequeo de la variable de respuesta
# -------------------------------------------------------------------------

# summary by cultivar
datos %>% 
  group_by(cultivar) %>% 
  summarize(mean    = mean(yield),
            std.dev = sd(yield),
            cv      = std.dev/mean )

# summary by block
datos %>% 
  group_by(block) %>% 
  summarize(mean    = mean(yield),
            std.dev = sd(yield),
            cv      = std.dev/mean )

# points 
g2 <- datos %>% 
  ggplot(aes(x = cultivar,  y = yield))+
  geom_point()+
  theme_bw()
g2

# points + color
g3 <- datos %>% 
  ggplot(aes(x = cultivar,  y = yield, color = block, group= block))+
  geom_point(size = 4)+
  # geom_line()+
  theme_bw()
g3

# boxplots by cultivar
g4 <- datos %>% 
  ggplot(aes(x = cultivar,  y = yield, fill = cultivar))+
  geom_boxplot()+
  theme_bw()+
  labs(title = "Yield distribution by cultivar")
g4

# boxplots by block
g5 <- datos %>% 
  ggplot(aes(x = block,  y = yield, fill = block))+
  geom_boxplot()+
  theme_bw()+
  labs(title = "Yield distribution by block")
g5

# interactive plots
ggplotly(g3)


# -------------------------------------------------------------------------
# modelling
# -------------------------------------------------------------------------


############# model without blocking

model_1 <- lm(formula = yield ~ cultivar, data = datos )

# summary
summary(model_1)

# ANOVA
ao_1 <- anova(model_1)
ao_1

# goodness of fit
gf_1 <- glance(model_1)
gf_1

# Mean comparisons
mean_comparisons_1 <- model_1 %>% 
  emmeans(pairwise ~ "cultivar", adjust="tukey") %>% 
  pluck("emmeans") %>% 
  cld(details=TRUE, Letters=letters) # add letter display

mean_comparisons_1$emmeans # adjusted cultivar means

############# model with block

model_2 <- lm(formula = yield ~ cultivar + block, data = datos )

# summary
summary(model_2)

# ANOVA
ao_2 <- anova(model_2)
ao_2

# goodness of fit
gf_2 <- glance(model_2)
gf_2

# Mean comparisons
mean_comparisons_2 <- model_2 %>% 
  emmeans(pairwise ~ "cultivar", adjust="tukey") %>% 
  pluck("emmeans") %>% 
  cld(details=TRUE, Letters=letters) # add letter display

mean_comparisons_2$emmeans # adjusted cultivar means


# -------------------------------------------------------------------------
# Comparing Models --------------------------------------------------------
# -------------------------------------------------------------------------

anova(model_1,model_2)

ao_1
ao_2

gf_1$model <- "no_block"
gf_2$model <- "si_block"
gf <- rbind(gf_1,gf_2)
gf


# -------------------------------------------------------------------------
# Plots -------------------------------------------------------------------
# -------------------------------------------------------------------------


ggplot() +
  # black dots representing the raw data
  geom_point(
    data = datos,
    aes(y = yield, x = cultivar, color = block)
  ) +
  # red dots representing the adjusted means
  geom_point(
    data = mean_comparisons_2$emmeans,
    aes(y = emmean, x = cultivar),
    color = "red",
    position = position_nudge(x = 0.1)
  ) +
  # red error bars representing the confidence limits of the adjusted means
  geom_errorbar(
    data = mean_comparisons_2$emmeans,
    aes(ymin = lower.CL, ymax = upper.CL, x = cultivar),
    color = "red",
    width = 0.1,
    position = position_nudge(x = 0.1)
  ) +
  # red letters 
  geom_text(
    data = mean_comparisons_2$emmeans,
    aes(y = emmean, x = cultivar, label = .group),
    color = "red",
    position = position_nudge(x = 0.2)
  ) + 
  ylim(0, NA) + # force y-axis to start at 0
  ylab("Yield in t/ha") + # label y-axis
  xlab("Cultivar") +      # label x-axis
  labs(caption = "Black dots represent raw data
       Red dots and error bars represent adjusted mean with 95% confidence limits per cultivar
       Means followed by a common letter are not significantly different according to the Tukey-test") +
  theme_bw() # clearer plot format 


# -------------------------------------------------------------------------
# LSD ---------------------------------------------------------------------
# -------------------------------------------------------------------------


test <- agricolae::LSD.test(model_2, "cultivar", alpha = 0.05)
test$statistics


LSD_0 <- qt(p = 0.975,df = test$statistics$Df) * sqrt( test$statistics$MSerror * ( 1/test$means$r[1] + 1/test$means$r[2] ) )
LSD_0
