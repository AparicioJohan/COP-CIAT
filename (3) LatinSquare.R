

####################################
# Latin square design
####################################

library(plotly)
library(tidyverse)
library(desplot)
library(broom)
library(emmeans)
library(multcomp)
library(multcompView) # mean comparisons

datos <- read.csv("data/latinsquare.csv")
head(datos)
str(datos)
summary(datos)

# La variedad las filas y las columnas deben ser definidas como factor
datos$gen   <- as.factor(datos$gen)
datos$col_f <- as.factor(datos$col)
datos$row_f <- as.factor(datos$row)

# Grafico usando desplot
desplot(data = datos, flip = FALSE,
        form = gen ~ col + row, 
        out1 = col, out1.gpar=list(col="black", lwd=3),
        out2 = col, out2.gpar=list(col="black", lwd=3),
        text = gen, cex = 1, shorten = "no",
        main = "Field layout", 
        show.key = F)

# Grafico usando ggdesplot
g0 <- ggdesplot(data = datos, flip = FALSE,
        form = gen ~ col + row, 
        out1 = col, out1.gpar=list(col="black", lwd=3),
        out2 = col, out2.gpar=list(col="black", lwd=3),
        text = gen, cex = 1, shorten = "no",
        main = "Field layout", 
        show.key = T)

ggplotly(g0)

# Grafico usando ggplot2
g1 <- datos %>% 
  ggplot(aes(x = col, y = row, fill = gen )) +
  geom_tile(color = "black")+
  geom_text(aes(label = gen))+
  theme_bw()+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
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
            cv      = std.dev/mean )

# summary by row
datos %>% 
  group_by(row) %>% 
  summarize(mean    = mean(yield),
            std.dev = sd(yield),
            cv      = std.dev/mean )

# summary by col
datos %>% 
  group_by(col) %>% 
  summarize(mean    = mean(yield),
            std.dev = sd(yield),
            cv      = std.dev/mean )

# points 
g2 <- datos %>% 
  ggplot(aes(x = gen,  y = yield))+
  geom_point()+
  theme_bw()
g2

# points + color
g3 <- datos %>% 
  ggplot(aes(x = gen,  y = yield, color = row_f, shape = col_f))+
  geom_point(size = 4)+
  theme_bw()
g3

# boxplots by gen
g4 <- datos %>% 
  ggplot(aes(x = gen,  y = yield, fill = gen))+
  geom_boxplot()+
  theme_bw()+
  labs(title = "Yield distribution by gen")
g4


# interactive plots
ggplotly(g3)


# -------------------------------------------------------------------------
# modelling
# -------------------------------------------------------------------------


############# model without col-row

model_1 <- lm(formula = yield ~ gen, data = datos )

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
  emmeans(pairwise ~ "gen", adjust="tukey") %>% 
  pluck("emmeans") %>% 
  cld(details=TRUE, Letters=letters) # add letter display

mean_comparisons_1$emmeans # adjusted gen means

############# model with col-row

model_2 <- lm(formula = yield ~ gen + col_f + row_f, data = datos )

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
  emmeans(pairwise ~ "gen", adjust="tukey") %>% 
  pluck("emmeans") %>% 
  cld(details=TRUE, Letters=letters) # add letter display

mean_comparisons_2$emmeans # adjusted gen means


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


model_1 %>%
  emmeans(pairwise ~ "gen",
          adjust = "tukey") %>%
  pluck("contrasts") %>% # extract diffs
  as_tibble %>% # format to table
  pull("SE") %>% # extract s.e.d. column
  mean() # get arithmetic mean

model_2 %>%
  emmeans(pairwise ~ "gen",
          adjust = "tukey") %>%
  pluck("contrasts") %>% # extract diffs
  as_tibble %>% # format to table
  pull("SE") %>% # extract s.e.d. column
  mean() # get arithmetic mean

# Plot --------------------------------------------------------------------



ggplot() +
  # black dots representing the raw data
  geom_point(
    data = datos,
    aes(y = yield, x = gen)
  ) +
  # red dots representing the adjusted means
  geom_point(
    data = mean_comparisons_2$emmeans,
    aes(y = emmean, x = gen),
    color = "red",
    position = position_nudge(x = 0.1)
  ) +
  # red error bars representing the confidence limits of the adjusted means
  geom_errorbar(
    data = mean_comparisons_2$emmeans,
    aes(ymin = lower.CL, ymax = upper.CL, x = gen),
    color = "red",
    width = 0.1,
    position = position_nudge(x = 0.1)
  ) +
  # red letters 
  geom_text(
    data = mean_comparisons_2$emmeans,
    aes(y = emmean, x = gen, label = .group),
    color = "red",
    position = position_nudge(x = 0.2)
  ) + 
  ylim(0, NA) + # force y-axis to start at 0
  ylab("Yield") + # label y-axis
  xlab("Cucumber genotype") + # label x-axis
  labs(caption = "Black dots represent raw data
       Red dots and error bars represent adjusted mean with 95% confidence limits per genotype
       Means followed by a common letter are not significantly different according to the Tukey-test") +
  theme_bw() # clearer plot format 



