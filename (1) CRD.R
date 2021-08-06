

################################
# Completely Randomized design
################################

library(tidyverse)
library(desplot)
library(plotly)
library(broom)
library(emmeans)
library(multcomp)
library(multcompView) # mean comparisons

datos <- read.csv("data/CRD.csv")
head(datos)
str(datos)
summary(datos)

# La variedad debe ser definida como factor
datos$variety <- as.factor(datos$variety)

# Grafico usando desplot
desplot(data = datos, flip = FALSE,
        form = variety ~ col + row,              # fill color per genotype
        text = variety, cex = 1, shorten = "no", # show genotype names per plot
        main = "Field layout", show.key = FALSE)     

# Grafico usando ggplot2
g1 <- datos %>% 
  ggplot(aes(x = col, y = row, fill = variety )) +
  geom_tile(color = "black")+
  geom_text(aes(label = variety))+
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

# summary
datos %>% 
  group_by(variety) %>% 
  summarize(mean    = mean(yield),
            std.dev = sd(yield),
            cv      = std.dev/mean )

# points 
g2 <- datos %>% 
  ggplot(aes(x = variety,  y = yield))+
  geom_point()+
  theme_bw()
g2

# boxplots
g3 <- datos %>% 
  ggplot(aes(x = variety,  y = yield, fill = variety))+
  geom_boxplot()+
  theme_bw()+
  labs(title = "Yield distribution by variety")
g3

# interactive plots
ggplotly(g3)


# -------------------------------------------------------------------------
# Modelling
# -------------------------------------------------------------------------

# model
model <- lm(formula = yield ~ variety, data = datos )

# summary
summary(model)

# ANOVA
anova(model)

# goodness of fit
glance(model)

# Mean comparisons
mean_comparisons <- model %>% 
  emmeans(pairwise ~ "variety", adjust="tukey") %>% 
  pluck("emmeans") %>% 
  cld(details=TRUE, Letters=letters) # add letter display

mean_comparisons$emmeans # adjusted variety means

# Matrix comparison
objtmeans <- emmeans(model, "variety")
pwpm(objtmeans)
pwpp(objtmeans)+theme_bw()
plot(objtmeans)

# -------------------------------------------------------------------------
# Plot --------------------------------------------------------------------
# -------------------------------------------------------------------------


# plotting Results
ggplot() +
  # black dots representing the raw data
  geom_point(
    data = datos,
    aes(y = yield, x = variety)
  ) +
  # red dots representing the adjusted means
  geom_point(
    data = mean_comparisons$emmeans,
    aes(y = emmean, x = variety),
    color = "red",
    position = position_nudge(x = 0.1)
  ) +
  # red error bars representing the confidence limits of the adjusted means
  geom_errorbar(
    data = mean_comparisons$emmeans,
    aes(ymin = lower.CL, ymax = upper.CL, x = variety),
    color = "red",
    width = 0.1,
    position = position_nudge(x = 0.1)
  ) +
  # red letters 
  geom_text(
    data = mean_comparisons$emmeans,
    aes(y = emmean, x = variety, label = .group),
    color = "red",
    position = position_nudge(x = 0.2)
  ) + 
  ylim(0, NA) + # force y-axis to start at 0
  ylab("Yield in t/ha") + # label y-axis
  xlab("Variety") +      # label x-axis
  labs(caption = "Black dots represent raw data
       Red dots and error bars represent adjusted mean with 95% confidence limits per variety
       Means followed by a common letter are not significantly different according to the Tukey-test") +
  theme_bw() # clearer plot format 



# -------------------------------------------------------------------------
# LSD ---------------------------------------------------------------------
# -------------------------------------------------------------------------


test <- agricolae::LSD.test(model, "variety", alpha = 0.05)
test$statistics


LSD_0 <- qt(p = 0.975,df = test$statistics$Df) * sqrt( test$statistics$MSerror * ( 1/test$means$r[1] + 1/test$means$r[2] ) )
LSD_0



# -------------------------------------------------------------------------
# Residuals ---------------------------------------------------------------
# -------------------------------------------------------------------------


res <- augment(model)

res %>% 
  ggplot(aes(x = .fitted, y = .resid,  color= variety))+
  geom_point(size = 3)+
  theme_bw()

res %>% 
ggpubr::ggqqplot(x=".resid",
                 ggtheme = theme_bw(), 
                 ylab = "Sample Quantile",
                 xlab = "Theoretical Quantile") +
  labs(title = "QQ-Plot")


res$.resid %>% hist()
