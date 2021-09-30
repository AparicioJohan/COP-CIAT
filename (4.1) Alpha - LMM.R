

####################################
# Alpha-lattice
####################################

library(lmerTest)
library(plotly)
library(tidyverse)
library(desplot)
library(broom.mixed)
library(emmeans)
library(multcomp)
library(multcompView) # mean comparisons
library(ggpubr)

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
raw_m <- datos %>% 
  group_by(gen) %>% 
  summarize(mean    = mean(yield),
            std.dev = sd(yield),
            cv      = std.dev/mean, 
            n = n() ,
            n_mis = sum(is.na(yield))) %>% 
  print(n = Inf)
raw_m

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


# gen as fixed 
mod.fg <- lmer(yield ~ gen + rep + (1|rep:inc.block),
             data = datos)

# gen as random 
mod.rg <- lmer(yield ~ (1|gen) + rep + (1 | rep:inc.block),
               data = datos)

# summary
summary(mod.fg, ddf = "Kenward-Roger")
summary(mod.rg, ddf = "Kenward-Roger")

# Anova (fixed effects)
mod.fg %>% anova(ddf="Kenward-Roger")
mod.rg %>% anova(ddf="Kenward-Roger")

# Ranova (random effects)
mod.fg %>% ranova()
mod.rg %>% ranova()

# Variance components
as.data.frame(VarCorr(mod.fg))
as.data.frame(VarCorr(mod.rg))

# Coefficients
coef(summary(mod.fg,ddf = "Kenward-Roger"))
coef(summary(mod.rg,ddf = "Kenward-Roger"))

# BLUEs -------------------------------------------------------------------

BLUEs <- emmeans::emmeans(mod.fg, ~ gen) %>%
           as.data.frame() %>%
           transmute(gen ,  BLUE = emmean, std.error_BLUE = SE ) 

# BLUPs -------------------------------------------------------------------

mu_manual <-  fixef(mod.rg)[1]   +   sum(fixef(mod.rg)[2:3])/3   # mod.rg@frame[[1]] %>% mean(., na.rm = T) 

BLUPs <- augment(ranef(mod.rg)) %>%
          filter(grp == "gen") %>% 
          transmute( gen = level, BLUP =  mu_manual + estimate , std.error_BLUP = std.error)

# -------------------------------------------------------------------------

pvals <- merge(BLUEs, BLUPs, by = "gen")
pvals <- merge(pvals, raw_m, by = "gen")


# -------------------------------------------------------------------------
# Heritability ------------------------------------------------------------
# -------------------------------------------------------------------------

# classic 
vcomps <- as.data.frame(VarCorr(mod.rg))
vc.g <- vcomps[vcomps$grp=="gen","vcov"]
vc.e <- vcomps[vcomps$grp=="Residual","vcov"]
nreps = 3
hc <- vc.g / ( vc.g + vc.e/nreps)
hc

# cullis 
aveped <- mean(attr(ranef(mod.rg,drop=T)$gen,"postVar"))
hcullis <- 1-aveped/vc.g 
hcullis

# regression 
pvals %>% 
  ggplot(aes(x = BLUE, y = BLUP))+
  geom_smooth(se = F, color = "red",  size = 0.8, method = "lm", formula = y ~ x )+
  geom_abline(slope = 1, intercept = 0, color = "black", size = 0.8, linetype = 2)+
  geom_point(size = 3, alpha = 0.5)+
  stat_regline_equation()+
  theme_bw()+
  coord_fixed()
  
mo_h2 <- lm(BLUP ~  BLUE, data = pvals)
hregress <- coef(mo_h2)[2]

# piepho
vc.g <- vcomps[vcomps$grp=="gen","vcov"]

vdBLUE.avg <- mod.fg %>% 
                 emmeans(pairwise ~ gen) %>% 
                 pluck("contrasts") %>% 
                 as_tibble %>% 
                 mutate(Var=SE^2) %>% 
                 pull(Var) %>% 
                 mean # mean variance of a difference between genotypes

H2.p <- vc.g/(vc.g + vdBLUE.avg/2)
H2.p 

hc
hcullis
hregress
H2.p
