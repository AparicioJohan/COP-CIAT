

####################################
# Row-Col - LMM
####################################

library(lmerTest)
library(plotly)
library(tidyverse)
library(desplot)
library(broom)
library(broom.mixed)
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
raw_m <- datos %>% 
  group_by(gen) %>% 
  summarize(mean    = mean(yield, na.rm = T),
            std.dev = sd(yield, na.rm = T),
            n  = n(),
            cv      = std.dev/mean,
            n_missing  = sum(is.na(yield))
  ) %>% 
  arrange(desc(mean)) %>% 
  print(n=Inf) # print full table

raw_m

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

# Genotype as fixed  
mod.fg <- lmer(yield ~ gen + rep + (1|rep:row_fct) + (1|rep:col_fct),
               data = datos)

# Genotype as Random
mod.rg <- lmer(yield ~ (1|gen) + rep + (1|rep:row_fct) + (1|rep:col_fct),
               data = datos)


mod.rg %>% 
  VarCorr() %>% 
  as.data.frame() %>% 
  dplyr::select(grp, vcov)


# ANOVA
mod.rg %>% anova(ddf="Kenward-Roger")

# RANOVA
ranova(mod.rg)

# BLUEs -------------------------------------------------------------------

BLUEs <- emmeans::emmeans(mod.fg, ~ gen) %>%
  as.data.frame() %>%
  transmute(gen ,  BLUE = emmean, std.error_BLUE = SE ) 

# BLUPs -------------------------------------------------------------------

mu_manual <-  fixef(mod.rg)[1]   +   sum(fixef(mod.rg)[2])/2   # mod.rg@frame[[1]] %>% mean(., na.rm = T) 

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
nreps = 2
hc <- vc.g / ( vc.g + vc.e/nreps)
hc

# cullis 
aveped <- mean(attr(ranef(mod.rg,drop=T)$gen,"postVar"))
hcullis <- 1-aveped/vc.g 
hcullis
