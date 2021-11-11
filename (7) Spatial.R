

####################################
# Spatial
####################################

library(hrbrthemes)
library(lmerTest)
library(plotly)
library(tidyverse)
library(desplot)
library(broom.mixed)
library(emmeans)
library(multcomp)
library(multcompView) 
library(SpATS)
source("utilities_7_.R")
source("https://raw.githubusercontent.com/AparicioJohan/lme4.plus/master/functions.R")
source("https://raw.githubusercontent.com/AparicioJohan/SpATS.plus/master/All_additional.R")

datos <- read.csv("data/Example_Beans.csv") %>% type.convert()
datos$R <- as.factor(datos$row)
datos$C <- as.factor(datos$col)


head(datos)
str(datos)
summary(datos)

# Grafico usando desplot
desplot(data = datos, flip = FALSE,
        form =  line ~ col + row | rep,            # fill color per linetype, headers per replicate
        text =  line, cex = 0.7, shorten = "no",   # show linetype names per plot
        out1 = row, out1.gpar=list(col="black"), # lines between rows
        out2 = col, out2.gpar=list(col="black"), # lines between columns
        main = "Field layout", show.key = F)     # formatting


# -------------------------------------------------------------------------
# Chequeo de la variable de respuesta
# -------------------------------------------------------------------------

# summary by gen
raw_mean <- datos %>% 
  group_by(line) %>% 
  summarize(mean    = mean(YdHa_clean, na.rm = T),
            std.dev = sd(YdHa_clean, na.rm = T),
            cv      = std.dev/mean,
            n_missing  = sum(is.na(YdHa_clean))
  ) %>% 
  arrange(desc(mean)) %>% 
  print(n=Inf) # print full table
raw_mean

# summary by rep
datos %>% 
  group_by(rep) %>% 
  summarize(mean    = mean(YdHa_clean, na.rm = T),
            std.dev = sd(YdHa_clean, na.rm = T),
            cv      = std.dev/mean )

# summary by col
datos %>% 
  group_by(col) %>% 
  summarize(mean    = mean(YdHa_clean, na.rm = T),
            std.dev = sd(YdHa_clean, na.rm = T),
            cv      = std.dev/mean ) %>% 
  print(n=Inf)

# summary by row
datos %>% 
  group_by(row) %>% 
  summarize(mean    = mean(YdHa_clean, na.rm = T),
            std.dev = sd(YdHa_clean, na.rm = T),
            cv      = std.dev/mean ) %>% 
  print(n=Inf)

# points 
g2 <- datos %>% 
  ggplot(aes(x = line,  y = YdHa_clean))+
  geom_point()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
g2

# points + mean 
plotdata <- datos %>% 
  group_by(line) %>% 
  mutate(mean_YdHa_clean = mean(YdHa_clean, na.rm = T)) %>% # add column with mean YdHa_clean per gen
  ungroup() %>% 
  mutate(gen = fct_reorder(.f = line, .x = mean_YdHa_clean)) # sort factor variable by mean YdHa_clean

g3 <- ggplot(data = plotdata, 
             aes(x = gen)) +
  geom_point(aes(y = YdHa_clean, shape = rep)) +  # scatter plot observed
  geom_point(aes(y = mean_YdHa_clean), color = "cornflowerblue") + # scatter plot mean
  ylim(0, NA) +   # force y-axis to start at 0
  labs(caption = "Blue dots represent arithmetic mean per linetype") +
  theme_bw(base_size = 6) + # clearer plot format 
  theme(axis.text.x = element_text(angle=90, vjust=1)) # rotate x-axis labels
g3

# boxplots by gen
g4 <- datos %>% 
  ggplot(aes(x = line,  y = YdHa_clean))+
  geom_boxplot()+
  theme_bw()+
  labs(title = "YdHa_clean distribution by gen")+
  theme(axis.text.x = element_text(angle=90, vjust=1))
g4

# boxplots by rep
g5 <- datos %>% 
  ggplot(aes(x = rep,  y = YdHa_clean, fill = rep))+
  geom_boxplot()+
  theme_bw()+
  labs(title = "YdHa_clean distribution by rep")
g5

# boxplots by col
g6 <- datos %>% 
  ggplot(aes(x = C,  y = YdHa_clean, fill = C))+
  geom_boxplot()+
  theme_bw()+
  labs(title = "YdHa_clean distribution by COl")
g6


# boxplots by row
g6 <- datos %>% 
  ggplot(aes(x = R,  y = YdHa_clean, fill = R))+
  geom_boxplot()+
  theme_bw()+
  labs(title = "YdHa_clean distribution by gen")
g6


# interactive plots
# ggplotly(g3)

ggplot(data = datos,
       aes(y = YdHa_clean, x = line, color=rep)) +
  geom_point() +  # scatter plot
  ylim(0, NA) +   # force y-axis to start at 0
  theme_classic() + # clearer plot format 
  theme(axis.text.x = element_text(angle=90, vjust=0.5), # rotate x-axis label
        panel.grid.major.x = element_line(), # add vertikal grid lines
        legend.position = "top") # legend on top 


# -------------------------------------------------------------------------
# Modelling ---------------------------------------------------------------
# -------------------------------------------------------------------------


# CRD ---------------------------------------------------------------------
mo_crd <- lmer(formula = YdHa_clean ~ (1|line), data = datos)

gf_crd <- glance(mo_crd) %>% mutate(model = "crd")
gf_crd

lme4.plot(mo_crd, datos, col="col", row="row", gen="line", main = "Yield (CRD)")

# RCBD --------------------------------------------------------------------
mo_rcbd <- lmer(formula = YdHa_clean ~ (1|line) + rep, data = datos)

gf_rcbd <- glance(mo_rcbd) %>% mutate(model = "rcbd")
gf_rcbd

lme4.plot(mo_rcbd, datos, col="col", row="row", gen="line", main = "Yield (RCBD)")

# alpha -------------------------------------------------------------------
mo_alpha <- lmer(formula = YdHa_clean ~ (1|line) + rep + (1|rep:block), data = datos)

gf_alpha <- glance(mo_alpha) %>% mutate(model = "alpha")
gf_alpha

lme4.plot(mo_alpha, datos, col="col", row="row", gen="line", main = "Yield (alpha)")

# row_col -----------------------------------------------------------------
mo_rowcol <- lmer(formula = YdHa_clean ~ (1|line) + rep + (1|rep:C) + (1|rep:R), data = datos)

gf_rowcol <- glance(mo_rowcol) %>% mutate(model = "row_col")
gf_rowcol

lme4.plot(mo_rowcol, datos, col="col", row="row", gen="line", main = "Yield (row-col)")

# row_col + inblock -------------------------------------------------------
mo_rowcolinb <- lmer(formula = YdHa_clean ~ (1|line) + rep + (1|rep:block) + (1|rep:C) + (1|rep:R), data = datos)

gf_rowcolinb <- glance(mo_rowcolinb) %>% mutate(model = "row_col_inb")
gf_rowcolinb

lme4.plot(mo_rowcolinb, datos, col="col", row="row", gen="line", main = "Yield (row-col + inblk)")

# comparison
rbind.data.frame(gf_crd, gf_rcbd, gf_alpha, gf_rowcol, gf_rowcolinb)

# spatial -----------------------------------------------------------------

ncols = length(unique(datos$C))
nrows = length(unique(datos$R))

mo_spatial <- SpATS(response='YdHa_clean',
                      genotype='line', 
                      genotype.as.random=T,
                      fixed=NULL ,
                      spatial = ~ PSANOVA(col,
                                          row, 
                                          nseg = c(ncols, nrows), 
                                          degree = c(3,3),
                                          nest.div=2),
                      random = ~ rep:R + rep:C , 
                      data=datos,
                      control = list(tolerance=1e-03, monitoring=1))
sqrt(mo_spatial$psi[1])
plot(mo_spatial)

# -------------------------------------------------------------------------
# Criterios ---------------------------------------------------------------
# -------------------------------------------------------------------------

# crd
crd_vcov <- data.frame(VarCorr(mo_crd))
crd_g <- crd_vcov[crd_vcov$grp == "line", "vcov"]
crd_e <- crd_vcov[crd_vcov$grp == "Residual", "vcov"]
crd_h2 <- h.cullis(mo_crd, "line")
crd_r2 <- r.lme4(mo_crd)
crd_aic <- AIC(mo_crd)
crd_bic <- BIC(mo_crd)

# rcbd
rcbd_vcov <- data.frame(VarCorr(mo_rcbd))
rcbd_g <- rcbd_vcov[rcbd_vcov$grp == "line", "vcov"]
rcbd_e <- rcbd_vcov[rcbd_vcov$grp == "Residual", "vcov"]
rcbd_h2 <- h.cullis(mo_rcbd, "line")
rcbd_r2 <- r.lme4(mo_rcbd)
rcbd_aic <- AIC(mo_rcbd)
rcbd_bic <- BIC(mo_rcbd)

# alpha
alpha_vcov <- data.frame(VarCorr(mo_alpha))
alpha_g <- alpha_vcov[alpha_vcov$grp == "line", "vcov"]
alpha_e <- alpha_vcov[alpha_vcov$grp == "Residual", "vcov"]
alpha_h2 <- h.cullis(mo_alpha, "line")
alpha_r2 <- r.lme4(mo_alpha)
alpha_aic <- AIC(mo_alpha)
alpha_bic <- BIC(mo_alpha)

# row_col 
rc_vcov <- data.frame(VarCorr(mo_rowcol))
rc_g <- rc_vcov[rc_vcov$grp == "line", "vcov"]
rc_e <- rc_vcov[rc_vcov$grp == "Residual", "vcov"]
rc_h2 <- h.cullis(mo_rowcol, "line")
rc_r2 <- r.lme4(mo_rowcol)
rc_aic <- AIC(mo_rowcol)
rc_bic <- BIC(mo_rowcol)

# row_col + alpha
rc_alpha_vcov <- data.frame(VarCorr(mo_rowcolinb))
rc_alpha_g <- rc_alpha_vcov[rc_alpha_vcov$grp == "line", "vcov"]
rc_alpha_e <- rc_alpha_vcov[rc_alpha_vcov$grp == "Residual", "vcov"]
rc_alpha_h2 <- h.cullis(mo_rowcolinb, "line")
rc_alpha_r2 <- r.lme4(mo_rowcolinb)
rc_alpha_aic <- AIC(mo_rowcolinb)
rc_alpha_bic <- BIC(mo_rowcolinb)

# spats
spat_g <- mo_spatial$var.comp["line"]
spat_e <- mo_spatial$psi[1]
spat_h2 <- getHeritability(mo_spatial)
spat_r2 <- R.square(mo_spatial)
spat_aic <- AIC(mo_spatial)
spat_bic <- BIC(mo_spatial)

stats_g <- c(crd_g, rcbd_g, alpha_g,rc_g, rc_alpha_g, spat_g)
stats_e <- c(crd_e, rcbd_e, alpha_e, rc_e, rc_alpha_e, spat_e)
stats_h2 <- c(crd_h2, rcbd_h2, alpha_h2, rc_alpha_h2, rc_alpha_h2, spat_h2)
stats_r2 <- c(crd_r2, rcbd_r2, alpha_r2, rc_r2, rc_alpha_r2, spat_r2)
stats_aic <- c(crd_aic, rcbd_aic, alpha_aic, rc_aic, rc_alpha_aic, spat_aic)
stats_bic <- c(crd_bic, rcbd_bic, alpha_bic, rc_bic, rc_alpha_bic, spat_bic)
stats_model <- c("crd", "rcbd", "alpha", "row_col" ,"row_col_alpha", "spats")


stack <- data.frame(stats_model, 
                    stats_g, stats_e,
                    stats_h2, stats_r2,
                    stats_aic, stats_bic)

v <- c("crd","rcbd", "alpha", "row_col", "row_col_alpha" ,"spats")
o1 <- stack %>%
  gather(key = "statistic", value = "value", -1) %>%
  ggplot(aes(x = stats_model , y = value, group = statistic))+
  geom_point()+
  geom_line()+
  facet_wrap(~ statistic, scales = "free")+
  scale_x_discrete(limits=v)+
  theme_bw()+
  labs(x = " ",
       y  = " ",
       title = paste0("Comparison"))+
  theme(axis.text.x = element_text(hjust = 1, angle = 75))
o1


# -------------------------------------------------------------------------
# BLUPs -------------------------------------------------------------------
# -------------------------------------------------------------------------

ng <- length(mo_spatial$terms$geno$geno_names)

gen = "line"
g.SpATS  <- mo_spatial$coeff[1:ng] %>% data.frame(BLUP = .) %>% rownames_to_column(gen) %>% mutate( type = "spats")
g.CRD    <- ranef(mo_crd)[[gen]] %>% as.data.frame() %>%  rownames_to_column(gen) %>% mutate(BLUP = `(Intercept)` , type = "crd") %>% dplyr::select(-`(Intercept)`)
g.RCBD   <- ranef(mo_rcbd)[[gen]] %>% as.data.frame() %>%  rownames_to_column(gen) %>% mutate(BLUP = `(Intercept)` , type = "rcbd") %>% dplyr::select(-`(Intercept)`)
g.alpha  <- ranef(mo_alpha)[[gen]] %>% as.data.frame() %>%  rownames_to_column(gen) %>% mutate(BLUP = `(Intercept)` , type = "alpha") %>% dplyr::select(-`(Intercept)`)
g.rowcol  <- ranef(mo_rowcol)[[gen]] %>% as.data.frame() %>%  rownames_to_column(gen) %>% mutate(BLUP = `(Intercept)` , type = "row_col") %>% dplyr::select(-`(Intercept)`)
g.rc_alpha  <- ranef(mo_rowcolinb)[[gen]] %>% as.data.frame() %>%  rownames_to_column(gen) %>% mutate(BLUP = `(Intercept)` , type = "row_col_alpha") %>% dplyr::select(-`(Intercept)`)

BLUPs <- rbind(g.CRD, g.RCBD, g.alpha, g.rowcol, g.rc_alpha, g.SpATS) %>% 
          spread(key = "type", value = "BLUP")

# Genetic Gain
names(stats_h2) <- v
gg_BLUPs <- Genetic_Gain(data = BLUPs, heritability = stats_h2, genotype = "line", exp = v, percentage= 0.15 )
gg_BLUPs %>% 
  ggplot(aes(x = Model, y = S, label = paste( round(S,2), "\n (kg/ha)" ),  fill = Model))+
  geom_bar(stat = "identity", color = "black")+
  scale_x_discrete(limits=v)+
  geom_text(nudge_y = 20)+
  scale_fill_brewer(palette="Set3")+
  theme_ipsum()+
  labs(title = "Ganancia Genética")+
  theme(axis.text.x = element_text(angle=75, hjust=1))


