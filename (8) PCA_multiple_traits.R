
library(hrbrthemes)
library(FactoMineR)
library(factoextra)
library(explor)
library(tidyverse)
library(desplot)
library(SpATS)
source("https://raw.githubusercontent.com/darizasu/work/master/scripts/ggCor.R")
source("https://raw.githubusercontent.com/AparicioJohan/SpATS.plus/master/All_additional.R")
source("utilities_8_.R")

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

traits <- c("DF", "DPM", "SW100", "TSW", "YdHa_clean")

table(datos$line, datos$rep)


# -------------------------------------------------------------------------
# Descriptive -------------------------------------------------------------
# -------------------------------------------------------------------------

summ <- datos[, traits] %>%
          gather(data = ., key = "traits", value ="value") %>% 
          group_by(traits) %>% 
          summarise( mean = mean(value, na.rm = T),
                     sd = sd(value, na.rm = T),
                     cv = sd / mean,
                     n = n(),
                     n_miss = sum(is.na(value)),
                     miss_perc = n_miss/n
          )
  
summ

# -------------------------------------------------------------------------
# Modelling ---------------------------------------------------------------
# -------------------------------------------------------------------------

spat_g <- spat_e <- spat_h2 <- spat_r2 <- c()
BLUPs <- BLUEs <- mo_rg <- mo_fg <- list()

ncols = length(unique(datos$C))
nrows = length(unique(datos$R))

for (var in traits) {
  
  message("|--------",var,"--------|")
  
  # Genotype as Fixed
  mo_fg[[var]] <- SpATS(response = var,
                             genotype = 'line', 
                             genotype.as.random = F,
                             fixed = NULL ,
                             spatial = ~ PSANOVA(col,
                                                 row, 
                                                 nseg = c(ncols, nrows), 
                                                 degree = c(3,3),
                                                 nest.div=2),
                             random = ~ rep:R + rep:C , 
                             data=datos,
                             control = list(tolerance=1e-03, monitoring=1))
  
  # Genotype as Random
  mo_rg[[var]] <- SpATS(response = var,
                      genotype = 'line', 
                      genotype.as.random = T,
                      fixed = NULL ,
                      spatial = ~ PSANOVA(col,
                                          row, 
                                          nseg = c(ncols, nrows), 
                                          degree = c(3,3),
                                          nest.div=2),
                      random = ~ rep:R + rep:C , 
                      data=datos,
                      control = list(tolerance=1e-03, monitoring=0))
  plot(mo_rg[[var]])
  
  # parameters
  spat_g[var] <- mo_rg[[var]]$var.comp["line"]
  spat_e[var] <- mo_rg[[var]]$psi[1]
  spat_h2[var] <- getHeritability(mo_rg[[var]])
  spat_r2[var] <- R.square(mo_rg[[var]])
  
  # BLUPs
  BLUPs[[var]] <- predict(mo_rg[[var]], which = "line",  predFixed = 'marginal') %>%
                   dplyr::select(line, predicted.values, standard.errors )
  
  # BLUEs
  BLUEs[[var]] <- predict(mo_fg[[var]], which = "line",  predFixed = 'marginal') %>%
                   dplyr::select(line, predicted.values, standard.errors )

}

BLUPs_table <- data.frame(plyr::ldply(BLUPs[], data.frame, .id = "trait")) %>% 
                dplyr::select(- standard.errors) %>% 
                rename(BLUPs = predicted.values) %>% 
                pivot_wider(names_from = c(trait), values_from = c(BLUPs))

BLUEs_table <- data.frame(plyr::ldply(BLUEs[], data.frame, .id = "trait")) %>% 
                dplyr::select(- standard.errors) %>% 
                rename(BLUEs = predicted.values) %>% 
                pivot_wider(names_from = c(trait), values_from = c(BLUEs))


# Save table --------------------------------------------------------------

names(BLUEs_table)[-1] <- paste0(names(BLUEs_table)[-1], "_BLUE")
names(BLUPs_table)[-1] <- paste0(names(BLUPs_table)[-1], "_BLUP")

all <- merge(BLUEs_table, BLUPs_table, by = "line", all = T)


# Summary by trait --------------------------------------------------------

stack <- data.frame(traits, spat_g, spat_e, spat_h2, spat_r2, row.names = NULL)
stack <- merge(stack, summ)


o1 <- stack %>%
  gather(key = "statistic", value = "value", -1) %>%
  filter(statistic %in% c("spat_h2")) %>% 
  ggplot(aes(x = reorder(traits, - value), y = value, label = paste( round(value,2)),  fill = traits))+
  geom_bar(stat = "identity", color = "black")+
  geom_text(nudge_y = 0.05)+
  scale_fill_brewer(palette="Set3")+
  theme_ipsum()+
  labs(title = "Broad-sense Heritability", y = "heritability", x = "traits")+
  theme(axis.text.x = element_text(angle=75, hjust=1))
o1

o2 <- stack %>%
  gather(key = "statistic", value = "value", -1) %>%
  filter(statistic %in% c("cv")) %>% 
  ggplot(aes(x = reorder(traits, - value), y = value, label = paste( round(value,2)),  fill = traits))+
  geom_bar(stat = "identity", color = "black")+
  geom_text(nudge_y = 0.02)+
  scale_fill_brewer(palette="Set3")+
  theme_ipsum()+
  labs(title = "Coefficient of Variation", y = "cv", x = "traits")+
  theme(axis.text.x = element_text(angle=75, hjust=1))
o2


# Correlation -------------------------------------------------------------

g1 <- ggCor(BLUPs_table, colours = c("#db4437","white","#4285f4")) +
        labs(title = "Genotypic Correlation", subtitle = "BLUPs")
g1

g2 <- ggCor(BLUEs_table, colours = c("#db4437","white","#4285f4")) +
  labs(title = "Genotypic Correlation", subtitle = "BLUEs")
g2

ggarrange(g1,g2)

# -------------------------------------------------------------------------
# Dendogram
# -------------------------------------------------------------------------

CC <- cor(BLUPs_table %>% dplyr::select_if(is.numeric), use = "pairwise.complete.obs")
res <- factoextra::hcut(CC, k =  3, stand = FALSE)
g3 <- factoextra::fviz_dend(res,
                            rect = F, 
                            cex = 0.9, 
                            # palette = c("#00AFBB","#2E9FDF", "red"),
                            palette = "jco", 
                            lwd = 0.5, main = "Cluster Dendogram ", 
                            horiz = T, legend = "top") + labs(subtitle = "Correlation")

g3
ggsave("images/correlation_dendogram.png", plot = g3, units = "in", dpi = 300, width = 10, height = 7)

# PCA
pcs <- BLUPs_table  %>% column_to_rownames("line") %>%  dplyr::select_if(is.numeric)
pca <- PCA(pcs)

# Interactiva
explor(pca)

# VAR
var_plot <- fviz_pca_var(pca, col.var = "black", repel = T) + ggtitle("PCA")
ggsave("images/PCA_var.png", plot = var_plot, units = "in", dpi = 300, width = 8, height = 8)

# IND
IND <- fviz_pca_ind(pca, repel = F, alpha.ind = 0.2, col.ind = "grey20", labelsize = 2)
IND

# BIPLOT
BI <- fviz_pca_biplot(pca, repel = T, alpha.ind = 0.8, col.ind = "grey60", labelsize = 4, col.var = "red") # geom.ind = "point"
BI
ggsave("images/biplot.png", plot = BI, units = "in", dpi = 300, width = 12, height = 10)


# -------------------------------------------------------------------------
# Index -------------------------------------------------------------------
# -------------------------------------------------------------------------

bweights <- c(-0.2,-0.2,0.4,0.2,0.3)

res.pca <- pca_index(data = BLUEs_table, id = "line", b = bweights, percentage = 0.15 )
res.pca$final
res.pca$selection

res.pca$results


g0 <- res.pca$results %>% 
  rownames_to_column("line") %>% 
  ggplot(aes(x = reorder(line, -index) , y = index))+
  geom_point()+
  theme(axis.text.x = element_text(angle=90, hjust=1, size = 5))+
  labs(x = "")


library(plotly)
ggplotly(g0)

