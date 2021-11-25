



library(explor)
library(FactoMineR)
library(factoextra)

selIndex <- function (Y, b, scale = FALSE) 
{
  if (scale) {
    return(scale(Y) %*% b)
  }
  return(Y %*% b)
}



pca_index <- function(data, id, variables = NULL, percentage = 0.15, b ){
  
  data <- as.data.frame(data)
  rownames(data) <- data[,id]
  
  if(is.null(variables)) variables <- names(data)[names(data)!=id]
  data <- data[,variables]
  
  index <- selIndex(Y = as.matrix(data), b = b, scale = T)
  index <- c(index)
  data$index <- index
  data <- data %>% arrange(desc(index))
  data$selected <- NA
  data$selected[1:(round(percentage * nrow(data)))] <- TRUE
  data$selected <- ifelse(is.na(data$selected), FALSE, data$selected)
  
  res.pca <- PCA(data,  graph = FALSE, scale.unit = T, quali.sup = ncol(data))
  
  
  p1 <- fviz_pca_var(res.pca, col.var = "black", repel = T)
  p2 <- fviz_pca_ind(res.pca, label = "none", habillage = data$selected,
                     palette = c("#00AFBB", "#FC4E07" ),  addEllipses = TRUE) 
  
  
  final <- ggarrange(p1,p2, legend = "bottom", common.legend = T) 
  final <- annotate_figure(final, 
                           top = text_grob( paste("Selection:", 
                                                  paste0(percentage*100,"%"), "\n",
                                                  paste("Weights:", "(", paste0(b, collapse = ', '),")","\n"  )), color = "black", face = "bold", size = 14))
  
  selection <- data %>% filter(selected == T)
  
  return(list(res.pca = res.pca, final = final, results = data , selection = selection))
}
