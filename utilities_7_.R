


# lme4 Spatial Plot -------------------------------------------------------


lme4.plot <- function(x,Datos,col="col",row="row",gen="line",all.in.one = TRUE, main = NULL, annotated = FALSE, depict.missing = FALSE, ...) {
  
  Datos <- as.data.frame(Datos,row.names = 1:dim(Datos)[1])
  xlab <- col
  ylab <- row
  x.coord <- Datos[,col]
  y.coord <- Datos[,row]
  response <- Datos[, names(x@frame)[1]]  # !is.na(Datos[,names(x@frame)[1]]) 
  Datos$rowname <- as.numeric(rownames(Datos))
  
  residuals <- data.frame(rowname=as.numeric(names(residuals(x))), residuals= as.numeric(residuals(x)))
  residuals <- merge(residuals,Datos,by="rowname", all = T, sort=T)
  residuals <- residuals$residuals
  
  fitted <- data.frame(rowname=as.numeric(names(fitted.values(x))), fitted=as.numeric(fitted.values(x)))
  fitted <- merge(fitted,Datos,by="rowname", all = T, sort=T)
  fitted <- fitted$fitted
  
  geno.pred <- data.frame(genotype=rownames(ranef(x)[[gen]])   , predicted.value=ranef(x)[[gen]][,1])
  names(geno.pred)[1] <- gen
  columns <- seq(min(x.coord), max(x.coord), by = min(diff(sort(unique(x.coord)))))
  rows <- seq(min(y.coord), max(y.coord), by = min(diff(sort(unique(y.coord)))))
  xy.coord <- data.table::data.table(expand.grid(columns = columns, rows = rows))
  v <- merge(geno.pred,Datos,by = gen,all.y = TRUE, sort=F)
  xx <- arrange(v,rowname)
  geno.pred <- as.vector(xx$predicted.value)
  
  environment <- fitted-geno.pred-fixef(x)["(Intercept)"]
  environment <- round(environment,6)
  
  data.table::setNumericRounding(2)
  
  if(is.null(main)) main = paste("Trait: ",names(x@frame)[1] , sep = "")
  
  data.table::setkeyv(xy.coord, c("rows", "columns"))
  ONE <- rep(1, length(x.coord))    
  df <- data.table::data.table(columns = x.coord, rows = y.coord, 
                   response = response, fitted = fitted,environment,
                   residuals = residuals, geno.pred = geno.pred, ONE = ONE)
  data.table::setkeyv(df, c("rows", "columns"))
  df <- df[xy.coord]
  df <- df[order(df$columns, df$rows),]
  
  colors = topo.colors(100)
  
  main.legends <- c('Raw data', 'Fitted data', 'Residuals',"Effect Design"  ,"Genotypic BLUPs", 'Histogram')
  if(all.in.one) {
    op <- par(mfrow = c(2,3), oma = c(ifelse(annotated, 12, 2), 1, 3, 2), mar = c(2.7, 4, 2.5, 2.5), mgp = c(1.7, 0.5, 0))                
  } else {
    if(!is.null(main))
      main.legends <- rep(main, length(main.legends))
  }
  
  range <- range(c(response, fitted), na.rm = TRUE)
  fields::image.plot(columns, rows, t(matrix(df$response, ncol = length(columns), nrow = length(rows))), main = main.legends[1], col = colors, xlab = xlab, ylab = ylab, zlim = range, graphics.reset = TRUE, ...)
  if(!all.in.one)
    readline("Press return for next page....")
  fields::image.plot(columns, rows, t(matrix(df$fitted, ncol = length(columns), nrow = length(rows))), main = main.legends[2], col = colors, xlab = xlab, ylab = ylab, zlim = range, graphics.reset = TRUE, ...)
  if(!all.in.one)
    readline("Press return for next page....")
  fields::image.plot(columns, rows, t(matrix(df$residuals, ncol = length(columns), nrow = length(rows))), main = main.legends[3], col = colors, xlab = xlab, ylab = ylab, graphics.reset = TRUE, ...)
  if(!all.in.one)
    readline("Press return for next page....")
  fields::image.plot(columns, rows, t(matrix(df$environment, ncol = length(columns), nrow = length(rows))), main = main.legends[4], col = colors, xlab = xlab, ylab = ylab, graphics.reset = TRUE, ...)
  if(!all.in.one)
    readline("Press return for next page....")
  fields::image.plot(columns, rows, t(matrix(df$geno.pred, ncol = length(columns), nrow = length(rows))), main = main.legends[5], col = colors, xlab = xlab, ylab = ylab, graphics.reset = TRUE, ...)
  if(!all.in.one)
    readline("Press return for next page....")
  # 
  suppressWarnings(hist(unique(geno.pred), main = main.legends[6], xlab = main.legends[6], ...))        
  title("")
  mtext(main, cex = 1.5, outer = TRUE, side = 3)
  invisible(df)
}


# Genetic Gain ------------------------------------------------------------



Genetic_Gain <- function(data, heritability, genotype, exp , trace = F, percentage= 0.15){
  
  data_results <- data.frame("Model" = as.character(), "Selected" = as.numeric(), "Total" = as.numeric(),
                             "U_selected" = as.numeric(), "U_pop" = as.numeric(), "S" = as.numeric(), 
                             "h2" = as.numeric())
  
  for (i in exp) {
    tmp <- na.omit(data[,c(genotype,i)])
    total_g <- length(unique(tmp[[genotype]]))
    ng <- round(total_g*percentage)
    
    u_selected <- tmp %>% arrange(desc(.data[[i]])) %>% .[1:ng,] %>% summarise(u=mean(.data[[i]])) %>% pull(u)
    u_populati <- tmp %>% summarise(u=mean(.data[[i]])) %>% pull(u)
    S <- u_selected-u_populati
    h2_tmp <- heritability[i]
    # R_tmp <- h2_tmp*(S)         
    # delta_R_tmp <- R_tmp/14       
    
    data_tmp <- data.frame(Model = i, Selected =  ng , Total= total_g,
                           U_selected = u_selected, U_pop = u_populati, S=S,  h2 = h2_tmp)
    data_results <- rbind(data_results, data_tmp)
    
    if(trace){
      if(i==exp[1]) cat("Model","\tSelected", "\tTotal", "\tU_selected", "\tU_pop", "\tS" , "\th2", "\tR")
      cat("\n",i, "\t" , ng , "\t" , total_g, "\t", u_selected, "\t", u_populati, "\t", S, "\t", h2_tmp)
    }
  }
  return(data_results)
}



# Cullis Heritability -----------------------------------------------------

h.cullis <- function(model, gen){
  aveped <- mean(attr(ranef(model,drop=T)[[gen]],"postVar"))
  vc.g <- as.data.frame(VarCorr(model))
  vc.g <- vc.g[vc.g$grp==gen,"vcov"]
  ifelse(vc.g==0, 0 , round(1-aveped/vc.g,2) )
}
