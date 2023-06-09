ntiles <- function(X , n = 4){
  cut( X , breaks = quantile( X , probs = (0:n)/n , na.rm = TRUE ) , include.lowest = TRUE )
}

univglms <- function(y,varnames=varnames,data=data,family=binomial,verbose=T){
  counter <- 1
  modelfits <- vector(length(varnames), mode = "list")
  names(modelfits) <- varnames
  for(i in varnames){
    if(verbose){print(counter)}
    modelformula <- paste(y,"~",i)
    modelfits[[i]] <- glm(as.formula(modelformula),data=data,family=family)
    fit <- modelfits[[i]]
    restab <- cbind(counter,exp(coef(summary(fit))[,1]),coef(summary(fit))[,-1],exp(confint.default(fit)),AIC(fit),fit$df.null,fit$df.residual)
    colnames(restab)[2] <- "OR"; colnames(restab)[8:10] <- c("AIC","DF.NULL","DF.RESID")
    if(verbose==T){print(summary(modelfits[[i]]))}
    if(counter==1){results <- restab}else{results <- rbind(results,restab)}
    counter=counter+1
  }
  return(results)
}

model.sel.tab <- function(list){
  n.var <- length(list)
  table <- data.frame(matrix(, nrow = n.var, ncol = 2))
  names(table) <- c("Formula","AIC")
  for(i in 1:n.var){
    table[i,1] <- as.character(formula(list[[i]]))[3]
    table[i,2] <- AIC(list[[i]])
  }
  table <- table[order(table$AIC),]
  return(table)
}

multiv.tab.glmer <- function(mod){
  tab <- as.data.frame(summary(mod)[10]$coefficients)[-1,]
  tab$OR <- exp(tab[,1])
  tab$CIlo <- exp(tab[,1]-1.96*tab[,2])
  tab$CIhi <- exp(tab[,1]+1.96*tab[,2])
  res <- tab[,c(5,6,7,4)]
  return(res)
}

multiv.tab.glm <- function(mod){
  tab <- as.data.frame(summary(mod)[12]$coefficients)
  tab$OR <- exp(tab[,1])
  tab$CIlo <- exp(tab[,1]-1.96*tab[,2])
  tab$CIhi <- exp(tab[,1]+1.96*tab[,2])
  res <- tab[,c(5,6,7,4)]
  return(res)
}

zero.to.min <- function(x){
  zeros <- which(x==0)
  ans <- x
  ans[zeros] <- min(x[-zeros],na.rm=T)
  return(ans)
}

ws_area <- function(site,atos,year){
  newsite <- rep("MBAY",length(site))
  newsite[which(atos>371.5)] <- "MPEN"
  newsite[which(atos>440.5)] <- "BIGS"
  newsite[which(atos>681.5)] <- "SNLO"
  newsite[which(atos>1040)] <- "SABA"
  newsite[which(newsite=="SABA" & year<2012)] <- "SNLO"
  newsite[which(atos>10000)] <- "SNIC"
  newsite[which(atos==321)] <- "ELKS"
  newsite[which(newsite=="ELKS" & year<2012)] <- "MBAY"
  newsite[which(is.na(site))] <- NA
  return(newsite)
}

prop.ci <- function(x,n){
  require(PropCIs)
  limlo <- c()
  limhi <- c()
  for (i in 1:length(x)){
    ans <- blakerci(x[i],n[i],0.95)
    limlo[i] <- ans$conf.int[1]
    limhi[i] <- ans$conf.int[2]
  }
  return(as.data.frame(cbind(limlo,limhi)))
}

# Non-parametric test for associations:
spearman.tab <- function(list,y){
  n.var <- length(list)
  results <- as.data.frame(matrix(, nrow = n.var, ncol = 3))
  for(i in 1:n.var){
    test <- cor.test(list[[i]],y,"two.sided","spearman")
    results[i,] <- with(test,c(names(list)[i],round(estimate,3),round(p.value,4)))
  }
  names(results) <- c("Var","rho","p-value")
  return(results)
}
