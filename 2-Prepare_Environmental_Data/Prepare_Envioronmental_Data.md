---
title: "Prepare_Environmental_Data"
author: "Enrico"
date: "24/1/2022"
output: html_document
editor_options:
  chunk_output_type: console
---

In this markdown I will describe the process of extracting and selecting environmental data.

## Extracting

## Selecting environmental variables

### correlation and VIF between variables
Prepare R and load data
```{R}
library(tidyverse)
library("Hmisc")
library(corrplot)
library(viridis)
library(psych)
library(VIF)
library(fmsb)
library(usdm)

vif_func<-function(in_frame,thresh=10,trace=T,...){
  
  require(fmsb)
  
  if(class(in_frame) != 'data.frame') in_frame<-data.frame(in_frame)
  
  #get initial vif value for all comparisons of variables
  vif_init<-NULL
  var_names <- names(in_frame)
  for(val in var_names){
    regressors <- var_names[-which(var_names == val)]
    form <- paste(regressors, collapse = '+')
    form_in <- formula(paste(val, '~', form))
    vif_init<-rbind(vif_init, c(val, VIF(lm(form_in, data = in_frame, ...))))
  }
  vif_max<-max(as.numeric(vif_init[,2]), na.rm = TRUE)
  
  if(vif_max < thresh){
    if(trace==T){ #print output of each iteration
      prmatrix(vif_init,collab=c('var','vif'),rowlab=rep('',nrow(vif_init)),quote=F)
      cat('\n')
      cat(paste('All variables have VIF < ', thresh,', max VIF ',round(vif_max,2), sep=''),'\n\n')
    }
    return(var_names)
  }
  else{
    
    in_dat<-in_frame
    
    #backwards selection of explanatory variables, stops when all VIF values are below 'thresh'
    while(vif_max >= thresh){
      
      vif_vals<-NULL
      var_names <- names(in_dat)
      
      for(val in var_names){
        regressors <- var_names[-which(var_names == val)]
        form <- paste(regressors, collapse = '+')
        form_in <- formula(paste(val, '~', form))
        vif_add<-VIF(lm(form_in, data = in_dat, ...))
        vif_vals<-rbind(vif_vals,c(val,vif_add))
      }
      max_row<-which(vif_vals[,2] == max(as.numeric(vif_vals[,2]), na.rm = TRUE))[1]
      
      vif_max<-as.numeric(vif_vals[max_row,2])
      
      if(vif_max<thresh) break
      
      if(trace==T){ #print output of each iteration
        prmatrix(vif_vals,collab=c('var','vif'),rowlab=rep('',nrow(vif_vals)),quote=F)
        cat('\n')
        cat('removed: ',vif_vals[max_row,1],vif_max,'\n\n')
        flush.console()
      }
      
      in_dat<-in_dat[,!names(in_dat) %in% vif_vals[max_row,1]]
      
    }
    
    return(names(in_dat))
    
  }
  
}

varmatrix <- read_tsv("2-Prepare_Environmental_Data/WorldClim_table_persample.tsv") %>% column_to_rownames(., var="sample")
snowdata <- read_tsv("2-Prepare_Environmental_Data/Snow_table_persample.tsv") %>% column_to_rownames(., var="sample")
varmatrix2 <- as.matrix(cbind(varmatrix, snowdata))
varmatrix3 <- t(as.matrix(cbind(varmatrix, snowdata)))
```
Extract Matrix with correlations
```{R}
pdf(file = paste0("2-Prepare_Environmental_Data/plots/allvars_pairs.panels.pdf"),
    width = 13,
    height = 8)
pairs.panels(varmatrix2, scale=T, cex.cor = 1.5)
dev.off()
```
Selection of 4 main categories based on correlation values (<0.7):
 - 1,2,3,4,6,7,9,11,12,14,15,17,19
 - 5,10
 - 8
 - 13,16,18
 - Snow days
 - Jan depth
Selecting 5, 8 and 13 and snow variables for easy interpretation, discarding 3 and 4 for harder interpretation, and 12 for correlation with 13, we need to choose between:
 - 1,2,4,6,7,9,11,14,15,17,19

Reitarative VIF method
```{R}
# all variables
varmatrix_short <- varmatrix2[,c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21)]
# all competing variables (removed 3,4,10,12,13,16,18 for reasons mentioned above)
varmatrix_short <- varmatrix2[,c(1,2,7,6,5,8,9,11,13,14,15,17,20,21)]

# manually do by removing highest VIF from c() or use function below (depends on order of columns) 
varmatrix_short <- varmatrix2[,c(1,2,7,6,5,8,9,11,13,14,15,17,20,21)]
usdm::vifcor(varmatrix_short, th=10)

# final set of variables - loose
vif_func(varmatrix_short, trace =T, thresh = 4)
varmatrix_final <- varmatrix2[,c(2,5,6,8,13,20,21)]

usdm::vifcor(varmatrix_final, th=10)
#   Variables   VIF
#      bio2     3.976499
#      bio5     2.158984
#      bio6     3.750249
#      bio8     1.594988
#     bio13     1.412701
# jan_depth     2.236679
# snow_days     3.329262


# final set of variables - strict
vif_func(varmatrix_short, trace =T, thresh = 2)
varmatrix_final <- varmatrix2[,c(5,6,8,13,20)]

usdm::vifcor(varmatrix_final, th=10)
# Variables      VIF
#      bio6 1.261347
#      bio5 1.467656
#      bio8 1.328688
#     bio13 1.218006
# jan_depth 1.241692

# Export table of selected variables only
export_final <- data.frame(varmatrix_final) %>% rownames_to_column(., var = "sample")

write_tsv(export_final, "2-Prepare_Environmental_Data/uncorrelated_variables_matrix.tsv") 
```
Correlation plot of uncorrelated set of variables
```{R}
varmatrix <- read_tsv("2-Prepare_Environmental_Data/uncorrelated_variables_matrix.tsv") %>% column_to_rownames(., var="sample")
pdf(file = paste0("2-Prepare_Environmental_Data/plots/uncorrelated_variables_pairs.panels.pdf"),
    width = 13,
    height = 8)
pairs.panels(varmatrix, scale=T, cex.cor = 1.5)
dev.off()
```
