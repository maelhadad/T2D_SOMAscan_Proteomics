################################################################################
# Example code for the proteome wide association analysis:
## First: load the OMICs and phenotype data:

################################################################################
## Second: use for loop to iterate the analysis:
result<-data.frame(matrix(nrow = length(iteration),ncol = 10))
colnames(result)<-c("Protein","Beta","SE","Test_statistic","Pvalue")


for(i in 2:1096){
  mod<- glm(t2d ~ colnames(dat)[i]+ age+ sex+ bmi + smoking+ hypertension ,
            family = binomial(),data = dat,na.action = na.omit)
  result[r,1]=colnames(dat)[i]
  result[r,2:5]=summary(mod)$coefficients[2,]
}

## Third: save results:
xlsx::write.xlsx(result,file = "result_file.xlsx",row.names = FALSE)
