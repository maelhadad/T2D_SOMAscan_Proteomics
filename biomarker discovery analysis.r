################################################################################
# Example code for the biomarker discovery analysis:
## First: load the OMICs and phenotype data for both cohorts.

################################################################################
## Second: clean phenotype data to calculate the GDRS score:

### Adaptation of the GDRS risk score was necessary as some of the variables were 
### missing from one or both cohorts, and was performed as follows: 
### We defined smoking status using only information on current and former smoking 
### per se without regard to the amount of cigarettes smoked. We used the average 
### of the original GDRS score weights for each smoking category to represent our 
### combined categories (former: (15+45)/2 = 30, current: (23+77)/2 = 50).

### Family history of diabetes was defined in KORA as having at least one parent 
### or sibling with diabetes and in HUNT as having at least one parent, sibling or 
### a child with diabetes. We calculated the risk of a positive family history by 
### averaging the original GDRS risk scores of having one parent, both parents or 
### at least one sibling with type 2 diabetes ((56+106+48)/3 = 70).

### Our final adapted GDRS score was calculated as follows: 
### 5.1×age in years + 7.6×waist circumference in cm - 2.7×height in cm + 
### 47×hypertension status - 2×physical activity (at least one hour per week) + 
### 30×former smoking + 50×current smoking + 70×family history of type 2 diabetes

koradata$gdrs<-
  with(koradata, I(5.1*age)+ I(7.6*waist_circumference)+ I(-2.7*height)+ 
         I(47*hypertension)+I(-2*physical_activity)+gdrs_smoking+ I(70*diabetes_family_history))

huntdata$gdrs<-
  with(huntdata, I(5.1*age)+ I(7.6*waist_circumference)+ I(-2.7*height)+ 
         I(47*hypertension)+I(-2*physical_activity)+gdrs_smoking+ I(70*diabetes_family_history))

################################################################################
## Third: test ROC-AUC:
mod<-glm(inci.t2d ~ gdrs,family = "binomial", data = koradata)

ypred.test<-predict(mod,koradata,type="response")
ypred.valid<-predict(mod,huntdata,type="response")

ROC_Test <- with(koradata,pROC::roc(response = inci.t2d, predictor = ypred.test))
plot(ROC_Test)
ROC_Valid <- with(huntdata,pROC::roc(response = inci.t2d, predictor = ypred.valid))
plot(ROC_Valid)

### Model with proteins:
#### add replicated proteins that are available in both cohorts: sig.res
mod2<-glm(inci.t2d ~ gdrs + IGFBP2 + Growthhormonereceptor + TGFbRIII + 
            TrATPase + Aminoacylase1 + PAPPA + Afamin + sCD163 + UNC5H4,
          family = "binomial", data = koradata)

ypred.test2<-predict(mod2,koradata,type="response")
ypred.valid2<-predict(mod2,huntdata,type="response")

ROC_Test2 <- with(koradata,pROC::roc(response = inci.t2d, predictor = ypred.test2))
plot(ROC_Test2)
ROC_Valid2 <- with(huntdata,pROC::roc(response = inci.t2d, predictor = ypred.valid2))
plot(ROC_Valid2)

delong_t<-pROC::roc.test(ROC_Test,ROC_Test2)
delong_v<-pROC::roc.test(ROC_Valid,ROC_Valid2)

################################################################################
## LASSO selected proteins:
################################################################################
dat <- na.omit(koradata[,c("inci.t2d","gdrs",sig.res$Protein)])
X <- as.matrix(dat[,c("gdrs",sig.res$Protein)])
Y <- dat$inci.t2d
na_index <- is.na(Y)
fit <- glmnet::glmnet(x = X[!na_index,], y = Y[!na_index],family = "binomial", 
                      alpha = 1, penalty.factor = c(0, rep(1, ncol(X) - 1)))

##Cross validation:
crossval <-  glmnet::cv.glmnet(x = X[!na_index,], y = Y[!na_index],family="binomial",
                               penalty.factor = c(0, rep(1, ncol(X) - 1)))

penalty <- crossval$lambda.min #optimal lambda #minimal shrinkage
fit1 <- glmnet::glmnet(x = X[!na_index,], y = Y[!na_index],family = "binomial", alpha = 1, 
                      lambda = penalty, penalty.factor = c(0, rep(1, ncol(X) - 1))) #estimate the model with that

### Model with LASSO selected proteins:
lasso_mod<-glm(inci.t2d ~ gdrs + TGFbRIII + TrATPase + PAPPA + Afamin + sCD163,
          family = "binomial", data = koradata)

lasso_ypred.test<-predict(lasso_mod,koradata,type="response")
lasso_ypred.valid<-predict(lasso_mod,huntdata,type="response")

lasso_ROC_Test <- with(koradata,pROC::roc(response = inci.t2d, predictor = lasso_ypred.test))
plot(lasso_ROC_Test)
lasso_ROC_Valid <- with(huntdata,pROC::roc(response = inci.t2d, predictor = lasso_ypred.valid))
plot(lasso_ROC_Valid)
