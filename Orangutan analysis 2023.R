library(fitdistrplus)
library(devtools)

source("https://raw.githubusercontent.com/sealavi/GAM-first-derivative-functions/main/gam%20first%20derivitave%20plot%20function_general.R")


data=read.csv("C:/Users/salavi/Documents/20211031 To Send to Shauhin for analyses/20211031 To Send to Shauhin for analyses/nutrition_urine_biomarkers_by_follow_RSAB_29oct2021.csv")
data$urea_time_collected[which(!is.na(data$urea_time_collected))]=paste(data$date[which(!is.na(data$urea_time_collected))], data$urea_time_collected[which(!is.na(data$urea_time_collected))], sep=" ")
data$ucp_time_collected[which(!is.na(data$ucp_time_collected))]=paste(data$date[which(!is.na(data$ucp_time_collected))], data$ucp_time_collected[which(!is.na(data$ucp_time_collected))], sep=" ")
data$isotope_time_collected[which(!is.na(data$isotope_time_collected))]=paste(data$date[which(!is.na(data$isotope_time_collected))], data$isotope_time_collected[which(!is.na(data$isotope_time_collected))], sep=" ")
data$urea_time_collected=as.POSIXct(data$urea_time_collected, format="%Y-%m-%d %I:%M:%S %p",origin="01-01-1900", tz="Asia/Jakarta")
data$ucp_time_collected=as.POSIXct(data$ucp_time_collected, format="%Y-%m-%d %I:%M:%S %p",origin="01-01-1900", tz="Asia/Jakarta")
data$isotope_time_collected=as.POSIXct(data$isotope_time_collected, format="%Y-%m-%d %I:%M:%S %p",origin="01-01-1900", tz="Asia/Jakarta")

data$urea_time_collected=as.numeric(hms::as_hms(data$urea_time_collected))
data$ucp_time_collected=as.numeric(hms::as_hms(data$ucp_time_collected))
data$isotope_time_collected=as.numeric(hms::as_hms(data$isotope_time_collected))

data$date=as.POSIXct(data$date, format="%Y-%m-%d",origin="01-01-1900", tz="Asia/Jakarta")
data$name_focal=as.factor(data$name_focal)
data$class_focal=as.factor(data$class_focal)



##assess distributional family of response variables
hist(na.omit(data$urea_sgcor), breaks=100) ##very skewed
ureahist=fitdist(as.numeric(na.omit(data$urea_sgcor)), distr = "gamma", method="mme")
summary(ureahist)
plot(ureahist)

ureahist2=fitdist(as.numeric(na.omit(data$urea_sgcor)), distr = "lnorm", method="mle")
plot(ureahist2)  ##gamma looks slightly better than lognormal, might have to try both

hist(na.omit(data$ucp_sgcor), breaks=100) ##very skewed
ucphist=fitdist(as.numeric(na.omit(data$ucp_sgcor)), distr = "gamma", method="mme")
plot(ucphist)

ucphist2=fitdist(as.numeric(na.omit(data$ucp_sgcor)), distr = "lnorm", method="mle")
plot(ucphist2)  ##lognormal looks way better than gamma!
summary(ucphist2)

hist(na.omit(data$dn15_result), breaks=100) ##pretty normal 

dn15hist=fitdist(as.numeric(na.omit(data$dn15_result)), distr = "norm", method="mle")
plot(dn15hist)  ##Gaussian fits well

##Will use winning families as model families for regressions 
detach("package:fitdistrplus", unload = TRUE)

library(ggplot2)
library(ggpubr)
library(brms)
options(mc.cores = parallel::detectCores())

##FAI only models###
get_prior(urea_sgcor~s(fai) + (1 +  fai  | name_focal), 
              data = data[complete.cases(data[,c("urea_sgcor","fai")]),],
              family=Gamma(link="log"))


urea_FAI_model<-brm(urea_sgcor~s(fai) + s(urea_time_collected) + (1 +  fai + urea_time_collected  | name_focal), 
                    data = data[complete.cases(data[,c("urea_sgcor","fai")]),],
                    family=Gamma(link="log"),iter = 8000, inits=0,
                    prior=c(prior(normal(0,10),class="Intercept"),
                            prior(normal(0,10),class="b", coef =sfai_1),
                            prior(normal(0,10),class="b", coef =surea_time_collected_1),
                            prior(gamma(.01,.01),class="shape"),
                            prior(normal(0,10),class="sd")),
                    control = list(max_treedepth = 11, adapt_delta = .999))

summary(urea_FAI_model, digits=9)
plot(urea_FAI_model)
plot(conditional_effects(urea_FAI_model,spaghetti=TRUE)) ##plots the regression lines/conditional effects
pp_check(urea_FAI_model)
loo_R2(urea_FAI_model)
bayes_R2(urea_FAI_model)
loo(urea_FAI_model)
urea_FAI_model=add_criterion(urea_FAI_model, c("kfold"))
deriv_plot(model=urea_FAI_model,term="s(fai)",main="fai",eps=0.001,output="urea_FAI_deriv")

urea_sgcor_FAI_model2<-brm(urea_sgcor~s(fai) + s(urea_time_collected) + (1 +  fai + urea_time_collected  | name_focal), 
                         data = data[complete.cases(data[,c("urea_sgcor","fai")]),],
                         family=lognormal(),iter = 2000, 
                         prior=c(prior(normal(0,10),class="Intercept"),
                                 prior(normal(0,10),class="b", coef =sfai_1),
                                 prior(normal(0,10),class="b", coef =surea_time_collected_1),
                                 prior(normal(0,10),class="sigma"),
                                 prior(normal(0,10),class="sd")),
                         control = list(max_treedepth = 11, adapt_delta = .999))

summary(urea_sgcor_FAI_model2, digits=9)
plot(urea_sgcor_FAI_model2)
plot(conditional_effects(urea_sgcor_FAI_model2,spaghetti=TRUE)) ##plots the regression lines/conditional effects
pp_check(urea_sgcor_FAI_model2) ###looks worse than gamma
loo_R2(urea_sgcor_FAI_model2)
bayes_R2(urea_sgcor_FAI_model2)
loo(urea_sgcor_FAI_model2)


get_prior(ucp_sgcor~s(fai) + s(urea_time_collected) + (1 +  fai + urea_time_collected  | name_focal), 
          data = data[complete.cases(data[,c("urea_sgcor","fai")]),],
          family=lognormal())

ucp_sgcor_FAI_model<-brm(ucp_sgcor~s(fai) + s(ucp_time_collected) + (1 +  fai + ucp_time_collected  | name_focal), 
                    data = data[complete.cases(data[,c("ucp_sgcor","ucp_time_collected","fai")]),],
                    family=lognormal(),iter = 2000, 
                    prior=c(prior(normal(0,10),class="Intercept"),
                            prior(normal(0,10),class="b", coef =sfai_1),
                            prior(normal(0,10),class="b", coef =sucp_time_collected_1),
                            prior(normal(0,10),class="sigma"),
                            prior(normal(0,10),class="sd")),
                    control = list(max_treedepth = 10, adapt_delta = .999))

summary(ucp_sgcor_FAI_model, digits=9)
plot(conditional_effects(ucp_sgcor_FAI_model,spaghetti=FALSE)) ##plots the regression lines/conditional effects
pp_check(ucp_sgcor_FAI_model)
loo_R2(ucp_sgcor_FAI_model)
bayes_R2(ucp_sgcor_FAI_model)
deriv_plot(model=ucp_sgcor_FAI_model,term="s(fai)",main="fai",eps=0.001,output="ucp_sgcor_FAI_deriv")


data_FAI_model<-brm(dn15_result~s(fai) + s(isotope_time_collected) + (1 +  fai + isotope_time_collected  | name_focal), 
                         data = data[complete.cases(data[,c("dn15_result","isotope_time_collected","fai")]),],
                         family=gaussian(),iter = 2000, 
                         prior=c(prior(normal(0,10),class="Intercept"),
                                 prior(normal(0,10),class="b", coef =sfai_1),
                                 prior(normal(0,10),class="b", coef =sisotope_time_collected_1),
                                 prior(normal(0,10),class="sigma"),
                                 prior(normal(0,10),class="sd")),
                         control = list(max_treedepth = 10, adapt_delta = .999))

summary(data_FAI_model, digits=9)
plot(conditional_effects(data_FAI_model,spaghetti=FALSE)) ##plots the regression lines/conditional effects
pp_check(data_FAI_model)
loo_R2(data_FAI_model)
bayes_R2(data_FAI_model)
deriv_plot(model=data_FAI_model,term="s(fai)",main="fai",eps=0.001,output="dn15_FAI_deriv")

get_prior(pos_neg_dn~s(fai) + (1 +  fai  | name_focal), 
          data = data[complete.cases(data[,c("pos_neg_dn","fai")]),],
          family=bernoulli())

pos_neg_dn_FAI_model<-brm(pos_neg_dn~s(fai) + (1 + fai | name_focal), 
                    data = data[complete.cases(data[,c("pos_neg_dn","fai")]),],
                    family=bernoulli(),iter = 3000, 
                    prior=c(prior(student_t(4,0,1.5),class="Intercept"),
                            prior(student_t(4,0,1.5),class="b"),
                            prior(student_t(4,0,1.5),class="sds"),
                            prior(student_t(4,0,1.5),class="sd")),
                    control = list(max_treedepth = 10, adapt_delta = .999))

summary(pos_neg_dn_FAI_model, digits=9)
plot(conditional_effects(pos_neg_dn_FAI_model,spaghetti=FALSE), points = TRUE) ##plots the regression lines/conditional effects
pp_check(pos_neg_dn_FAI_model)
loo_R2(pos_neg_dn_FAI_model)
bayes_R2(pos_neg_dn_FAI_model)
deriv_plot(model=pos_neg_dn_FAI_model,term="s(fai)",main="fai",eps=0.001,output="ket_FAI_deriv")


#####Ketone MODELS######

Ketone_KCAL<-brm(bf(pos_neg_dn|trials(1)~s(fai) + class_focal + s(total_kcal_using_ap_low_fermentation) + (1 +  fai + total_kcal_using_ap_low_fermentation   | name_focal)), 
                 data = data[complete.cases(data[,c("pos_neg_dn","fai")]),],
                 family=zero_inflated_binomial(link = "logit"),iter = 3000, 
                 prior=c(prior(normal(0,1),class="Intercept"),
                         prior(normal(0,1),class="b"),
                         prior(exponential(2),class="sd")),
                 control = list(max_treedepth = 11, adapt_delta = .9999))

summary(Ketone_KCAL, digits=9)
plot(conditional_effects(Ketone_KCAL,spaghetti=FALSE),points=TRUE) ##plots the regression lines/conditional effects
pp_check(Ketone_KCAL)
loo_R2(Ketone_KCAL)
bayes_R2(Ketone_KCAL)
gc()
Ketone_KCAL=add_criterion(Ketone_KCAL, c("loo", "waic"))



Ketone_NPE<-brm(pos_neg_dn~s(fai) + class_focal + s(total_kcal_npe_low_fermentation) + (1 +  fai + total_kcal_npe_low_fermentation   | name_focal), 
                data = data[complete.cases(data[,c("pos_neg_dn","fai","total_kcal_npe_low_fermentation")]),],
                family=bernoulli(),iter = 3000, 
                prior=c(prior(normal(0,10),class="Intercept"),
                        prior(normal(0,10),class="b"),
                        
                        prior(normal(0,10),class="sd")),
                control = list(max_treedepth = 11, adapt_delta = .999))

summary(Ketone_NPE, digits=9)
plot(conditional_effects(Ketone_NPE,spaghetti=FALSE),points=TRUE) ##plots the regression lines/conditional effects
pp_check(Ketone_NPE)
loo_R2(Ketone_NPE)
bayes_R2(Ketone_NPE)
gc()



Ketone_AP<-brm(pos_neg_dn~s(fai) + class_focal + s(ap_kcal) + (1 +  fai + ap_kcal   | name_focal), 
               data = data[complete.cases(data[,c("pos_neg_dn","fai","ap_kcal")]),],
               family=bernoulli(),iter = 3000, 
               prior=c(prior(normal(0,10),class="Intercept"),
                       prior(normal(0,10),class="b"),
                       
                       prior(normal(0,10),class="sd")),
               control = list(max_treedepth = 11, adapt_delta = .9999))

summary(Ketone_AP, digits=9)
plot(conditional_effects(Ketone_AP,spaghetti=FALSE),points=TRUE) ##plots the regression lines/conditional effects
pp_check(Ketone_AP)
loo_R2(Ketone_AP)
bayes_R2(Ketone_AP)
gc()



Ketone_Lipid<-brm(pos_neg_dn~s(fai) + class_focal + s(lipid_kcal) + (1 +  fai + lipid_kcal   | name_focal), 
                  data = data[complete.cases(data[,c("pos_neg_dn","fai","lipid_kcal")]),],
                  family=bernoulli(),iter = 3000, 
                  prior=c(prior(normal(0,10),class="Intercept"),
                          prior(normal(0,10),class="b"),
                          
                          prior(normal(0,10),class="sd")),
                  control = list(max_treedepth = 11, adapt_delta = .99999))

summary(Ketone_Lipid, digits=9)
plot(conditional_effects(Ketone_Lipid,spaghetti=FALSE),points=TRUE) ##plots the regression lines/conditional effects
pp_check(Ketone_Lipid)
loo_R2(Ketone_Lipid)
bayes_R2(Ketone_Lipid)
gc()




Ketone_TNC<-brm(pos_neg_dn~s(fai) + class_focal + s(tnc_kcal) + (1 +  fai + tnc_kcal   | name_focal), 
                data = data[complete.cases(data[,c("pos_neg_dn","fai","tnc_kcal")]),],
                family=bernoulli(),iter = 3000, 
                prior=c(prior(normal(0,10),class="Intercept"),
                        prior(normal(0,10),class="b"),
                        
                        prior(normal(0,10),class="sd")),
                control = list(max_treedepth = 11, adapt_delta = .999))

summary(Ketone_TNC, digits=9)
plot(conditional_effects(Ketone_TNC,spaghetti=FALSE),points=TRUE) ##plots the regression lines/conditional effects
pp_check(Ketone_TNC)
loo_R2(Ketone_TNC)
bayes_R2(Ketone_TNC)
gc()


Ketone_NDF<-brm(pos_neg_dn~s(fai) + class_focal + s(ndf_kcal_low) + (1 +  fai + ndf_kcal_low   | name_focal), 
                data = data[complete.cases(data[,c("pos_neg_dn","fai","ndf_kcal_low")]),],
                family=bernoulli(),iter = 3000, 
                prior=c(prior(normal(0,10),class="Intercept"),
                        prior(normal(0,10),class="b"),
                        
                        prior(normal(0,10),class="sd")),
                control = list(max_treedepth = 11, adapt_delta = .999))

summary(Ketone_NDF, digits=9)
plot(conditional_effects(Ketone_NDF,spaghetti=FALSE),points=TRUE) ##plots the regression lines/conditional effects
pp_check(Ketone_NDF)
loo_R2(Ketone_NDF)
bayes_R2(Ketone_NDF)
gc()



Ketone_npe_ap<-brm(pos_neg_dn~s(fai) + class_focal + s(npe_to_ap) + (1 +  fai + npe_to_ap   | name_focal), 
                   data = data[complete.cases(data[,c("pos_neg_dn","fai","npe_to_ap")]),],
                   family=bernoulli(),iter = 3000, 
                   prior=c(prior(normal(0,10),class="Intercept"),
                           prior(normal(0,10),class="b"),
                           
                           prior(normal(0,10),class="sd")),
                   control = list(max_treedepth = 11, adapt_delta = .999))

summary(Ketone_npe_ap, digits=9)
plot(conditional_effects(Ketone_npe_ap,spaghetti=FALSE),points=TRUE) ##plots the regression lines/conditional effects
pp_check(Ketone_npe_ap)
loo_R2(Ketone_npe_ap)
bayes_R2(Ketone_npe_ap)
gc()



Ketone_lipid_ap<-brm(pos_neg_dn~s(fai) + class_focal + s(lipid_to_ap) + (1 +  fai + lipid_to_ap   | name_focal), 
                     data = data[complete.cases(data[,c("pos_neg_dn","fai","lipid_to_ap")]),],
                     family=bernoulli(),iter = 3000, 
                     prior=c(prior(normal(0,10),class="Intercept"),
                             prior(normal(0,10),class="b"),
                             
                             prior(normal(0,10),class="sd")),
                     control = list(max_treedepth = 11, adapt_delta = .99999))

summary(Ketone_lipid_ap, digits=9)
plot(conditional_effects(Ketone_lipid_ap,spaghetti=FALSE),points=TRUE) ##plots the regression lines/conditional effects
pp_check(Ketone_lipid_ap)
loo_R2(Ketone_lipid_ap)
bayes_R2(Ketone_lipid_ap)
gc()






Ketone_lipid_tnc<-brm(pos_neg_dn~s(fai) + class_focal + s(lipid_to_tnc) + (1 +  fai + lipid_to_tnc   | name_focal), 
                      data = data[complete.cases(data[,c("pos_neg_dn","fai","lipid_to_tnc")]),],
                      family=bernoulli(),iter = 3000, 
                      prior=c(prior(normal(0,10),class="Intercept"),
                              prior(normal(0,10),class="b"),
                              
                              prior(normal(0,10),class="sd")),
                      control = list(max_treedepth = 11, adapt_delta = .999))

summary(Ketone_lipid_tnc, digits=9)
plot(conditional_effects(Ketone_lipid_tnc,spaghetti=FALSE),points=TRUE) ##plots the regression lines/conditional effects
pp_check(Ketone_lipid_tnc)
loo_R2(Ketone_lipid_tnc)
bayes_R2(Ketone_lipid_tnc)
gc()





Ketone_lipid_ndf<-brm(pos_neg_dn~s(fai) + class_focal + s(lipid_to_ndf_low) + (1 +  fai + lipid_to_ndf_low   | name_focal), 
                      data = data[complete.cases(data[,c("pos_neg_dn","fai","lipid_to_ndf_low")]),],
                      family=bernoulli(),iter = 3000, 
                      prior=c(prior(normal(0,10),class="Intercept"),
                              prior(normal(0,10),class="b"),
                              
                              prior(normal(0,10),class="sd")),
                      control = list(max_treedepth = 11, adapt_delta = .99999))

summary(Ketone_lipid_ndf, digits=9)
plot(conditional_effects(Ketone_lipid_ndf,spaghetti=FALSE),points=TRUE) ##plots the regression lines/conditional effects
pp_check(Ketone_lipid_ndf)
loo_R2(Ketone_lipid_ndf)
bayes_R2(Ketone_lipid_ndf)
gc()


Ketone_tnc_ap<-brm(pos_neg_dn~s(fai) + class_focal + s(tnc_to_ap) + (1 +  fai + tnc_to_ap   | name_focal), 
                   data = data[complete.cases(data[,c("pos_neg_dn","fai","tnc_to_ap")]),],
                   family=bernoulli(),iter = 3000, 
                   prior=c(prior(normal(0,10),class="Intercept"),
                           prior(normal(0,10),class="b"),
                           
                           prior(normal(0,10),class="sd")),
                   control = list(max_treedepth = 11, adapt_delta = .9999))

summary(Ketone_tnc_ap, digits=9)
plot(conditional_effects(Ketone_tnc_ap,spaghetti=FALSE),points=TRUE) ##plots the regression lines/conditional effects
pp_check(Ketone_tnc_ap)
loo_R2(Ketone_tnc_ap)
bayes_R2(Ketone_tnc_ap)
gc()



Ketone_tnc_ndf<-brm(pos_neg_dn~s(fai) + class_focal + s(tnc_to_ndf_low) + (1 +  fai + tnc_to_ndf_low   | name_focal), 
                    data = data[complete.cases(data[,c("pos_neg_dn","fai","tnc_to_ndf_low")]),],
                    family=bernoulli(),iter = 3000, 
                    prior=c(prior(normal(0,10),class="Intercept"),
                            prior(normal(0,10),class="b"),
                            
                            prior(normal(0,10),class="sd")),
                    control = list(max_treedepth = 11, adapt_delta = .9999))

summary(Ketone_tnc_ndf, digits=9)
plot(conditional_effects(Ketone_tnc_ndf,spaghetti=FALSE),points=TRUE) ##plots the regression lines/conditional effects
pp_check(Ketone_tnc_ndf)
loo_R2(Ketone_tnc_ndf)
bayes_R2(Ketone_tnc_ndf)
gc()




Ketone_ap_ndf<-brm(pos_neg_dn~s(fai) + class_focal + s(ap_to_ndf_low) + (1 +  fai + ap_to_ndf_low   | name_focal), 
                   data = data[complete.cases(data[,c("pos_neg_dn","fai","ap_to_ndf_low")]),],
                   family=bernoulli(),iter = 3000, 
                   prior=c(prior(normal(0,10),class="Intercept"),
                           prior(normal(0,10),class="b"),
                           
                           prior(normal(0,10),class="sd")),
                   control = list(max_treedepth = 11, adapt_delta = .999))

summary(Ketone_ap_ndf, digits=9)
plot(conditional_effects(Ketone_ap_ndf,spaghetti=FALSE),points=TRUE) ##plots the regression lines/conditional effects
pp_check(Ketone_ap_ndf)
loo_R2(Ketone_ap_ndf)
bayes_R2(Ketone_ap_ndf)
gc()






Ketone_PREV_AP<-brm(bf(pos_neg_dn ~s(fai) + class_focal + s(prev_day_ap_kcal) + (1 +  fai + prev_day_ap_kcal   | name_focal)), 
                    data = data[complete.cases(data[,c("pos_neg_dn","fai","prev_day_ap_kcal")]),],
                    family=bernoulli(),iter = 3000, 
                    prior=c(prior(student_t(4,0,1.5),class="Intercept"),
                            prior(student_t(4,0,1.5),class="b"),
                            prior(student_t(4,0,1.5),class="sds"),
                            prior(student_t(4,0,1.5),class="sd")),
                    control = list(max_treedepth = 10, adapt_delta = .9))

summary(Ketone_PREV_AP, digits=9)
plot(conditional_effects(Ketone_PREV_AP,spaghetti=FALSE)) ##plots the regression lines/conditional effects
pp_check(Ketone_PREV_AP)
loo_R2(Ketone_PREV_AP)
bayes_R2(Ketone_PREV_AP)
gc()





testKetone_PREV_AP<-brm(bf(pos_neg_dn|trials(1) ~s(fai) + class_focal + s(prev_day_ap_kcal) + (1 +  fai + prev_day_ap_kcal   | name_focal)), 
                    data = data[complete.cases(data[,c("pos_neg_dn","fai","prev_day_ap_kcal")]),],
                    sample_prior="only",
                    family=zero_inflated_binomial(link = "logit"),iter = 3000, 
                    prior=c(prior(normal(0,.1),class="Intercept"),
                            prior(normal(0,.1),class="b"),
                            prior(normal(0,1.5),class="sds"),
                            prior(normal(0,.1),class="sd")),
                    control = list(max_treedepth = 12, adapt_delta = .9999999999999999))

summary(Ketone_PREV_AP, digits=9)
ytest <- data$pos_neg_dn[complete.cases(data[,c("pos_neg_dn","fai","prev_day_ap_kcal")])]
yreptest <- posterior_predict(testKetone_PREV_AP ,ndraws = 500)

bayesplot::ppc_dens_overlay(ytest, 
                 yreptest) + coord_cartesian(ylim = c(0, 50))

pp_check(testKetone_PREV_AP)





Ketone_PREV_lipid<-brm(bf(pos_neg_dn~s(fai) + class_focal + s(prev_day_lipid_kcal) + (1 +  fai + prev_day_lipid_kcal   | name_focal)), 
                       data = data[complete.cases(data[,c("pos_neg_dn","fai","prev_day_lipid_kcal")]),],
                       family=bernoulli(),iter = 3000, 
                       prior=c(prior(student_t(4,0,1.5),class="Intercept"),
                               prior(student_t(4,0,1.5),class="b"),
                               prior(student_t(4,0,1.5),class="sds"),
                               prior(student_t(4,0,1.5),class="sd")),
                       control = list(max_treedepth = 10, adapt_delta = .9))

summary(Ketone_PREV_lipid, digits=9)
plot(conditional_effects(Ketone_PREV_lipid,spaghetti=FALSE)) ##plots the regression lines/conditional effects
pp_check(Ketone_PREV_lipid)
loo_R2(Ketone_PREV_lipid)
bayes_R2(Ketone_PREV_lipid)
gc()

fitBinom=fitdistrplus::fitdist(data=data$pos_neg_dn[complete.cases(data[,c("pos_neg_dn","fai","prev_day_lipid_kcal")])], dist="binom", fix.arg=list(size=1), start=list(prob=0.1), method="mle")
plot(fitBinom)
testKetone_PREV_lipid<-brm(bf(pos_neg_dn|trials(1)~s(fai) + class_focal + s(prev_day_lipid_kcal) + (1 +  fai + prev_day_lipid_kcal   | name_focal)), 
                           data = data[complete.cases(data[,c("pos_neg_dn","fai","prev_day_lipid_kcal")]),],
                        sample_prior="only",
                        family=zero_inflated_binomial(link = "logit"),iter = 3000, 
                        prior=c(prior(normal(0,.1),class="Intercept"),
                                prior(normal(0,.1),class="b"),
                                prior(normal(0,1.5),class="sds"),
                                prior(normal(0,.1),class="sd")),
                        control = list(max_treedepth = 12, adapt_delta = .9999999999999999))

ytest <- data$pos_neg_dn[complete.cases(data[,c("pos_neg_dn","fai","prev_day_lipid_kcal")])]
yreptest <- posterior_predict(testKetone_PREV_lipid ,ndraws = 500)

bayesplot::ppc_dens_overlay(ytest, 
                            yreptest) + coord_cartesian(ylim = c(0, 50))

pp_check(testKetone_PREV_lipid)



Ketone_PREV_TNC<-brm(bf(pos_neg_dn~s(fai) + class_focal + s(prev_day_tnc_kcal) + (1 +  fai + prev_day_tnc_kcal   | name_focal)), 
                     data = data[complete.cases(data[,c("pos_neg_dn","fai","prev_day_tnc_kcal")]),],
                     family=bernoulli(),iter = 3000, 
                     prior=c(prior(student_t(4,0,1.5),class="Intercept"),
                             prior(student_t(4,0,1.5),class="b"),
                             prior(student_t(4,0,1.5),class="sds"),
                             prior(student_t(4,0,1.5),class="sd")),
                     control = list(max_treedepth = 10, adapt_delta = .99))

summary(Ketone_PREV_TNC, digits=9)
plot(conditional_effects(Ketone_PREV_TNC,spaghetti=FALSE)) ##plots the regression lines/conditional effects
pp_check(Ketone_PREV_TNC)
loo_R2(Ketone_PREV_TNC)
bayes_R2(Ketone_PREV_TNC)
gc()




Ketone_PREV_NDF<-brm(pos_neg_dn~s(fai) + class_focal + s(prev_day_ndf_kcal_low) + (1 +  fai + prev_day_ndf_kcal_low   | name_focal), 
                     data = data[complete.cases(data[,c("pos_neg_dn","fai","prev_day_ndf_kcal_low")]),],
                     family=bernoulli(),iter = 3000, 
                     prior=c(prior(student_t(4,0,1.5),class="Intercept"),
                             prior(student_t(4,0,1.5),class="b"),
                             prior(student_t(4,0,1.5),class="sds"),
                             prior(student_t(4,0,1.5),class="sd")),
                     control = list(max_treedepth = 10, adapt_delta = .99))


summary(Ketone_PREV_NDF, digits=9)
plot(conditional_effects(Ketone_PREV_NDF,spaghetti=FALSE),points=TRUE) ##plots the regression lines/conditional effects
pp_check(Ketone_PREV_NDF)
loo_R2(Ketone_PREV_NDF)
bayes_R2(Ketone_PREV_NDF)
gc()


Ketone_PREV_NDF<-brm(pos_neg_dn~s(fai) + class_focal + s(prev_day_ndf_kcal_low) + (1 +  fai + prev_day_ndf_kcal_low   | name_focal), 
                     data = data[complete.cases(data[,c("pos_neg_dn","fai","prev_day_ndf_kcal_low")]),],
                     family=bernoulli(),iter = 3000, 
                     prior=c(prior(student_t(4,0,1.5),class="Intercept"),
                             prior(student_t(4,0,1.5),class="b"),
                             prior(student_t(4,0,1.5),class="sds"),
                             prior(student_t(4,0,1.5),class="sd")),
                     control = list(max_treedepth = 10, adapt_delta = .99))

summary(Ketone_PREV_NDF, digits=9)
plot(conditional_effects(Ketone_PREV_NDF,spaghetti=FALSE),points=TRUE) ##plots the regression lines/conditional effects
pp_check(Ketone_PREV_NDF)
loo_R2(Ketone_PREV_NDF)
bayes_R2(Ketone_PREV_NDF)
gc()



Ketone_PREV_KCAL<-brm(pos_neg_dn~s(fai) + class_focal + s(prev_day_total_kcal_using_ap_low_fermentation) + (1 +  fai + prev_day_total_kcal_using_ap_low_fermentation   | name_focal), 
                      data = data[complete.cases(data[,c("pos_neg_dn","fai","prev_day_total_kcal_using_ap_low_fermentation")]),],
                      family=bernoulli(),iter = 3000, 
                      prior=c(prior(student_t(4,0,1.5),class="Intercept"),
                              prior(student_t(4,0,1.5),class="b"),
                              prior(student_t(4,0,1.5),class="sds"),
                              prior(student_t(4,0,1.5),class="sd")),
                      control = list(max_treedepth = 10, adapt_delta = .99))

summary(Ketone_PREV_KCAL, digits=9)
plot(conditional_effects(Ketone_PREV_KCAL,spaghetti=FALSE),points=TRUE) ##plots the regression lines/conditional effects
pp_check(Ketone_PREV_KCAL)
loo_R2(Ketone_PREV_KCAL)
bayes_R2(Ketone_PREV_KCAL)
gc()


Ketone_PREV_NPE<-brm(pos_neg_dn~s(fai) + class_focal + s(prev_day_total_kcal_npe_low_fermentation) + (1 +  fai + prev_day_total_kcal_npe_low_fermentation   | name_focal), 
                     data = data[complete.cases(data[,c("pos_neg_dn","fai","prev_day_total_kcal_npe_low_fermentation")]),],
                     family=bernoulli(),iter = 3000, 
                     prior=c(prior(student_t(4,0,1.5),class="Intercept"),
                             prior(student_t(4,0,1.5),class="b"),
                             prior(student_t(4,0,1.5),class="sds"),
                             prior(student_t(4,0,1.5),class="sd")),
                     control = list(max_treedepth = 10, adapt_delta = .99))


summary(Ketone_PREV_NPE, digits=9)
plot(conditional_effects(Ketone_PREV_NPE,spaghetti=FALSE),points=TRUE) ##plots the regression lines/conditional effects
pp_check(Ketone_PREV_NPE)
loo_R2(Ketone_PREV_NPE)
bayes_R2(Ketone_PREV_NPE)
gc()


data$prev_day_npe_ab=data$prev_day_total_kcal_npe_low_fermentation/data$prev_day_ap_kcal
Ketone_PREV_npe_ab<-brm(pos_neg_dn~s(fai) + class_focal + s(prev_day_npe_ab) + (1 +  fai + prev_day_total_kcal_npe_low_fermentation   | name_focal), 
                     data = data[complete.cases(data[,c("pos_neg_dn","fai","prev_day_npe_ab")]),],
                     family=bernoulli(),iter = 3000, 
                     prior=c(prior(student_t(4,0,1.5),class="Intercept"),
                             prior(student_t(4,0,1.5),class="b"),
                             prior(student_t(4,0,1.5),class="sds"),
                             prior(student_t(4,0,1.5),class="sd")),
                     control = list(max_treedepth = 10, adapt_delta = .99))


summary(Ketone_PREV_npe_ab, digits=9)
plot(conditional_effects(Ketone_PREV_npe_ab,spaghetti=FALSE),points=TRUE) ##plots the regression lines/conditional effects
pp_check(Ketone_PREV_npe_ab)
loo_R2(Ketone_PREV_npe_ab)
bayes_R2(Ketone_PREV_npe_ab)
gc()




data$PREV_lipid_ap=data$prev_day_lipid_kcal/data$prev_day_ap_kcal
Ketone_PREV_lipid_ap<-brm(pos_neg_dn~s(fai) + class_focal + s((PREV_lipid_ap)) + (1 +  fai + (PREV_lipid_ap)   | name_focal), 
                          data = data[complete.cases(data[,c("pos_neg_dn","fai","PREV_lipid_ap")]),],
                          family=bernoulli(),iter = 3000, 
                          prior=c(prior(student_t(4,0,1.5),class="Intercept"),
                                  prior(student_t(4,0,1.5),class="b"),
                                  prior(student_t(4,0,1.5),class="sds"),
                                  prior(student_t(4,0,1.5),class="sd")),
                          control = list(max_treedepth = 10, adapt_delta = .99))


summary(Ketone_PREV_lipid_ap, digits=9)
plot(conditional_effects(Ketone_PREV_lipid_ap,spaghetti=TRUE),points=TRUE) ##plots the regression lines/conditional effects
pp_check(Ketone_PREV_lipid_ap)
loo_R2(Ketone_PREV_lipid_ap)
bayes_R2(Ketone_PREV_lipid_ap)
gc()



data$PREV_lipid_tnc=data$prev_day_lipid_kcal/data$prev_day_tnc_kcal
Ketone_PREV_lipid_tnc<-brm(pos_neg_dn~s(fai) + class_focal + s((PREV_lipid_tnc)) + (1 +  fai + (PREV_lipid_tnc)   | name_focal), 
                           data = data[complete.cases(data[,c("pos_neg_dn","fai","PREV_lipid_tnc")]),],
                           family=bernoulli(),iter = 3000, 
                           prior=c(prior(student_t(4,0,1.5),class="Intercept"),
                                   prior(student_t(4,0,1.5),class="b"),
                                   prior(student_t(4,0,1.5),class="sds"),
                                   prior(student_t(4,0,1.5),class="sd")),
                           control = list(max_treedepth = 10, adapt_delta = .99))


summary(Ketone_PREV_lipid_tnc, digits=9)
plot(conditional_effects(Ketone_PREV_lipid_tnc,spaghetti=FALSE),points=TRUE) ##plots the regression lines/conditional effects
pp_check(Ketone_PREV_lipid_tnc)
loo_R2(Ketone_PREV_lipid_tnc)
bayes_R2(Ketone_PREV_lipid_tnc)
gc()



data$PREV_lipid_ndf=data$prev_day_lipid_kcal/data$prev_day_ndf_kcal_low
Ketone_PREV_lipid_ndf<-brm(pos_neg_dn~s(fai) + class_focal + s((PREV_lipid_ndf)) + (1 +  fai + (PREV_lipid_ndf)   | name_focal), 
                           data = data[complete.cases(data[,c("pos_neg_dn","fai","PREV_lipid_ndf")]),],
                           family=bernoulli(),iter = 3000, 
                           prior=c(prior(student_t(4,0,1.5),class="Intercept"),
                                   prior(student_t(4,0,1.5),class="b"),
                                   prior(student_t(4,0,1.5),class="sds"),
                                   prior(student_t(4,0,1.5),class="sd")),
                           control = list(max_treedepth = 10, adapt_delta = .99))


summary(Ketone_PREV_lipid_ndf, digits=9)
plot(conditional_effects(Ketone_PREV_lipid_ndf,spaghetti=FALSE),points=TRUE) ##plots the regression lines/conditional effects
pp_check(Ketone_PREV_lipid_ndf)
loo_R2(Ketone_PREV_lipid_ndf)
bayes_R2(Ketone_PREV_lipid_ndf)
gc()





data$PREV_tnc_ap=data$prev_day_tnc_kcal/data$prev_day_ap_kcal
Ketone_PREV_tnc_ap<-brm(pos_neg_dn~s(fai) + class_focal + s((PREV_tnc_ap)) + (1 +  fai + (PREV_tnc_ap)   | name_focal), 
                        data = data[complete.cases(data[,c("pos_neg_dn","fai","PREV_tnc_ap")]),],
                        family=bernoulli(),iter = 3000, 
                        prior=c(prior(student_t(4,0,1.5),class="Intercept"),
                                prior(student_t(4,0,1.5),class="b"),
                                prior(student_t(4,0,1.5),class="sds"),
                                prior(student_t(4,0,1.5),class="sd")),
                        control = list(max_treedepth = 10, adapt_delta = .99))


summary(Ketone_PREV_tnc_ap, digits=9)
plot(conditional_effects(Ketone_PREV_tnc_ap,spaghetti=FALSE),points=TRUE) ##plots the regression lines/conditional effects
pp_check(Ketone_PREV_tnc_ap)
loo_R2(Ketone_PREV_tnc_ap)
bayes_R2(Ketone_PREV_tnc_ap)
gc()



data$PREV_tnc_ndf=data$prev_day_tnc_kcal/data$prev_day_ndf_kcal_low
Ketone_PREV_tnc_ndf<-brm(pos_neg_dn~s(fai) + class_focal + s((PREV_tnc_ndf)) + (1 +  fai + (PREV_tnc_ndf)   | name_focal), 
                         data = data[complete.cases(data[,c("pos_neg_dn","fai","PREV_tnc_ndf")]),],
                         family=bernoulli(),iter = 3000, 
                         prior=c(prior(student_t(4,0,1.5),class="Intercept"),
                                 prior(student_t(4,0,1.5),class="b"),
                                 prior(student_t(4,0,1.5),class="sds"),
                                 prior(student_t(4,0,1.5),class="sd")),
                         control = list(max_treedepth = 10, adapt_delta = .99))


summary(Ketone_PREV_tnc_ndf, digits=9)
plot(conditional_effects(Ketone_PREV_tnc_ndf,spaghetti=FALSE),points=TRUE) ##plots the regression lines/conditional effects
pp_check(Ketone_PREV_tnc_ndf)
loo_R2(Ketone_PREV_tnc_ndf)
bayes_R2(Ketone_PREV_tnc_ndf)
gc()



data$PREV_ap_ndf=data$prev_day_ap_kcal/data$prev_day_ndf_kcal_low
Ketone_PREV_ap_ndf<-brm(pos_neg_dn~s(fai) + class_focal + s((PREV_ap_ndf)) + (1 +  fai + (PREV_ap_ndf)   | name_focal), 
                        data = data[complete.cases(data[,c("pos_neg_dn","fai","PREV_ap_ndf")]),],
                        family=bernoulli(),iter = 3000, 
                        prior=c(prior(student_t(4,0,1.5),class="Intercept"),
                                prior(student_t(4,0,1.5),class="b"),
                                prior(student_t(4,0,1.5),class="sds"),
                                prior(student_t(4,0,1.5),class="sd")),
                        control = list(max_treedepth = 10, adapt_delta = .99))


summary(Ketone_PREV_ap_ndf, digits=9)
plot(conditional_effects(Ketone_PREV_ap_ndf,spaghetti=FALSE)) ##plots the regression lines/conditional effects
pp_check(Ketone_PREV_ap_ndf)
loo_R2(Ketone_PREV_ap_ndf)
bayes_R2(Ketone_PREV_ap_ndf)
gc()




#####UCP MODELS######

ucp_sgcor_KCAL<-brm(ucp_sgcor~s(fai) + class_focal + s(total_kcal_using_ap_low_fermentation) + s(ucp_time_collected) + (1 +  fai + total_kcal_using_ap_low_fermentation + ucp_time_collected  | name_focal), 
                         data = data[complete.cases(data[,c("ucp_sgcor","ucp_time_collected","fai")]),],
                         family=lognormal(),iter = 2000, 
                         prior=c(prior(normal(0,10),class="Intercept"),
                                 prior(normal(0,10),class="b"),
                                 prior(normal(0,10),class="sigma"),
                                 prior(normal(0,10),class="sd")),
                         control = list(max_treedepth = 11, adapt_delta = .999))

summary(ucp_sgcor_KCAL, digits=9)
plot(conditional_effects(ucp_sgcor_KCAL,spaghetti=FALSE)) ##plots the regression lines/conditional effects
pp_check(ucp_sgcor_KCAL)
loo_R2(ucp_sgcor_KCAL)
bayes_R2(ucp_sgcor_KCAL)




ucp_sgcor_NPE<-brm(ucp_sgcor~s(fai) + class_focal + s(total_kcal_npe_low_fermentation) + s(ucp_time_collected) + (1 +  fai + total_kcal_npe_low_fermentation + ucp_time_collected  | name_focal), 
                    data = data[complete.cases(data[,c("ucp_sgcor","ucp_time_collected","fai","total_kcal_npe_low_fermentation")]),],
                    family=lognormal(),iter = 2000, 
                    prior=c(prior(normal(0,10),class="Intercept"),
                            prior(normal(0,10),class="b"),
                            prior(normal(0,10),class="sigma"),
                            prior(normal(0,10),class="sd")),
                    control = list(max_treedepth = 11, adapt_delta = .999))

summary(ucp_sgcor_NPE, digits=9)
plot(conditional_effects(ucp_sgcor_NPE,spaghetti=FALSE)) ##plots the regression lines/conditional effects
pp_check(ucp_sgcor_NPE)
loo_R2(ucp_sgcor_NPE)
bayes_R2(ucp_sgcor_NPE)





ucp_sgcor_AP<-brm(ucp_sgcor~s(fai) + class_focal + s(ap_kcal) + s(ucp_time_collected) + (1 +  fai + ap_kcal + ucp_time_collected  | name_focal), 
                   data = data[complete.cases(data[,c("ucp_sgcor","ucp_time_collected","fai","ap_kcal")]),],
                   family=lognormal(),iter = 2000, 
                   prior=c(prior(normal(0,10),class="Intercept"),
                           prior(normal(0,10),class="b"),
                           prior(normal(0,10),class="sigma"),
                           prior(normal(0,10),class="sd")),
                   control = list(max_treedepth = 11, adapt_delta = .999))

summary(ucp_sgcor_AP, digits=9)
plot(conditional_effects(ucp_sgcor_AP,spaghetti=FALSE)) ##plots the regression lines/conditional effects
pp_check(ucp_sgcor_AP)
loo_R2(ucp_sgcor_AP)
bayes_R2(ucp_sgcor_AP)



ucp_sgcor_Lipid<-brm(ucp_sgcor~s(fai) + class_focal + s(lipid_kcal) + s(ucp_time_collected) + (1 +  fai + lipid_kcal + ucp_time_collected  | name_focal), 
                  data = data[complete.cases(data[,c("ucp_sgcor","ucp_time_collected","fai","lipid_kcal")]),],
                  family=lognormal(),iter = 2000, 
                  prior=c(prior(normal(0,10),class="Intercept"),
                          prior(normal(0,10),class="b"),
                          prior(normal(0,10),class="sigma"),
                          prior(normal(0,10),class="sd")),
                  control = list(max_treedepth = 11, adapt_delta = .999))

summary(ucp_sgcor_Lipid, digits=9)
plot(conditional_effects(ucp_sgcor_Lipid,spaghetti=FALSE)) ##plots the regression lines/conditional effects
pp_check(ucp_sgcor_Lipid)
loo_R2(ucp_sgcor_Lipid)
bayes_R2(ucp_sgcor_Lipid)




ucp_sgcor_TNC<-brm(ucp_sgcor~s(fai) + class_focal + s(tnc_kcal) + s(ucp_time_collected) + (1 +  fai + tnc_kcal + ucp_time_collected  | name_focal), 
                     data = data[complete.cases(data[,c("ucp_sgcor","ucp_time_collected","fai","tnc_kcal")]),],
                     family=lognormal(),iter = 2000, 
                     prior=c(prior(normal(0,10),class="Intercept"),
                             prior(normal(0,10),class="b"),
                             prior(normal(0,10),class="sigma"),
                             prior(normal(0,10),class="sd")),
                     control = list(max_treedepth = 11, adapt_delta = .999))

summary(ucp_sgcor_TNC, digits=9)
plot(conditional_effects(ucp_sgcor_TNC,spaghetti=FALSE)) ##plots the regression lines/conditional effects
pp_check(ucp_sgcor_TNC)
loo_R2(ucp_sgcor_TNC)
bayes_R2(ucp_sgcor_TNC)






ucp_sgcor_NDF<-brm(ucp_sgcor~s(fai) + class_focal + s(ndf_kcal_low) + s(ucp_time_collected) + (1 +  fai + ndf_kcal_low + ucp_time_collected  | name_focal), 
                   data = data[complete.cases(data[,c("ucp_sgcor","ucp_time_collected","fai","ndf_kcal_low")]),],
                   family=lognormal(),iter = 2000, 
                   prior=c(prior(normal(0,10),class="Intercept"),
                           prior(normal(0,10),class="b"),
                           prior(normal(0,10),class="sigma"),
                           prior(normal(0,10),class="sd")),
                   control = list(max_treedepth = 11, adapt_delta = .999))

summary(ucp_sgcor_NDF, digits=9)
plot(conditional_effects(ucp_sgcor_NDF,spaghetti=FALSE)) ##plots the regression lines/conditional effects
pp_check(ucp_sgcor_NDF)
loo_R2(ucp_sgcor_NDF)
bayes_R2(ucp_sgcor_NDF)




ucp_sgcor_npe_ap<-brm(ucp_sgcor~s(fai) + class_focal + s(npe_to_ap) + s(ucp_time_collected) + (1 +  fai + npe_to_ap + ucp_time_collected  | name_focal), 
                   data = data[complete.cases(data[,c("ucp_sgcor","ucp_time_collected","fai","npe_to_ap")]),],
                   family=lognormal(),iter = 2000, 
                   prior=c(prior(normal(0,10),class="Intercept"),
                           prior(normal(0,10),class="b"),
                           prior(normal(0,10),class="sigma"),
                           prior(normal(0,10),class="sd")),
                   control = list(max_treedepth = 11, adapt_delta = .999))

summary(ucp_sgcor_npe_ap, digits=9)
plot(conditional_effects(ucp_sgcor_npe_ap,spaghetti=FALSE)) ##plots the regression lines/conditional effects
pp_check(ucp_sgcor_npe_ap)
loo_R2(ucp_sgcor_npe_ap)
bayes_R2(ucp_sgcor_npe_ap)



ucp_sgcor_lipid_ap<-brm(ucp_sgcor~s(fai) + class_focal + s(lipid_to_ap) + s(ucp_time_collected) + (1 +  fai + lipid_to_ap + ucp_time_collected  | name_focal), 
                      data = data[complete.cases(data[,c("ucp_sgcor","ucp_time_collected","fai","lipid_to_ap")]),],
                      family=lognormal(),iter = 2000, 
                      prior=c(prior(normal(0,10),class="Intercept"),
                              prior(normal(0,10),class="b"),
                              prior(normal(0,10),class="sigma"),
                              prior(normal(0,10),class="sd")),
                      control = list(max_treedepth = 11, adapt_delta = .999))

summary(ucp_sgcor_lipid_ap, digits=9)
plot(conditional_effects(ucp_sgcor_lipid_ap,spaghetti=FALSE)) ##plots the regression lines/conditional effects
pp_check(ucp_sgcor_lipid_ap)
loo_R2(ucp_sgcor_lipid_ap)
bayes_R2(ucp_sgcor_lipid_ap)



ucp_sgcor_lipid_tnc<-brm(ucp_sgcor~s(fai) + class_focal + s(lipid_to_tnc) + s(ucp_time_collected) + (1 +  fai + lipid_to_tnc + ucp_time_collected  | name_focal), 
                        data = data[complete.cases(data[,c("ucp_sgcor","ucp_time_collected","fai","lipid_to_tnc")]),],
                        family=lognormal(),iter = 2000, 
                        prior=c(prior(normal(0,10),class="Intercept"),
                                prior(normal(0,10),class="b"),
                                prior(normal(0,10),class="sigma"),
                                prior(normal(0,10),class="sd")),
                        control = list(max_treedepth = 11, adapt_delta = .999))

summary(ucp_sgcor_lipid_tnc, digits=9)
plot(conditional_effects(ucp_sgcor_lipid_tnc,spaghetti=FALSE)) ##plots the regression lines/conditional effects
pp_check(ucp_sgcor_lipid_tnc)
loo_R2(ucp_sgcor_lipid_tnc)
bayes_R2(ucp_sgcor_lipid_tnc)




ucp_sgcor_lipid_ndf<-brm(ucp_sgcor~s(fai) + class_focal + s(lipid_to_ndf_low) + s(ucp_time_collected) + (1 +  fai + lipid_to_ndf_low + ucp_time_collected  | name_focal), 
                         data = data[complete.cases(data[,c("ucp_sgcor","ucp_time_collected","fai","lipid_to_ndf_low")]),],
                         family=lognormal(),iter = 2000, 
                         prior=c(prior(normal(0,10),class="Intercept"),
                                 prior(normal(0,10),class="b"),
                                 prior(normal(0,10),class="sigma"),
                                 prior(normal(0,10),class="sd")),
                         control = list(max_treedepth = 11, adapt_delta = .999))

summary(ucp_sgcor_lipid_ndf, digits=9)
plot(conditional_effects(ucp_sgcor_lipid_ndf,spaghetti=FALSE)) ##plots the regression lines/conditional effects
pp_check(ucp_sgcor_lipid_ndf)
loo_R2(ucp_sgcor_lipid_ndf)
bayes_R2(ucp_sgcor_lipid_ndf)





ucp_sgcor_tnc_ap<-brm(ucp_sgcor~s(fai) + class_focal + s(tnc_to_ap) + s(ucp_time_collected) + (1 +  fai + tnc_to_ap + ucp_time_collected  | name_focal), 
                         data = data[complete.cases(data[,c("ucp_sgcor","ucp_time_collected","fai","tnc_to_ap")]),],
                         family=lognormal(),iter = 2000, 
                         prior=c(prior(normal(0,10),class="Intercept"),
                                 prior(normal(0,10),class="b"),
                                 prior(normal(0,10),class="sigma"),
                                 prior(normal(0,10),class="sd")),
                         control = list(max_treedepth = 11, adapt_delta = .999))

summary(ucp_sgcor_tnc_ap, digits=9)
plot(conditional_effects(ucp_sgcor_tnc_ap,spaghetti=FALSE)) ##plots the regression lines/conditional effects
pp_check(ucp_sgcor_tnc_ap)
loo_R2(ucp_sgcor_tnc_ap)
bayes_R2(ucp_sgcor_tnc_ap)






ucp_sgcor_tnc_ndf<-brm(ucp_sgcor~s(fai) + class_focal + s(tnc_to_ndf_low) + s(ucp_time_collected) + (1 +  fai + tnc_to_ndf_low + ucp_time_collected  | name_focal), 
                      data = data[complete.cases(data[,c("ucp_sgcor","ucp_time_collected","fai","tnc_to_ndf_low")]),],
                      family=lognormal(),iter = 2000, 
                      prior=c(prior(normal(0,10),class="Intercept"),
                              prior(normal(0,10),class="b"),
                              prior(normal(0,10),class="sigma"),
                              prior(normal(0,10),class="sd")),
                      control = list(max_treedepth = 11, adapt_delta = .999))

summary(ucp_sgcor_tnc_ndf, digits=9)
plot(conditional_effects(ucp_sgcor_tnc_ndf,spaghetti=FALSE)) ##plots the regression lines/conditional effects
pp_check(ucp_sgcor_tnc_ndf)
loo_R2(ucp_sgcor_tnc_ndf)
bayes_R2(ucp_sgcor_tnc_ndf)




ucp_sgcor_ap_ndf<-brm(ucp_sgcor~s(fai) + class_focal + s(ap_to_ndf_low) + s(ucp_time_collected) + (1 +  fai + ap_to_ndf_low + ucp_time_collected  | name_focal), 
                       data = data[complete.cases(data[,c("ucp_sgcor","ucp_time_collected","fai","ap_to_ndf_low")]),],
                       family=lognormal(),iter = 2000, 
                       prior=c(prior(normal(0,10),class="Intercept"),
                               prior(normal(0,10),class="b"),
                               prior(normal(0,10),class="sigma"),
                               prior(normal(0,10),class="sd")),
                       control = list(max_treedepth = 11, adapt_delta = .999))

summary(ucp_sgcor_ap_ndf, digits=9)
plot(conditional_effects(ucp_sgcor_ap_ndf,spaghetti=FALSE)) ##plots the regression lines/conditional effects
pp_check(ucp_sgcor_ap_ndf)
loo_R2(ucp_sgcor_ap_ndf)
bayes_R2(ucp_sgcor_ap_ndf)





ucp_sgcor_PREV_AP<-brm(ucp_sgcor~s(fai) + class_focal + s(prev_day_ap_kcal) + s(ucp_time_collected) + (1 +  fai + prev_day_lipid_kcal + ucp_time_collected  | name_focal), 
                      data = data[complete.cases(data[,c("ucp_sgcor","ucp_time_collected","fai","prev_day_ap_kcal")]),],
                      family=lognormal(),iter = 2000, 
                      prior=c(prior(normal(0,10),class="Intercept"),
                              prior(normal(0,10),class="b"),
                              prior(normal(0,10),class="sigma"),
                              prior(normal(0,10),class="sd")),
                      control = list(max_treedepth = 11, adapt_delta = .999))

summary(ucp_sgcor_PREV_AP, digits=9)
plot(conditional_effects(ucp_sgcor_PREV_AP,spaghetti=FALSE)) ##plots the regression lines/conditional effects
pp_check(ucp_sgcor_PREV_AP)
loo_R2(ucp_sgcor_PREV_AP)
bayes_R2(ucp_sgcor_PREV_AP)





ucp_sgcor_PREV_lipid<-brm(ucp_sgcor~s(fai) + class_focal + s(prev_day_lipid_kcal) + s(ucp_time_collected) + (1 +  fai + prev_day_lipid_kcal + ucp_time_collected  | name_focal), 
                       data = data[complete.cases(data[,c("ucp_sgcor","ucp_time_collected","fai","prev_day_lipid_kcal")]),],
                       family=lognormal(),iter = 2000, 
                       prior=c(prior(normal(0,10),class="Intercept"),
                               prior(normal(0,10),class="b"),
                               prior(normal(0,10),class="sigma"),
                               prior(normal(0,10),class="sd")),
                       control = list(max_treedepth = 11, adapt_delta = .999))

summary(ucp_sgcor_PREV_lipid, digits=9)
plot(conditional_effects(ucp_sgcor_PREV_lipid,spaghetti=FALSE)) ##plots the regression lines/conditional effects
pp_check(ucp_sgcor_PREV_lipid)
loo_R2(ucp_sgcor_PREV_lipid)
bayes_R2(ucp_sgcor_PREV_lipid)



ucp_sgcor_PREV_TNC<-brm(ucp_sgcor~s(fai) + class_focal + s(prev_day_tnc_kcal) + s(ucp_time_collected) + (1 +  fai + prev_day_tnc_kcal + ucp_time_collected  | name_focal), 
                          data = data[complete.cases(data[,c("ucp_sgcor","ucp_time_collected","fai","prev_day_tnc_kcal")]),],
                          family=lognormal(),iter = 2000, 
                          prior=c(prior(normal(0,10),class="Intercept"),
                                  prior(normal(0,10),class="b"),
                                  prior(normal(0,10),class="sigma"),
                                  prior(normal(0,10),class="sd")),
                          control = list(max_treedepth = 11, adapt_delta = .999))

summary(ucp_sgcor_PREV_TNC, digits=9)
plot(conditional_effects(ucp_sgcor_PREV_TNC,spaghetti=FALSE)) ##plots the regression lines/conditional effects
pp_check(ucp_sgcor_PREV_TNC)
loo_R2(ucp_sgcor_PREV_TNC)
bayes_R2(ucp_sgcor_PREV_TNC)




ucp_sgcor_PREV_NDF<-brm(ucp_sgcor~s(fai) + class_focal + s(prev_day_ndf_kcal_low) + s(ucp_time_collected) + (1 +  fai + prev_day_ndf_kcal_low + ucp_time_collected  | name_focal), 
                        data = data[complete.cases(data[,c("ucp_sgcor","ucp_time_collected","fai","prev_day_ndf_kcal_low")]),],
                        family=lognormal(),iter = 2000, 
                        prior=c(prior(normal(0,10),class="Intercept"),
                                prior(normal(0,10),class="b"),
                                prior(normal(0,10),class="sigma"),
                                prior(normal(0,10),class="sd")),
                        control = list(max_treedepth = 11, adapt_delta = .999))

summary(ucp_sgcor_PREV_NDF, digits=9)
plot(conditional_effects(ucp_sgcor_PREV_NDF,spaghetti=FALSE)) ##plots the regression lines/conditional effects
pp_check(ucp_sgcor_PREV_NDF)
loo_R2(ucp_sgcor_PREV_NDF)
bayes_R2(ucp_sgcor_PREV_NDF)





ucp_sgcor_PREV_NDF<-brm(ucp_sgcor~s(fai) + class_focal + s(prev_day_ndf_kcal_low) + s(ucp_time_collected) + (1 +  fai + prev_day_ndf_kcal_low + ucp_time_collected  | name_focal), 
                        data = data[complete.cases(data[,c("ucp_sgcor","ucp_time_collected","fai","prev_day_ndf_kcal_low")]),],
                        family=lognormal(),iter = 2000, 
                        prior=c(prior(normal(0,10),class="Intercept"),
                                prior(normal(0,10),class="b"),
                                prior(normal(0,10),class="sigma"),
                                prior(normal(0,10),class="sd")),
                        control = list(max_treedepth = 11, adapt_delta = .999))

summary(ucp_sgcor_PREV_NDF, digits=9)
plot(conditional_effects(ucp_sgcor_PREV_NDF,spaghetti=FALSE)) ##plots the regression lines/conditional effects
pp_check(ucp_sgcor_PREV_NDF)
loo_R2(ucp_sgcor_PREV_NDF)
bayes_R2(ucp_sgcor_PREV_NDF)



ucp_sgcor_PREV_KCAL<-brm(ucp_sgcor~s(fai) + class_focal + s(prev_day_total_kcal_using_ap_low_fermentation) + s(ucp_time_collected) + (1 +  fai + prev_day_total_kcal_using_ap_low_fermentation + ucp_time_collected  | name_focal), 
                        data = data[complete.cases(data[,c("ucp_sgcor","ucp_time_collected","fai","prev_day_total_kcal_using_ap_low_fermentation")]),],
                        family=lognormal(),iter = 2000, 
                        prior=c(prior(normal(0,10),class="Intercept"),
                                prior(normal(0,10),class="b"),
                                prior(normal(0,10),class="sigma"),
                                prior(normal(0,10),class="sd")),
                        control = list(max_treedepth = 11, adapt_delta = .999))

summary(ucp_sgcor_PREV_KCAL, digits=9)
plot(conditional_effects(ucp_sgcor_PREV_KCAL,spaghetti=FALSE)) ##plots the regression lines/conditional effects
pp_check(ucp_sgcor_PREV_KCAL)
loo_R2(ucp_sgcor_PREV_KCAL)
bayes_R2(ucp_sgcor_PREV_KCAL)





ucp_sgcor_PREV_NPE<-brm(ucp_sgcor~s(fai) + class_focal + s(prev_day_total_kcal_npe_low_fermentation) + s(ucp_time_collected) + (1 +  fai + prev_day_total_kcal_npe_low_fermentation + ucp_time_collected  | name_focal), 
                         data = data[complete.cases(data[,c("ucp_sgcor","ucp_time_collected","fai","prev_day_total_kcal_npe_low_fermentation")]),],
                         family=lognormal(),iter = 2000, 
                         prior=c(prior(normal(0,10),class="Intercept"),
                                 prior(normal(0,10),class="b"),
                                 prior(normal(0,10),class="sigma"),
                                 prior(normal(0,10),class="sd")),
                         control = list(max_treedepth = 11, adapt_delta = .999))

summary(ucp_sgcor_PREV_NPE, digits=9)
plot(conditional_effects(ucp_sgcor_PREV_NPE,spaghetti=FALSE)) ##plots the regression lines/conditional effects
pp_check(ucp_sgcor_PREV_NPE)
loo_R2(ucp_sgcor_PREV_NPE)
bayes_R2(ucp_sgcor_PREV_NPE)



ucp_sgcor_PREV_lipid_ap<-brm(ucp_sgcor~s(fai) + class_focal + s((PREV_lipid_ap)) + s(ucp_time_collected) + (1 +  fai + PREV_lipid_ap + ucp_time_collected  | name_focal), 
                        data = data[complete.cases(data[,c("ucp_sgcor","ucp_time_collected","fai","PREV_lipid_ap")]),],
                        family=lognormal(),iter = 2000, 
                        prior=c(prior(normal(0,10),class="Intercept"),
                                prior(normal(0,10),class="b"),
                                prior(normal(0,10),class="sigma"),
                                prior(normal(0,10),class="sd")),
                        control = list(max_treedepth = 11, adapt_delta = .999))

summary(ucp_sgcor_PREV_lipid_ap, digits=9)
plot(conditional_effects(ucp_sgcor_PREV_lipid_ap,spaghetti=FALSE)) ##plots the regression lines/conditional effects
pp_check(ucp_sgcor_PREV_lipid_ap)
loo_R2(ucp_sgcor_PREV_lipid_ap)
bayes_R2(ucp_sgcor_PREV_lipid_ap)





ucp_sgcor_PREV_lipid_tnc<-brm(ucp_sgcor~s(fai) + class_focal + s((PREV_lipid_tnc)) + s(ucp_time_collected) + (1 +  fai + (PREV_lipid_tnc) + ucp_time_collected  | name_focal), 
                             data = data[complete.cases(data[,c("ucp_sgcor","ucp_time_collected","fai","PREV_lipid_tnc")]),],
                             family=lognormal(),iter = 2000, 
                             prior=c(prior(normal(0,10),class="Intercept"),
                                     prior(normal(0,10),class="b"),
                                     prior(normal(0,10),class="sigma"),
                                     prior(normal(0,10),class="sd")),
                             control = list(max_treedepth = 11, adapt_delta = .999))

summary(ucp_sgcor_PREV_lipid_tnc, digits=9)
plot(conditional_effects(ucp_sgcor_PREV_lipid_tnc,spaghetti=FALSE)) ##plots the regression lines/conditional effects
pp_check(ucp_sgcor_PREV_lipid_tnc)
loo_R2(ucp_sgcor_PREV_lipid_tnc)
bayes_R2(ucp_sgcor_PREV_lipid_tnc)




data$PREV_lipid_ndf=data$prev_day_lipid_kcal/data$prev_day_ndf_kcal_low
ucp_sgcor_PREV_lipid_ndf<-brm(ucp_sgcor~s(fai) + class_focal + s((PREV_lipid_ndf)) + s(ucp_time_collected) + (1 +  fai + PREV_lipid_ndf + ucp_time_collected  | name_focal), 
                              data = data[complete.cases(data[,c("ucp_sgcor","ucp_time_collected","fai","PREV_lipid_ndf")]),],
                              family=lognormal(),iter = 2000, 
                              prior=c(prior(normal(0,10),class="Intercept"),
                                      prior(normal(0,10),class="b"),
                                      prior(normal(0,10),class="sigma"),
                                      prior(normal(0,10),class="sd")),
                              control = list(max_treedepth = 11, adapt_delta = .999))

summary(ucp_sgcor_PREV_lipid_ndf, digits=9)
plot(conditional_effects(ucp_sgcor_PREV_lipid_ndf,spaghetti=FALSE)) ##plots the regression lines/conditional effects
pp_check(ucp_sgcor_PREV_lipid_ndf)
loo_R2(ucp_sgcor_PREV_lipid_ndf)
bayes_R2(ucp_sgcor_PREV_lipid_ndf)



ucp_sgcor_PREV_tnc_ap<-brm(ucp_sgcor~s(fai) + class_focal + s(PREV_tnc_ap) + s(ucp_time_collected) + (1 +  fai + PREV_tnc_ap + ucp_time_collected  | name_focal), 
                              data = data[complete.cases(data[,c("ucp_sgcor","ucp_time_collected","fai","PREV_tnc_ap")]),],
                              family=lognormal(),iter = 2000, 
                              prior=c(prior(normal(0,10),class="Intercept"),
                                      prior(normal(0,10),class="b"),
                                      prior(normal(0,10),class="sigma"),
                                      prior(normal(0,10),class="sd")),
                              control = list(max_treedepth = 11, adapt_delta = .999))

summary(ucp_sgcor_PREV_tnc_ap, digits=9)
plot(conditional_effects(ucp_sgcor_PREV_tnc_ap,spaghetti=FALSE)) ##plots the regression lines/conditional effects
pp_check(ucp_sgcor_PREV_tnc_ap)
loo_R2(ucp_sgcor_PREV_tnc_ap)
bayes_R2(ucp_sgcor_PREV_tnc_ap)





ucp_sgcor_PREV_tnc_ndf<-brm(ucp_sgcor~s(fai) + class_focal + s(PREV_tnc_ndf) + s(ucp_time_collected) + (1 +  fai + PREV_tnc_ndf + ucp_time_collected  | name_focal), 
                           data = data[complete.cases(data[,c("ucp_sgcor","ucp_time_collected","fai","PREV_tnc_ndf")]),],
                           family=lognormal(),iter = 2000, 
                           prior=c(prior(normal(0,10),class="Intercept"),
                                   prior(normal(0,10),class="b"),
                                   prior(normal(0,10),class="sigma"),
                                   prior(normal(0,10),class="sd")),
                           control = list(max_treedepth = 11, adapt_delta = .999))

summary(ucp_sgcor_PREV_tnc_ndf, digits=9)
plot(conditional_effects(ucp_sgcor_PREV_tnc_ndf,spaghetti=FALSE)) ##plots the regression lines/conditional effects
pp_check(ucp_sgcor_PREV_tnc_ndf)
loo_R2(ucp_sgcor_PREV_tnc_ndf)
bayes_R2(ucp_sgcor_PREV_tnc_ndf)




ucp_sgcor_PREV_ap_ndf<-brm(ucp_sgcor~s(fai) + class_focal + s(PREV_ap_ndf) + s(ucp_time_collected) + (1 +  fai + PREV_ap_ndf + ucp_time_collected  | name_focal), 
                            data = data[complete.cases(data[,c("ucp_sgcor","ucp_time_collected","fai","PREV_ap_ndf")]),],
                            family=lognormal(),iter = 2000, 
                            prior=c(prior(normal(0,10),class="Intercept"),
                                    prior(normal(0,10),class="b"),
                                    prior(normal(0,10),class="sigma"),
                                    prior(normal(0,10),class="sd")),
                            control = list(max_treedepth = 11, adapt_delta = .999))

summary(ucp_sgcor_PREV_ap_ndf, digits=9)
plot(conditional_effects(ucp_sgcor_PREV_ap_ndf,spaghetti=FALSE)) ##plots the regression lines/conditional effects
pp_check(ucp_sgcor_PREV_ap_ndf)
loo_R2(ucp_sgcor_PREV_ap_ndf)
bayes_R2(ucp_sgcor_PREV_ap_ndf)



ucp_sgcor_PREV_npe_ab<-brm(ucp_sgcor~s(fai) + class_focal + s(prev_day_npe_ab) + s(ucp_time_collected) + (1 +  fai + prev_day_npe_ab + ucp_time_collected  | name_focal), 
                           data = data[complete.cases(data[,c("ucp_sgcor","ucp_time_collected","fai","prev_day_npe_ab")]),],
                           family=lognormal(),iter = 2000, 
                           prior=c(prior(normal(0,10),class="Intercept"),
                                   prior(normal(0,10),class="b"),
                                   prior(normal(0,10),class="sigma"),
                                   prior(normal(0,10),class="sd")),
                           control = list(max_treedepth = 11, adapt_delta = .999))

summary(ucp_sgcor_PREV_npe_ab, digits=9)
plot(conditional_effects(ucp_sgcor_PREV_npe_ab,spaghetti=FALSE)) ##plots the regression lines/conditional effects
pp_check(ucp_sgcor_PREV_npe_ab)
loo_R2(ucp_sgcor_PREV_npe_ab)
bayes_R2(ucp_sgcor_PREV_npe_ab)


###Isotope models####

dn15_result_FAI_model<-brm(dn15_result~s(fai) + s(isotope_time_collected) + (1 +  fai + isotope_time_collected  | name_focal), 
                           data = data[complete.cases(data[,c("dn15_result","isotope_time_collected","fai")]),],
                           family=gaussian(),iter = 2000, 
                           prior=c(prior(normal(0,10),class="Intercept"),
                                   prior(normal(0,10),class="b", coef =sfai_1),
                                   prior(normal(0,10),class="b", coef =sisotope_time_collected_1),
                                   prior(normal(0,10),class="sigma"),
                                   prior(normal(0,10),class="sd")),
                           control = list(max_treedepth = 10, adapt_delta = .999))

summary(dn15_result_FAI_model, digits=9)
plot(conditional_effects(dn15_result_FAI_model,spaghetti=FALSE)) ##plots the regression lines/conditional effects
pp_check(dn15_result_FAI_model)
loo_R2(dn15_result_FAI_model)
bayes_R2(dn15_result_FAI_model)




dn15_result_KCAL<-brm(dn15_result~s(fai) + class_focal + s(total_kcal_using_ap_low_fermentation) + s(isotope_time_collected) + (1 +  fai + total_kcal_using_ap_low_fermentation + isotope_time_collected  | name_focal), 
                      data = data[complete.cases(data[,c("dn15_result","isotope_time_collected","fai")]),],
                      family=gaussian(),iter = 2000, 
                      prior=c(prior(normal(0,10),class="Intercept"),
                              prior(normal(0,10),class="b"),
                              prior(normal(0,10),class="sigma"),
                              prior(normal(0,10),class="sd")),
                      control = list(max_treedepth = 12, adapt_delta = .999))

summary(dn15_result_KCAL, digits=9)
plot(conditional_effects(dn15_result_KCAL,spaghetti=FALSE)) ##plots the regression lines/conditional effects
pp_check(dn15_result_KCAL)
loo_R2(dn15_result_KCAL)
bayes_R2(dn15_result_KCAL)




dn15_result_NPE<-brm(dn15_result~s(fai) + class_focal + s(total_kcal_npe_low_fermentation) + s(isotope_time_collected) + (1 +  fai + total_kcal_npe_low_fermentation + isotope_time_collected  | name_focal), 
                     data = data[complete.cases(data[,c("dn15_result","isotope_time_collected","fai","total_kcal_npe_low_fermentation")]),],
                     family=gaussian(),iter = 2000, 
                     prior=c(prior(normal(0,10),class="Intercept"),
                             prior(normal(0,10),class="b"),
                             prior(normal(0,10),class="sigma"),
                             prior(normal(0,10),class="sd")),
                     control = list(max_treedepth = 11, adapt_delta = .999))

summary(dn15_result_NPE, digits=9)
plot(conditional_effects(dn15_result_NPE,spaghetti=FALSE)) ##plots the regression lines/conditional effects
pp_check(dn15_result_NPE)
loo_R2(dn15_result_NPE)
bayes_R2(dn15_result_NPE)



dn15_result_AP<-brm(dn15_result~s(fai) + class_focal + s(ap_kcal) + s(isotope_time_collected) + (1 +  fai + ap_kcal + isotope_time_collected  | name_focal), 
                    data = data[complete.cases(data[,c("dn15_result","isotope_time_collected","fai","ap_kcal")]),],
                    family=gaussian(),iter = 2000, 
                    prior=c(prior(normal(0,10),class="Intercept"),
                            prior(normal(0,10),class="b"),
                            prior(normal(0,10),class="sigma"),
                            prior(normal(0,10),class="sd")),
                    control = list(max_treedepth = 12, adapt_delta = .999))

summary(dn15_result_AP, digits=9)
plot(conditional_effects(dn15_result_AP,spaghetti=FALSE)) ##plots the regression lines/conditional effects
pp_check(dn15_result_AP)
loo_R2(dn15_result_AP)
bayes_R2(dn15_result_AP)





dn15_result_Lipid<-brm(dn15_result~s(fai) + class_focal + s(lipid_kcal) + s(isotope_time_collected) + (1 +  fai + lipid_kcal + isotope_time_collected  | name_focal), 
                       data = data[complete.cases(data[,c("dn15_result","isotope_time_collected","fai","lipid_kcal")]),],
                       family=gaussian(),iter = 2000, 
                       prior=c(prior(normal(0,10),class="Intercept"),
                               prior(normal(0,10),class="b"),
                               prior(normal(0,10),class="sigma"),
                               prior(normal(0,10),class="sd")),
                       control = list(max_treedepth = 11, adapt_delta = .999))

summary(dn15_result_Lipid, digits=9)
plot(conditional_effects(dn15_result_Lipid,spaghetti=FALSE)) ##plots the regression lines/conditional effects
pp_check(dn15_result_Lipid)
loo_R2(dn15_result_Lipid)
bayes_R2(dn15_result_Lipid)





dn15_result_TNC<-brm(dn15_result~s(fai) + class_focal + s(tnc_kcal) + s(isotope_time_collected) + (1 +  fai + tnc_kcal + isotope_time_collected  | name_focal), 
                     data = data[complete.cases(data[,c("dn15_result","isotope_time_collected","fai","tnc_kcal")]),],
                     family=gaussian(),iter = 2000, 
                     prior=c(prior(normal(0,10),class="Intercept"),
                             prior(normal(0,10),class="b"),
                             prior(normal(0,10),class="sigma"),
                             prior(normal(0,10),class="sd")),
                     control = list(max_treedepth = 11, adapt_delta = .999))

summary(dn15_result_TNC, digits=9)
plot(conditional_effects(dn15_result_TNC,spaghetti=FALSE)) ##plots the regression lines/conditional effects
pp_check(dn15_result_TNC)
loo_R2(dn15_result_TNC)
bayes_R2(dn15_result_TNC)



dn15_result_NDF<-brm(dn15_result~s(fai) + class_focal + s(ndf_kcal_low) + s(isotope_time_collected) + (1 +  fai + ndf_kcal_low + isotope_time_collected  | name_focal), 
                     data = data[complete.cases(data[,c("dn15_result","isotope_time_collected","fai","ndf_kcal_low")]),],
                     family=gaussian(),iter = 2000, 
                     prior=c(prior(normal(0,10),class="Intercept"),
                             prior(normal(0,10),class="b"),
                             prior(normal(0,10),class="sigma"),
                             prior(normal(0,10),class="sd")),
                     control = list(max_treedepth = 11, adapt_delta = .999))

summary(dn15_result_NDF, digits=9)
plot(conditional_effects(dn15_result_NDF,spaghetti=FALSE)) ##plots the regression lines/conditional effects
pp_check(dn15_result_NDF)
loo_R2(dn15_result_NDF)
bayes_R2(dn15_result_NDF)




dn15_result_npe_ap<-brm(dn15_result~s(fai) + class_focal + s(npe_to_ap) + s(isotope_time_collected) + (1 +  fai + npe_to_ap + isotope_time_collected  | name_focal), 
                        data = data[complete.cases(data[,c("dn15_result","isotope_time_collected","fai","npe_to_ap")]),],
                        family=gaussian(),iter = 2000, 
                        prior=c(prior(normal(0,10),class="Intercept"),
                                prior(normal(0,10),class="b"),
                                prior(normal(0,10),class="sigma"),
                                prior(normal(0,10),class="sd")),
                        control = list(max_treedepth = 11, adapt_delta = .999))

summary(dn15_result_npe_ap, digits=9)
plot(conditional_effects(dn15_result_npe_ap,spaghetti=FALSE)) ##plots the regression lines/conditional effects
pp_check(dn15_result_npe_ap)
loo_R2(dn15_result_npe_ap)
bayes_R2(dn15_result_npe_ap)




dn15_result_lipid_ap<-brm(dn15_result~s(fai) + class_focal + s(lipid_to_ap) + s(isotope_time_collected) + (1 +  fai + lipid_to_ap + isotope_time_collected  | name_focal), 
                          data = data[complete.cases(data[,c("dn15_result","isotope_time_collected","fai","lipid_to_ap")]),],
                          family=gaussian(),iter = 2000, 
                          prior=c(prior(normal(0,10),class="Intercept"),
                                  prior(normal(0,10),class="b"),
                                  prior(normal(0,10),class="sigma"),
                                  prior(normal(0,10),class="sd")),
                          control = list(max_treedepth = 11, adapt_delta = .999))

summary(dn15_result_lipid_ap, digits=9)
plot(conditional_effects(dn15_result_lipid_ap,spaghetti=FALSE)) ##plots the regression lines/conditional effects
pp_check(dn15_result_lipid_ap)
loo_R2(dn15_result_lipid_ap)
bayes_R2(dn15_result_lipid_ap)



dn15_result_lipid_tnc<-brm(dn15_result~s(fai) + class_focal + s(lipid_to_tnc) + s(isotope_time_collected) + (1 +  fai + lipid_to_tnc + isotope_time_collected  | name_focal), 
                           data = data[complete.cases(data[,c("dn15_result","isotope_time_collected","fai","lipid_to_tnc")]),],
                           family=gaussian(),iter = 3000, 
                           prior=c(prior(normal(0,10),class="Intercept"),
                                   prior(normal(0,10),class="b"),
                                   prior(normal(0,10),class="sigma"),
                                   prior(normal(0,10),class="sd")),
                           control = list(max_treedepth = 11, adapt_delta = .999))

summary(dn15_result_lipid_tnc, digits=9)
plot(conditional_effects(dn15_result_lipid_tnc,spaghetti=FALSE)) ##plots the regression lines/conditional effects
pp_check(dn15_result_lipid_tnc)
loo_R2(dn15_result_lipid_tnc)
bayes_R2(dn15_result_lipid_tnc)





dn15_result_lipid_ndf<-brm(dn15_result~s(fai) + class_focal + s(lipid_to_ndf_low) + s(isotope_time_collected) + (1 +  fai + lipid_to_ndf_low + isotope_time_collected  | name_focal), 
                           data = data[complete.cases(data[,c("dn15_result","isotope_time_collected","fai","lipid_to_ndf_low")]),],
                           family=gaussian(),iter = 2000, 
                           prior=c(prior(normal(0,10),class="Intercept"),
                                   prior(normal(0,10),class="b"),
                                   prior(normal(0,10),class="sigma"),
                                   prior(normal(0,10),class="sd")),
                           control = list(max_treedepth = 11, adapt_delta = .999))

summary(dn15_result_lipid_ndf, digits=9)
plot(conditional_effects(dn15_result_lipid_ndf,spaghetti=FALSE)) ##plots the regression lines/conditional effects
pp_check(dn15_result_lipid_ndf)
loo_R2(dn15_result_lipid_ndf)
bayes_R2(dn15_result_lipid_ndf)


dn15_result_tnc_ap<-brm(dn15_result~s(fai) + class_focal + s(tnc_to_ap) + s(isotope_time_collected) + (1 +  fai + tnc_to_ap + isotope_time_collected  | name_focal), 
                        data = data[complete.cases(data[,c("dn15_result","isotope_time_collected","fai","tnc_to_ap")]),],
                        family=gaussian(),iter = 3000, 
                        prior=c(prior(normal(0,10),class="Intercept"),
                                prior(normal(0,10),class="b"),
                                prior(normal(0,10),class="sigma"),
                                prior(normal(0,10),class="sd")),
                        control = list(max_treedepth = 11, adapt_delta = .999))

summary(dn15_result_tnc_ap, digits=9)
plot(conditional_effects(dn15_result_tnc_ap,spaghetti=FALSE)) ##plots the regression lines/conditional effects
pp_check(dn15_result_tnc_ap)
loo_R2(dn15_result_tnc_ap)
bayes_R2(dn15_result_tnc_ap)



dn15_result_tnc_ndf<-brm(dn15_result~s(fai) + class_focal + s(tnc_to_ndf_low) + s(isotope_time_collected) + (1 +  fai + tnc_to_ndf_low + isotope_time_collected  | name_focal), 
                         data = data[complete.cases(data[,c("dn15_result","isotope_time_collected","fai","tnc_to_ndf_low")]),],
                         family=gaussian(),iter = 2000, 
                         prior=c(prior(normal(0,10),class="Intercept"),
                                 prior(normal(0,10),class="b"),
                                 prior(normal(0,10),class="sigma"),
                                 prior(normal(0,10),class="sd")),
                         control = list(max_treedepth = 11, adapt_delta = .999))

summary(dn15_result_tnc_ndf, digits=9)
plot(conditional_effects(dn15_result_tnc_ndf,spaghetti=FALSE)) ##plots the regression lines/conditional effects
pp_check(dn15_result_tnc_ndf)
loo_R2(dn15_result_tnc_ndf)
bayes_R2(dn15_result_tnc_ndf)




dn15_result_ap_ndf<-brm(dn15_result~s(fai) + class_focal + s(ap_to_ndf_low) + s(isotope_time_collected) + (1 +  fai + ap_to_ndf_low + isotope_time_collected  | name_focal), 
                        data = data[complete.cases(data[,c("dn15_result","isotope_time_collected","fai","ap_to_ndf_low")]),],
                        family=gaussian(),iter = 3000, 
                        prior=c(prior(normal(0,10),class="Intercept"),
                                prior(normal(0,10),class="b"),
                                prior(normal(0,10),class="sigma"),
                                prior(normal(0,10),class="sd")),
                        control = list(max_treedepth = 11, adapt_delta = .999))

summary(dn15_result_ap_ndf, digits=9)
plot(conditional_effects(dn15_result_ap_ndf,spaghetti=FALSE)) ##plots the regression lines/conditional effects
pp_check(dn15_result_ap_ndf)
loo_R2(dn15_result_ap_ndf)
bayes_R2(dn15_result_ap_ndf)







dn15_result_ap_ndf<-brm(dn15_result~s(fai) + class_focal + s(ap_to_ndf_low) + s(isotope_time_collected) + (1 +  fai + ap_to_ndf_low + isotope_time_collected  | name_focal), 
                        data = data[complete.cases(data[,c("dn15_result","isotope_time_collected","fai","ap_to_ndf_low")]),],
                        family=gaussian(),iter = 2000, 
                        prior=c(prior(normal(0,10),class="Intercept"),
                                prior(normal(0,10),class="b"),
                                prior(normal(0,10),class="sigma"),
                                prior(normal(0,10),class="sd")),
                        control = list(max_treedepth = 11, adapt_delta = .999))

summary(dn15_result_ap_ndf, digits=9)
plot(conditional_effects(dn15_result_ap_ndf,spaghetti=FALSE)) ##plots the regression lines/conditional effects
pp_check(dn15_result_ap_ndf)
loo_R2(dn15_result_ap_ndf)
bayes_R2(dn15_result_ap_ndf)








dn15_result_PREV_AP<-brm(dn15_result~s(fai) + class_focal + s(prev_day_ap_kcal) + s(isotope_time_collected) + (1 +  fai + prev_day_lipid_kcal + isotope_time_collected  | name_focal), 
                         data = data[complete.cases(data[,c("dn15_result","isotope_time_collected","fai","prev_day_lipid_kcal")]),],
                         family=gaussian(),iter = 2000, 
                         prior=c(prior(normal(0,10),class="Intercept"),
                                 prior(normal(0,10),class="b"),
                                 prior(normal(0,10),class="sigma"),
                                 prior(normal(0,10),class="sd")),
                         control = list(max_treedepth = 11, adapt_delta = .999))

summary(dn15_result_PREV_AP, digits=9)
plot(conditional_effects(dn15_result_PREV_AP,spaghetti=FALSE)) ##plots the regression lines/conditional effects
pp_check(dn15_result_PREV_AP)
loo_R2(dn15_result_PREV_AP)
bayes_R2(dn15_result_PREV_AP)







dn15_result_PREV_lipid<-brm(dn15_result~s(fai) + class_focal + s(prev_day_lipid_kcal) + s(isotope_time_collected) + (1 +  fai + prev_day_lipid_kcal + isotope_time_collected  | name_focal), 
                            data = data[complete.cases(data[,c("dn15_result","isotope_time_collected","fai","prev_day_lipid_kcal")]),],
                            family=gaussian(),iter = 2000, 
                            prior=c(prior(normal(0,10),class="Intercept"),
                                    prior(normal(0,10),class="b"),
                                    prior(normal(0,10),class="sigma"),
                                    prior(normal(0,10),class="sd")),
                            control = list(max_treedepth = 11, adapt_delta = .999))

summary(dn15_result_PREV_lipid, digits=9)
plot(conditional_effects(dn15_result_PREV_lipid,spaghetti=FALSE)) ##plots the regression lines/conditional effects
pp_check(dn15_result_PREV_lipid)
loo_R2(dn15_result_PREV_lipid)
bayes_R2(dn15_result_PREV_lipid)





dn15_result_PREV_TNC<-brm(dn15_result~s(fai) + class_focal + s(prev_day_tnc_kcal) + s(isotope_time_collected) + (1 +  fai + prev_day_tnc_kcal + isotope_time_collected  | name_focal), 
                          data = data[complete.cases(data[,c("dn15_result","isotope_time_collected","fai","prev_day_tnc_kcal")]),],
                          family=gaussian(),iter = 2000, 
                          prior=c(prior(normal(0,10),class="Intercept"),
                                  prior(normal(0,10),class="b"),
                                  prior(normal(0,10),class="sigma"),
                                  prior(normal(0,10),class="sd")),
                          control = list(max_treedepth = 11, adapt_delta = .999))

summary(dn15_result_PREV_TNC, digits=9)
plot(conditional_effects(dn15_result_PREV_TNC,spaghetti=FALSE)) ##plots the regression lines/conditional effects
pp_check(dn15_result_PREV_TNC)
loo_R2(dn15_result_PREV_TNC)
bayes_R2(dn15_result_PREV_TNC)





dn15_result_PREV_NDF<-brm(dn15_result~s(fai) + class_focal + s(prev_day_ndf_kcal_low) + s(isotope_time_collected) + (1 +  fai + prev_day_ndf_kcal_low + isotope_time_collected  | name_focal), 
                          data = data[complete.cases(data[,c("dn15_result","isotope_time_collected","fai","prev_day_ndf_kcal_low")]),],
                          family=gaussian(),iter = 2000, 
                          prior=c(prior(normal(0,10),class="Intercept"),
                                  prior(normal(0,10),class="b"),
                                  prior(normal(0,10),class="sigma"),
                                  prior(normal(0,10),class="sd")),
                          control = list(max_treedepth = 11, adapt_delta = .999))

summary(dn15_result_PREV_NDF, digits=9)
plot(conditional_effects(dn15_result_PREV_NDF,spaghetti=FALSE)) ##plots the regression lines/conditional effects
pp_check(dn15_result_PREV_NDF)
loo_R2(dn15_result_PREV_NDF)
bayes_R2(dn15_result_PREV_NDF)






dn15_result_PREV_NDF<-brm(dn15_result~s(fai) + class_focal + s(prev_day_ndf_kcal_low) + s(isotope_time_collected) + (1 +  fai + prev_day_ndf_kcal_low + isotope_time_collected  | name_focal), 
                          data = data[complete.cases(data[,c("dn15_result","isotope_time_collected","fai","prev_day_ndf_kcal_low")]),],
                          family=gaussian(),iter = 2000, 
                          prior=c(prior(normal(0,10),class="Intercept"),
                                  prior(normal(0,10),class="b"),
                                  prior(normal(0,10),class="sigma"),
                                  prior(normal(0,10),class="sd")),
                          control = list(max_treedepth = 11, adapt_delta = .999))

summary(dn15_result_PREV_NDF, digits=9)
plot(conditional_effects(dn15_result_PREV_NDF,spaghetti=FALSE)) ##plots the regression lines/conditional effects
pp_check(dn15_result_PREV_NDF)
loo_R2(dn15_result_PREV_NDF)
bayes_R2(dn15_result_PREV_NDF)

dn15_result_PREV_KCAL<-brm(dn15_result~s(fai) + class_focal + s(prev_day_total_kcal_using_ap_low_fermentation) + s(isotope_time_collected) + (1 +  fai + prev_day_total_kcal_using_ap_low_fermentation + isotope_time_collected  | name_focal), 
                           data = data[complete.cases(data[,c("dn15_result","isotope_time_collected","fai","prev_day_total_kcal_using_ap_low_fermentation")]),],
                           family=gaussian(),iter = 2000, 
                           prior=c(prior(normal(0,10),class="Intercept"),
                                   prior(normal(0,10),class="b"),
                                   prior(normal(0,10),class="sigma"),
                                   prior(normal(0,10),class="sd")),
                           control = list(max_treedepth = 11, adapt_delta = .999))

summary(dn15_result_PREV_KCAL, digits=9)
plot(conditional_effects(dn15_result_PREV_KCAL,spaghetti=FALSE)) ##plots the regression lines/conditional effects
pp_check(dn15_result_PREV_KCAL)
loo_R2(dn15_result_PREV_KCAL)
bayes_R2(dn15_result_PREV_KCAL)





dn15_result_PREV_NPE<-brm(dn15_result~s(fai) + class_focal + s(prev_day_total_kcal_npe_low_fermentation) + s(isotope_time_collected) + (1 +  fai + prev_day_total_kcal_npe_low_fermentation + isotope_time_collected  | name_focal), 
                          data = data[complete.cases(data[,c("dn15_result","isotope_time_collected","fai","prev_day_total_kcal_npe_low_fermentation")]),],
                          family=gaussian(),iter = 2000, 
                          prior=c(prior(normal(0,10),class="Intercept"),
                                  prior(normal(0,10),class="b"),
                                  prior(normal(0,10),class="sigma"),
                                  prior(normal(0,10),class="sd")),
                          control = list(max_treedepth = 12, adapt_delta = .999))

summary(dn15_result_PREV_NPE, digits=9)
plot(conditional_effects(dn15_result_PREV_NPE,spaghetti=FALSE)) ##plots the regression lines/conditional effects
pp_check(dn15_result_PREV_NPE)
loo_R2(dn15_result_PREV_NPE)
bayes_R2(dn15_result_PREV_NPE)






dn15_result_PREV_lipid_ap<-brm(dn15_result~s(fai) + class_focal + s((PREV_lipid_ap)) + s(isotope_time_collected) + (1 +  fai + PREV_lipid_ap + isotope_time_collected  | name_focal), 
                               data = data[complete.cases(data[,c("dn15_result","isotope_time_collected","fai","PREV_lipid_ap")]),],
                               family=gaussian(),iter = 3000, 
                               prior=c(prior(normal(0,10),class="Intercept"),
                                       prior(normal(0,10),class="b"),
                                       prior(normal(0,10),class="sigma"),
                                       prior(normal(0,10),class="sd")),
                               control = list(max_treedepth = 12, adapt_delta = .999))

summary(dn15_result_PREV_lipid_ap, digits=9)
plot(conditional_effects(dn15_result_PREV_lipid_ap,spaghetti=FALSE)) ##plots the regression lines/conditional effects
pp_check(dn15_result_PREV_lipid_ap)
loo_R2(dn15_result_PREV_lipid_ap)
bayes_R2(dn15_result_PREV_lipid_ap)






dn15_result_PREV_lipid_tnc<-brm(dn15_result~s(fai) + class_focal + s((PREV_lipid_tnc)) + s(isotope_time_collected) + (1 +  fai + (PREV_lipid_tnc) + isotope_time_collected  | name_focal), 
                                data = data[complete.cases(data[,c("dn15_result","isotope_time_collected","fai","PREV_lipid_tnc")]),],
                                family=gaussian(),iter = 3000, 
                                prior=c(prior(normal(0,10),class="Intercept"),
                                        prior(normal(0,10),class="b"),
                                        prior(normal(0,10),class="sigma"),
                                        prior(normal(0,10),class="sd")),
                                control = list(max_treedepth = 12, adapt_delta = .999))

summary(dn15_result_PREV_lipid_tnc, digits=9)
plot(conditional_effects(dn15_result_PREV_lipid_tnc,spaghetti=FALSE)) ##plots the regression lines/conditional effects
pp_check(dn15_result_PREV_lipid_tnc)
loo_R2(dn15_result_PREV_lipid_tnc)
bayes_R2(dn15_result_PREV_lipid_tnc)





dn15_result_PREV_lipid_ndf<-brm(dn15_result~s(fai) + class_focal + s((PREV_lipid_ndf)) + s(isotope_time_collected) + (1 +  fai + PREV_lipid_ndf + isotope_time_collected  | name_focal), 
                                data = data[complete.cases(data[,c("dn15_result","isotope_time_collected","fai","PREV_lipid_ndf")]),],
                                family=gaussian(),iter = 3000, 
                                prior=c(prior(normal(0,10),class="Intercept"),
                                        prior(normal(0,10),class="b"),
                                        prior(normal(0,10),class="sigma"),
                                        prior(normal(0,10),class="sd")),
                                control = list(max_treedepth = 12, adapt_delta = .999))

summary(dn15_result_PREV_lipid_ndf, digits=9)
plot(conditional_effects(dn15_result_PREV_lipid_ndf,spaghetti=FALSE)) ##plots the regression lines/conditional effects
pp_check(dn15_result_PREV_lipid_ndf)
loo_R2(dn15_result_PREV_lipid_ndf)
bayes_R2(dn15_result_PREV_lipid_ndf)




dn15_result_PREV_tnc_ap<-brm(dn15_result~s(fai) + class_focal + s(PREV_tnc_ap) + s(isotope_time_collected) + (1 +  fai + PREV_tnc_ap + isotope_time_collected  | name_focal), 
                             data = data[complete.cases(data[,c("dn15_result","isotope_time_collected","fai","PREV_tnc_ap")]),],
                             family=gaussian(),iter = 3000, 
                             prior=c(prior(normal(0,10),class="Intercept"),
                                     prior(normal(0,10),class="b"),
                                     prior(normal(0,10),class="sigma"),
                                     prior(normal(0,10),class="sd")),
                             control = list(max_treedepth = 12, adapt_delta = .999))

summary(dn15_result_PREV_tnc_ap, digits=9)
plot(conditional_effects(dn15_result_PREV_tnc_ap,spaghetti=FALSE)) ##plots the regression lines/conditional effects
pp_check(dn15_result_PREV_tnc_ap)
loo_R2(dn15_result_PREV_tnc_ap)
bayes_R2(dn15_result_PREV_tnc_ap)




dn15_result_PREV_tnc_ndf<-brm(dn15_result~s(fai) + class_focal + s(PREV_tnc_ndf) + s(isotope_time_collected) + (1 +  fai + PREV_tnc_ndf + isotope_time_collected  | name_focal), 
                              data = data[complete.cases(data[,c("dn15_result","isotope_time_collected","fai","PREV_tnc_ndf")]),],
                              family=gaussian(),iter = 3000, 
                              prior=c(prior(normal(0,10),class="Intercept"),
                                      prior(normal(0,10),class="b"),
                                      prior(normal(0,10),class="sigma"),
                                      prior(normal(0,10),class="sd")),
                              control = list(max_treedepth = 12, adapt_delta = .999))

summary(dn15_result_PREV_tnc_ndf, digits=9)
plot(conditional_effects(dn15_result_PREV_tnc_ndf,spaghetti=FALSE)) ##plots the regression lines/conditional effects
pp_check(dn15_result_PREV_tnc_ndf)
loo_R2(dn15_result_PREV_tnc_ndf)
bayes_R2(dn15_result_PREV_tnc_ndf)



dn15_result_PREV_ap_ndf<-brm(dn15_result~s(fai) + class_focal + s(PREV_ap_ndf) + s(isotope_time_collected) + (1 +  fai + PREV_ap_ndf + isotope_time_collected  | name_focal), 
                             data = data[complete.cases(data[,c("dn15_result","isotope_time_collected","fai","PREV_ap_ndf")]),],
                             family=gaussian(),iter = 3000, 
                             prior=c(prior(normal(0,10),class="Intercept"),
                                     prior(normal(0,10),class="b"),
                                     prior(normal(0,10),class="sigma"),
                                     prior(normal(0,10),class="sd")),
                             control = list(max_treedepth = 12, adapt_delta = .999))

summary(dn15_result_PREV_ap_ndf, digits=9)
plot(conditional_effects(dn15_result_PREV_ap_ndf,spaghetti=FALSE)) ##plots the regression lines/conditional effects
pp_check(dn15_result_PREV_ap_ndf)
loo_R2(dn15_result_PREV_ap_ndf)
bayes_R2(dn15_result_PREV_ap_ndf)



dn15_result_PREV_npe_ab<-brm(dn15_result~s(fai) + class_focal + s(prev_day_npe_ab) + s(isotope_time_collected) + (1 +  fai + prev_day_npe_ab + isotope_time_collected  | name_focal), 
                             data = data[complete.cases(data[,c("dn15_result","isotope_time_collected","fai","prev_day_npe_ab")]),],
                             family=gaussian(),iter = 3000, 
                             prior=c(prior(normal(0,10),class="Intercept"),
                                     prior(normal(0,10),class="b"),
                                     prior(normal(0,10),class="sigma"),
                                     prior(normal(0,10),class="sd")),
                             control = list(max_treedepth = 12, adapt_delta = .999))

summary(dn15_result_PREV_npe_ab, digits=9)
plot(conditional_effects(dn15_result_PREV_npe_ab,spaghetti=FALSE)) ##plots the regression lines/conditional effects
pp_check(dn15_result_PREV_npe_ab)
loo_R2(dn15_result_PREV_npe_ab)
bayes_R2(dn15_result_PREV_npe_ab)




#####Urea MODELS######

urea_FAI_model<-brm(urea_sgcor~s(fai) + s(urea_time_collected) + (1 +  fai + urea_time_collected  | name_focal), 
                    data = data[complete.cases(data[,c("urea_sgcor","fai")]),],
                    family=Gamma(link="log"),iter = 8000, inits=0,
                    prior=c(prior(normal(0,10),class="Intercept"),
                            prior(normal(0,10),class="b"),
                            prior(gamma(.01,.01),class="shape"),
                            prior(normal(0,10),class="sd")),
                    control = list(max_treedepth = 11, adapt_delta = .999))

summary(urea_FAI_model, digits=9)
plot(urea_FAI_model)
plot(conditional_effects(urea_FAI_model,spaghetti=TRUE)) ##plots the regression lines/conditional effects
pp_check(urea_FAI_model)
loo_R2(urea_FAI_model)
bayes_R2(urea_FAI_model)
loo(urea_FAI_model)



urea_sgcor_KCAL<-brm(urea_sgcor~s(fai) + class_focal + s(total_kcal_using_ap_low_fermentation) + s(urea_time_collected) + (1 +  fai + total_kcal_using_ap_low_fermentation + urea_time_collected  | name_focal), 
                     data = data[complete.cases(data[,c("urea_sgcor","urea_time_collected","fai")]),],
                     family=Gamma(link="log"),iter = 8000, inits=0, 
                     prior=c(prior(normal(0,10),class="Intercept"),
                             prior(normal(0,10),class="b"),
                             prior(gamma(.01,.01),class="shape"),
                             prior(normal(0,10),class="sd")),
                     control = list(max_treedepth = 12, adapt_delta = .999))

summary(urea_sgcor_KCAL, digits=9)
plot(conditional_effects(urea_sgcor_KCAL,spaghetti=FALSE)) ##plots the regression lines/conditional effects
pp_check(urea_sgcor_KCAL)
loo_R2(urea_sgcor_KCAL)
bayes_R2(urea_sgcor_KCAL)



urea_sgcor_NPE<-brm(urea_sgcor~s(fai) + class_focal + s(total_kcal_npe_low_fermentation) + s(urea_time_collected) + (1 +  fai + total_kcal_npe_low_fermentation + urea_time_collected  | name_focal), 
                    data = data[complete.cases(data[,c("urea_sgcor","urea_time_collected","fai","total_kcal_npe_low_fermentation")]),],
                    family=Gamma(link="log"),iter = 8000, inits=0, 
                    prior=c(prior(normal(0,10),class="Intercept"),
                            prior(normal(0,10),class="b"),
                            prior(gamma(.01,.01),class="shape"),
                            prior(normal(0,10),class="sd")),
                    control = list(max_treedepth = 11, adapt_delta = .999))

summary(urea_sgcor_NPE, digits=9)
plot(conditional_effects(urea_sgcor_NPE,spaghetti=FALSE)) ##plots the regression lines/conditional effects
pp_check(urea_sgcor_NPE)
loo_R2(urea_sgcor_NPE)
bayes_R2(urea_sgcor_NPE)



urea_sgcor_AP<-brm(urea_sgcor~s(fai) + class_focal + s(ap_kcal) + s(urea_time_collected) + (1 +  fai + ap_kcal + urea_time_collected  | name_focal), 
                   data = data[complete.cases(data[,c("urea_sgcor","urea_time_collected","fai","ap_kcal")]),],
                   family=Gamma(link="log"),iter = 8000, inits=0, 
                   prior=c(prior(normal(0,10),class="Intercept"),
                           prior(normal(0,10),class="b"),
                           prior(gamma(.01,.01),class="shape"),
                           prior(normal(0,10),class="sd")),
                   control = list(max_treedepth = 12, adapt_delta = .999))

summary(urea_sgcor_AP, digits=9)
plot(conditional_effects(urea_sgcor_AP,spaghetti=FALSE)) ##plots the regression lines/conditional effects
pp_check(urea_sgcor_AP)
loo_R2(urea_sgcor_AP)
bayes_R2(urea_sgcor_AP)







urea_sgcor_Lipid<-brm(urea_sgcor~s(fai) + class_focal + s(lipid_kcal) + s(urea_time_collected) + (1 +  fai + lipid_kcal + urea_time_collected  | name_focal), 
                      data = data[complete.cases(data[,c("urea_sgcor","urea_time_collected","fai","lipid_kcal")]),],
                      family=Gamma(link="log"),iter = 8000, inits=0, 
                      prior=c(prior(normal(0,10),class="Intercept"),
                              prior(normal(0,10),class="b"),
                              prior(gamma(.01,.01),class="shape"),
                              prior(normal(0,10),class="sd")),
                      control = list(max_treedepth = 11, adapt_delta = .999))

summary(urea_sgcor_Lipid, digits=9)
plot(conditional_effects(urea_sgcor_Lipid,spaghetti=FALSE)) ##plots the regression lines/conditional effects
pp_check(urea_sgcor_Lipid)
loo_R2(urea_sgcor_Lipid)
bayes_R2(urea_sgcor_Lipid)






urea_sgcor_TNC<-brm(urea_sgcor~s(fai) + class_focal + s(tnc_kcal) + s(urea_time_collected) + (1 +  fai + tnc_kcal + urea_time_collected  | name_focal), 
                    data = data[complete.cases(data[,c("urea_sgcor","urea_time_collected","fai","tnc_kcal")]),],
                    family=Gamma(link="log"),iter = 8000, inits=0, 
                    prior=c(prior(normal(0,10),class="Intercept"),
                            prior(normal(0,10),class="b"),
                            prior(gamma(.01,.01),class="shape"),
                            prior(normal(0,10),class="sd")),
                    control = list(max_treedepth = 12, adapt_delta = .999))

summary(urea_sgcor_TNC, digits=9)
plot(conditional_effects(urea_sgcor_TNC,spaghetti=FALSE)) ##plots the regression lines/conditional effects
pp_check(urea_sgcor_TNC)
loo_R2(urea_sgcor_TNC)
bayes_R2(urea_sgcor_TNC)





urea_sgcor_NDF<-brm(urea_sgcor~s(fai) + class_focal + s(ndf_kcal_low) + s(urea_time_collected) + (1 +  fai + ndf_kcal_low + urea_time_collected  | name_focal), 
                    data = data[complete.cases(data[,c("urea_sgcor","urea_time_collected","fai","ndf_kcal_low")]),],
                    family=Gamma(link="log"),iter = 8000, inits=0, 
                    prior=c(prior(normal(0,10),class="Intercept"),
                            prior(normal(0,10),class="b"),
                            prior(gamma(.01,.01),class="shape"),
                            prior(normal(0,10),class="sd")),
                    control = list(max_treedepth = 11, adapt_delta = .999))

summary(urea_sgcor_NDF, digits=9)
plot(conditional_effects(urea_sgcor_NDF,spaghetti=FALSE)) ##plots the regression lines/conditional effects
pp_check(urea_sgcor_NDF)
loo_R2(urea_sgcor_NDF)
bayes_R2(urea_sgcor_NDF)





urea_sgcor_npe_ap<-brm(urea_sgcor~s(fai) + class_focal + s(npe_to_ap) + s(urea_time_collected) + (1 +  fai + npe_to_ap + urea_time_collected  | name_focal), 
                       data = data[complete.cases(data[,c("urea_sgcor","urea_time_collected","fai","npe_to_ap")]),],
                       family=Gamma(link="log"),iter = 8000, inits=0, 
                       prior=c(prior(normal(0,10),class="Intercept"),
                               prior(normal(0,10),class="b"),
                               prior(gamma(.01,.01),class="shape"),
                               prior(normal(0,10),class="sd")),
                       control = list(max_treedepth = 12, adapt_delta = .999))

summary(urea_sgcor_npe_ap, digits=9)
plot(conditional_effects(urea_sgcor_npe_ap,spaghetti=FALSE)) ##plots the regression lines/conditional effects
pp_check(urea_sgcor_npe_ap)
loo_R2(urea_sgcor_npe_ap)
bayes_R2(urea_sgcor_npe_ap)






urea_sgcor_lipid_ap<-brm(urea_sgcor~s(fai) + class_focal + s(lipid_to_ap) + s(urea_time_collected) + (1 +  fai + lipid_to_ap + urea_time_collected  | name_focal), 
                         data = data[complete.cases(data[,c("urea_sgcor","urea_time_collected","fai","lipid_to_ap")]),],
                         family=Gamma(link="log"),iter = 8000, inits=0, 
                         prior=c(prior(normal(0,10),class="Intercept"),
                                 prior(normal(0,10),class="b"),
                                 prior(gamma(.01,.01),class="shape"),
                                 prior(normal(0,10),class="sd")),
                         control = list(max_treedepth = 11, adapt_delta = .999))

summary(urea_sgcor_lipid_ap, digits=9)
plot(conditional_effects(urea_sgcor_lipid_ap,spaghetti=FALSE)) ##plots the regression lines/conditional effects
pp_check(urea_sgcor_lipid_ap)
loo_R2(urea_sgcor_lipid_ap)
bayes_R2(urea_sgcor_lipid_ap)




urea_sgcor_lipid_tnc<-brm(urea_sgcor~s(fai) + class_focal + s(lipid_to_tnc) + s(urea_time_collected) + (1 +  fai + lipid_to_tnc + urea_time_collected  | name_focal), 
                          data = data[complete.cases(data[,c("urea_sgcor","urea_time_collected","fai","lipid_to_tnc")]),],
                          family=Gamma(link="log"),iter = 8000, inits=0, 
                          prior=c(prior(normal(0,10),class="Intercept"),
                                  prior(normal(0,10),class="b"),
                                  prior(gamma(.01,.01),class="shape"),
                                  prior(normal(0,10),class="sd")),
                          control = list(max_treedepth = 12, adapt_delta = .999))

summary(urea_sgcor_lipid_tnc, digits=9)
plot(conditional_effects(urea_sgcor_lipid_tnc,spaghetti=FALSE)) ##plots the regression lines/conditional effects
pp_check(urea_sgcor_lipid_tnc)
loo_R2(urea_sgcor_lipid_tnc)
bayes_R2(urea_sgcor_lipid_tnc)




urea_sgcor_lipid_ndf<-brm(urea_sgcor~s(fai) + class_focal + s(lipid_to_ndf_low) + s(urea_time_collected) + (1 +  fai + lipid_to_ndf_low + urea_time_collected  | name_focal), 
                          data = data[complete.cases(data[,c("urea_sgcor","urea_time_collected","fai","lipid_to_ndf_low")]),],
                          family=Gamma(link="log"),iter = 8000, inits=0, 
                          prior=c(prior(normal(0,10),class="Intercept"),
                                  prior(normal(0,10),class="b"),
                                  prior(gamma(.01,.01),class="shape"),
                                  prior(normal(0,10),class="sd")),
                          control = list(max_treedepth = 12, adapt_delta = .999))

summary(urea_sgcor_lipid_ndf, digits=9)
plot(conditional_effects(urea_sgcor_lipid_ndf,spaghetti=FALSE)) ##plots the regression lines/conditional effects
pp_check(urea_sgcor_lipid_ndf)
loo_R2(urea_sgcor_lipid_ndf)
bayes_R2(urea_sgcor_lipid_ndf)









urea_sgcor_tnc_ap<-brm(urea_sgcor~s(fai) + class_focal + s(tnc_to_ap) + s(urea_time_collected) + (1 +  fai + tnc_to_ap + urea_time_collected  | name_focal), 
                       data = data[complete.cases(data[,c("urea_sgcor","urea_time_collected","fai","tnc_to_ap")]),],
                       family=Gamma(link="log"),iter = 8000, inits=0, 
                       prior=c(prior(normal(0,10),class="Intercept"),
                               prior(normal(0,10),class="b"),
                               prior(gamma(.01,.01),class="shape"),
                               prior(normal(0,10),class="sd")),
                       control = list(max_treedepth = 11, adapt_delta = .999))

summary(urea_sgcor_tnc_ap, digits=9)
plot(conditional_effects(urea_sgcor_tnc_ap,spaghetti=FALSE)) ##plots the regression lines/conditional effects
pp_check(urea_sgcor_tnc_ap)
loo_R2(urea_sgcor_tnc_ap)
bayes_R2(urea_sgcor_tnc_ap)





urea_sgcor_tnc_ndf<-brm(urea_sgcor~s(fai) + class_focal + s(tnc_to_ndf_low) + s(urea_time_collected) + (1 +  fai + tnc_to_ndf_low + urea_time_collected  | name_focal), 
                        data = data[complete.cases(data[,c("urea_sgcor","urea_time_collected","fai","tnc_to_ndf_low")]),],
                        family=Gamma(link="log"),iter = 8000, inits=0, 
                        prior=c(prior(normal(0,10),class="Intercept"),
                                prior(normal(0,10),class="b"),
                                prior(gamma(.01,.01),class="shape"),
                                prior(normal(0,10),class="sd")),
                        control = list(max_treedepth = 11, adapt_delta = .999))

summary(urea_sgcor_tnc_ndf, digits=9)
plot(conditional_effects(urea_sgcor_tnc_ndf,spaghetti=FALSE)) ##plots the regression lines/conditional effects
pp_check(urea_sgcor_tnc_ndf)
loo_R2(urea_sgcor_tnc_ndf)
bayes_R2(urea_sgcor_tnc_ndf)







urea_sgcor_ap_ndf<-brm(urea_sgcor~s(fai) + class_focal + s(ap_to_ndf_low) + s(urea_time_collected) + (1 +  fai + ap_to_ndf_low + urea_time_collected  | name_focal), 
                       data = data[complete.cases(data[,c("urea_sgcor","urea_time_collected","fai","ap_to_ndf_low")]),],
                       family=Gamma(link="log"),iter = 8000, inits=0, 
                       prior=c(prior(normal(0,10),class="Intercept"),
                               prior(normal(0,10),class="b"),
                               prior(gamma(.01,.01),class="shape"),
                               prior(normal(0,10),class="sd")),
                       control = list(max_treedepth = 12, adapt_delta = .999))

summary(urea_sgcor_ap_ndf, digits=9)
plot(conditional_effects(urea_sgcor_ap_ndf,spaghetti=FALSE)) ##plots the regression lines/conditional effects
pp_check(urea_sgcor_ap_ndf)
loo_R2(urea_sgcor_ap_ndf)
bayes_R2(urea_sgcor_ap_ndf)






urea_sgcor_PREV_AP<-brm(urea_sgcor~s(fai) + class_focal + s(prev_day_ap_kcal) + s(urea_time_collected) + (1 +  fai + prev_day_ap_kcal + urea_time_collected  | name_focal), 
                        data = data[complete.cases(data[,c("urea_sgcor","urea_time_collected","fai","prev_day_ap_kcal")]),],
                        family=Gamma(link="log"),iter = 8000, inits=0, 
                        prior=c(prior(normal(0,10),class="Intercept"),
                                prior(normal(0,10),class="b"),
                                prior(gamma(.01,.01),class="shape"),
                                prior(normal(0,10),class="sd")),
                        control = list(max_treedepth = 11, adapt_delta = .999))

summary(urea_sgcor_PREV_AP, digits=9)
plot(conditional_effects(urea_sgcor_PREV_AP,spaghetti=FALSE)) ##plots the regression lines/conditional effects
pp_check(urea_sgcor_PREV_AP)
loo_R2(urea_sgcor_PREV_AP)
bayes_R2(urea_sgcor_PREV_AP)





urea_sgcor_PREV_lipid<-brm(urea_sgcor~s(fai) + class_focal + s(prev_day_lipid_kcal) + s(urea_time_collected) + (1 +  fai + prev_day_lipid_kcal + urea_time_collected  | name_focal), 
                           data = data[complete.cases(data[,c("urea_sgcor","urea_time_collected","fai","prev_day_lipid_kcal")]),],
                           family=Gamma(link="log"),iter = 8000, inits=0, 
                           prior=c(prior(normal(0,10),class="Intercept"),
                                   prior(normal(0,10),class="b"),
                                   prior(gamma(.01,.01),class="shape"),
                                   prior(normal(0,10),class="sd")),
                           control = list(max_treedepth = 11, adapt_delta = .999))

summary(urea_sgcor_PREV_lipid, digits=9)
plot(conditional_effects(urea_sgcor_PREV_lipid,spaghetti=FALSE)) ##plots the regression lines/conditional effects
pp_check(urea_sgcor_PREV_lipid)
loo_R2(urea_sgcor_PREV_lipid)
bayes_R2(urea_sgcor_PREV_lipid)







urea_sgcor_PREV_TNC<-brm(urea_sgcor~s(fai) + class_focal + s(prev_day_tnc_kcal) + s(urea_time_collected) + (1 +  fai + prev_day_tnc_kcal + urea_time_collected  | name_focal), 
                         data = data[complete.cases(data[,c("urea_sgcor","urea_time_collected","fai","prev_day_tnc_kcal")]),],
                         family=Gamma(link="log"),iter = 8000, inits=0, 
                         prior=c(prior(normal(0,10),class="Intercept"),
                                 prior(normal(0,10),class="b"),
                                 prior(gamma(.01,.01),class="shape"),
                                 prior(normal(0,10),class="sd")),
                         control = list(max_treedepth = 12, adapt_delta = .999))

summary(urea_sgcor_PREV_TNC, digits=9)
plot(conditional_effects(urea_sgcor_PREV_TNC,spaghetti=FALSE)) ##plots the regression lines/conditional effects
pp_check(urea_sgcor_PREV_TNC)
loo_R2(urea_sgcor_PREV_TNC)
bayes_R2(urea_sgcor_PREV_TNC)




urea_sgcor_PREV_NDF<-brm(urea_sgcor~s(fai) + class_focal + s(prev_day_ndf_kcal_low) + s(urea_time_collected) + (1 +  fai + prev_day_ndf_kcal_low + urea_time_collected  | name_focal), 
                         data = data[complete.cases(data[,c("urea_sgcor","urea_time_collected","fai","prev_day_ndf_kcal_low")]),],
                         family=Gamma(link="log"),iter = 8000, inits=0, 
                         prior=c(prior(normal(0,10),class="Intercept"),
                                 prior(normal(0,10),class="b"),
                                 prior(gamma(.01,.01),class="shape"),
                                 prior(normal(0,10),class="sd")),
                         control = list(max_treedepth = 11, adapt_delta = .999))

summary(urea_sgcor_PREV_NDF, digits=9)
plot(conditional_effects(urea_sgcor_PREV_NDF,spaghetti=FALSE)) ##plots the regression lines/conditional effects
pp_check(urea_sgcor_PREV_NDF)
loo_R2(urea_sgcor_PREV_NDF)
bayes_R2(urea_sgcor_PREV_NDF)





urea_sgcor_PREV_KCAL<-brm(urea_sgcor~s(fai) + class_focal + s(prev_day_total_kcal_using_ap_low_fermentation) + s(urea_time_collected) + (1 +  fai + prev_day_total_kcal_using_ap_low_fermentation + urea_time_collected  | name_focal), 
                          data = data[complete.cases(data[,c("urea_sgcor","urea_time_collected","fai","prev_day_total_kcal_using_ap_low_fermentation")]),],
                          family=Gamma(link="log"),iter = 8000, inits=0, 
                          prior=c(prior(normal(0,10),class="Intercept"),
                                  prior(normal(0,10),class="b"),
                                  prior(gamma(.01,.01),class="shape"),
                                  prior(normal(0,10),class="sd")),
                          control = list(max_treedepth = 11, adapt_delta = .999))

summary(urea_sgcor_PREV_KCAL, digits=9)
plot(conditional_effects(urea_sgcor_PREV_KCAL,spaghetti=FALSE)) ##plots the regression lines/conditional effects
pp_check(urea_sgcor_PREV_KCAL)
loo_R2(urea_sgcor_PREV_KCAL)
bayes_R2(urea_sgcor_PREV_KCAL)





urea_sgcor_PREV_NPE<-brm(urea_sgcor~s(fai) + class_focal + s(prev_day_total_kcal_npe_low_fermentation) + s(urea_time_collected) + (1 +  fai + prev_day_total_kcal_npe_low_fermentation + urea_time_collected  | name_focal), 
                         data = data[complete.cases(data[,c("urea_sgcor","urea_time_collected","fai","prev_day_total_kcal_npe_low_fermentation")]),],
                         family=Gamma(link="log"),iter = 8000, inits=0, 
                         prior=c(prior(normal(0,10),class="Intercept"),
                                 prior(normal(0,10),class="b"),
                                 prior(gamma(.01,.01),class="shape"),
                                 prior(normal(0,10),class="sd")),
                         control = list(max_treedepth = 12, adapt_delta = .999))

summary(urea_sgcor_PREV_NPE, digits=9)
plot(conditional_effects(urea_sgcor_PREV_NPE,spaghetti=FALSE)) ##plots the regression lines/conditional effects
pp_check(urea_sgcor_PREV_NPE)
loo_R2(urea_sgcor_PREV_NPE)
bayes_R2(urea_sgcor_PREV_NPE)




urea_sgcor_PREV_lipid_ap<-brm(urea_sgcor~s(fai) + class_focal + s((PREV_lipid_ap)) + s(urea_time_collected) + (1 +  fai + PREV_lipid_ap + urea_time_collected  | name_focal), 
                              data = data[complete.cases(data[,c("urea_sgcor","urea_time_collected","fai","PREV_lipid_ap")]),],
                              family=Gamma(link="log"),iter = 8000, inits=0, 
                              prior=c(prior(normal(0,10),class="Intercept"),
                                      prior(normal(0,10),class="b"),
                                      prior(gamma(.01,.01),class="shape"),
                                      prior(normal(0,10),class="sd")),
                              control = list(max_treedepth = 12, adapt_delta = .999))

summary(urea_sgcor_PREV_lipid_ap, digits=9)
plot(conditional_effects(urea_sgcor_PREV_lipid_ap,spaghetti=FALSE)) ##plots the regression lines/conditional effects
pp_check(urea_sgcor_PREV_lipid_ap)
loo_R2(urea_sgcor_PREV_lipid_ap)
bayes_R2(urea_sgcor_PREV_lipid_ap)






urea_sgcor_PREV_lipid_tnc<-brm(urea_sgcor~s(fai) + class_focal + s((PREV_lipid_tnc)) + s(urea_time_collected) + (1 +  fai + (PREV_lipid_tnc) + urea_time_collected  | name_focal), 
                               data = data[complete.cases(data[,c("urea_sgcor","urea_time_collected","fai","PREV_lipid_tnc")]),],
                               family=Gamma(link="log"),iter = 8000, inits=0, 
                               prior=c(prior(normal(0,10),class="Intercept"),
                                       prior(normal(0,10),class="b"),
                                       prior(gamma(.01,.01),class="shape"),
                                       prior(normal(0,10),class="sd")),
                               control = list(max_treedepth = 12, adapt_delta = .999))

summary(urea_sgcor_PREV_lipid_tnc, digits=9)
plot(conditional_effects(urea_sgcor_PREV_lipid_tnc,spaghetti=FALSE)) ##plots the regression lines/conditional effects
pp_check(urea_sgcor_PREV_lipid_tnc)
loo_R2(urea_sgcor_PREV_lipid_tnc)
bayes_R2(urea_sgcor_PREV_lipid_tnc)



urea_sgcor_PREV_lipid_ndf<-brm(urea_sgcor~s(fai) + class_focal + s((PREV_lipid_ndf)) + s(urea_time_collected) + (1 +  fai + PREV_lipid_ndf + urea_time_collected  | name_focal), 
                               data = data[complete.cases(data[,c("urea_sgcor","urea_time_collected","fai","PREV_lipid_ndf")]),],
                               family=Gamma(link="log"),iter = 8000, inits=0, 
                               prior=c(prior(normal(0,10),class="Intercept"),
                                       prior(normal(0,10),class="b"),
                                       prior(gamma(.01,.01),class="shape"),
                                       prior(normal(0,10),class="sd")),
                               control = list(max_treedepth = 12, adapt_delta = .999))

summary(urea_sgcor_PREV_lipid_ndf, digits=9)
plot(conditional_effects(urea_sgcor_PREV_lipid_ndf,spaghetti=FALSE)) ##plots the regression lines/conditional effects
pp_check(urea_sgcor_PREV_lipid_ndf)
loo_R2(urea_sgcor_PREV_lipid_ndf)
bayes_R2(urea_sgcor_PREV_lipid_ndf)





urea_sgcor_PREV_tnc_ap<-brm(urea_sgcor~s(fai) + class_focal + s(PREV_tnc_ap) + s(urea_time_collected) + (1 +  fai + PREV_tnc_ap + urea_time_collected  | name_focal), 
                            data = data[complete.cases(data[,c("urea_sgcor","urea_time_collected","fai","PREV_tnc_ap")]),],
                            family=Gamma(link="log"),iter = 8000, inits=0, 
                            prior=c(prior(normal(0,10),class="Intercept"),
                                    prior(normal(0,10),class="b"),
                                    prior(gamma(.01,.01),class="shape"),
                                    prior(normal(0,10),class="sd")),
                            control = list(max_treedepth = 12, adapt_delta = .999))

summary(urea_sgcor_PREV_tnc_ap, digits=9)
plot(conditional_effects(urea_sgcor_PREV_tnc_ap,spaghetti=FALSE)) ##plots the regression lines/conditional effects
pp_check(urea_sgcor_PREV_tnc_ap)
loo_R2(urea_sgcor_PREV_tnc_ap)
bayes_R2(urea_sgcor_PREV_tnc_ap)




urea_sgcor_PREV_tnc_ndf<-brm(urea_sgcor~s(fai) + class_focal + s(PREV_tnc_ndf) + s(urea_time_collected) + (1 +  fai + PREV_tnc_ndf + urea_time_collected  | name_focal), 
                             data = data[complete.cases(data[,c("urea_sgcor","urea_time_collected","fai","PREV_tnc_ndf")]),],
                             family=Gamma(link="log"),iter = 8000, inits=0, 
                             prior=c(prior(normal(0,10),class="Intercept"),
                                     prior(normal(0,10),class="b"),
                                     prior(gamma(.01,.01),class="shape"),
                                     prior(normal(0,10),class="sd")),
                             control = list(max_treedepth = 12, adapt_delta = .999))

summary(urea_sgcor_PREV_tnc_ndf, digits=9)
plot(conditional_effects(urea_sgcor_PREV_tnc_ndf,spaghetti=FALSE)) ##plots the regression lines/conditional effects
pp_check(urea_sgcor_PREV_tnc_ndf)
loo_R2(urea_sgcor_PREV_tnc_ndf)
bayes_R2(urea_sgcor_PREV_tnc_ndf)





urea_sgcor_PREV_ap_ndf<-brm(urea_sgcor~s(fai) + class_focal + s(PREV_ap_ndf) + s(urea_time_collected) + (1 +  fai + PREV_ap_ndf + urea_time_collected  | name_focal), 
                            data = data[complete.cases(data[,c("urea_sgcor","urea_time_collected","fai","PREV_ap_ndf")]),],
                            family=Gamma(link="log"),iter = 8000, inits=0, 
                            prior=c(prior(normal(0,10),class="Intercept"),
                                    prior(normal(0,10),class="b"),
                                    prior(gamma(.01,.01),class="shape"),
                                    prior(normal(0,10),class="sd")),
                            control = list(max_treedepth = 12, adapt_delta = .999))

summary(urea_sgcor_PREV_ap_ndf, digits=9)
plot(conditional_effects(urea_sgcor_PREV_ap_ndf,spaghetti=FALSE)) ##plots the regression lines/conditional effects
pp_check(urea_sgcor_PREV_ap_ndf)
loo_R2(urea_sgcor_PREV_ap_ndf)
bayes_R2(urea_sgcor_PREV_ap_ndf)






urea_sgcor_PREV_npe_ab<-brm(urea_sgcor~s(fai) + class_focal + s(prev_day_npe_ab) + s(urea_time_collected) + (1 +  fai + prev_day_npe_ab + urea_time_collected  | name_focal), 
                            data = data[complete.cases(data[,c("urea_sgcor","urea_time_collected","fai","prev_day_npe_ab")]),],
                            family=Gamma(link="log"),iter = 8000, inits=0, 
                            prior=c(prior(normal(0,10),class="Intercept"),
                                    prior(normal(0,10),class="b"),
                                    prior(gamma(.01,.01),class="shape"),
                                    prior(normal(0,10),class="sd")),
                            control = list(max_treedepth = 12, adapt_delta = .999))

summary(urea_sgcor_PREV_npe_ab, digits=9)
plot(conditional_effects(urea_sgcor_PREV_npe_ab,spaghetti=FALSE)) ##plots the regression lines/conditional effects
pp_check(urea_sgcor_PREV_npe_ab)
loo_R2(urea_sgcor_PREV_npe_ab)
bayes_R2(urea_sgcor_PREV_npe_ab)







############################## PLOTS ######################


deriv_plot(model=Ketone_PREV_KCAL,term="s(prev_day_total_kcal_using_ap_low_fermentation)",main="prev_day_total_kcal_using_ap_low_fermentation",eps=0.001,output="Ketone_PREV_KCAL_plot")
deriv_plot(model=Ketone_PREV_AP,term="s(prev_day_ap_kcal)",main="prev_day_ap_kcal",eps=0.001,output="Ketone_PREV_AP_plot")
deriv_plot(model=Ketone_PREV_lipid,term="s(prev_day_lipid_kcal)",main="prev_day_lipid_kcal",eps=0.001,output="Ketone_PREV_lipid_plot")
deriv_plot(model=Ketone_PREV_TNC,term="s(prev_day_tnc_kcal)",main="prev_day_tnc_kcal",eps=0.001,output="Ketone_PREV_TNC_plot")
deriv_plot(model=Ketone_PREV_NDF,term="s(prev_day_ndf_kcal_low)",main="prev_day_ndf_kcal_low",eps=0.001,output="Ketone_PREV_NDF_plot")
deriv_plot(model=Ketone_PREV_NPE,term="s(prev_day_total_kcal_npe_low_fermentation)",main="prev_day_total_kcal_npe_low_fermentation",eps=0.001,output="Ketone_PREV_NPE_plot")
deriv_plot(model=Ketone_PREV_npe_ab,term="s(prev_day_npe_ab)",main="prev_day_npe_ab",eps=0.001,output="Ketone_PREV_npe_ab_plot")
deriv_plot(model=Ketone_PREV_lipid_ap,term="s((PREV_lipid_ap))",main="(PREV_lipid_ap)",eps=0.001,output="Ketone_PREV_lipid_ap_plot")
deriv_plot(model=Ketone_PREV_lipid_tnc,term="s((PREV_lipid_tnc))",main="(PREV_lipid_tnc)",eps=0.001,output="Ketone_PREV_lipid_tnc_plot")
deriv_plot(model=Ketone_PREV_lipid_ndf,term="s((PREV_lipid_ndf))",main="(PREV_lipid_ndf)",eps=0.001,output="Ketone_PREV_lipid_ndf_plot")
deriv_plot(model=Ketone_PREV_tnc_ap,term="s((PREV_tnc_ap))",main="(PREV_tnc_ap)",eps=0.001,output="Ketone_PREV_tnc_ap_plot")
deriv_plot(model=Ketone_PREV_tnc_ndf,term="s((PREV_tnc_ndf))",main="(PREV_tnc_ndf)",eps=0.001,output="Ketone_PREV_tnc_ndf_plot")
deriv_plot(model=Ketone_PREV_ap_ndf,term="s((PREV_ap_ndf))",main="(PREV_ap_ndf)",eps=0.001,output="Ketone_PREV_ap_ndf_plot")

Ketone_PREV_KCAL_plot=Ketone_PREV_KCAL_plot+ylab("Ketone body presence")+xlab("Total caloric intake")
Ketone_PREV_AP_plot=Ketone_PREV_AP_plot+ylab("Ketone pres/abs")+xlab("Previous day protein (kcal)")
Ketone_PREV_lipid_plot=Ketone_PREV_lipid_plot+ylab("Ketone pres/abs")+xlab("Previous day lipid (kcal)")
Ketone_PREV_TNC_plot=Ketone_PREV_TNC_plot+ylab("Ketone pres/abs")+xlab("Previous day TNC (kcal)")
Ketone_PREV_NDF_plot=Ketone_PREV_NDF_plot+ylab("Ketone pres/abs")+xlab("Previous NDF protein (kcal)")
Ketone_PREV_NPE_plot=Ketone_PREV_NPE_plot+ylab("Ketone pres/abs")+xlab("Previous day NPE")
Ketone_PREV_lipid_ap_plot=Ketone_PREV_lipid_ap_plot+ylab("Ketone pres/abs")+xlab("Previous day Lipid/Ap")

Ketone_PREV_KCAL_plot=Ketone_PREV_KCAL_plot+coord_cartesian(ylim = c(0,.15))+ylab("Ketone body presence")+xlab("Total caloric intake")+ggtitle("A")
Ketone_PREV_AP_plot=Ketone_PREV_AP_plot+coord_cartesian(ylim = c(0,.15))+ylab("Ketone body presence")+xlab("P")+ggtitle("B")
Ketone_PREV_lipid_plot=Ketone_PREV_lipid_plot+coord_cartesian(ylim = c(0,.15))+ylab("Ketone body presence")+xlab("Lipid")+ggtitle("C")
Ketone_PREV_TNC_plot=Ketone_PREV_TNC_plot+coord_cartesian(ylim = c(0,.15))+ylab("Ketone body presence")+xlab("TNC")+ggtitle("D")
Ketone_PREV_NDF_plot=Ketone_PREV_NDF_plot+coord_cartesian(ylim = c(0,.15))+ylab("Ketone body presence")+xlab("NDF")+ggtitle("E")
Ketone_PREV_NPE_plot=Ketone_PREV_NPE_plot+coord_cartesian(ylim = c(0,.15))+ylab("Ketone body presence")+xlab("NPe")+ggtitle("F")
Ketone_PREV_npe_ab_plot=Ketone_PREV_npe_ab_plot+coord_cartesian(ylim = c(0,.15))+ylab("Ketone body presence")+xlab("NPe:P")+ggtitle("G")
Ketone_PREV_lipid_ap_plot=Ketone_PREV_lipid_ap_plot+coord_cartesian(ylim = c(0,.15))+ylab("Ketone body presence")+xlab("Lipid:P")+ggtitle("H")
Ketone_PREV_lipid_tnc_plot=Ketone_PREV_lipid_tnc_plot+coord_cartesian(ylim = c(0,.15))+ylab("Ketone body presence")+xlab("Lipid:TNC")+ggtitle("I")  + geom_vline(xintercept = 0.32, color = "gray")+ geom_vline(xintercept = c(0.18), color = "gray", linetype = "dashed")
Ketone_PREV_lipid_ndf_plot=Ketone_PREV_lipid_ndf_plot+coord_cartesian(ylim = c(0,.15))+ylab("Ketone body presence")+xlab("Lipid:NDF")+ggtitle("J")
Ketone_PREV_tnc_ap_plot=Ketone_PREV_tnc_ap_plot+coord_cartesian(ylim = c(0,.15))+ylab("Ketone body presence")+xlab("TNC:P")+ggtitle("K")
Ketone_PREV_tnc_ndf_plot=Ketone_PREV_tnc_ndf_plot+coord_cartesian(ylim = c(0,.15))+ylab("Ketone body presence")+xlab("TNC:NDF")+ggtitle("L")
Ketone_PREV_ap_ndf_plot=Ketone_PREV_ap_ndf_plot+coord_cartesian(ylim = c(0,.15))+ylab("Ketone body presence")+xlab("P:NDF")+ggtitle("M")

Ketone_PREV_lipid_ndf_plot=Ketone_PREV_lipid_ndf_plot+xlab("Previous Lipid/NDF")

ggarrange(Ketone_PREV_KCAL_plot,Ketone_PREV_AP_plot,Ketone_PREV_lipid_plot,
          Ketone_PREV_TNC_plot, Ketone_PREV_NDF_plot,Ketone_PREV_NPE_plot, Ketone_PREV_npe_ab_plot,
          Ketone_PREV_lipid_ap_plot,
          Ketone_PREV_lipid_tnc_plot,Ketone_PREV_lipid_ndf_plot,Ketone_PREV_tnc_ap_plot,
          Ketone_PREV_tnc_ndf_plot,Ketone_PREV_ap_ndf_plot)


ggsave("Ketone previous day results zoomed 95CI MODES.pdf", units="in", width=10.366, height=7.730,dpi=300, device = "pdf")



deriv_plot(model=ucp_sgcor_PREV_KCAL,term="s(prev_day_total_kcal_using_ap_low_fermentation)",main="prev_day_total_kcal_using_ap_low_fermentation",eps=0.001,output="ucp_sgcor_PREV_KCAL_plot")
deriv_plot(model=ucp_sgcor_PREV_AP,term="s(prev_day_ap_kcal)",main="prev_day_ap_kcal",eps=0.001,output="ucp_sgcor_PREV_AP_plot")
deriv_plot(model=ucp_sgcor_PREV_lipid,term="s(prev_day_lipid_kcal)",main="prev_day_lipid_kcal",eps=0.001,output="ucp_sgcor_PREV_lipid_plot")
deriv_plot(model=ucp_sgcor_PREV_TNC,term="s(prev_day_tnc_kcal)",main="prev_day_tnc_kcal",eps=0.001,output="ucp_sgcor_PREV_TNC_plot")
deriv_plot(model=ucp_sgcor_PREV_NDF,term="s(prev_day_ndf_kcal_low)",main="prev_day_ndf_kcal_low",eps=0.001,output="ucp_sgcor_PREV_NDF_plot")
deriv_plot(model=ucp_sgcor_PREV_NPE,term="s(prev_day_total_kcal_npe_low_fermentation)",main="prev_day_total_kcal_npe_low_fermentation",eps=0.001,output="ucp_sgcor_PREV_NPE_plot")
deriv_plot(model=ucp_sgcor_PREV_npe_ab,term="s(prev_day_npe_ab)",main="prev_day_npe_ab",eps=0.001,output="ucp_sgcor_PREV_npe_ab_plot")
deriv_plot(model=ucp_sgcor_PREV_lipid_ap,term="s((PREV_lipid_ap))",main="(PREV_lipid_ap)",eps=0.001,output="ucp_sgcor_PREV_lipid_ap_plot")
deriv_plot(model=ucp_sgcor_PREV_lipid_tnc,term="s((PREV_lipid_tnc))",main="(PREV_lipid_tnc)",eps=0.001,output="ucp_sgcor_PREV_lipid_tnc_plot")
deriv_plot(model=ucp_sgcor_PREV_lipid_ndf,term="s((PREV_lipid_ndf))",main="(PREV_lipid_ndf)",eps=0.001,output="ucp_sgcor_PREV_lipid_ndf_plot")
deriv_plot(model=ucp_sgcor_PREV_tnc_ap,term="s(PREV_tnc_ap)",main="PREV_tnc_ap",eps=0.001,output="ucp_sgcor_PREV_tnc_ap_plot")
deriv_plot(model=ucp_sgcor_PREV_tnc_ndf,term="s(PREV_tnc_ndf)",main="PREV_tnc_ndf",eps=0.001,output="ucp_sgcor_PREV_tnc_ndf_plot")
deriv_plot(model=ucp_sgcor_PREV_ap_ndf,term="s(PREV_ap_ndf)",main="PREV_ap_ndf",eps=0.001,output="ucp_sgcor_PREV_ap_ndf_plot")

ucp_sgcor_PREV_KCAL_plot=ucp_sgcor_PREV_KCAL_plot+coord_cartesian(ylim = c(0,8000))+ylab("C-peptide of insulin\n(pg/ml)")+xlab("Total caloric intake")+ggtitle("A")  + geom_vline(xintercept = 2521.21, color = "gray")+ geom_vline(xintercept = c(842, 1773, 2213), color = "gray", linetype = "dashed")
ucp_sgcor_PREV_AP_plot=ucp_sgcor_PREV_AP_plot+coord_cartesian(ylim = c(0,8000))+ylab("C-peptide of insulin\n(pg/ml)")+xlab("P")+ggtitle("B")
ucp_sgcor_PREV_lipid_plot=ucp_sgcor_PREV_lipid_plot+coord_cartesian(ylim = c(0,8000))+ylab("C-peptide of insulin\n(pg/ml)")+xlab("Lipid")+ggtitle("C")
ucp_sgcor_PREV_TNC_plot=ucp_sgcor_PREV_TNC_plot+coord_cartesian(ylim = c(0,8000))+ylab("C-peptide of insulin\n(pg/ml)")+xlab("TNC")+ggtitle("D")
ucp_sgcor_PREV_NDF_plot=ucp_sgcor_PREV_NDF_plot+coord_cartesian(ylim = c(0,8000))+ylab("C-peptide of insulin\n(pg/ml)")+xlab("NDF")+ggtitle("E")
ucp_sgcor_PREV_NPE_plot=ucp_sgcor_PREV_NPE_plot+coord_cartesian(ylim = c(0,8000))+ylab("C-peptide of insulin\n(pg/ml)")+xlab("NPe")+ggtitle("F") + geom_vline(xintercept = 2265.66, color = "gray")+ geom_vline(xintercept = c(631, 795, 1004, 1488, 1494, 1863, 1988, 2844), color = "gray", linetype = "dashed")
ucp_sgcor_PREV_npe_ab_plot=ucp_sgcor_PREV_npe_ab_plot+coord_cartesian(ylim = c(0,8000))+ylab("C-peptide of insulin\n(pg/ml)")+xlab("NPe:P")+ggtitle("G") + geom_vline(xintercept = 9.71, color = "gray")+ geom_vline(xintercept = c(3.21), color = "gray", linetype = "dashed")
ucp_sgcor_PREV_lipid_ap_plot=ucp_sgcor_PREV_lipid_ap_plot+coord_cartesian(ylim = c(0,8000))+ylab("C-peptide of insulin\n(pg/ml)")+xlab("Lipid:P")+ggtitle("H")
ucp_sgcor_PREV_lipid_tnc_plot=ucp_sgcor_PREV_lipid_tnc_plot+coord_cartesian(ylim = c(0,8000))+ylab("C-peptide of insulin\n(pg/ml)")+xlab("Lipid:TNC")+ggtitle("I")
ucp_sgcor_PREV_lipid_ndf_plot=ucp_sgcor_PREV_lipid_ndf_plot+coord_cartesian(ylim = c(0,8000))+ylab("C-peptide of insulin\n(pg/ml)")+xlab("Lipid:NDF")+ggtitle("J")
ucp_sgcor_PREV_tnc_ap_plot=ucp_sgcor_PREV_tnc_ap_plot+coord_cartesian(ylim = c(0,8000))+ylab("C-peptide of insulin\n(pg/ml)")+xlab("TNC:P")+ggtitle("K") + geom_vline(xintercept = 6.69, color = "gray")+ geom_vline(xintercept = c(2.24, 2.3), color = "gray", linetype = "dashed")
ucp_sgcor_PREV_tnc_ndf_plot=ucp_sgcor_PREV_tnc_ndf_plot+coord_cartesian(ylim = c(0,8000))+ylab("C-peptide of insulin\n(pg/ml)")+xlab("TNC:NDF")+ggtitle("L")
ucp_sgcor_PREV_ap_ndf_plot=ucp_sgcor_PREV_ap_ndf_plot+coord_cartesian(ylim = c(0,8000))+ylab("C-peptide of insulin\n(pg/ml)")+xlab("P:NDF")+ggtitle("M")


ucp_sgcor_PREV_npe_ab_plot=ucp_sgcor_PREV_npe_ab_plot+ylab("UCP")+xlab("Previous NPE/Ap")




ggarrange(ucp_sgcor_PREV_KCAL_plot,ucp_sgcor_PREV_AP_plot,ucp_sgcor_PREV_lipid_plot,
          ucp_sgcor_PREV_TNC_plot, ucp_sgcor_PREV_NDF_plot,ucp_sgcor_PREV_NPE_plot, ucp_sgcor_PREV_npe_ab_plot,
          ucp_sgcor_PREV_lipid_ap_plot,
          ucp_sgcor_PREV_lipid_tnc_plot,ucp_sgcor_PREV_lipid_ndf_plot,ucp_sgcor_PREV_tnc_ap_plot,
          ucp_sgcor_PREV_tnc_ndf_plot,ucp_sgcor_PREV_ap_ndf_plot)
ggsave("UCP previous day results zoomed 95CI MODES.pdf", units="in", width=10.366, height=7.730,dpi=300, device = "pdf")

deriv_plot(model=dn15_result_PREV_KCAL,term="s(prev_day_total_kcal_using_ap_low_fermentation)",main="prev_day_total_kcal_using_ap_low_fermentation",eps=0.001,output="dn15_result_PREV_KCAL_plot")
deriv_plot(model=dn15_result_PREV_AP,term="s(prev_day_ap_kcal)",main="prev_day_ap_kcal",eps=0.001,output="dn15_result_PREV_AP_plot")
deriv_plot(model=dn15_result_PREV_lipid,term="s(prev_day_lipid_kcal)",main="prev_day_lipid_kcal",eps=0.001,output="dn15_result_PREV_lipid_plot")
deriv_plot(model=dn15_result_PREV_TNC,term="s(prev_day_tnc_kcal)",main="prev_day_tnc_kcal",eps=0.001,output="dn15_result_PREV_TNC_plot")
deriv_plot(model=dn15_result_PREV_NDF,term="s(prev_day_ndf_kcal_low)",main="prev_day_ndf_kcal_low",eps=0.001,output="dn15_result_PREV_NDF_plot")
deriv_plot(model=dn15_result_PREV_NPE,term="s(prev_day_total_kcal_npe_low_fermentation)",main="prev_day_total_kcal_npe_low_fermentation",eps=0.001,output="dn15_result_PREV_NPE_plot")
deriv_plot(model=dn15_result_PREV_npe_ab,term="s(prev_day_npe_ab)",main="prev_day_npe_ab",eps=0.001,output="dn15_result_PREV_npe_ab_plot")
deriv_plot(model=dn15_result_PREV_lipid_ap,term="s((PREV_lipid_ap))",main="(PREV_lipid_ap)",eps=0.001,output="dn15_result_PREV_lipid_ap_plot")
deriv_plot(model=dn15_result_PREV_lipid_tnc,term="s((PREV_lipid_tnc))",main="(PREV_lipid_tnc)",eps=0.001,output="dn15_result_PREV_lipid_tnc_plot")
deriv_plot(model=dn15_result_PREV_lipid_ndf,term="s((PREV_lipid_ndf))",main="(PREV_lipid_ndf)",eps=0.001,output="dn15_result_PREV_lipid_ndf_plot")
deriv_plot(model=dn15_result_PREV_tnc_ap,term="s(PREV_tnc_ap)",main="PREV_tnc_ap",eps=0.001,output="dn15_result_PREV_tnc_ap_plot")
deriv_plot(model=dn15_result_PREV_tnc_ndf,term="s(PREV_tnc_ndf)",main="PREV_tnc_ndf",eps=0.001,output="dn15_result_PREV_tnc_ndf_plot")
deriv_plot(model=dn15_result_PREV_ap_ndf,term="s(PREV_ap_ndf)",main="PREV_ap_ndf",eps=0.001,output="dn15_result_PREV_ap_ndf_plot")

dn15_result_PREV_KCAL_plot=dn15_result_PREV_KCAL_plot+coord_cartesian(ylim = c(-5,5))+ylab(expression(paste(delta^15,"N(‰)")))+xlab("Total caloric intake")+ggtitle("A")
dn15_result_PREV_AP_plot=dn15_result_PREV_AP_plot+coord_cartesian(ylim = c(-5,5))+ylab(expression(paste(delta^15,"N(‰)")))+xlab("P")+ggtitle("B")
dn15_result_PREV_lipid_plot=dn15_result_PREV_lipid_plot+coord_cartesian(ylim = c(-5,5))+ylab(expression(paste(delta^15,"N(‰)")))+xlab("Lipid")+ggtitle("C")
dn15_result_PREV_TNC_plot=dn15_result_PREV_TNC_plot+coord_cartesian(ylim = c(-5,5))+ylab(expression(paste(delta^15,"N(‰)")))+xlab("TNC")+ggtitle("D")
dn15_result_PREV_NDF_plot=dn15_result_PREV_NDF_plot+coord_cartesian(ylim = c(-5,5))+ylab(expression(paste(delta^15,"N(‰)")))+xlab("NDF")+ggtitle("E")
dn15_result_PREV_NPE_plot=dn15_result_PREV_NPE_plot+coord_cartesian(ylim = c(-5,5))+ylab(expression(paste(delta^15,"N(‰)")))+xlab("NPe")+ggtitle("F")
dn15_result_PREV_npe_ab_plot=dn15_result_PREV_npe_ab_plot+coord_cartesian(ylim = c(-5,5))+ylab(expression(paste(delta^15,"N(‰)")))+xlab("NPe:P")+ggtitle("G")
dn15_result_PREV_lipid_ap_plot=dn15_result_PREV_lipid_ap_plot+coord_cartesian(ylim = c(-5,5))+ylab(expression(paste(delta^15,"N(‰)")))+xlab("Lipid:P")+ggtitle("H")
dn15_result_PREV_lipid_tnc_plot=dn15_result_PREV_lipid_tnc_plot+coord_cartesian(ylim = c(-5,5))+ylab(expression(paste(delta^15,"N(‰)")))+xlab("Lipid:TNC")+ggtitle("I")
dn15_result_PREV_lipid_ndf_plot=dn15_result_PREV_lipid_ndf_plot+coord_cartesian(ylim = c(-5,5))+ylab(expression(paste(delta^15,"N(‰)")))+xlab("Lipid:NDF")+ggtitle("J")
dn15_result_PREV_tnc_ap_plot=dn15_result_PREV_tnc_ap_plot+coord_cartesian(ylim = c(-5,5))+ylab(expression(paste(delta^15,"N(‰)")))+xlab("TNC:P")+ggtitle("K")
dn15_result_PREV_tnc_ndf_plot=dn15_result_PREV_tnc_ndf_plot+coord_cartesian(ylim = c(-5,5))+ylab(expression(paste(delta^15,"N(‰)")))+xlab("TNC:NDF")+ggtitle("L")
dn15_result_PREV_ap_ndf_plot=dn15_result_PREV_ap_ndf_plot+coord_cartesian(ylim = c(-5,5))+ylab(expression(paste(delta^15,"N(‰)")))+xlab("P:NDF")+ggtitle("M")

dn15_result_PREV_KCAL_plot=dn15_result_PREV_KCAL_plot+coord_cartesian(ylim = c(-5,5)) + geom_vline(xintercept = 2521.21, color = "gray")+ geom_vline(xintercept = c(842, 1773, 2213), color = "gray", linetype = "dashed")
dn15_result_PREV_AP_plot=dn15_result_PREV_AP_plot+coord_cartesian(ylim = c(-5,5))
dn15_result_PREV_lipid_plot=dn15_result_PREV_lipid_plot+coord_cartesian(ylim = c(-5,5))
dn15_result_PREV_TNC_plot=dn15_result_PREV_TNC_plot+coord_cartesian(ylim = c(-5,5))+ geom_vline(xintercept = 1556.79, color = "gray")+ geom_vline(xintercept = 887, color = "gray", linetype = "dashed")
dn15_result_PREV_NDF_plot=dn15_result_PREV_NDF_plot+coord_cartesian(ylim = c(-5,5))
dn15_result_PREV_NPE_plot=dn15_result_PREV_NPE_plot+coord_cartesian(ylim = c(-5,5))+ geom_vline(xintercept = 2265.66, color = "gray")+ geom_vline(xintercept = c(631, 795, 1004, 1488, 1494, 1863, 1988, 2844), color = "gray", linetype = "dashed")
dn15_result_PREV_npe_ab_plot=dn15_result_PREV_npe_ab_plot+coord_cartesian(ylim = c(-5,5))
dn15_result_PREV_lipid_ap_plot=dn15_result_PREV_lipid_ap_plot+coord_cartesian(ylim = c(-5,5))
dn15_result_PREV_lipid_tnc_plot=dn15_result_PREV_lipid_tnc_plot+coord_cartesian(ylim = c(-5,5))
dn15_result_PREV_lipid_ndf_plot=dn15_result_PREV_lipid_ndf_plot+coord_cartesian(ylim = c(-5,5))
dn15_result_PREV_tnc_ap_plot=dn15_result_PREV_tnc_ap_plot+coord_cartesian(ylim = c(-5,5))
dn15_result_PREV_tnc_ndf_plot=dn15_result_PREV_tnc_ndf_plot+coord_cartesian(ylim = c(-5,5))
dn15_result_PREV_ap_ndf_plot=dn15_result_PREV_ap_ndf_plot+coord_cartesian(ylim = c(-5,5))


dn15_result_PREV_npe_ab_plot=dn15_result_PREV_npe_ab_plot+xlab("Previous NpE/Ap")
ggarrange(dn15_result_PREV_KCAL_plot,dn15_result_PREV_AP_plot,dn15_result_PREV_lipid_plot,
          dn15_result_PREV_TNC_plot, dn15_result_PREV_NDF_plot,dn15_result_PREV_NPE_plot, dn15_result_PREV_npe_ab_plot,
          dn15_result_PREV_lipid_ap_plot,
          dn15_result_PREV_lipid_tnc_plot,dn15_result_PREV_lipid_ndf_plot,dn15_result_PREV_tnc_ap_plot,
          dn15_result_PREV_tnc_ndf_plot,dn15_result_PREV_ap_ndf_plot)
ggsave("dn15 previous day results zoomed 95CI axis MODES.pdf", units="in", width=10.366, height=7.730,dpi=300, device = "pdf")



deriv_plot(model=urea_sgcor_PREV_KCAL,term="s(prev_day_total_kcal_using_ap_low_fermentation)",main="prev_day_total_kcal_using_ap_low_fermentation",eps=0.001,output="urea_sgcor_PREV_KCAL_plot")
deriv_plot(model=urea_sgcor_PREV_AP,term="s(prev_day_ap_kcal)",main="prev_day_ap_kcal",eps=0.001,output="urea_sgcor_PREV_AP_plot")
deriv_plot(model=urea_sgcor_PREV_lipid,term="s(prev_day_lipid_kcal)",main="prev_day_lipid_kcal",eps=0.001,output="urea_sgcor_PREV_lipid_plot")
deriv_plot(model=urea_sgcor_PREV_TNC,term="s(prev_day_tnc_kcal)",main="prev_day_tnc_kcal",eps=0.001,output="urea_sgcor_PREV_TNC_plot")
deriv_plot(model=urea_sgcor_PREV_NDF,term="s(prev_day_ndf_kcal_low)",main="prev_day_ndf_kcal_low",eps=0.001,output="urea_sgcor_PREV_NDF_plot")
deriv_plot(model=urea_sgcor_PREV_NPE,term="s(prev_day_total_kcal_npe_low_fermentation)",main="prev_day_total_kcal_npe_low_fermentation",eps=0.001,output="urea_sgcor_PREV_NPE_plot")
deriv_plot(model=urea_sgcor_PREV_npe_ab,term="s(prev_day_npe_ab)",main="prev_day_npe_ab",eps=0.001,output="urea_sgcor_PREV_npe_ab_plot")
deriv_plot(model=urea_sgcor_PREV_lipid_ap,term="s((PREV_lipid_ap))",main="(PREV_lipid_ap)",eps=0.001,output="urea_sgcor_PREV_lipid_ap_plot")
deriv_plot(model=urea_sgcor_PREV_lipid_tnc,term="s((PREV_lipid_tnc))",main="(PREV_lipid_tnc)",eps=0.001,output="urea_sgcor_PREV_lipid_tnc_plot")
deriv_plot(model=urea_sgcor_PREV_lipid_ndf,term="s((PREV_lipid_ndf))",main="(PREV_lipid_ndf)",eps=0.001,output="urea_sgcor_PREV_lipid_ndf_plot")
deriv_plot(model=urea_sgcor_PREV_tnc_ap,term="s(PREV_tnc_ap)",main="PREV_tnc_ap",eps=0.001,output="urea_sgcor_PREV_tnc_ap_plot")
deriv_plot(model=urea_sgcor_PREV_tnc_ndf,term="s(PREV_tnc_ndf)",main="PREV_tnc_ndf",eps=0.001,output="urea_sgcor_PREV_tnc_ndf_plot")
deriv_plot(model=urea_sgcor_PREV_ap_ndf,term="s(PREV_ap_ndf)",main="PREV_ap_ndf",eps=0.001,output="urea_sgcor_PREV_ap_ndf_plot")

urea_sgcor_PREV_KCAL_plot=urea_sgcor_PREV_KCAL_plot+coord_cartesian(ylim = c(0,2000))+ylab("Urea concentration\n(mg/ml)")+xlab("Total caloric intake")+ggtitle("A")
urea_sgcor_PREV_AP_plot=urea_sgcor_PREV_AP_plot+coord_cartesian(ylim = c(0,2000))+ylab("Urea concentration\n(mg/ml)")+xlab("P")+ggtitle("B")
urea_sgcor_PREV_lipid_plot=urea_sgcor_PREV_lipid_plot+coord_cartesian(ylim = c(0,2000))+ylab("Urea concentration\n(mg/ml)")+xlab("Lipid")+ggtitle("C")
urea_sgcor_PREV_TNC_plot=urea_sgcor_PREV_TNC_plot+coord_cartesian(ylim = c(0,2000))+ylab("Urea concentration\n(mg/ml)")+xlab("TNC")+ggtitle("D")
urea_sgcor_PREV_NDF_plot=urea_sgcor_PREV_NDF_plot+coord_cartesian(ylim = c(0,2000))+ylab("Urea concentration\n(mg/ml)")+xlab("NDF")+ggtitle("E")
urea_sgcor_PREV_NPE_plot=urea_sgcor_PREV_NPE_plot+coord_cartesian(ylim = c(0,2000))+ylab("Urea concentration\n(mg/ml)")+xlab("NPe")+ggtitle("F")
urea_sgcor_PREV_npe_ab_plot=urea_sgcor_PREV_npe_ab_plot+coord_cartesian(ylim = c(0,2000))+ylab("Urea concentration\n(mg/ml)")+xlab("NPe:P")+ggtitle("G")
urea_sgcor_PREV_lipid_ap_plot=urea_sgcor_PREV_lipid_ap_plot+coord_cartesian(ylim = c(0,2000))+ylab("Urea concentration\n(mg/ml)")+xlab("Lipid:P")+ggtitle("H")
urea_sgcor_PREV_lipid_tnc_plot=urea_sgcor_PREV_lipid_tnc_plot+coord_cartesian(ylim = c(0,2000))+ylab("Urea concentration\n(mg/ml)")+xlab("Lipid:TNC")+ggtitle("I")
urea_sgcor_PREV_lipid_ndf_plot=urea_sgcor_PREV_lipid_ndf_plot+coord_cartesian(ylim = c(0,2000))+ylab("Urea concentration\n(mg/ml)")+xlab("Lipid:NDF")+ggtitle("J")
urea_sgcor_PREV_tnc_ap_plot=urea_sgcor_PREV_tnc_ap_plot+coord_cartesian(ylim = c(0,2000))+ylab("Urea concentration\n(mg/ml)")+xlab("TNC:P")+ggtitle("K")
urea_sgcor_PREV_tnc_ndf_plot=urea_sgcor_PREV_tnc_ndf_plot+coord_cartesian(ylim = c(0,2000))+ylab("Urea concentration\n(mg/ml)")+xlab("TNC:NDF")+ggtitle("L")
urea_sgcor_PREV_ap_ndf_plot=urea_sgcor_PREV_ap_ndf_plot+coord_cartesian(ylim = c(0,2000))+ylab("Urea concentration\n(mg/ml)")+xlab("P:NDF")+ggtitle("M")

urea_sgcor_PREV_KCAL_plot=urea_sgcor_PREV_KCAL_plot+coord_cartesian(ylim = c(0,2000)) + geom_vline(xintercept = 2521.21, color = "gray")+ geom_vline(xintercept = 842, color = "gray", linetype = "dashed")+ geom_vline(xintercept = 1773, color = "gray", linetype = "dashed")+ geom_vline(xintercept = 2213, color = "gray", linetype = "dashed")
urea_sgcor_PREV_AP_plot=urea_sgcor_PREV_AP_plot+coord_cartesian(ylim = c(0,2000))
urea_sgcor_PREV_lipid_plot=urea_sgcor_PREV_lipid_plot+coord_cartesian(ylim = c(0,2000))  + geom_vline(xintercept = 453.09, color = "gray")+ geom_vline(xintercept = 143, color = "gray", linetype = "dashed")
urea_sgcor_PREV_TNC_plot=urea_sgcor_PREV_TNC_plot+coord_cartesian(ylim = c(0,2000)) + geom_vline(xintercept = 1556.79, color = "gray")+ geom_vline(xintercept = 887, color = "gray", linetype = "dashed")
urea_sgcor_PREV_NDF_plot=urea_sgcor_PREV_NDF_plot+coord_cartesian(ylim = c(0,2000))
urea_sgcor_PREV_NPE_plot=urea_sgcor_PREV_NPE_plot+coord_cartesian(ylim = c(0,2000)) + geom_vline(xintercept = 2265.66, color = "gray")+ geom_vline(xintercept = c(631, 795, 1004, 1488, 1494, 1863, 1988, 2844), color = "gray", linetype = "dashed")
urea_sgcor_PREV_npe_ab_plot=urea_sgcor_PREV_npe_ab_plot+coord_cartesian(ylim = c(0,2000))  + geom_vline(xintercept = 9.71, color = "gray")+ geom_vline(xintercept = c(3.21), color = "gray", linetype = "dashed")
urea_sgcor_PREV_lipid_ap_plot=urea_sgcor_PREV_lipid_ap_plot+coord_cartesian(ylim = c(0,2000))  + geom_vline(xintercept = 1.9, color = "gray")+ geom_vline(xintercept = c(0.39), color = "gray", linetype = "dashed")
urea_sgcor_PREV_lipid_tnc_plot=urea_sgcor_PREV_lipid_tnc_plot+coord_cartesian(ylim = c(0,2000))  + geom_vline(xintercept = 0.32, color = "gray")+ geom_vline(xintercept = c(0.18), color = "gray", linetype = "dashed")
urea_sgcor_PREV_lipid_ndf_plot=urea_sgcor_PREV_lipid_ndf_plot+coord_cartesian(ylim = c(0,2000))   + geom_vline(xintercept = 1.85, color = "gray")+ geom_vline(xintercept = c(0.91), color = "gray", linetype = "dashed")
urea_sgcor_PREV_tnc_ap_plot=urea_sgcor_PREV_tnc_ap_plot+coord_cartesian(ylim = c(0,2000))  + geom_vline(xintercept = 6.69, color = "gray")+ geom_vline(xintercept = c(2.24, 2.3), color = "gray", linetype = "dashed")
urea_sgcor_PREV_tnc_ndf_plot=urea_sgcor_PREV_tnc_ndf_plot+coord_cartesian(ylim = c(0,2000))  + geom_vline(xintercept = 6.82, color = "gray")+ geom_vline(xintercept = c(5.38), color = "gray", linetype = "dashed")
urea_sgcor_PREV_ap_ndf_plot=urea_sgcor_PREV_ap_ndf_plot+coord_cartesian(ylim = c(0,2000))  + geom_vline(xintercept = 1.43, color = "gray")+ geom_vline(xintercept = c(0.7), color = "gray", linetype = "dashed")


urea_sgcor_PREV_npe_ab_plot=urea_sgcor_PREV_npe_ab_plot+xlab("Previous day NPE/AP")
ggarrange(urea_sgcor_PREV_KCAL_plot,urea_sgcor_PREV_AP_plot,urea_sgcor_PREV_lipid_plot,
          urea_sgcor_PREV_TNC_plot, urea_sgcor_PREV_NDF_plot,urea_sgcor_PREV_NPE_plot, urea_sgcor_PREV_npe_ab_plot,
          urea_sgcor_PREV_lipid_ap_plot,
          urea_sgcor_PREV_lipid_tnc_plot,urea_sgcor_PREV_lipid_ndf_plot,urea_sgcor_PREV_tnc_ap_plot,
          urea_sgcor_PREV_tnc_ndf_plot,urea_sgcor_PREV_ap_ndf_plot)
ggsave("Urea previous day results zoomed axis 95CI MODES.pdf", units="in", width=10.366, height=7.730,dpi=300, device = "pdf")



deriv_plot(model=ucp_sgcor_FAI_model,term="s(fai)",main="fai",eps=0.001,output="ucp_sgcor_FAI_deriv")
ucp_sgcor_FAI_deriv=ucp_sgcor_FAI_deriv+ylab("C-peptide of insulin (pg/ml)")+xlab("FAI")+ggtitle("A")

deriv_plot(model=data_FAI_model,term="s(fai)",main="fai",eps=0.001,output="dn15_FAI_deriv")
dn15_FAI_deriv=dn15_FAI_deriv+ylab(expression(paste(delta^15,"N(‰)")))+xlab("FAI")+ggtitle("D")

deriv_plot(model=pos_neg_dn_FAI_model,term="s(fai)",main="fai",eps=0.001,output="ket_FAI_deriv")
ket_FAI_deriv=ket_FAI_deriv+ylab("Ketone body presence")+xlab("FAI")+ggtitle("B")

deriv_plot(model=urea_FAI_model,term="s(fai)",main="fai",eps=0.001,output="urea_sgcore_KCAL_FAI_deriv")

urea_sgcore_KCAL_FAI_deriv=urea_sgcore_KCAL_FAI_deriv+ylab("Urea concentration (mg/ml)")+xlab("FAI")+ggtitle("C")


ggarrange(ucp_sgcor_FAI_deriv,ket_FAI_deriv,urea_sgcore_KCAL_FAI_deriv,dn15_FAI_deriv)
ggsave("FAI and metabolites results zoomed axis 95CI.pdf", dpi=300, device = "pdf")

valuesDF=list(                       
VOI_dn15_result_PREV_KCAL,                     
VOI_dn15_result_PREV_NPE,                      
VOI_dn15_result_PREV_TNC,                      
VOI_Ketone_PREV_lipid_tnc,                     
VOI_ucp_sgcor_PREV_KCAL,                       
VOI_ucp_sgcor_PREV_NPE,                        
VOI_ucp_sgcor_PREV_npe_ab,                     
VOI_ucp_sgcor_PREV_tnc_ap,                     
VOI_urea_sgcor_PREV_ap_ndf,                    
VOI_urea_sgcor_PREV_KCAL,                      
VOI_urea_sgcor_PREV_lipid,                     
VOI_urea_sgcor_PREV_lipid_ap,                  
VOI_urea_sgcor_PREV_lipid_ndf,                 
VOI_urea_sgcor_PREV_lipid_tnc,                 
VOI_urea_sgcor_PREV_NPE,                       
VOI_urea_sgcor_PREV_npe_ab,                    
VOI_urea_sgcor_PREV_TNC,                       
VOI_urea_sgcor_PREV_tnc_ap,                    
VOI_urea_sgcor_PREV_tnc_ndf)  

valuesDF <- setNames(valuesDF,c("VOI_dn15_result_PREV_KCAL",                     
                                "VOI_dn15_result_PREV_NPE",                      
                                "VOI_dn15_result_PREV_TNC",                      
                                "VOI_Ketone_PREV_lipid_tnc",                     
                                "VOI_ucp_sgcor_PREV_KCAL",                       
                                "VOI_ucp_sgcor_PREV_NPE",                        
                                "VOI_ucp_sgcor_PREV_npe_ab",                     
                                "VOI_ucp_sgcor_PREV_tnc_ap",                     
                                "VOI_urea_sgcor_PREV_ap_ndf",                    
                                "VOI_urea_sgcor_PREV_KCAL",                      
                                "VOI_urea_sgcor_PREV_lipid",                     
                                "VOI_urea_sgcor_PREV_lipid_ap",                  
                                "VOI_urea_sgcor_PREV_lipid_ndf",                 
                                "VOI_urea_sgcor_PREV_lipid_tnc",                 
                                "VOI_urea_sgcor_PREV_NPE",                       
                                "VOI_urea_sgcor_PREV_npe_ab",                    
                                "VOI_urea_sgcor_PREV_TNC",                       
                                "VOI_urea_sgcor_PREV_tnc_ap",                    
                                "VOI_urea_sgcor_PREV_tnc_ndf"))

save(valuesDF,file="valuesDF_95.RData")

deriv_plot2(model=dn15_result_PREV_KCAL,term="s(prev_day_total_kcal_using_ap_low_fermentation)",main="prev_day_total_kcal_using_ap_low_fermentation",eps=0.001,output="dn15_result_PREV_KCAL_deriv")
deriv_plot(model=dn15_result_PREV_KCAL,term="s(prev_day_total_kcal_using_ap_low_fermentation)",main="prev_day_total_kcal_using_ap_low_fermentation",eps=0.001,output="dn15_result_PREV_KCAL_plot")

dn15_KCAL_deriv_plot=dn15_result_PREV_KCAL_deriv+
  ylab(expression(paste(italic("f(x)'"))))+xlab("Total caloric intake")+ggtitle("A")
dn15_KCAL_deriv_reg_plot=dn15_result_PREV_KCAL_plot+
  ylab(expression(paste(delta^15,"N(‰)")))+xlab("Total caloric intake")+ggtitle("B")

ggarrange(dn15_KCAL_deriv_plot,dn15_KCAL_deriv_reg_plot,nrow=2, legend = "top")
ggsave("Supplementary derivative example plot.pdf",dpi=300, device = "pdf")


deriv_plot(model=dn15_result_PREV_TNC,term="s(prev_day_tnc_kcal)",main="prev_day_tnc_kcal",eps=0.001,output="dn15_result_PREV_TNC_plot", confidence = 95)
dn15_TNC_cutoff_plot=dn15_result_PREV_TNC_plot+
  coord_cartesian(ylim = c(-5,5), xlim = c(0,5000))+ylab(expression(paste(delta^15,"N(‰)")))+xlab("TNC")+ggtitle("A") + geom_vline(xintercept = 870, color = "dodgerblue", linetype = "dashed") + geom_vline(xintercept = 1556.79, color = "gray")+ geom_vline(xintercept = 887, color = "gray", linetype = "dashed")



deriv_plot(model=urea_sgcor_PREV_TNC,term="s(prev_day_tnc_kcal)",main="prev_day_tnc_kcal",eps=0.001,output="urea_sgcor_PREV_TNC_plot", confidence = 95)

urea_TNC_cutoff_plot=urea_sgcor_PREV_TNC_plot+
  coord_cartesian( xlim = c(0,5000)) + ylab("Urea concentration (mg/ml)")+xlab("TNC")+ggtitle("B") + geom_vline(xintercept = 4250, color = "dodgerblue", linetype = "dashed") + geom_vline(xintercept = 1556.79, color = "gray")+ geom_vline(xintercept = 887, color = "gray", linetype = "dashed")




ggarrange(dn15_TNC_cutoff_plot, urea_TNC_cutoff_plot,nrow=2)
ggsave("Urea against dn15 TNC 95CI 2023.pdf",dpi=300, device = "pdf")


all_models=list(Ketone_PREV_KCAL,Ketone_PREV_AP,Ketone_PREV_lipid,
          Ketone_PREV_TNC, Ketone_PREV_NDF,Ketone_PREV_NPE, Ketone_PREV_npe_ab,
          Ketone_PREV_lipid_ap,
          Ketone_PREV_lipid_tnc,Ketone_PREV_lipid_ndf,Ketone_PREV_tnc_ap,
          Ketone_PREV_tnc_ndf,Ketone_PREV_ap_ndf, ucp_sgcor_PREV_KCAL,ucp_sgcor_PREV_AP,ucp_sgcor_PREV_lipid,
          ucp_sgcor_PREV_TNC, ucp_sgcor_PREV_NDF,ucp_sgcor_PREV_NPE, ucp_sgcor_PREV_npe_ab,
          ucp_sgcor_PREV_lipid_ap,
          ucp_sgcor_PREV_lipid_tnc,ucp_sgcor_PREV_lipid_ndf,ucp_sgcor_PREV_tnc_ap,
          ucp_sgcor_PREV_tnc_ndf,ucp_sgcor_PREV_ap_ndf, dn15_result_PREV_KCAL,dn15_result_PREV_AP,dn15_result_PREV_lipid,
          dn15_result_PREV_TNC, dn15_result_PREV_NDF,dn15_result_PREV_NPE, dn15_result_PREV_npe_ab,
          dn15_result_PREV_lipid_ap,
          dn15_result_PREV_lipid_tnc,dn15_result_PREV_lipid_ndf,dn15_result_PREV_tnc_ap,
          dn15_result_PREV_tnc_ndf,dn15_result_PREV_ap_ndf, urea_sgcor_PREV_KCAL,urea_sgcor_PREV_AP,urea_sgcor_PREV_lipid,
          urea_sgcor_PREV_TNC, urea_sgcor_PREV_NDF,urea_sgcor_PREV_NPE, urea_sgcor_PREV_npe_ab,
          urea_sgcor_PREV_lipid_ap,
          urea_sgcor_PREV_lipid_tnc,urea_sgcor_PREV_lipid_ndf,urea_sgcor_PREV_tnc_ap,
          urea_sgcor_PREV_tnc_ndf,urea_sgcor_PREV_ap_ndf)

sumtable=function(model){
  
  temp=summary(model)
  temp_fix=data.frame(temp$fixed)
  temp_spline=data.frame(temp$splines)
  coeffs=rbind(temp_fix,temp_spline)
  Mod_name=paste(temp$formula$resp,trimws(strsplit(as.character(temp$formula$formula[3]), "+", fixed = TRUE)[[1]][3]), sep = " ~ ")
  coeffs=cbind(rownames(coeffs),coeffs)
  coeffs=cbind(rep(Mod_name,nrow(coeffs)),coeffs)
  colnames(coeffs)[1:2]=c("Model_ID","Covariate")
  rownames(coeffs)=NULL
  
  return(coeffs)
}

mytables=lapply(all_models,sumtable)

mytables=do.call(rbind,mytables)
mytables$Model_ID=gsub("ureasgcor", "Urea",mytables$Model_ID)
mytables$Model_ID=gsub("posnegdn ","Ketone presence" ,mytables$Model_ID)
mytables$Model_ID=gsub("ucpsgcor","C-peptide" ,mytables$Model_ID)
mytables$Model_ID=gsub("dn15result","d15N" ,mytables$Model_ID)
mytables$Model_ID=gsub("prev_day_total_kcal_using_ap_low_fermentation","Total kcal" ,mytables$Model_ID)
mytables$Model_ID=gsub("prev_day_ap_kcal","Protein" ,mytables$Model_ID)
mytables$Model_ID=gsub("prev_day_lipid_kcal","Lipid" ,mytables$Model_ID)
mytables$Model_ID=gsub("prev_day_tnc_kcal","TNC" ,mytables$Model_ID)
mytables$Model_ID=gsub("prev_day_ndf_kcal_low","NDF" ,mytables$Model_ID)
mytables$Model_ID=gsub("prev_day_total_kcal_npe_low_fermentation","Npe" ,mytables$Model_ID)
mytables$Model_ID=gsub("prev_day_npe_ab","Npe:P" ,mytables$Model_ID)
mytables$Model_ID=gsub("PREV_lipid_ap","Lipid:P" ,mytables$Model_ID)
mytables$Model_ID=gsub("PREV_lipid_tnc","Lipid:TNC" ,mytables$Model_ID)
mytables$Model_ID=gsub("PREV_lipid_ndf","Lipid:NDF" ,mytables$Model_ID)
mytables$Model_ID=gsub("PREV_tnc_ap","TNC:P" ,mytables$Model_ID)
mytables$Model_ID=gsub("PREV_tnc_ndf","TNC:NDF" ,mytables$Model_ID)
mytables$Model_ID=gsub("PREV_ap_ndf","P:NDF" ,mytables$Model_ID)
mytables$Model_ID=gsub("s(", "" ,mytables$Model_ID, fixed = TRUE)
mytables$Model_ID=gsub("(","" ,mytables$Model_ID, fixed = TRUE)
mytables$Model_ID=gsub(")","" ,mytables$Model_ID, fixed = TRUE)

write.csv(mytables,file="PNAS model summary tables.csv")

sumtable=function(model){
  
  temp=summary(model)
  temp_fix=data.frame(temp$fixed)
  temp_spline=data.frame(temp$splines)
  coeffs=rbind(temp_fix,temp_spline)
  Mod_name=paste(temp$formula$resp,trimws(strsplit(as.character(temp$formula$formula[3]), "+", fixed = TRUE)[[1]][1]), sep = " ~ ")
  coeffs=cbind(rownames(coeffs),coeffs)
  coeffs=cbind(rep(Mod_name,nrow(coeffs)),coeffs)
  colnames(coeffs)[1:2]=c("Model_ID","Covariate")
  rownames(coeffs)=NULL
  
  return(coeffs)
}

fai_models=list(ucp_sgcor_FAI_model,data_FAI_model,pos_neg_dn_FAI_model,urea_FAI_model)
fai_tables=lapply(fai_models,sumtable)

fai_tables=do.call(rbind,fai_tables)

fai_tables$Model_ID=gsub("fai", "FAI",fai_tables$Model_ID)
fai_tables$Model_ID=gsub("ureasgcor", "Urea",fai_tables$Model_ID)
fai_tables$Model_ID=gsub("posnegdn ","Ketone presence" ,fai_tables$Model_ID)
fai_tables$Model_ID=gsub("ucpsgcor","C-peptide" ,fai_tables$Model_ID)
fai_tables$Model_ID=gsub("dn15result","d15N" ,fai_tables$Model_ID)
fai_tables$Model_ID=gsub("s(", "" ,fai_tables$Model_ID, fixed = TRUE)
fai_tables$Model_ID=gsub("(","" ,fai_tables$Model_ID, fixed = TRUE)
fai_tables$Model_ID=gsub(")","" ,fai_tables$Model_ID, fixed = TRUE)
write.csv(fai_tables,file="PNAS FAI metabolite tables.csv")


Modes <- function(x) {
  ux <- unique(x)
  tab <- tabulate(match(x, ux))
  ux[tab == max(tab)]
}
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

aphist=ggplot(data, aes(x=ap_kcal)) + 
  geom_histogram(binwidth=1)
aphist=data.frame(ggplot_build(aphist)$data[[1]])
max(aphist$xmax[which(aphist$y==max(aphist$y))])
ggplot(data, aes(x=ap_kcal)) + 
  geom_histogram(binwidth=10)+geom_vline(aes(xintercept = max(aphist$xmax[which(aphist$y==max(aphist$y))]),color="red"))+
  theme_classic()+ geom_text(x=max(aphist$xmax[which(aphist$y==max(aphist$y))]), y=-3, label=paste(max(aphist$xmax[which(aphist$y==max(aphist$y))])), color="red")

lipidhist=ggplot(data, aes(x=lipid_kcal)) + 
  geom_histogram(binwidth=1)
lipidhist=data.frame(ggplot_build(lipidhist)$data[[1]])
max(lipidhist$xmax[which(lipidhist$y==max(lipidhist$y))])
ggplot(data, aes(x=lipid_kcal)) + 
  geom_histogram(binwidth=10)+geom_vline(aes(xintercept = max(lipidhist$xmax[which(lipidhist$y==max(lipidhist$y))]),color="red"))+
  theme_classic()+ geom_text(x=max(lipidhist$xmax[which(lipidhist$y==max(lipidhist$y))]), y=-3, label=paste(max(lipidhist$xmax[which(lipidhist$y==max(lipidhist$y))])), color="red")


tnchist=ggplot(data, aes(x=tnc_kcal)) + 
  geom_histogram(binwidth=1)
tnchist=data.frame(ggplot_build(tnchist)$data[[1]])
max(tnchist$xmax[which(tnchist$y==max(tnchist$y))])
ggplot(data, aes(x=tnc_kcal)) + 
  geom_histogram(binwidth=10)+geom_vline(aes(xintercept = max(tnchist$xmax[which(tnchist$y==max(tnchist$y))]),color="red"))+
  theme_classic()+ geom_text(x=max(tnchist$xmax[which(tnchist$y==max(tnchist$y))]), y=-1, label=paste(max(tnchist$xmax[which(tnchist$y==max(tnchist$y))])), color="red")


nfdhist=ggplot(data, aes(x=ndf_kcal_low)) + 
  geom_histogram(binwidth=1)#+geom_vline(aes(xintercept = getmode(round(data$ndf_kcal_low, digits = 2)),color="red"))+
  #theme_classic()+ geom_text(x=getmode(round(data$ndf_kcal_low, digits = 1)), y=-3, label=paste(getmode(round(data$ndf_kcal_low, digits = 2))), color="red")
nfdhist
nfdhist=data.frame(ggplot_build(nfdhist)$data[[1]])
min(nfdhist$xmax[which(nfdhist$count==max(nfdhist$count))])
ggplot(data, aes(x=ndf_kcal_low)) + 
  geom_histogram(binwidth=10)+geom_vline(aes(xintercept = max(nfdhist$xmax[which(nfdhist$y==max(nfdhist$y))]),color="red"))+
theme_classic()+ geom_text(x=max(nfdhist$xmax[which(nfdhist$y==max(nfdhist$y))]), y=-3, label=paste(max(nfdhist$xmax[which(nfdhist$y==max(nfdhist$y))])), color="red")


kcalhist=ggplot(data, aes(x=total_kcal_using_ap_low_fermentation)) + 
  geom_histogram(binwidth=1)
kcalhist=data.frame(ggplot_build(kcalhist)$data[[1]])
min(kcalhist$xmax[which(kcalhist$count==max(kcalhist$count))])
ggplot(data, aes(x=total_kcal_using_ap_low_fermentation)) + 
  geom_histogram(binwidth=10)+geom_vline(aes(xintercept = min(kcalhist$xmax[which(kcalhist$count==max(kcalhist$count))]),color="red"))+
  theme_classic()+ geom_text(x=min(kcalhist$xmax[which(kcalhist$count==max(kcalhist$count))]), y=-1, label=paste(min(kcalhist$xmax[which(kcalhist$count==max(kcalhist$count))])), color="red")


npehist=ggplot(data, aes(x=total_kcal_npe_low_fermentation)) + 
  geom_histogram(binwidth=1)
npehist=data.frame(ggplot_build(npehist)$data[[1]])
(npehist$xmax[which(npehist$count==max(npehist$count))])
ggplot(data, aes(x=total_kcal_npe_low_fermentation)) + 
  geom_histogram(binwidth=10)+geom_vline(aes(xintercept = min(npehist$xmax[which(npehist$y==max(npehist$y))]),color="red"))+
  theme_classic()+ geom_text(x=min(npehist$xmax[which(npehist$y==max(npehist$y))]), y=-.5, label=paste(min(npehist$xmax[which(npehist$y==max(npehist$y))])), color="red")


npeAPhist=ggplot(data, aes(x=npe_to_ap)) + 
  geom_histogram(binwidth=1)
npeAPhist=data.frame(ggplot_build(npeAPhist)$data[[1]])
(npeAPhist$xmax[which(npeAPhist$count==max(npeAPhist$count))])
ggplot(data, aes(x=npe_to_ap)) + 
  geom_histogram(binwidth=1)+geom_vline(aes(xintercept = min(npeAPhist$xmax[which(npeAPhist$y==max(npeAPhist$y))]),color="red"))+
  theme_classic()+ geom_text(x=min(npeAPhist$xmax[which(npeAPhist$y==max(npeAPhist$y))]), y=-.5, label=paste(min(npeAPhist$xmax[which(npeAPhist$y==max(npeAPhist$y))])), color="red")


LIPIDAPhist=ggplot(data, aes(x=lipid_to_ap)) + 
  geom_histogram(binwidth=1)
LIPIDAPhist=data.frame(ggplot_build(LIPIDAPhist)$data[[1]])
(LIPIDAPhist$xmax[which(LIPIDAPhist$count==max(LIPIDAPhist$count))])
ggplot(data, aes(x=lipid_to_ap)) + 
  geom_histogram(binwidth=1)+geom_vline(aes(xintercept = min(LIPIDAPhist$xmax[which(LIPIDAPhist$y==max(LIPIDAPhist$y))]),color="red"))+
  theme_classic()+ geom_text(x=min(LIPIDAPhist$xmax[which(LIPIDAPhist$y==max(LIPIDAPhist$y))]), y=-.5, label=paste(min(LIPIDAPhist$xmax[which(LIPIDAPhist$y==max(LIPIDAPhist$y))])), color="red")


LIPIDNDFhist=ggplot(data, aes(x=lipid_to_tnc)) + 
  geom_histogram(binwidth=1)
LIPIDNDFhist=data.frame(ggplot_build(LIPIDNDFhist)$data[[1]])
(LIPIDNDFhist$xmax[which(LIPIDNDFhist$count==max(LIPIDNDFhist$count))])
ggplot(data, aes(x=lipid_to_tnc)) + 
  geom_histogram(binwidth=1)+geom_vline(aes(xintercept = min(LIPIDNDFhist$xmax[which(LIPIDNDFhist$y==max(LIPIDNDFhist$y))]),color="red"))+
  theme_classic()+ geom_text(x=min(LIPIDNDFhist$xmax[which(LIPIDNDFhist$y==max(LIPIDNDFhist$y))]), y=-10, label=paste(min(LIPIDNDFhist$xmax[which(LIPIDNDFhist$y==max(LIPIDNDFhist$y))])), color="red")


LIPIDNDFhist=ggplot(data, aes(x=lipid_to_ndf_low)) + 
  geom_histogram(binwidth=1)
LIPIDNDFhist=data.frame(ggplot_build(LIPIDNDFhist)$data[[1]])
(LIPIDNDFhist$xmax[which(LIPIDNDFhist$count==max(LIPIDNDFhist$count))])
ggplot(data, aes(x=lipid_to_ndf_low)) + 
  geom_histogram(binwidth=1)+geom_vline(aes(xintercept = min(LIPIDNDFhist$xmax[which(LIPIDNDFhist$y==max(LIPIDNDFhist$y))]),color="red"))+
  theme_classic()+ geom_text(x=min(LIPIDNDFhist$xmax[which(LIPIDNDFhist$y==max(LIPIDNDFhist$y))]), y=-10, label=paste(min(LIPIDNDFhist$xmax[which(LIPIDNDFhist$y==max(LIPIDNDFhist$y))])), color="red")


tncAPhist=ggplot(data, aes(x=tnc_to_ap)) + 
  geom_histogram(binwidth=1)
tncAPhist=data.frame(ggplot_build(tncAPhist)$data[[1]])
(tncAPhist$xmax[which(tncAPhist$count==max(tncAPhist$count))])
ggplot(data, aes(x=tnc_to_ap)) + 
  geom_histogram(binwidth=1)+geom_vline(aes(xintercept = min(tncAPhist$xmax[which(tncAPhist$y==max(tncAPhist$y))]),color="red"))+
  theme_classic()+ geom_text(x=min(tncAPhist$xmax[which(tncAPhist$y==max(tncAPhist$y))]), y=-10, label=paste(min(tncAPhist$xmax[which(tncAPhist$y==max(tncAPhist$y))])), color="red")


tncNDFhist=ggplot(data, aes(x=tnc_to_ndf_low)) + 
  geom_histogram(binwidth=1)
tncNDFhist=data.frame(ggplot_build(tncNDFhist)$data[[1]])
(tncNDFhist$xmax[which(tncNDFhist$count==max(tncNDFhist$count))])
ggplot(data, aes(x=tnc_to_ndf_low)) + 
  geom_histogram(binwidth=1)+geom_vline(aes(xintercept = min(tncNDFhist$xmax[which(tncNDFhist$y==max(tncNDFhist$y))]),color="red"))+
  theme_classic()+ geom_text(x=min(tncNDFhist$xmax[which(tncNDFhist$y==max(tncNDFhist$y))]), y=-10, label=paste(min(tncNDFhist$xmax[which(tncNDFhist$y==max(tncNDFhist$y))])), color="red")


apNDFhist=ggplot(data, aes(x=ap_to_ndf_low)) + 
  geom_histogram(binwidth=1)
apNDFhist=data.frame(ggplot_build(apNDFhist)$data[[1]])
(apNDFhist$xmax[which(apNDFhist$count==max(apNDFhist$count))])
ggplot(data, aes(x=ap_to_ndf_low)) + 
  geom_histogram(binwidth=1)+geom_vline(aes(xintercept = min(apNDFhist$xmax[which(apNDFhist$y==max(apNDFhist$y))]),color="red"))+
  theme_classic()+ geom_text(x=min(apNDFhist$xmax[which(apNDFhist$y==max(apNDFhist$y))]), y=-10, label=paste(min(apNDFhist$xmax[which(apNDFhist$y==max(apNDFhist$y))])), color="red")


ggplot()+geom_point(data=data,aes(x=tnc_kcal, y=lipid_kcal,color=I("black")))+theme_classic()+geom_point(data=datasmall, aes(x=tnc_kcal, y=lipid_kcal, color="blue"))
