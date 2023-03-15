library(ggplot2)
library(ggpubr)
library(brms)

source("https://raw.githubusercontent.com/sealavi/GAM-first-derivative-functions/main/gam%20first%20derivitave%20plot%20function_general.R")

options(mc.cores = parallel::detectCores())


data$total_kcal_using_ap_low_fermentation=data$total_kcal_using_ap_low_fermentation+.001
data2=read.csv("C:/Users/salavi/Documents/20211031 To Send to Shauhin for analyses/20211031 To Send to Shauhin for analyses/nutrition_urine_biomarkers_by_follow_RSAB_29oct2021.csv")

KCAL<-brm(total_kcal_using_ap_low_fermentation~s(fai) + class_focal  + (1 +  fai + class_focal  | name_focal), 
          data = data[complete.cases(data[,c("total_kcal_using_ap_low_fermentation","fai")]),],
          family=Gamma(link="log"),iter = 4000, inits=0,chains = 2, 
          prior=c(prior(normal(0,10),class="Intercept"),
                  prior(normal(0,10),class="b"),
                  prior(gamma(.01,.01),class="shape"),
                  prior(normal(0,10),class="sd")),
          control = list(max_treedepth = 11, adapt_delta = .999))

summary(KCAL, digits=9)
plot(conditional_effects(KCAL,spaghetti=FALSE)) ##plots the regression lines/conditional effects
pp_check(KCAL)
loo_R2(KCAL)
bayes_R2(KCAL)

deriv_plot(model=KCAL,term="s(fai)",main="fai",eps=0.001,output="KCAL_plot", confidence = 95)

data$total_kcal_npe_low_fermentation=data$total_kcal_npe_low_fermentation+.001

NPE<-brm(total_kcal_npe_low_fermentation~s(fai) + class_focal  + (1 +  fai + class_focal  | name_focal), 
          data = data[complete.cases(data[,c("total_kcal_npe_low_fermentation","fai")]),],
          family=Gamma(link="log"),iter = 4000, inits=0,chains = 2, 
          prior=c(prior(normal(0,10),class="Intercept"),
                  prior(normal(0,10),class="b"),
                  prior(gamma(.01,.01),class="shape"),
                  prior(normal(0,10),class="sd")),
          control = list(max_treedepth = 11, adapt_delta = .999))

summary(NPE, digits=9)
plot(conditional_effects(NPE,spaghetti=FALSE)) ##plots the regression lines/conditional effects
pp_check(NPE)
loo_R2(NPE)
bayes_R2(NPE)

deriv_plot(model=NPE,term="s(fai)",main="fai",eps=0.001,output="NPEL_plot", confidence = 95)

data$ap_kcal=data$ap_kcal+.001

AP<-brm(ap_kcal~s(fai) + class_focal  + (1 +  fai + class_focal  | name_focal), 
         data = data[complete.cases(data[,c("ap_kcal","fai")]),],
         family=Gamma(link="log"),iter = 4000, inits=0,chains = 2, 
         prior=c(prior(normal(0,10),class="Intercept"),
                 prior(normal(0,10),class="b"),
                 prior(gamma(.01,.01),class="shape"),
                 prior(normal(0,10),class="sd")),
         control = list(max_treedepth = 11, adapt_delta = .999))

summary(AP, digits=9)
plot(conditional_effects(AP,spaghetti=FALSE)) ##plots the regression lines/conditional effects
pp_check(AP)
loo_R2(AP)
bayes_R2(AP)

deriv_plot(model=AP,term="s(fai)",main="fai",eps=0.001,output="AP_plot", confidence = 95)






data$lipid_kcal=data$lipid_kcal+.001

Lipid<-brm(lipid_kcal~s(fai) + class_focal  + (1 +  fai + class_focal  | name_focal), 
        data = data[complete.cases(data[,c("lipid_kcal","fai")]),],
        family=Gamma(link="log"),iter = 4000, inits=0,chains = 2, 
        prior=c(prior(normal(0,10),class="Intercept"),
                prior(normal(0,10),class="b"),
                prior(gamma(.01,.01),class="shape"),
                prior(normal(0,10),class="sd")),
        control = list(max_treedepth = 11, adapt_delta = .999))

summary(Lipid, digits=9)
plot(conditional_effects(Lipid,spaghetti=FALSE)) ##plots the regression lines/conditional effects
pp_check(Lipid)
loo_R2(Lipid)
bayes_R2(Lipid)

deriv_plot(model=Lipid,term="s(fai)",main="fai",eps=0.001,output="Lipid_plot", confidence = 95)


data$tnc_kcal=data$tnc_kcal+.001

TNC<-brm(tnc_kcal~s(fai) + class_focal  + (1 +  fai + class_focal  | name_focal), 
           data = data[complete.cases(data[,c("tnc_kcal","fai")]),],
           family=Gamma(link="log"),iter = 4000, inits=0,chains = 2, 
           prior=c(prior(normal(0,10),class="Intercept"),
                   prior(normal(0,10),class="b"),
                   prior(gamma(.01,.01),class="shape"),
                   prior(normal(0,10),class="sd")),
           control = list(max_treedepth = 11, adapt_delta = .999))

summary(TNC, digits=9)
plot(conditional_effects(TNC,spaghetti=FALSE)) ##plots the regression lines/conditional effects
pp_check(TNC)
loo_R2(TNC)
bayes_R2(TNC)

deriv_plot(model=TNC,term="s(fai)",main="fai",eps=0.001,output="TNC_plot", confidence = 95)


data$ndf_kcal_low=data$ndf_kcal_low+.001

NDF<-brm(ndf_kcal_low~s(fai) + class_focal  + (1 +  fai + class_focal  | name_focal), 
         data = data[complete.cases(data[,c("ndf_kcal_low","fai")]),],
         family=Gamma(link="log"),iter = 4000, inits=0,chains = 2, 
         prior=c(prior(normal(0,10),class="Intercept"),
                 prior(normal(0,10),class="b"),
                 prior(gamma(.01,.01),class="shape"),
                 prior(normal(0,10),class="sd")),
         control = list(max_treedepth = 11, adapt_delta = .999))

summary(NDF, digits=9)
plot(conditional_effects(NDF,spaghetti=FALSE)) ##plots the regression lines/conditional effects
pp_check(NDF)
loo_R2(NDF)
bayes_R2(NDF)

deriv_plot(model=NDF,term="s(fai)",main="fai",eps=0.001,output="NDF_plot", confidence = 95)



data$npe_to_ap=data$npe_to_ap+.001

NPeAp<-brm(npe_to_ap~s(fai) + class_focal  + (1 +  fai + class_focal  | name_focal), 
         data = data[complete.cases(data[,c("npe_to_ap","fai")]),],
         family=Gamma(link="log"),iter = 4000, inits=0,chains = 2, 
         prior=c(prior(normal(0,10),class="Intercept"),
                 prior(normal(0,10),class="b"),
                 prior(gamma(.01,.01),class="shape"),
                 prior(normal(0,10),class="sd")),
         control = list(max_treedepth = 11, adapt_delta = .999))

summary(NPeAp, digits=9)
plot(conditional_effects(NPeAp,spaghetti=FALSE)) ##plots the regression lines/conditional effects
pp_check(NPeAp)
loo_R2(NPeAp)
bayes_R2(NPeAp)

deriv_plot(model=NPeAp,term="s(fai)",main="fai",eps=0.001,output="NPeAp_plot", confidence = 95)


data$lipid_to_ap=data$lipid_to_ap+.001

LipidAp<-brm(lipid_to_ap~s(fai) + class_focal  + (1 +  fai + class_focal  | name_focal), 
           data = data[complete.cases(data[,c("lipid_to_ap","fai")]),],
           family=Gamma(link="log"),iter = 4000, inits=0,chains = 2, 
           prior=c(prior(normal(0,10),class="Intercept"),
                   prior(normal(0,10),class="b"),
                   prior(gamma(.01,.01),class="shape"),
                   prior(normal(0,10),class="sd")),
           control = list(max_treedepth = 11, adapt_delta = .999))

summary(LipidAp, digits=9)
plot(conditional_effects(LipidAp,spaghetti=FALSE)) ##plots the regression lines/conditional effects
pp_check(LipidAp)
loo_R2(LipidAp)
bayes_R2(LipidAp)

deriv_plot(model=LipidAp,term="s(fai)",main="fai",eps=0.001,output="LipidAp_plot", confidence = 95)



data$lipid_to_tnc=data$lipid_to_tnc+.001

LipidTNC<-brm(lipid_to_tnc~s(fai) + class_focal  + (1 +  fai + class_focal  | name_focal), 
             data = data[complete.cases(data[,c("lipid_to_tnc","fai")]),],
             family=Gamma(link="log"),iter = 4000, inits=0,chains = 2, 
             prior=c(prior(normal(0,10),class="Intercept"),
                     prior(normal(0,10),class="b"),
                     prior(gamma(.01,.01),class="shape"),
                     prior(normal(0,10),class="sd")),
             control = list(max_treedepth = 11, adapt_delta = .999))

summary(LipidTNC, digits=9)
plot(conditional_effects(LipidTNC,spaghetti=FALSE)) ##plots the regression lines/conditional effects
pp_check(LipidTNC)
loo_R2(LipidTNC)
bayes_R2(LipidTNC)

deriv_plot(model=LipidTNC,term="s(fai)",main="fai",eps=0.001,output="LipidTNC_plot", confidence = 95)






data$lipid_to_ndf_low=data$lipid_to_ndf_low+.001

LipidNDF<-brm(lipid_to_ndf_low~s(fai) + class_focal  + (1 +  fai + class_focal  | name_focal), 
              data = data[complete.cases(data[,c("lipid_to_ndf_low","fai")]),],
              family=Gamma(link="log"),iter = 4000, inits=0,chains = 2, 
              prior=c(prior(normal(0,10),class="Intercept"),
                      prior(normal(0,10),class="b"),
                      prior(gamma(.01,.01),class="shape"),
                      prior(normal(0,10),class="sd")),
              control = list(max_treedepth = 11, adapt_delta = .999))

summary(LipidNDF, digits=9)
plot(conditional_effects(LipidNDF,spaghetti=FALSE)) ##plots the regression lines/conditional effects
pp_check(LipidNDF)
loo_R2(LipidNDF)
bayes_R2(LipidNDF)

deriv_plot(model=LipidNDF,term="s(fai)",main="fai",eps=0.001,output="LipidNDF_plot", confidence = 95)








data$tnc_to_ap=data$tnc_to_ap+.001

TNC_AP<-brm(tnc_to_ap~s(fai) + class_focal  + (1 +  fai + class_focal  | name_focal), 
              data = data[complete.cases(data[,c("tnc_to_ap","fai")]),],
              family=Gamma(link="log"),iter = 4000, inits=0,chains = 2, 
              prior=c(prior(normal(0,10),class="Intercept"),
                      prior(normal(0,10),class="b"),
                      prior(gamma(.01,.01),class="shape"),
                      prior(normal(0,10),class="sd")),
              control = list(max_treedepth = 11, adapt_delta = .999))

summary(TNC_AP, digits=9)
plot(conditional_effects(TNC_AP,spaghetti=FALSE)) ##plots the regression lines/conditional effects
pp_check(TNC_AP)
loo_R2(TNC_AP)
bayes_R2(TNC_AP)

deriv_plot(model=TNC_AP,term="s(fai)",main="fai",eps=0.001,output="TNC_AP_plot", confidence = 95)





data$tnc_to_ndf_low=data$tnc_to_ndf_low+.001

TNC_NDF<-brm(tnc_to_ndf_low~s(fai) + class_focal  + (1 +  fai + class_focal  | name_focal), 
            data = data[complete.cases(data[,c("tnc_to_ndf_low","fai")]),],
            family=Gamma(link="log"),iter = 4000, inits=0,chains = 2, 
            prior=c(prior(normal(0,10),class="Intercept"),
                    prior(normal(0,10),class="b"),
                    prior(gamma(.01,.01),class="shape"),
                    prior(normal(0,10),class="sd")),
            control = list(max_treedepth = 11, adapt_delta = .999))

summary(TNC_NDF, digits=9)
plot(conditional_effects(TNC_NDF,spaghetti=FALSE)) ##plots the regression lines/conditional effects
pp_check(TNC_NDF)
loo_R2(TNC_NDF)
bayes_R2(TNC_NDF)

deriv_plot(model=TNC_NDF,term="s(fai)",main="fai",eps=0.001,output="TNC_NDF_plot", confidence = 95)






data$ap_to_ndf_low=data$ap_to_ndf_low+.001

AP_NDF<-brm(ap_to_ndf_low~s(fai) + class_focal  + (1 +  fai + class_focal  | name_focal), 
             data = data[complete.cases(data[,c("ap_to_ndf_low","fai")]),],
             family=Gamma(link="log"),iter = 4000, inits=0,chains = 2, 
             prior=c(prior(normal(0,10),class="Intercept"),
                     prior(normal(0,10),class="b"),
                     prior(gamma(.01,.01),class="shape"),
                     prior(normal(0,10),class="sd")),
             control = list(max_treedepth = 11, adapt_delta = .999))

summary(AP_NDF, digits=9)
plot(conditional_effects(AP_NDF,spaghetti=FALSE)) ##plots the regression lines/conditional effects
pp_check(AP_NDF)
loo_R2(AP_NDF)
bayes_R2(AP_NDF)

deriv_plot(model=AP_NDF,term="s(fai)",main="fai",eps=0.001,output="AP_NDF_plot", confidence = 95)


KCAL_plot=KCAL_plot+ylab("Total kcal")+xlab("FAI") + geom_vline(xintercept = 4.5, color = "gray") + geom_vline(xintercept = 4, color = "gray", linetype = "dashed")
NPEL_plot=NPEL_plot+ylab("Npe")+xlab("FAI")+ geom_vline(xintercept = 4.5, color = "gray") + geom_vline(xintercept = 4, color = "gray", linetype = "dashed") + geom_vline(xintercept = 4.5, color = "gray") + geom_vline(xintercept = 4, color = "gray", linetype = "dashed")
AP_plot=AP_plot+ylab("Protein")+xlab("FAI") + geom_vline(xintercept = 4.5, color = "gray") + geom_vline(xintercept = 4, color = "gray", linetype = "dashed")
Lipid_plot=Lipid_plot+ylab("Lipid")+xlab("FAI") + geom_vline(xintercept = 4.5, color = "gray") + geom_vline(xintercept = 4, color = "gray", linetype = "dashed")
TNC_plot=TNC_plot+ylab("TNC")+xlab("FAI") + geom_vline(xintercept = 4.5, color = "gray") + geom_vline(xintercept = 4, color = "gray", linetype = "dashed")
NDF_plot=NDF_plot+ylab("NDF")+xlab("FAI") + geom_vline(xintercept = 4.5, color = "gray") + geom_vline(xintercept = 4, color = "gray", linetype = "dashed")
NPeAp_plot=NPeAp_plot+ylab("Npe:P")+xlab("FAI") + geom_vline(xintercept = 4.5, color = "gray") + geom_vline(xintercept = 4, color = "gray", linetype = "dashed")
LipidAp_plot=LipidAp_plot+ylab("Lipid:P")+xlab("FAI") + geom_vline(xintercept = 4.5, color = "gray") + geom_vline(xintercept = 4, color = "gray", linetype = "dashed")
LipidTNC_plot=LipidTNC_plot+ylab("Lipid:TNC")+xlab("FAI") + geom_vline(xintercept = 4.5, color = "gray") + geom_vline(xintercept = 4, color = "gray", linetype = "dashed")
LipidNDF_plot=LipidNDF_plot+ylab("Lipid:NDF")+xlab("FAI") + geom_vline(xintercept = 4.5, color = "gray") + geom_vline(xintercept = 4, color = "gray", linetype = "dashed")
TNC_AP_plot=TNC_AP_plot+ylab("TNC:P")+xlab("FAI") + geom_vline(xintercept = 4.5, color = "gray") + geom_vline(xintercept = 4, color = "gray", linetype = "dashed")
TNC_NDF_plot=TNC_NDF_plot+ylab("TNC:NDF")+xlab("FAI") + geom_vline(xintercept = 4.5, color = "gray") + geom_vline(xintercept = 4, color = "gray", linetype = "dashed")
AP_NDF_plot=AP_NDF_plot+ylab("P:NDF")+xlab("FAI") + geom_vline(xintercept = 4.5, color = "gray") + geom_vline(xintercept = 4, color = "gray", linetype = "dashed")

ggarrange(KCAL_plot,NPEL_plot,AP_plot,Lipid_plot,TNC_plot,NDF_plot,NPeAp_plot,
          LipidAp_plot,LipidTNC_plot,LipidNDF_plot,TNC_AP_plot,TNC_NDF_plot,AP_NDF_plot)
ggsave("Intake and FAI 95CI 2023.pdf",dpi=300, device = "pdf")


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

all_models=list(KCAL,NPE,AP,Lipid,TNC,NDF,NPeAp,
          LipidAp,LipidTNC,LipidNDF,TNC_AP,TNC_NDF,AP_NDF)

mytables=lapply(all_models,sumtable)

mytables=do.call(rbind,mytables)
mytables$Model_ID=gsub("fai","FAI" ,mytables$Model_ID)
mytables$Model_ID=gsub("totalkcalusingaplowfermentation","Total kcal" ,mytables$Model_ID)
mytables$Model_ID=gsub("apkcal","Protein" ,mytables$Model_ID)
mytables$Model_ID=gsub("lipidkcal","Lipid" ,mytables$Model_ID)
mytables$Model_ID=gsub("tnckcal","TNC" ,mytables$Model_ID)
mytables$Model_ID=gsub("ndfkcallow","NDF" ,mytables$Model_ID)
mytables$Model_ID=gsub("totalkcalnpelowfermentation","Npe" ,mytables$Model_ID)
mytables$Model_ID=gsub("npetoap","Npe:P" ,mytables$Model_ID)
mytables$Model_ID=gsub("lipidtoap","Lipid:P" ,mytables$Model_ID)
mytables$Model_ID=gsub("lipidtotnc","Lipid:TNC" ,mytables$Model_ID)
mytables$Model_ID=gsub("lipidtondflow","Lipid:NDF" ,mytables$Model_ID)
mytables$Model_ID=gsub("tnctoap","TNC:P" ,mytables$Model_ID)
mytables$Model_ID=gsub("tnctondflow","TNC:NDF" ,mytables$Model_ID)
mytables$Model_ID=gsub("aptondflow","P:NDF" ,mytables$Model_ID)
mytables$Model_ID=gsub("s(", "" ,mytables$Model_ID, fixed = TRUE)
mytables$Model_ID=gsub("(","" ,mytables$Model_ID, fixed = TRUE)
mytables$Model_ID=gsub(")","" ,mytables$Model_ID, fixed = TRUE)

write.csv(mytables,file="PNAS fai intake tables.csv")




