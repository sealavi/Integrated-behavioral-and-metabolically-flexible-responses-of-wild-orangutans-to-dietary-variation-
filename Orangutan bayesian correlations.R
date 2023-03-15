data=read.csv("C:/Users/salavi/Documents/documents 2020/Nutritional_DATA_28Feb2019_NOutliers_FAIketones.csv")
data$Name.to.be.used=as.factor(data$Name.to.be.used)
data=na.omit(data)
library(rstan)
library(brms)

options(mc.cores = parallel::detectCores()) ##parallelize the chains

get_prior(bf(mvbind(kcal_AvailProt_DM, kcal_TNC_DM) ~ 1 +(1|Name.to.be.used))+set_rescor(rescor = TRUE),
          data = data, 
          family = student)
prot_tnc <- brm(bf(mvbind(kcal_AvailProt_DM, kcal_TNC_DM) ~ 1 +(1|Name.to.be.used))+set_rescor(rescor = TRUE),
                data = data, family = student)
      #prior = c(prior(gamma(2, .1), class = nu),
      #          prior(normal(0, 100), class = Intercept),
      #          prior(normal(0, 100), class = sigma, resp = kcal_AvailProt_DM),
      #          prior(normal(0, 100), class = sigma, resp = kcal_TNC_DM),
      #          prior(lkj(1), class = rescor)),
      
sum1=summary(prot_tnc)

prot_lipid <- brm(bf(mvbind(kcal_AvailProt_DM, kcal_Lipid_DM) ~ 1 +(1|Name.to.be.used))+set_rescor(rescor = TRUE),
                data = data, family = student)
#prior = c(prior(gamma(2, .1), class = nu),
#          prior(normal(0, 100), class = Intercept),
#          prior(normal(0, 100), class = sigma, resp = kcal_AvailProt_DM),
#          prior(normal(0, 100), class = sigma, resp = kcal_TNC_DM),
#          prior(lkj(1), class = rescor)),

sum2=summary(prot_lipid)

prot_NDF <- brm(bf(mvbind(kcal_AvailProt_DM, kcal_NDF_DM) ~ 1 +(1|Name.to.be.used))+set_rescor(rescor = TRUE),
                  data = data, family = student)
#prior = c(prior(gamma(2, .1), class = nu),
#          prior(normal(0, 100), class = Intercept),
#          prior(normal(0, 100), class = sigma, resp = kcal_AvailProt_DM),
#          prior(normal(0, 100), class = sigma, resp = kcal_TNC_DM),
#          prior(lkj(1), class = rescor)),

sum3=summary(prot_NDF)

TNC_NDF <- brm(bf(mvbind(kcal_TNC_DM, kcal_NDF_DM) ~ 1 +(1|Name.to.be.used))+set_rescor(rescor = TRUE),
                data = data, family = student)
#prior = c(prior(gamma(2, .1), class = nu),
#          prior(normal(0, 100), class = Intercept),
#          prior(normal(0, 100), class = sigma, resp = kcal_AvailProt_DM),
#          prior(normal(0, 100), class = sigma, resp = kcal_TNC_DM),
#          prior(lkj(1), class = rescor)),

sum4=summary(TNC_NDF)


NDF_Lipid <- brm(bf(mvbind(kcal_NDF_DM, kcal_Lipid_DM) ~ 1 +(1|Name.to.be.used))+set_rescor(rescor = TRUE),
                 data = data, family = student)
#prior = c(prior(gamma(2, .1), class = nu),
#          prior(normal(0, 100), class = Intercept),
#          prior(normal(0, 100), class = sigma, resp = kcal_AvailProt_DM),
#          prior(normal(0, 100), class = sigma, resp = kcal_TNC_DM),
#          prior(lkj(1), class = rescor)),

sum5=summary(NDF_Lipid)

TNC_Lipid <- brm(bf(mvbind(kcal_TNC_DM, kcal_Lipid_DM) ~ 1 +(1|Name.to.be.used))+set_rescor(rescor = TRUE),
                 data = data, family = student)
#prior = c(prior(gamma(2, .1), class = nu),
#          prior(normal(0, 100), class = Intercept),
#          prior(normal(0, 100), class = sigma, resp = kcal_AvailProt_DM),
#          prior(normal(0, 100), class = sigma, resp = kcal_TNC_DM),
#          prior(lkj(1), class = rescor)),

sum6=summary(TNC_Lipid)
plot(TNC_Lipid) ##plots the regression lines/conditional effects


bayescorr = matrix(numeric(4*4), nrow = 4, ncol = 4) # empty matrix
colnames(bayescorr)=c("Protein", "Carbohydrates", "Lipid", "Fiber")
rownames(bayescorr)=c("Protein", "Carbohydrates", "Lipid", "Fiber")

bayescorr[1,2]=sum1$rescor_pars[1]
bayescorr[2,1]=sum1$rescor_pars[1]

bayescorr[1,3]=sum2$rescor_pars[1]
bayescorr[3,1]=sum2$rescor_pars[1]

bayescorr[1,4]=sum3$rescor_pars[1]
bayescorr[4,1]=sum3$rescor_pars[1]

bayescorr[2,4]=sum4$rescor_pars[1]
bayescorr[4,2]=sum4$rescor_pars[1]

bayescorr[4,3]=sum5$rescor_pars[1]
bayescorr[3,4]=sum5$rescor_pars[1]

bayescorr[2,3]=sum6$rescor_pars[1]
bayescorr[3,2]=sum6$rescor_pars[1]

bayescorr[1,1]=1
bayescorr[2,2]=1
bayescorr[3,3]=1
bayescorr[4,4]=1


bayescorr2 = matrix(numeric(4*4), nrow = 4, ncol = 4) # empty matrix
colnames(bayescorr2)=c("Protein", "Carbohydrates", "Lipid", "Fiber")
rownames(bayescorr2)=c("Protein", "Carbohydrates", "Lipid", "Fiber")

bayescorr2[1,2]=sum1$rescor_pars[4]-sum1$rescor_pars[3]
bayescorr2[2,1]=sum1$rescor_pars[4]-sum1$rescor_pars[3]

bayescorr2[1,3]=sum2$rescor_pars[4]-sum2$rescor_pars[3]
bayescorr2[3,1]=sum2$rescor_pars[4]-sum2$rescor_pars[3]

bayescorr2[1,4]=sum3$rescor_pars[4]-sum3$rescor_pars[3]
bayescorr2[4,1]=sum3$rescor_pars[4]-sum3$rescor_pars[3]

bayescorr2[2,4]=sum4$rescor_pars[4]-sum4$rescor_pars[3]
bayescorr2[4,2]=sum4$rescor_pars[4]-sum4$rescor_pars[3]

bayescorr2[4,3]=sum5$rescor_pars[4]-sum5$rescor_pars[3]
bayescorr2[3,4]=sum5$rescor_pars[4]-sum5$rescor_pars[3]

bayescorr2[2,3]=sum6$rescor_pars[4]-sum6$rescor_pars[3]
bayescorr2[3,2]=sum6$rescor_pars[4]-sum6$rescor_pars[3]

bayescorr2[1,1]=0
bayescorr2[2,2]=0
bayescorr2[3,3]=0
bayescorr2[4,4]=0

bayescorrdf=data.frame(row=rownames(bayescorr)[row(bayescorr)], col=colnames(bayescorr)[col(bayescorr)], Corr=c(bayescorr))

bayescorr2df=data.frame(row=rownames(bayescorr2)[row(bayescorr2)], col=colnames(bayescorr2)[col(bayescorr2)], Uncertainty=c(bayescorr2))
bayescorrdf$Uncertainty=bayescorr2df$Uncertainty
bayescorrdf$row=as.factor(bayescorrdf$row)
bayescorrdf$col=as.factor(bayescorrdf$col)

library(ggplot2)
ggplot(bayescorrdf, aes(x=col, y=row, size = Corr, color=Uncertainty)) +
  geom_point()+scale_colour_viridis_c()+coord_equal()+theme_classic()+
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank())

ggplot(data, aes(x=kcal_Prot_DM, y=kcal_TNC_DM)) +
  geom_point()+theme_classic()

