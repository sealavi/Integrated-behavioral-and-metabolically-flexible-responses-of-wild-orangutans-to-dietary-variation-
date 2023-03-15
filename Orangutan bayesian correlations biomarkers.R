library(brms)
library(cmdstanr)

options(mc.cores = parallel::detectCores()) ##parallelize the chains

urea_sgcor_mod=bf(log(urea_sgcor) ~ 1 + (1 | name_focal)) + student()
ucp_sgcor_mod=bf(log(ucp_sgcor) ~ 1 + (1 | name_focal)) + student()
dn15_result_mod=bf(dn15_result ~ 1 + (1 | name_focal)) + student()


biomarkercorrelation<-brm(urea_sgcor_mod + ucp_sgcor_mod + dn15_result_mod + set_rescor(rescor = TRUE), 
                      data = data[complete.cases(data[,c("urea_sgcor", "ucp_sgcor", "dn15_result")]),],
                      iter = 2000, 
                      prior=c(prior(student_t(3,0,10000),class="Intercept", resp = dn15result),
                              prior(student_t(3,0,10000),class="Intercept", resp = logucpsgcor),
                              prior(student_t(3,0,10000),class="Intercept", resp = logureasgcor),
                              prior(student_t(3,0,10000),class="sd", resp = dn15result),
                              prior(student_t(3,0,10000),class="sd", resp = logucpsgcor),
                              prior(student_t(3,0,10000),class="sd", resp = logureasgcor)
                              ),
                      save_pars = save_pars(all = TRUE),
                      backend = "cmdstanr",
                      threads = threading(2),
                      control = list(max_treedepth = 10, adapt_delta = .99))

table_data = summary(biomarkercorrelation)
table_data=table_data$rescor_pars
table_data=data.frame(table_data)

bayescorr = matrix(numeric(3*3), nrow = 3, ncol = 3) # empty matrix
colnames(bayescorr)=c("Urea", "UCP", "DN15")
rownames(bayescorr)=c("Urea", "UCP", "DN15")

bayescorr[1,2]=as.numeric(table_data$Estimate[1])
bayescorr[2,1]=as.numeric(table_data$Estimate[1])

bayescorr[1,3]=as.numeric(table_data$Estimate[2])
bayescorr[3,1]=as.numeric(table_data$Estimate[2])

bayescorr[2,3]=as.numeric(table_data$Estimate[3])
bayescorr[3,2]=as.numeric(table_data$Estimate[3])

bayescorr[1,1]=1
bayescorr[2,2]=1
bayescorr[3,3]=1


bayescorr2 = matrix(numeric(3*3), nrow = 3, ncol = 3) # empty matrix
colnames(bayescorr2)=c("Urea", "UCP", "DN15")
rownames(bayescorr2)=c("Urea", "UCP", "DN15")

bayescorr2[1,2]=as.numeric(table_data$u.95..CI[1])-as.numeric(table_data$l.95..CI[1])
bayescorr2[2,1]=as.numeric(table_data$u.95..CI[1])-as.numeric(table_data$l.95..CI[1])

bayescorr2[1,3]=as.numeric(table_data$u.95..CI[2])-as.numeric(table_data$l.95..CI[2])
bayescorr2[3,1]=as.numeric(table_data$u.95..CI[2])-as.numeric(table_data$l.95..CI[2])

bayescorr2[2,3]=as.numeric(table_data$u.95..CI[3])-as.numeric(table_data$l.95..CI[3])
bayescorr2[3,2]=as.numeric(table_data$u.95..CI[3])-as.numeric(table_data$l.95..CI[3])

bayescorr2[1,1]=0
bayescorr2[2,2]=0
bayescorr2[3,3]=0

bayescorrdf=data.frame(row=rownames(bayescorr)[row(bayescorr)], col=colnames(bayescorr)[col(bayescorr)], Corr=c(bayescorr))

bayescorr2df=data.frame(row=rownames(bayescorr2)[row(bayescorr2)], col=colnames(bayescorr2)[col(bayescorr2)], Uncertainty=c(bayescorr2))
bayescorrdf$Uncertainty=bayescorr2df$Uncertainty
bayescorrdf$row=as.factor(bayescorrdf$row)
bayescorrdf$col=as.factor(bayescorrdf$col)

library(ggplot2)
correlogram=ggplot(bayescorrdf, aes(x=col, y=row, size = Corr, color=Uncertainty)) +
  geom_point()+scale_colour_viridis_c()+coord_equal()+theme_classic()+ggtitle("Biomarker correlogram")+
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank())
correlogram
ggsave("Biomarker correlogram.pdf",dpi=300, device = "pdf")


Ketone_biomarkers<-brm(bf(pos_neg_dn~urea_sgcor + ucp_sgcor + dn15_result + fai +(1 +  fai + urea_sgcor + ucp_sgcor + dn15_result  | name_focal),decomp = "QR"), 
                      data = data[complete.cases(data[,c("pos_neg_dn","urea_sgcor", "ucp_sgcor", "dn15_result")]),],
                      family=bernoulli(),iter = 3000, 
                      prior=c(prior(student_t(4,0,1.5),class="Intercept"),
                              prior(student_t(4,0,1.5),class="b"),
                              prior(student_t(4,0,1.5),class="sd")),
                      save_pars = save_pars(all = TRUE),
                      backend = "cmdstanr",
                      threads = threading(2),
                      control = list(max_treedepth = 10, adapt_delta = .99))

summary(Ketone_biomarkers)
conditional_effects(Ketone_biomarkers)

results=summary(Ketone_biomarkers, digits=9)
results=data.frame(results$fixed)
results$covariate=rownames(results)
results=results[-which(results$covariate=="fai"),]
results=results[-which(results$covariate=="Intercept"),]


posteriors=ggplot(results,aes(x=Estimate,y=covariate,))+geom_point()+geom_linerange(aes(xmin=l.95..CI,xmax=u.95..CI))+
  geom_vline(xintercept = 0)+theme_classic()
posteriors

ggsave("Biomarker correlogram.pdf",dpi=300, device = "pdf")

