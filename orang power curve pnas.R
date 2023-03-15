
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

library(brms)
options(mc.cores = parallel::detectCores())
library(ggplot2)
prior1 <- prior(normal(0, 99999), nlpar = "b1") +
  prior(normal(0, 99999), nlpar = "a")
fit1 <- brm(bf(total_kcal_using_ap_low_fermentation ~ b1 * percent_ap^a, 
               b1 ~ 1 + (1 | name_focal),
               a ~ 1 + (1 | name_focal), nl = TRUE),
            family=Gamma(link="identity"),iter = 4000,chains = 2,
            data = data[complete.cases(data[,c("total_kcal_using_ap_low_fermentation","percent_ap")]),],
            prior = prior1,
            control = list(max_treedepth = 11,adapt_delta = .99999999))
fitplot=plot(conditional_effects(fit1,spaghetti=TRUE), points = TRUE, plot = FALSE)[[1]] ##plots the regression lines/conditional effects
fitplot+xlab("Percent protein (kcal)")+ylab("Total kcal")+theme_classic()
summary(fit1)
summary(fit1, prob = .9)

bayes_R2(fit1)
