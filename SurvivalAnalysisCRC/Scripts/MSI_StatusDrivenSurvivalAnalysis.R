## SURVIVAL ANALYSIS
#Simpliest Survival Analysis possible on the sample.
#It will be used to compare with more accurate and specific future 
# survival model.

#K-M estimation of Survival Function for both MSS and MSI
km_tot = survfit(formula = survival.data$SurvObj ~ 1,
                 data = survival.data)

#Kaplan-Meier estimation of Survival Function for SurvObj, splitting it
# into MSS and MSI case

#K-M estimation with survfit function from 'survival' package
km_mss_msi = survfit(formula = survival.data$SurvObj ~ survival.data$MSI_status,
             data = survival.data)

#good format representation for K-M
summary(km_mss_msi)

#Log-Rank Test to compare each group
log_rank = survdiff(formula = survival.data$SurvObj ~ survival.data$MSI_status,
                    data = clinical.data)
#if (chisq > alpha)
#H0 not verified --> groups doesn't follow the same survival distribution
# else
#H0 verified --> each group follow the same survival distrbution


#plot K-M survival step function
#km for all the samples:
plot(km_tot, 
     mark.time = TRUE, conf.int = TRUE,
     col = "forestgreen", lwd = 1.5,
     xlab = "Time(Days)", ylab = "Survival%",
     main = "Survival Function (full dataset)")
grid(NULL, NULL, lty = 5, lwd = 1.2)
legend('bottomleft',
       lty = 2, lwd = 1.5, col = "forestgreen", 
       legend = "conf. int.")
#km for mss and msi
plot(km_mss_msi, 
     xlab='Time(Days)', ylab='Survival %', 
     col= c('red','blue'),
     main='Survival Function for different MSI_status',
     mark.time = TRUE)
grid(NULL, NULL, lty = 5, lwd = 1.2)
legend("bottomleft", 
       col = c('red','blue'), 
       lwd=c(1, 1), 
       cex=0.75,
       pch=1, 
       pt.cex = 1,
       legend = c('MSS','MSI-H'))

#compare survival functions
plot(km_mss_msi, 
     xlab='Time(Days)', ylab='Survival %', 
     col= c('red','blue'),
     main='Survival Functions',
     mark.time = TRUE)
grid(NULL, NULL, lty = 5, lwd = 1.2)
lines(km_tot, col = 'green', conf.int = FALSE, lwd = 1.5)
legend("bottomleft", 
       col = c('red','blue','green'), 
       lwd = c(1, 1, 1), 
       cex = 0.75,
       pch = 1, 
       pt.cex = 1,
       legend = c('MSS','MSI-H','Both'))

#hazard function for full dataset
hazard_tot = kphaz.fit(survival.data$time, survival.data$status)
kphaz.plot(hazard_tot,  
           main = "Hazard Function (full dataset)",
           col="orange",
           lwd=1.5) #line width
grid(NULL, NULL, lty = 5, lwd = 1.2)

#hazard function for MSI and MSS
hazard_mss = kphaz.fit(MSSsurvival.data$time, MSSsurvival.data$status)
smoothed_hazard_mss = muhaz(MSSsurvival.data$time, MSSsurvival.data$status)
hazard_msi = kphaz.fit(MSIsurvival.data$time, MSIsurvival.data$status)
smoothed_hazard_msi = muhaz(MSIsurvival.data$time, MSIsurvival.data$status)
kphaz.plot(hazard_mss,  
           main = "Hazard Function for MSS patients", 
           col = "orange",
           lwd = 1.5)
grid(NULL, NULL, lty = 5, lwd = 1.2)
lines(smoothed_hazard_mss, col = "blue", lwd = 1.2)
legend('topright',
       lwd = c(1.5, 1.5), col = c("orange", "blue"), cex = 0.75,
       legend = c("step","smoothed"))
kphaz.plot(hazard_msi,  
           main = "Hazard Function for MSI patients", 
           col = "royalblue3",
           lwd = 1.5)
grid(NULL, NULL, lty = 5, lwd = 1.2)
lines(smoothed_hazard_msi, col = "red", lwd = 1.2)
legend('top',
       lwd = c(1.5, 1.5), col = c("red", "royalblue3"), cex = 0.75,
       legend = c("step","smoothed"))

#cumulative hazard function with nelson aalen estimation:
#cumulative hazard function for full dataset
cumulative_hazard_temp = km_tot$n.event / km_tot$n.risk
cumulative_hazard_tot = cumsum(cumulative_hazard_temp)
#cumulative hazard function with his definition:
#cumulative hazard function for full dataset
cumulative_hazard_temp = summary(km_tot)
cumulative_hazard_tot2 = -log(cumulative_hazard_temp$surv)
#plots full dataset
plot(cumulative_hazard_tot ~ km_tot$time,
     main='Cumulative Hazard Function (full dataset)',
     xlab='Time (Time Points)', ylab='Cumulative Hazard',
     type = 's', #type STEP function
     col = 'red')
grid(NULL, NULL, lty = 5, lwd = 1.2)
lines(cumulative_hazard_tot2 ~ cumulative_hazard_temp$time, col = 'blue', type = 's')
legend('topleft', 
       lwd = c(1.5, 1.5), col = c('red', 'blue'), cex = 0.75,
       legend = c('Nelson-Aalen', '-log(S)'))

#Life Table:
#lifetab(tis: time points, ninit = # patients,
#        nlost = # censored data, nevent = # events --> deaths )
lifetab_tot = lifetab(tis = append(km_tot$time, NA), ninit = km_tot$n, 
                      nlost = km_tot$n.censor, nevent = km_tot$n.event)
#write output in a CSV file after checking the existence of the directory
#check directory existence or create it
ifelse(!dir.exists(file.path("CsvLifeTables")), 
       dir.create(file.path("CsvLifeTables")), FALSE)
#write output
write.csv(lifetab_tot, file = "CsvLifeTables/LifeTableFullDataset.csv")

#Cox Model using MSI_status as covariate
#Suppose that MSI_status is seen as a covariate, so for each group
# every subject has 1-0 variable (1 = MSI-H, 0 = MSS).
#Not stratyfing Cox Model we assume that every covariate has the same 
# influence of the others.
cox_mss_msi = coxph(formula = SurvObj ~ MSI_status, data = survival.data)
plot(survfit(cox_mss_msi), 
     main = "Survival Function of Cox Model\nwith MSI_status as covariate",
     xlab = "Time(Days)", ylab = "Survival %",
     col = 'blue')
grid(NULL, NULL, lty = 5, lwd = 1.2)
#baseline hazard function estimation
baseline_cumulative_hazard = basehaz(cox_mss_msi)
#now we know that: cum_hazard = integral[0,inf] hazard_function
plot(baseline_cumulative_hazard$hazard ~ baseline_cumulative_hazard$time, 
     col = "darkgreen", type = 's',
     xlab = "Time (Days)", ylab = "Cumulative Hazard",
     main = "Baseline Cumulative Hazard Function for\nCox Model(full dataset) with MSI_status as covariate")
grid(NULL, NULL, lty = 5, lwd = 1.2)

#Stratified Cox Model

#Stratifying by MSI_status
strata_cox_mss_msi = coxph(formula = SurvObj ~ strata(MSI_status),
                           data = survival.data)
#OBS : Null Model cause it creates 2 strata (MSI-H/MSS) but it doesn't have
# any covariate. --> no effect parameters
#plot
plot(survfit(strata_cox_mss_msi), 
     main = "Survival Function of Cox Model\nstratified by MSI_status and without covariates",
     xlab = "Time(Days)", ylab = "Survival %",
     col = c('blue','red'), mark.time = TRUE)
grid(NULL, NULL, lty = 5, lwd = 1.2)
legend('bottomleft', lwd = c(1.5, 1.5), col = c('red', 'blue'), cex = 0.75,
       legend = c('MSI-H strata', 'MSS strata'))
#baseline cumulative hazard function estimation
baseline_cumulative_hazard_strata = basehaz(strata_cox_mss_msi)
baseline_cumulative_hazard_mss = subset(baseline_cumulative_hazard_strata,
                                        strata == "MSS")
baseline_cumulative_hazard_msi = subset(baseline_cumulative_hazard_strata,
                                        strata == "MSI-H")
haz_mss = baseline_cumulative_hazard_mss$hazard
time_mss = baseline_cumulative_hazard_mss$time
haz_msi = baseline_cumulative_hazard_msi$hazard
time_msi = baseline_cumulative_hazard_msi$time
plot(haz_mss ~ time_mss, 
     col = "darkgreen", type = 's',
     xlab = "Time (Days)", ylab = "Cumulative Hazard",
     main = "Baseline Cumulative Hazard Function\nfor different strata (MSS, MSI-H)")
grid(NULL, NULL, lty = 5, lwd = 1.2)
lines(haz_msi ~ time_msi, type = 's', col = 'red')
legend('topleft', col = c('darkgreen', 'red'), lwd = c(1.5, 1.5),
       cex = 0.85, legend = c("mss", "msi"))

#Stratifying by MSI_status AND modeling by two covariates
#Now we assume that MSS and MSI-H have a different baseline hazard function, 
# that leads to a stratification. --> Make 2 strata, 1 for MSS and 1 for MSI
#use fake covariates only for example:
cov1 = sample(0:1, length(survival.data$time), replace = TRUE)
cov2 = sample(0:1, length(survival.data$time), replace = TRUE)
survival.data_temp = survival.data
survival.data_temp$cov1 = cov1
survival.data_temp$cov2 = cov2
strata_cox_mss_msi_fake = coxph(formula = SurvObj ~ cov1 + cov2 + strata(MSI_status),
                                data = survival.data)
#print info about strata_cox_mss_msi_fake
# it shows the influence parameters of the covariates
summary(strata_cox_mss_msi_fake)
#plots survival function for each strata
plot(survfit(strata_cox_mss_msi_fake), 
     main = "Survival Function of Cox Model\nstratified by MSI_status and with 2 covariates",
     xlab = "Time(Days)", ylab = "Survival %",
     col = c('blue','red'), mark.time = TRUE)
grid(NULL, NULL, lty = 5, lwd = 1.2)
legend('bottomleft', lwd = c(1.5, 1.5), col = c('red', 'blue'), cex = 0.75,
       legend = c('MSI-H strata', 'MSS strata'))
legend('bottomright', legend = c('OBS: random covariates','generated by sample'), 
       bty = "n", cex = 0.75)
#Check model validity : Cox-Snell Residuals
martingale_residuals = residuals(strata_cox_mss_msi_fake, type = 'martingale')
cox_snell_residuals = survival.data$status - martingale_residuals
res_surv = survfit(Surv(cox_snell_residuals,survival.data$status) ~ 1)
#plot cumulative hazard of residuals
plot(res_surv$time, -log(res_surv$surv), 
     type = 'S', lwd = 1.5, col = 'cornflowerblue',
     main = "Cox Snell Residuals of stratified Cox Model:\n check overall fit of the model",
     xlab = "residuals", ylab = "Cumulative Hazard of residuals")
abline(0,1, col = 'red', lwd = 1.5, lty = 2)
legend('bottomright',
       cex = 0.75, lty = c(1,2), lwd = c(1.5, 1.5), 
       col = c('cornflowerblue', 'red'),
       legend = c('cumulative hazard', 'reference line'))
legend('topleft', cex = 0.65, bty = "n",
       legend = c('Cumulative Hazard has to follow reference line', 
                  '(if it deviate too much some covariates doesn\'t ',
                  'respect PH assumption and we need Scaled',
                  'Schoenfeld Residuals to check covariates)'))
#Check proportional hazard assumption : Scaled Schoenfeld Residuals
ssr_check_ph = cox.zph(strata_cox_mss_msi_fake)
plot(ssr_check_ph, df = 2, main = "Scaled Schoenfeld Residual test",
     col = "darkblue", lwd = 1.5)
#df = 2 --> degree of freedom --> linear test
# so a covariate respect the proportional hazard assumption if is 
# scaled Schoenfeld residual plots it is a horizontal line (more or less) at 0
