library(dplyr)
library(data.table)
library(kableExtra)
library(survival)
library(ggfortify)
library(ggpubr)
library(survminer)
library(GGally)
library(gridExtra)
library(forestplot)
library(mice)
library(ggpattern)

#for each method, calculate Kaplan-Meier curve, Shoenfeld residuals, Cox model, RMSTD

#read in data
#aca data
aca <- read.csv("C:\\Users\\smala\\OneDrive\\Documents\\TRUST\\imputeTBculture\\data\\aca.csv")

#locf data
locf <- read.csv("C:\\Users\\smala\\OneDrive\\Documents\\TRUST\\imputeTBculture\\data\\locf.csv")

#fcs data
fcs <- read.csv("C:\\Users\\smala\\OneDrive\\Documents\\TRUST\\imputeTBculture\\data\\fcs.csv")


#read in covariate data and subset to relevant columns
covars <- read.csv("C:\\Users\\smala\\OneDrive\\Documents\\TRUST\\imputeTBculture\\data\\TCC_imputation_data.csv")
covars <- covars %>% dplyr::select(pid, screen_years, screen_sex, bl_hiv, cxr_cavity_chest_radiograph_1, smear_pos_TRUST_or_TB)
covars$smear_pos_TRUST_or_TB <- factor(covars$smear_pos_TRUST_or_TB, levels = c(1,0), labels = c("Smear negative", "Smear positive"))
#summary stats and plots for three methods
#freq converted for ACA
table(aca$status)

#freq converted for LOCF
table(locf$status)

#freq converted for FCS
table(fcs$status)

covars$cxr <- covars$cxr_cavity_chest_radiograph_1
covars$cxr[covars$cxr_cavity_chest_radiograph_1 == "Unknown"] <- NA
covars$age <- ifelse(covars$screen_years < 20, 1, ifelse(covars$screen_years < 30, 2, ifelse(covars$screen_years < 40, 3, ifelse(covars$screen_years < 50, 4, ifelse(covars$screen_years < 60, 5,6)))))
covars$age <- factor(covars$age, levels = c(1,2,3,4,5,6), labels = c("<20", "20-29", "30-39", "40-49", "50-59", "60+"))

####################################
#Available Case Analysis############
####################################

#merge ACA data with covariates
aca_covars <- dplyr::left_join(aca, covars, by= "pid")

#Cox models
#unadjusted (smear only)
aca_model <- coxph(Surv(time, status) ~ smear_pos_TRUST_or_TB, data = aca_covars)

#adjusted(age, sex, HIV, cavitaion)
aca_model_adj <- coxph(Surv(time, status) ~ smear_pos_TRUST_or_TB + screen_sex + age  + bl_hiv + cxr, data = aca_covars)

#summarize both models and format into table
aca_model_summary <- summary(aca_model)
aca_model_final <- data.frame(aca_model_summary$coefficients, aca_model_summary$conf.int) %>% dplyr::mutate(var = "Smear negative",  LB = round(lower..95, digits = 2), UB = round(upper..95, digits = 2),  pval = round(Pr...z.., digits=3), HR = round(exp.coef., digits = 2) ) %>% dplyr::select(var, HR, LB, UB, pval)

rownames(aca_model_final) <- ""

aca_model_adj_summary <- summary(aca_model_adj)
aca_model_adj_final <- data.frame(aca_model_adj_summary$coefficients, aca_model_adj_summary$conf.int) %>% dplyr::mutate(var = c("Smear negative", "Male", "Age, 20-29", "30-39", "40-49", "50-59", "60+", "HIV Positive", "Cavitation"), LB = round(lower..95, digits = 2), UB = round(upper..95, digits = 2),  pval = round(Pr...z.., digits=3), HR = round(exp.coef., digits = 2) ) %>% dplyr::select(var, HR, LB, UB, pval)
rownames(aca_model_adj_final) <- NULL
print(aca_model_final)
print(aca_model_adj_final)

#Kaplan-Meier curve
surv_aca <- Surv(time = aca_covars$time, event = aca_covars$status)
surv_aca_smear2 <- survfit(surv_aca ~ smear_pos_TRUST_or_TB , data = aca_covars,type='kaplan-meier')

surv_plot_aca2<-ggsurvplot(surv_aca_smear2, fun = "event", conf.int = T, linetype = c("solid", "dashed"), palette =c("black", "black"), censor = F, ylim = c(0,1), legend = "none", title = "A)", xlab= "Treatment Week", ylab = "Cumulative event probability") 
print(surv_plot_aca2)


#Shoenfeld residuals
testph_aca <- cox.zph(aca_model_adj)
smear <- ggcoxzph(testph_aca, var = "smear_pos_TRUST_or_TB") + ggtitle("") + ylab("Beta(t) smear status")
age <- ggcoxzph(testph_aca, var = "age") + ggtitle("") + ylab("Beta(t) age")
sex <- ggcoxzph(testph_aca, var = "screen_sex") + ggtitle("") + ylab("Beta(t) sex")
hiv <- ggcoxzph(testph_aca, var = "bl_hiv") + ggtitle("") + ylab("Beta(t) HIV status")
cxr <- ggcoxzph(testph_aca, var = "cxr") + ggtitle("") + ylab("Beta(t) cavitation")
png("aca_resids.png")
ggarrange(smear$`1`,sex$`2`,age$`3`, hiv$`4`, cxr$`5`,  nrow=3, ncol=2)
dev.off()


#RMSTD
x <- summary(surv_aca_smear2)
surv_aca_upper <- data.frame(strata=x$strata, time=x$time, surv1=x$surv) %>% dplyr::filter(strata == "smear_pos_TRUST_or_TB=Smear negative") %>% dplyr::select( surv1, time)
surv_aca_lower <- data.frame(strata=x$strata, time=x$time, surv2=x$surv) %>% dplyr::filter(strata == "smear_pos_TRUST_or_TB=Smear positive") %>% dplyr::select(surv2, time)
surv_aca_area <- dplyr::left_join(surv_aca_lower, surv_aca_upper, id = "time")
surv_aca_area$surv1[surv_aca_area$time==5] <- 0.6111111


surv_aca_summary<-summary(surv_aca_smear2)
aca_rmstd_summary <- data.frame(time=surv_aca_summary$time, strata=surv_aca_summary$strata, surv=surv_aca_summary$surv)

aca_smearpos <- aca_rmstd_summary %>% dplyr::filter(strata=="smear_pos_TRUST_or_TB=Smear positive") %>% dplyr::mutate(survcomp = 1-surv)
aca_smearneg <- aca_rmstd_summary %>% dplyr::filter(strata=="smear_pos_TRUST_or_TB=Smear negative") %>% dplyr::mutate(survcomp = 1-surv)

aca_rmstd <- dplyr::full_join(aca_smearpos, aca_smearneg, by="time")

aca_rmstd$surv.x[aca_rmstd$time==5] <- 0.6111111
aca_rmstd$strata.x[aca_rmstd$time==5] <-"smear_pos_TRUST_or_TB=Smear positive"
aca_rmstd$survcomp.x[aca_rmstd$time==5] <-1- 0.6111111
fix <- data.frame(time = c(0,1), strata.x = c("smear_pos_TRUST_or_TB=Smear positive", "smear_pos_TRUST_or_TB=Smear positive"), surv.x = c(1,1), survcomp.x = c(0,0), strata.y = c("smear_pos_TRUST_or_TB=Smear negative", "smear_pos_TRUST_or_TB=Smear negative"), surv.y = c(1,1) ,survcomp.y = c(0,0))
aca_rmstd <- rbind(fix, aca_rmstd)
aca_rmstd_fig <- ggplot(data= aca_rmstd)  + geom_rect(aes(xmin = time, xmax = time+1, ymin = survcomp.y, ymax = survcomp.x), fill = "#cccccc") + geom_rect(xmin=7, xmax=12, ymin=0, ymax=1, fill="white")+xlim(0,12) + theme_pubr()  +  geom_vline(xintercept=7, size=1, linetype="solid") + xlab("Treatment Week") + ylab("") +  geom_step(data=,aes(x=time, y=survcomp.x,linetype="solid"), size=1)+ 
  geom_step(aes(x=time, y=survcomp.y, linetype="dashed"), size=1) + ggtitle("A)")  + ylim(0,1) + ylab("Culture conversion(%)") + geom_text(x=5.25, y=.25, aes(label="\n weeks", family = "serif"))+theme(text=element_text(size=14, family="serif")) + scale_linetype_manual(name = "",values = c("dashed"=1, "solid" = 2),labels = c("Smear Negative", "Smear Positive")) + theme(legend.position = "none") +  scale_x_continuous(breaks = c(0, 4, 7,8, 12)) 
                                                                                                                                                                                                                                                                                                          
print(aca_rmstd_fig)                                                                                                                                                                                                                                                                                                          
                                                                                                                                                                                                                                                                                                          
                                                         
                                                                                                                                                                                                                                                                                                           
########################################
#Last Observation Carried Forward ######
########################################

#merge LOCF data with covariates
locf_covars <- dplyr::left_join(locf, covars, by= "pid")

#Cox models
#unadjusted (smear only)
locf_model <- coxph(Surv(time, status) ~ smear_pos_TRUST_or_TB, data = locf_covars)

#adjusted(age, sex, HIV, cavitaion)
locf_model_adj <- coxph(Surv(time, status) ~ smear_pos_TRUST_or_TB + screen_sex + age  + bl_hiv + cxr, data = locf_covars)

#summarize both models and format into table
locf_model_summary <- summary(locf_model)
locf_model_final <- data.frame(locf_model_summary$coefficients, locf_model_summary$conf.int) %>% dplyr::mutate(var = "Smear negative",  LB = round(lower..95, digits = 2), UB = round(upper..95, digits = 2),  pval = round(Pr...z.., digits=3), HR = round(exp.coef., digits = 2) ) %>% dplyr::select(var, HR, LB, UB, pval)

rownames(locf_model_final) <- ""

locf_model_adj_summary <- summary(locf_model_adj)
locf_model_adj_final <- data.frame(locf_model_adj_summary$coefficients, locf_model_adj_summary$conf.int) %>% dplyr::mutate(var = c("Smear negative", "Male", "Age, 20-29", "30-39", "40-49", "50-59", "60+", "HIV Positive", "Cavitation"), LB = round(lower..95, digits = 2), UB = round(upper..95, digits = 2),  pval = round(Pr...z.., digits=3), HR = round(exp.coef., digits = 2) ) %>% dplyr::select(var, HR, LB, UB, pval)
rownames(locf_model_adj_final) <- NULL
print(locf_model_final)
print(locf_model_adj_final)

#Kaplan-Meier curve

surv_locf <- Surv(time = locf_covars$time, event = locf_covars$status)
surv_locf_smear2 <- survfit(surv_locf ~ smear_pos_TRUST_or_TB , data = locf_covars,type='kaplan-meier')


surv_plot_locf2 <-ggsurvplot(surv_locf_smear2, fun = "event", conf.int = T, linetype = c("solid", "dashed"), palette =c("black", "black"), censor = F, ylim = c(0,1), legend = "none", title = "B)", xlab= "Treatment Week", ylab = "Cumulative event probability") 
print(surv_plot_locf2 )


#Shoenfeld residuals
testph_locf <- cox.zph(locf_model_adj)
smear <- ggcoxzph(testph_locf, var = "smear_pos_TRUST_or_TB") + ggtitle("") + ylab("Beta(t) smear status")
age <- ggcoxzph(testph_locf, var = "age") + ggtitle("") + ylab("Beta(t) age")
sex <- ggcoxzph(testph_locf, var = "screen_sex") + ggtitle("") + ylab("Beta(t) sex")
hiv <- ggcoxzph(testph_locf, var = "bl_hiv") + ggtitle("") + ylab("Beta(t) HIV status")
cxr <- ggcoxzph(testph_locf, var = "cxr") + ggtitle("") + ylab("Beta(t) cavitation")
png("locf_resids.png")
ggarrange(smear$`1`,sex$`2`,age$`3`, hiv$`4`, cxr$`5`,  nrow=3, ncol=2)
dev.off()


#RMSTD
x <- summary(surv_locf_smear2)
surv_locf_upper <- data.frame(strata=x$strata, time=x$time, surv1=x$surv) %>% dplyr::filter(strata == "smear_pos_TRUST_or_TB=Smear negative") %>% dplyr::select( surv1, time)
surv_locf_lower <- data.frame(strata=x$strata, time=x$time, surv2=x$surv) %>% dplyr::filter(strata == "smear_pos_TRUST_or_TB=Smear positive") %>% dplyr::select(surv2, time)
surv_locf_area <- dplyr::left_join(surv_locf_lower, surv_locf_upper, id = "time")


surv_locf_summary<-summary(surv_locf_smear2)
locf_rmstd_summary <- data.frame(time=surv_locf_summary$time, strata=surv_locf_summary$strata, surv=surv_locf_summary$surv)

locf_smearpos <- locf_rmstd_summary %>% dplyr::filter(strata=="smear_pos_TRUST_or_TB=Smear positive") %>% dplyr::mutate(survcomp = 1-surv)
locf_smearneg <- locf_rmstd_summary %>% dplyr::filter(strata=="smear_pos_TRUST_or_TB=Smear negative") %>% dplyr::mutate(survcomp = 1-surv)



locf_rmstd <- dplyr::full_join(locf_smearpos, locf_smearneg, by="time")

locf_rmstd$surv.x[locf_rmstd$time==9] <- 0.02732794
locf_rmstd$strata.x[locf_rmstd$time==9] <-"smear_pos_TRUST_or_TB=Smear negative"
locf_rmstd$survcomp.x[locf_rmstd$time==9] <-1- 0.02732794
locf_rmstd$surv.y[locf_rmstd$time==13] <- 0.2406699
locf_rmstd$strata.y[locf_rmstd$time==13] <-"smear_pos_TRUST_or_TB=Smear negative"
locf_rmstd$survcomp.y[locf_rmstd$time==13] <-1- 0.2406699
fix <- data.frame(time = c(0,1), strata.x = c("smear_pos_TRUST_or_TB=Smear positive", "smear_pos_TRUST_or_TB=Smear positive"), surv.x = c(1,1), survcomp.x = c(0,0), strata.y = c("smear_pos_TRUST_or_TB=Smear negative", "smear_pos_TRUST_or_TB=Smear negative"), surv.y = c(1,1) ,survcomp.y = c(0,0))
locf_rmstd <- rbind(fix, locf_rmstd)
locf_rmstd_fig <- ggplot(data= locf_rmstd)  + geom_rect(aes(xmin = time, xmax = time+1, ymin = survcomp.y, ymax = survcomp.x), fill = "#cccccc") + geom_rect(xmin=8, xmax=12, ymin=0, ymax=1, fill="white")+xlim(0,12) + theme_pubr()  +  geom_vline(xintercept=8, size=1, linetype="solid") + xlab("Treatment Week") + ylab("") +  geom_step(data=,aes(x=time, y=survcomp.x,linetype="solid"), size=1)+ 
  geom_step(aes(x=time, y=survcomp.y, linetype="dashed"), size=1) + ggtitle("B)")  + ylim(0,1) + ylab("Culture conversion(%)") + geom_text(x=5, y=.4, aes(label="\n weeks", family = "serif"))+theme(text=element_text(size=14, family="serif")) + scale_linetype_manual(name = "",values = c("dashed"=1, "solid" = 2),labels = c("Smear Negative", "Smear Positive")) + theme(legend.position = "none") +  scale_x_continuous(breaks = c(0, 4,8, 12)) 

print(locf_rmstd_fig)             

###########################
#Multiple Imputation#######
###########################
fcs_covars <- dplyr::left_join(fcs, covars, by = "pid")
imputed2 <- as.mids(fcs_covars)
#Cox models
x <- coxph(Surv(time, status) ~ smear_pos_TRUST_or_TB + screen_sex + age  + bl_hiv + cxr, data = fcs_covars[fcs_covars$.imp==1,])
y <- cox.zph(x)
  
smear <- ggcoxzph(y, var = "smear_pos_TRUST_or_TB") + ggtitle("") + ylab("Beta(t) smear status")
age <- ggcoxzph(y, var = "age") + ggtitle("") + ylab("Beta(t) age")
sex <- ggcoxzph(y, var = "screen_sex") + ggtitle("") + ylab("Beta(t) sex")
hiv <- ggcoxzph(y, var = "bl_hiv") + ggtitle("") + ylab("Beta(t) HIV status")
cxr <- ggcoxzph(y, var = "cxr") + ggtitle("") + ylab("Beta(t) cavitation")
png("mi_resids.png")
ggarrange(smear$`1`,sex$`2`,age$`3`, hiv$`4`, cxr$`5`,  nrow=3, ncol=2)
dev.off()

#unadjusted (smear only)
mi_model <- with(imputed2,coxph(Surv(time, status) ~ smear_pos_TRUST_or_TB))
mi_model_summary <- data.frame(summary(pool(mi_model)))

mi_model_summary$lb <-  mi_model_summary$estimate - 1.96*mi_model_summary$std.error
mi_model_summary$ub <-  mi_model_summary$estimate + 1.96*mi_model_summary$std.error
mi_model_summary <- mi_model_summary %>% dplyr::mutate(hr = round(exp(estimate), digits = 2), lb_exp = round(exp(lb), digits = 2), ub_exp = round(exp(ub), digits = 2), pval = round(p.value, digits=3), var = "Smear negative")
mi_model_final <- mi_model_summary %>% dplyr::select(var, hr, lb_exp, ub_exp, pval)


mi_model_adj <-  with(imputed2,coxph(Surv(time, status) ~ smear_pos_TRUST_or_TB + screen_sex + age + bl_hiv + cxr))

mi_model_adj_summary <- data.frame(summary(pool(mi_model_adj)))
mi_model_adj_summary$lb <- mi_model_adj_summary$estimate - 1.96*mi_model_adj_summary$std.error
mi_model_adj_summary$ub <- mi_model_adj_summary$estimate + 1.96*mi_model_adj_summary$std.error
mi_model_adj_summary <- mi_model_adj_summary %>% dplyr::mutate(var = c("Smear negative", "Male", "Age, 20-29", "30-39", "40-49", "50-59", "60+", "HIV Positive", "Cavitation"), hr = round(exp(estimate), digits = 2), lb_exp = round(exp(lb), digits = 2), ub_exp = round(exp(ub), digits = 2), pval = round(p.value, digits=3))

mi_model_adj_final <- mi_model_adj_summary %>% dplyr::select(var, hr, lb_exp, ub_exp, pval)



survprob_mi <- function(alldata = data.frame(), imp_number = numeric()){
  survest1 <- alldata %>% dplyr::filter(.imp == imp_number)
  surv1 <- survfit(Surv(time, status) ~ smear_pos_TRUST_or_TB , data = survest1)
  prob <- data.frame(summary(surv1, times=c(1,2,3,4, 5, 6, 7, 8, 9, 10, 11, 12))$surv)
  strata <- data.frame(summary(surv1, times=c(1,2,3,4, 5, 6, 7, 8, 9, 10, 11, 12))$strata)
  stderr <- data.frame(summary(surv1, times=c(1,2,3,4, 5, 6, 7, 8, 9, 10, 11, 12))$std.err)
  time <- data.frame(summary(surv1, times=c(1,2,3,4, 5, 6, 7, 8, 9, 10, 11, 12))$time)
  surv1_prob <- data.frame(prob, stderr, strata, time, imp = imp_number)
  colnames(surv1_prob) <- c("prob", "stderr", "strata","time", "imp")
  return(surv1_prob)
}



survest_1 <- survprob_mi(alldata = fcs_covars, imp_number=1)
survest_2 <- survprob_mi(alldata = fcs_covars, imp_number=2)
survest_3 <- survprob_mi(alldata = fcs_covars, imp_number=3)
survest_4 <- survprob_mi(alldata = fcs_covars, imp_number=4)
survest_5 <- survprob_mi(alldata = fcs_covars, imp_number=5)
survest_6 <- survprob_mi(alldata = fcs_covars, imp_number=6)
survest_7 <- survprob_mi(alldata = fcs_covars, imp_number=7)
survest_8 <- survprob_mi(alldata = fcs_covars, imp_number=8)
survest_9 <- survprob_mi(alldata = fcs_covars, imp_number=9)
survest_10 <- survprob_mi(alldata = fcs_covars, imp_number=10)
survest_11 <- survprob_mi(alldata = fcs_covars, imp_number=11)
survest_12 <- survprob_mi(alldata = fcs_covars, imp_number=12)
survest_13 <- survprob_mi(alldata = fcs_covars, imp_number=13)
survest_14 <- survprob_mi(alldata = fcs_covars, imp_number=14)
survest_15 <- survprob_mi(alldata = fcs_covars, imp_number=15)
survest_16 <- survprob_mi(alldata = fcs_covars, imp_number=16)
survest_17 <- survprob_mi(alldata = fcs_covars, imp_number=17)
survest_18 <- survprob_mi(alldata = fcs_covars, imp_number=18)
survest_19 <- survprob_mi(alldata = fcs_covars, imp_number=19)
survest_20 <- survprob_mi(alldata = fcs_covars, imp_number=20)
surv_all <- rbind(survest_1, survest_2, survest_3, survest_4, survest_5, survest_6, survest_7, survest_8, survest_9, survest_10, survest_11, survest_12, survest_13, survest_14, survest_15, survest_16, survest_17, survest_18, survest_19, survest_20)

surv_all$prob[surv_all$prob == 1 | surv_all$prob == 0] <- NA
surv_all$s_transform <- log(-log(1-surv_all$prob))
surv_all$var_transform <- ((surv_all$stderr)*(surv_all$prob^2))/((log(1-surv_all$prob)*(1-surv_all$prob))^2)
surv_all$var_transform[is.nan(surv_all$var_transform)] <- NA

surv_all2 <- surv_all %>% dplyr::group_by(time, strata) %>% summarise(s_bar = mean(s_transform, na.rm = T, ), vbar = mean(var_transform, na.rm = T), big_var = var(s_transform, na.rm=T))

surv_all2$prob_pool <- 1-exp(-exp(surv_all2$s_bar))
surv_all2$prob_cumulative <- 1- surv_all2$prob_pool

surv_all2$total_var <- surv_all2$vbar + (1+(1/20))*surv_all2$big_var
surv_all2$prob_pool <- 1-exp(-exp(surv_all2$s_bar))
surv_all2$prob_cumulative <- 1-surv_all2$prob_pool

surv_all2$lb <- surv_all2$s_bar - 1.96*sqrt(surv_all2$total_var)
surv_all2$ub <- surv_all2$s_bar + 1.96*sqrt(surv_all2$total_var)


surv_all2$lb_pool <- 1-exp(-exp(surv_all2$lb))
surv_all2$ub_pool <- 1-exp(-exp(surv_all2$ub))


surv_all2$ub_cumulative <- 1-surv_all2$lb_pool
surv_all2$lb_cumulative <- 1-surv_all2$ub_pool

surv_all2$lb[surv_all2$time==1]<-0
surv_all2$ub[surv_all2$time==1]<-0



surv_all2$prob_cumulative[surv_all2$time==10 & surv_all2$strata == 'smear_pos_TRUST_or_TB=Smear positive']<-NA
surv_all2$prob_cumulative[surv_all2$time==11 & surv_all2$strata == 'smear_pos_TRUST_or_TB=Smear positive']<-NA
surv_all2$prob_cumulative[surv_all2$time==12 & surv_all2$strata == 'smear_pos_TRUST_or_TB=Smear positive']<-NA
surv_all2$ub[surv_all2$time==7 & surv_all2$strata == 'smear_pos_TRUST_or_TB=Smear positive']<-1
surv_all2$ub[surv_all2$time==8 & surv_all2$strata == 'smear_pos_TRUST_or_TB=Smear positive']<-1
surv_all2$ub[surv_all2$time==9 & surv_all2$strata == 'smear_pos_TRUST_or_TB=Smear positive']<-1
surv_all2$ub[surv_all2$time==6 & surv_all2$strata == 'smear_pos_TRUST_or_TB=Smear positive']<-1

surv_all2 <- surv_all2 %>% dplyr::select(time, strata, prob_cumulative, lb_cumulative, ub_cumulative)
colnames(surv_all2) <- c("time", "strata", "prob_cumulative", "lb", "ub")

add <- data.frame(time = c(0,0), strata = c('smear_pos_TRUST_or_TB=Smear negative', "smear_pos_TRUST_or_TB=Smear positive"),   prob_cumulative = c(0,0),lb = c(0,0), ub = c(0,0))

surv_all2$prob_cumulative[is.nan(surv_all2$prob_cumulative)] <- 0
surv_all3 <- rbind(add, surv_all2)
rect <- surv_all3 %>% dplyr::select(time, lb, ub, strata) 
rect$timemax <- rect$time+1
rect$lb[rect$time==2 & rect$strata== "smear_pos_TRUST_or_TB=Smear positive"] <- 0
rect$ub[rect$time==9 & rect$strata== "smear_pos_TRUST_or_TB=Smear positive"] <- NA
rect$lb[rect$time==9 & rect$strata== "smear_pos_TRUST_or_TB=Smear positive"] <- NA

rect <- rect %>% filter(timemax != 13)

fcs_km <- ggplot() +   geom_rect(data = rect,
                                 aes(xmin = time, xmax = timemax, ymin = lb, ymax = ub), fill = "grey")+
  geom_step(data=surv_all3, mapping=aes(x=time, y=prob_cumulative, linetype = strata), size=1)+
  xlab("Treatment Week") + ylab("Cumulative Event Probability") + 
  ggtitle("C)") + 
  scale_x_continuous(breaks = c(0,3,6,9,12)) + theme_pubr() +ylim(0,1)  + 
  theme( legend.position = "none", axis.title  = element_text(size=14)) +
  scale_linetype_manual(values = c("solid", "longdash"))

fcs_rmstd <- ggplot() +  
  geom_step(data=surv_all3, mapping=aes(x=time, y=prob_cumulative, linetype = strata), size=1)+
  xlab("Treatment Week") + ylab("Culture Conversion(%)") + geom_vline(xintercept=8, size=1) +
  ggtitle("C)") + 
  scale_x_continuous(breaks = c(0,4,8,12)) + theme_pubr() +ylim(0,1)  + 
  theme( legend.position = "none", axis.title  = element_text(size=14)) +
  scale_linetype_manual(values = c("solid", "longdash"))  + geom_text(x=5, y=.4, aes(label="\n weeks", family = "serif"))+theme(text=element_text(size=14, family="serif"))


#plots for manuscript 
#Figure 1
aca2 <- aca
aca2 <- aca2[,2:13]
aca2$id <- rownames(aca2)
aca_long <- reshape(aca2, varying = list(sputum_specimen=(1:12)),
                           v.names=c("sputum_specimen"), 
                           direction="long",  
                           times=1:12,        
                           timevar="week")

aca_long$type <- "aca"


#set up carry forward data for table
locf2 <- locf[,2:13]
locf2$id <- rownames(locf2)
locf_long <- reshape(locf2, varying = list(sputum_specimen=(1:12)),
                   v.names=c("sputum_specimen"), 
                   # that was needed after changed 'varying' arg to a list to allow 'times' 
                   direction="long",  
                   times=1:12,        # substitutes number for T1 and T2
                   timevar="week")
locf_long$type <- "locf"


#set up fcs data for table
fcs <- fcs_covars %>% dplyr::filter(.imp != 0) %>% dplyr::select(paste("culture_conversion_sputum_specimen_", 1:12, sep = ""))
fcs_long <- reshape(fcs, varying = list(sputum_specimen=(1:12)),
                    v.names=c("sputum_specimen"), 
                    # that was needed after changed 'varying' arg to a list to allow 'times' 
                    direction="long",  
                    times=1:12,        # substitutes number for T1 and T2
                    timevar="week")
fcs_long$type <- "FCS"


alldata_long <- rbind(aca_long, locf_long, fcs_long)

summary_N <- alldata_long %>% dplyr::group_by(type, week, sputum_specimen) %>% dplyr::summarise(N = n())

alldata_wide <- dcast(data.table(summary_N), week + sputum_specimen ~ type , value.var = "N" )


alldata_wide$aca_perc <- round((alldata_wide$aca/dim(aca)[1])*100, digits = 1)

alldata_wide$locf_perc <- round((alldata_wide$locf/dim(locf)[1])*100, digits = 1)

alldata_wide$fcs_avg <- round(alldata_wide$FCS / 20, digits = 1)

alldata_wide$fcs_perc <- round((alldata_wide$FCS/dim(fcs)[1])*100, digits = 1)

alldata_wide$sputum_specimen[is.na(alldata_wide$sputum_specimen)] <- "Missing"
alldata_wide <- alldata_wide %>% complete(week, sputum_specimen)
alldata_wide$locf[is.na(alldata_wide$locf)] <- 0
alldata_wide$aca[is.na(alldata_wide$aca)] <- 0
alldata_wide$fcs_avg[is.na(alldata_wide$fcs_avg)] <- 0.0
alldata_wide$fcs_perc[is.na(alldata_wide$fcs_perc)] <- 0.0
alldata_wide$locf_perc[is.na(alldata_wide$locf_perc)] <- 0.0
alldata_wide$aca_perc[is.na(alldata_wide$aca_perc)] <- 0.0
alldata_wide$sputum_specimen <- factor(alldata_wide$sputum_specimen, levels = c("tb_positive", "tb_negative", "Censored", "Missing"), labels = c("TB Positive", "TB Negative", "Censored", "Missing"))

x1 <- alldata_wide %>% dplyr::select(week, aca_perc, sputum_specimen) %>% dplyr::mutate(perc = aca_perc, method = "aca") %>% dplyr::select(!aca_perc)

x5 <- aca  %>% 
  dplyr::group_by(time, status) %>% dplyr::summarize(count2 = n()) %>% dplyr::mutate(perc = 100*(count2/ 238) , method = "aca") %>% dplyr::filter(status == 1) %>% dplyr::select(time, perc) 

x2 <- alldata_wide %>% dplyr::select(week, locf_perc, sputum_specimen) %>% dplyr::mutate(perc = locf_perc, method = "locf") %>% dplyr::select(!locf_perc)

x3 <- alldata_wide %>% dplyr::select(week, fcs_perc, sputum_specimen) %>% dplyr::mutate(perc = fcs_perc, method = "fcs") %>% dplyr::select(!fcs_perc)

x <- rbind(x1, x2, x3)
x$method <- factor(x$method, levels = c("aca", "locf", "fcs"), labels = c("A)", "B)", "C)"))
x$sputum_specimen[is.na(x$sputum_specimen)] <- "Missing"
x$week <- factor(x$week, levels = c(1:12), labels= c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11"))
summary_plot <- ggplot(data=x, aes(x=week, y=perc, fill =sputum_specimen)) + geom_bar_pattern(aes(pattern = sputum_specimen, x=week, y=perc), stat = "identity", color = "black", pattern_density = 0.15, 
                                                                                                      pattern_fill    = 'black',
                                                                                                      pattern_colour  = 'black') + scale_pattern_manual(values = c( "stripe", "circle", "none", "none"))   + facet_wrap(~method, scales = "free")  + ylab("Partcipants(%)") + xlab("Treatment Week") + scale_fill_manual(name = "", values = c("#FFFFFF","#FFFFFF", "#FFFFFF", "#000000"))+theme_pubclean()+  theme(axis.text.x = element_text( color="black"), axis.text.y = element_text( color="black"), strip.text = element_text(size = 15, hjust = 0),title = element_text(size = 12)) + theme(text = element_text(size = 10,family = "serif"), legend.position = "none", strip.background = element_blank()) + ggtitle("")

jpeg("C:\\Users\\smala\\Documents\\TRUST\\imputeTBculture\\output\\figure1.jpeg", width = 1650, height =1200,res=300)
summary_plot
dev.off()


#Figure 2 Kaplan-Meier plots
jpeg("C:\\Users\\smala\\Documents\\TRUST\\imputeTBculture\\output\\kaplanmeier.jpeg", width = 2000, height =1200, res=300)
ggarrange(surv_plot_aca2$plot,surv_plot_locf2$plot , fcs_km, ncol=3)
dev.off()

#Figure 3 Hazard ratios from Cox Model
mi_model_final$y <- 1
mi_model_adj_final$y <- 2
fcs_hr <- rbind(mi_model_final, mi_model_adj_final) %>% dplyr::filter(var == "Smear negative") %>% dplyr::mutate(method = "FCS")
colnames(fcs_hr) <- c("Var","HR", "LB", "UB", "pval","y", "method")

locf_model_final$y <- 1
locf_model_adj_final$y <- 2
locf_hr <- rbind(locf_model_final, locf_model_adj_final) %>% dplyr::filter(var == "Smear negative") %>% dplyr::mutate(method = "LOCF")
colnames(locf_hr) <- c("Var","HR", "LB", "UB", "pval", "y", "method")

aca_model_final$y <- 1
aca_model_adj_final$y <- 2
aca_hr <- rbind(aca_model_final, aca_model_adj_final) %>% dplyr::filter(var == "Smear negative") %>% dplyr::mutate(method = "aca")
colnames(aca_hr) <- c("Var","HR", "LB", "UB", "pval","y", "method")

hrs <- rbind(aca_hr, locf_hr, fcs_hr)

hrs$method <- factor(paste0(hrs$method,hrs$y), levels = c("aca1", "aca2", "LOCF1", "LOCF2" ,"FCS1", "FCS2"))


hr_plot<-  ggplot(hrs, aes(x = HR, y = method, xmin = LB, xmax = UB))  + 
  geom_point(aes(x = HR, y = method, shape = factor(y)), fill = "black", size=3) +
  geom_errorbar(aes(y = method, xmin = LB, xmax = UB), width = 0)+
  scale_shape_manual(values = c(21, 22)) + 
  geom_vline(xintercept = 1, linetype = 3) +
  xlab("Hazard Ratio") +
  theme_pubr() +
  scale_colour_identity()   + xlim(0,15) + scale_x_continuous(breaks = c(0,1,5,10,15)) + 
  scale_y_discrete(limits = rev(hrs$method), labels = c("", "MI", "", "LOCF", "", "ACA"))  +
  theme(axis.title.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "none",text = element_text(size = 10,family = "serif"))+ geom_text(aes(label = paste0(HR, " (", LB, ", ", UB, ")" ), x = as.numeric(UB)), hjust=1, vjust=-1.1, size=3) 


ggsave("hrs.png", hr_plot, width = 5, height=6)

#Figure 4 RMSTD 





library(survRM2)
#function: calculate RMSTD(adjusted or unadjusted) and pool across imputed data sets
rmstd_mi <- function(longdata = data.frame(), imp = numeric(), tau = numeric(), alpha = 0.05, covars = NULL){
  rmstd_result <- data.frame(imp = NA, est = NA, var = NA)   
  
  if(is.null(covars)){
    for(i in 1:imp){
      
      data <- longdata %>% dplyr::filter(.imp == i )
      data$smear_pos_TRUST_or_TB <- as.numeric(data$smear_pos_TRUST_or_TB)
      data$smear_pos_TRUST_or_TB[data$smear_pos_TRUST_or_TB==2] <- 1
      rmst <- rmst2(data$time, data$status, data$smear_pos_TRUST_or_TB, tau = tau)
      rmstd_est <- rmst$RMST.arm1$result[1,] - rmst$RMST.arm0$result[1,]
      rmstd_var <- sqrt(rmst$RMST.arm1$rmst.var + rmst$RMST.arm0$rmst.var)
      imp_result <- data.frame(imp = i, est = rmstd_est, var = rmstd_var)[1,]
      rmstd_result[i,] <- imp_result
    }
    rownames(rmstd_result) <- NULL
    
    est_pooled <- mean(rmstd_result$est)
    var_pooled <- mean(rmstd_result$var)
    lb <- est_pooled - qnorm(1 - alpha/2)*sqrt(var_pooled)
    ub <- est_pooled + qnorm(1 - alpha/2)*sqrt(var_pooled) 
    pval <- pnorm(-abs(est_pooled)/sqrt(var_pooled)) * 2
    rmstd_pooled <- data.frame(est_pooled, var_pooled, lb, ub, pval)
    
    rmstd_mi_result <- list(rmstd_result, rmstd_pooled)
    
  }
  
  if(is.null(covars) == F){
    
    rmstd_result <- data.frame(imp = NA, est = NA, var = NA)   
    
    for(i in 1:imp){
      covars = y %>% dplyr::filter(.imp == i) %>% dplyr::select(-.imp)
      data <- longdata %>% dplyr::filter(.imp == i )
      rmst_adj <- rmst2(data$time, data$status, data$smear_pos_TRUST_or_TB, tau = tau, covariates = covars)
      adj_summary <- data.frame(rmst_adj$RMST.difference.adjusted)[2,1:2]
      imp_result <- data.frame(imp = i, est = adj_summary$coef, var = adj_summary$se.coef.)
      rmstd_result[i,] <- imp_result
    }
    rownames(rmstd_result) <- NULL
    
    est_pooled <- mean(rmstd_result$est)
    var_pooled <- mean(rmstd_result$var)
    lb <- est_pooled - qnorm(1 - alpha/2)*sqrt(var_pooled)
    ub <- est_pooled + qnorm(1 - alpha/2)*sqrt(var_pooled) 
    pval <- pnorm(-abs(est_pooled)/sqrt(var_pooled)) * 2
    rmstd_pooled <- data.frame(est_pooled, var_pooled, lb, ub, pval)
    
    rmstd_mi_result <- list(rmstd_result, rmstd_pooled)
  }
  
  return(rmstd_mi_result)
  
}


library(survRM2)
#CCA unadjusted (maximum tau that can be used is 7 because the all participants in smear negative arm convert or censored at week 7)
cca8_final <- data.frame(method = "CCA", type = c("Unadjusted", "Adjusted"), Estimate = NA, LB = NA, UB = NA, pval = NA)
#locf unadjusted
locf_8 <- rmst2(aca_rmstd_summary$time, aca_rmstd_summary$status, aca_rmstd_summary$smear_pos_TRUST_or_TB, tau = 7)
locf_8_summary <- data.frame(locf_8$unadjusted.result)
locf_8_summary <- locf_8_summary[1,]
rownames(locf_8_summary) <- NULL
locf8_unadj <- locf_8_summary %>% dplyr::mutate(method = "LOCF", type = "Unadjusted", Estimate = round(Est., digits = 2), LB = round(lower..95, digits = 2), UB = round(upper..95, digits=2), pval = round(p, digits=4) )%>% dplyr::select(method, type, Estimate, LB, UB, pval  ) 
plot(locf_8)
locf_12 <- rmst2(carry_forward$time, carry_forward$status, carry_forward$smear_pos_TRUST_or_TB, tau = 12)
locf_12_summary <- data.frame(locf_12$unadjusted.result)
locf_12_summary <- locf_12_summary[1,]
rownames(locf_12_summary) <- NULL
locf_12_final <- locf_12_summary %>% dplyr::mutate(RMSTD = "smear positive-smear negative", Estimate = round(Est., digits = 2), LB = round(lower..95, digits = 2), UB = round(upper..95, digits=2), pval = round(p, digits=4) )%>% dplyr::select(RMSTD,Estimate, LB, UB, pval  ) %>% kable(caption = "LOCF Unadjusted, tau = 12") %>% kable_styling()

#locf adjusted
D <- carry_forward %>% dplyr::filter(is.na(bl_hiv)==F & is.na(cxr_cavity_chest_radiograph_1) == F & cxr_cavity_chest_radiograph_1 != 2)
time=D$time
status=D$status
arm=as.double(as.character(D$smear_pos_TRUST_or_TB))
tau=NULL
x=D[,c("bl_hiv", "screen_sex", "screen_years", "cxr_cavity_chest_radiograph_1")]
x$bl_hiv2 <- ifelse(x$bl_hiv == "HIV Positive", 1, 0)
x$screen_sex2 <- ifelse(x$screen_sex == "Male", 1, 0)

y <- x[,c("bl_hiv2", "screen_sex2", "screen_years", "cxr_cavity_chest_radiograph_1")]

locf_8_adj <- rmst2(time, status, arm, tau = 8, covariates = y)
locf_8_adj_summary <- data.frame(locf_8_adj$RMST.difference.adjusted)
locf_8_adj_summary <- locf_8_adj_summary[1,]
rownames(locf_8_adj_summary) <- NULL
locf8_adj <- locf_8_adj_summary %>% dplyr::mutate(method = "LOCF", type = "Adjusted", Estimate = round(coef, digits = 2), LB = round(lower..95, digits = 2), UB = round(upper..95, digits=2), pval = round(p, digits=3) )%>% dplyr::select(method, type, Estimate, LB, UB, pval  )

locf8_final <- rbind(locf8_unadj, locf8_adj)

locf_12_adj <- rmst2(time, status, arm, tau = 12, covariates = y)
locf_12_adj_summary <- data.frame(locf_12_adj$RMST.difference.adjusted)
locf_12_adj_summary <- locf_12_adj_summary[1,]
rownames(locf_12_adj_summary) <- NULL
locf_12_adj_final <- locf_12_adj_summary %>% dplyr::mutate(RMSTD = "smear positive-smear negative", Estimate = round(coef, digits = 2), LB = round(lower..95, digits = 2), UB = round(upper..95, digits=2), pval = round(p, digits=3) )%>% dplyr::select(RMSTD,Estimate, LB, UB, pval  ) %>% kable(caption = "LOCF Adjusted, tau = 12") %>% kable_styling()

#mifcs unadjusted
data2 <- long1 %>% dplyr::filter(is.na(bl_hiv)==F & is.na(cxr_cavity_chest_radiograph_1) == F & cxr_cavity_chest_radiograph_1 != 2 )

x=data2[,c(".imp","bl_hiv", "screen_sex", "screen_years", "cxr_cavity_chest_radiograph_1")]
x$bl_hiv2 <- ifelse(x$bl_hiv == "HIV Positive", 1, 0)
x$screen_sex <- ifelse(x$screen_sex == "Male", 1, 0)

y <- x[,c(".imp", "bl_hiv2", "screen_sex", "screen_years", "cxr_cavity_chest_radiograph_1")]

mi8_unadj <-  data.frame(rmstd_mi(longdata = data2, imp = 20,  tau = 8)[2])
mi8_adj <- data.frame(rmstd_mi(longdata = data2, imp = 20,  tau = 8, covars = y)[2])

mi8_final <- rbind(mi8_unadj, mi8_adj)
mi8_final <- mi8_final %>% dplyr::mutate( method = "MIFCS",type = c("Unadjusted", "Adjusted"), Estimate = round(est_pooled, digits=2), LB = round(lb, digits=2), UB = round(ub, digits=2), pval = round(pval, digits=4)) %>% dplyr::select(method, type, Estimate, LB, UB, pval)

#RMSTD, tau = 8
rmstd <- rbind(cca8_final, locf8_final,mi8_final)
rmstd %>% kable(caption = "RMSTD(smear positive - smear negative), tau=8") %>% kable_styling()

data <- fcs_covars %>% dplyr::filter(.imp==20)
data$smear_pos_TRUST_or_TB <- as.numeric(data$smear_pos_TRUST_or_TB)
data$smear_pos_TRUST_or_TB[data$smear_pos_TRUST_or_TB==2] <- 0
rmst2(data$time, data$status, data$smear_pos, tau = 8)
2.76 + 2.83 + 2.74 + 2.62 + 2.57 + 2.82 +  3.03 + 2.89 + 2.62 + 2.69 +2.87 + 2.44 + 2.70 + 2.31 + 2.82 + 2.70 + 2.72 + 2.97 + 2.87 +2.75 
####################################################################################
# setEPS()
# postscript(paste0(directory, "FIGUREFILENAME.eps"))
tps <- c(0,4, 8, 12) # this is where the numbers at risk will be printed along x axis
fix <- data.frame(surv2 = c(1,1), time = c(0,1), surv1 = c(1,1))
# figure area parameters (just use these defaults)
dat <- rbind(fix,surv_cca_area)
dat$time <- dat$time-1
par(fig=c(0,1,0.2,1), new=TRUE)
par(mar = c(2, 4, 2, 4))
par(oma=c(0,0,0,0))

# set your xlim and ylim as appropriate then update axis commands accordingly
plot(NULL,xlim=c(0,12),ylim=c(0,1),frame=F, axes=F,ann=F, xaxs="i", yaxs="i")
axis(side=1, at=c(0,4,8,12), cex.axis=0.8) # same as above
axis(side=2, at=seq(0,1,0.2), cex.axis=0.8) # same as above
mtext("Time since randomization (weeks)",side=1,line=2, cex=1.2)
mtext("Without culture conversion (%)", side=2,line=2, cex=1.2)
#polygon(c(dat[dat$time<=12,]$time,12,12,rev(dat[dat$time<=12,]$time)),c(1-dat[dat$time<=12,]$surv1,0.3, 0.3,rev(1-dat[dat$time<=12,]$surv2)),col=c("#F7CB44FF"),  border=F)


# here you plot the shaded area (polygon) then the KM curve lines
polygon(c(dat[dat$time<=7,]$time,7, 7, rev(dat[dat$time<=7,]$time)),c(1-dat[dat$time<=7,]$surv2,0.3, 0.3,rev(1-dat[dat$time<=7,]$surv1)),col=c("#F7CB44FF"),  border=F)
lines(dat$time, 1-dat$surv1,type="s",col="#B8627dff",lwd=3,lty=1)
lines(dat$time, 1-dat$surv2,type="s",col="#39568CFF",lwd=3,lty=1)

# adds vertical line for time horizon
# if your area between curves is big enough, move the text() to be inside the curves
segments(7, 0.8, 7, 0, lty=2, col="#F7CB44FF", lwd=2)
text(10.75, 0.4, "RMSTD = 8 days", col="#F7CB44FF", cex=0.8)

# legend
legend("topright", legend=c("Smear Negative", "Smear Positive"), col=c("#B8627dff","#39568CFF"), lty=c(1,1), lwd=3, box.lty=0, cex=0.8)

# prepares numbers at risk part of figure
par(fig=c(0,1,0,0.3), new=TRUE)
par(mar = c(2, 2.5, 2, 3))
par(oma=c(0,0,0,0))
plot(NULL,xlim=c(0,12),ylim=c(0,0.1),ylab="",xlab="", axes=F,frame=F) #replace 12 by the max value you want on x axis


rmstd_plot <- ggarrange(aca_rmstd_fig, locf_rmstd_fig, fcs_rmstd, nrow=1)
ggsave(rmstd_plot, height = 4, width=6)

cca8_final <- data.frame(method = "CCA", type = c("Unadjusted", "Adjusted"), Estimate = NA, LB = NA, UB = NA, pval = NA)
#locf unadjusted
locf_8 <- rmst2(aca_covars$time, aca_covars$status, aca_covars$smear_pos_TRUST_or_TB, tau = 7)
locf_8_summary <- data.frame(locf_8$unadjusted.result)
locf_8_summary <- locf_8_summary[1,]
rownames(locf_8_summary) <- NULL
rmst

#--- sample data ---#

D=aca_covars
time=as.numeric(D$time)
status=as.numeric(D$status)
arm=as.numeric(D$smear_pos_TRUST_or_TB)
arm[arm==2] <- 0
tau=7
x=D[,c(4,6,7)]
#--- without covariates ----
a=rmst2(time = time, status = status, arm =arm , tau = tau)
print(a)

x=D[,c("bl_hiv", "screen_sex", "screen_years", "cxr_cavity_chest_radiograph_1")]
x$bl_hiv2 <- ifelse(x$bl_hiv == "HIV Positive", 1, 0)
x$screen_sex2 <- ifelse(x$screen_sex == "Male", 1, 0)
x$cxr2 <- ifelse(x$cxr_cavity_chest_radiograph_1 == "Yes", 1, ifelse(x$cxr_cavity_chest_radiograph_1=="No", 0, NA))

y <- x[,c("bl_hiv2", "screen_sex2", "screen_years", "cxr2")]
a=rmst2(time = time, status = status, arm =arm , tau = tau, covariates=y)
print(a)

D=locf_covars
time=as.numeric(D$time)
status=as.numeric(D$status)
arm=as.numeric(D$smear_pos_TRUST_or_TB)
arm[arm==2] <- 0
tau=8
x=D[,c(4,6,7)]
#--- without covariates ----
a=rmst2(time = time, status = status, arm =arm , tau = tau)
print(a)