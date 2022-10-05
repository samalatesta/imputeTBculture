library(dplyr)
library(data.table)
library(kableExtra)
library(mice)

data <- read.csv("C:\\Users\\smala\\Documents\\RA\\TRUST\\AIM 1\\culture_conversion\\data\\culture_complete.csv", stringsAsFactors = F)

sub <- data  %>% dplyr::filter( culture_conversion_sputum_specimen_1 %in% c("tb_positive", "tb_positive_contaminated")) 

#covariates needed for FCS
covars <- read.csv("C:\\Users\\smala\\Documents\\RA\\TRUST\\AIM 1\\culture_conversion\\data\\TCC_imputation_data.csv")
covars <- covars %>% dplyr::select(pid, paste0("s_concafb_sputum_specimen_", 1:12), screen_sex, screen_years, bl_hiv, bl_prevtb, fstrom1_baseline, bmi, smoked_substance_use, cxr_cavity_chest_radiograph_1,  problem_alcohol)

########################
#Available case analysis
########################

aca <- sub %>% dplyr::select(pid, paste("culture_conversion_sputum_specimen_", 1:12, sep = ""))


for (i in 2:13) {
  aca[[i]][aca[[i]] %in% c("Enrolled prior to COVID", "ESP'd before 12 wks", "Newly enrolled")] <- "Censored"
  aca[[i]][aca[[i]] %in% c("tb_positive_contaminated")] <- "tb_positive"
  aca[[i]][aca[[i]] %in% c("tb_negative_contaminated", "Lab closed", "Couldn't find participant", "Hospitalized", "Poor Quality")] <- NA
  aca[[i]][aca[[i]] %in% c("No sputum produced")] <- "tb_negative"
}



aca$missing <- 0
aca$num_miss <- 0
for(j in 1:dim(aca)[1]){
  for (i in 2:13) {
    if(is.na(aca[j,i]) == T){
      aca$missing[j] <- 1
      aca$num_miss[j]+1
    }
    
  }}


for(j in 1:dim(aca)[1]){
  for (i in 3:13) {
    
    #two samples missing in a row
    if(is.na(aca[j,i]) == T & is.na(aca[j,i-1])== T){
      for(k in i:14){
        aca[j,k-1] <- "Censored"
      }
    }
    #current sample missing and previous sample negative
    if(is.na(aca[j,i]) == T & is.na(aca[j,i-1])== F  & (aca[j,i-1] == "tb_negative")){
      for(k in i:14){
        aca[j,k] <- "Censored"
      }
    }
    #current sample negative, previous sample missing
    if(is.na(aca[j,i]) == F & is.na(aca[j,i-1])== T  & (aca[j,i] == "tb_negative")){
      for(k in i:14){
        aca[j,k-1] <- "Censored"
      }
    }
  }
}

for(j in 1:dim(aca)[1]){
  if((is.na(aca[j,2]) == T & (is.na(aca[j,3])== T) | (is.na(aca[j,2]) == T & aca[j,3] == "tb_negative") | (is.na(aca[j,2]) == T & aca[j,3] == "Censored"))){
    for(k in 2:13){
      aca[j,k] <- "Censored"
    }
  }
}

for(j in 1:dim(aca)[1]){
  if(is.na(aca[j,13]) == T & (is.na(aca[j,12])== T  |(is.na(aca[j,12])== F & aca[j,12] == "tb_negative" & is.na(aca[j,13])) | (is.na(aca[j,12])== F & aca[j,12] == "Censored"))){
    aca[j,13] <- "Censored"
  }
}






aca$status <- 0
aca$conversion_week <- apply(aca[,2:13], MARGIN = 1, FUN = function(x) {
  count = 0
  flag = FALSE
  for(i in 1:length(x)) { 
    if(is.na(x[i]) == F & x[i] == "tb_negative") {
      count = count + 1
      if(count == 2) {
        return(i-1)
      }
    }else{
      count = 0 #reset count
    }
  }
  return(-1) #No conversion
})




aca$status <- ifelse(aca$conversion_week != -1, 1,0)
aca$time <- 0

aca$time <- ifelse(aca$status == 1, aca$conversion_week, 0)



aca$weeks <- 0
for(j in 1:dim(aca)[1]){
  for (i in 2:13) {
    if(is.na(aca[j,i]) == T | (is.na(aca[j,i]) == F & aca[j,i] != "Censored" )){
      aca$weeks[j] <- aca$weeks[j] +1
    }
  }
}

aca$time <- ifelse(aca$status == 1,aca$conversion_week, aca$weeks )

aca_final <- aca %>% select(pid,paste0("culture_conversion_sputum_specimen_", 1:12), time, status )

write.csv(aca_final, "C:\\Users\\smala\\Documents\\RA\\TRUST\\AIM 1\\culture_conversion\\data/aca.csv", row.names = F)


#################################
#Last Observation Carried Forward
#################################

locf <- sub %>% dplyr::select(pid, paste("culture_conversion_sputum_specimen_", 1:12, sep = ""))

for (i in 2:13) {
  locf[[i]][locf[[i]] %in% c("Enrolled prior to COVID", "ESP'd before 12 wks", "Newly enrolled")] <- "Censored"
  locf[[i]][locf[[i]] %in% c("tb_positive_contaminated")] <- "tb_positive"
  locf[[i]][locf[[i]] %in% c("tb_negative_contaminated", "Lab closed", "Couldn't find participant", "Hospitalized", "Poor Quality")] <- NA
  locf[[i]][locf[[i]] %in% c("No sputum produced")] <- "tb_negative"
}

#implement carry forward imputation
count <- 0
for(j in 1:dim(locf)[1]){
  for(i in 2:13){
   
    if( is.na(locf[j,i])){
      locf[j,i] <- "Missing"
    }
    
    if(i > 1 & locf[j,i] == "Missing" ) {
      count = count + 1
    } else{
      count = 0
    }
    if(count == 1){
      locf[j,i] <- locf[j,i-1]
    }
    
  }
}

locf[locf =="Missing"] <- NA

for(j in 1:dim(locf)[1]){
  for (i in 3:13) {
    
    #two samples missing in a row
    if(is.na(locf[j,i]) == T & is.na(locf[j,i-1])== T){
      for(k in i:14){
        locf[j,k-1] <- "Censored"
      }
    }
    #current sample missing and previous sample negative
    if(is.na(locf[j,i]) == T & is.na(locf[j,i-1])== F  & (locf[j,i-1] == "tb_negative")){
      for(k in i:14){
        locf[j,k] <- "Censored"
      }
    }
    #current sample negatve, previous sample missing
    if(is.na(locf[j,i]) == F & is.na(locf[j,i-1])== T  & (locf[j,i] == "tb_negative")){
      for(k in i:14){
        locf[j,k-1] <- "Censored"
      }
    }
  }
}

for(j in 1:dim(locf)[1]){
  if((is.na(locf[j,2]) == T & (is.na(locf[j,3])== T) | (is.na(locf[j,2]) == T & locf[j,3] == "tb_negative") | (is.na(locf[j,2]) == T & locf[j,3] == "Censored"))){
    for(k in 2:13){
      locf[j,k] <- "Censored"
    }
  }
}

for(j in 1:dim(locf)[1]){
  if(is.na(locf[j,13]) == T & (is.na(locf[j,12])== T  |(is.na(locf[j,12])== F & locf[j,12] == "tb_negative" & is.na(locf[j,13])) | (is.na(locf[j,12])== F & locf[j,12] == "Censored"))){
    locf[j,13] <- "Censored"
  }
}


locf$status <- 0
locf$conversion_week <- apply(locf[,2:13], MARGIN = 1, FUN = function(x) {
  count = 0
  flag = FALSE
  for(i in 1:length(x)) { 
    if(is.na(x[i]) == F & x[i] == "tb_negative") {
      count = count + 1
      if(count == 2) {
        return(i-1)
      }
    }else{
      count = 0 #reset count
    }
  }
  return(-1) #No conversion
})




locf$status <- ifelse(locf$conversion_week != -1, 1,0)
locf$time <- 0


locf$time <- ifelse(locf$status == 1, locf$conversion_week, 0)

locf$weeks <- 0
for(j in 1:dim(locf)[1]){
  for (i in 2:13) {
    if(is.na(locf[j,i]) == T | (is.na(locf[j,i]) == F & locf[j,i] != "Censored" )){
      locf$weeks[j] <- locf$weeks[j] +1
    }
  }
}

locf$time <- ifelse(locf$status == 1,locf$conversion_week, locf$weeks )

locf_final <- locf %>% select(pid,paste0("culture_conversion_sputum_specimen_", 1:12), time, status )

write.csv(locf_final, "C:\\Users\\smala\\Documents\\RA\\TRUST\\AIM 1\\culture_conversion\\data/locf.csv", row.names = F)



#######################################################
#Multiple Imputation by Fully Conditional Specification
#######################################################

fcs <- sub
fcs <- dplyr::left_join(fcs, covars, by ="pid")
for (i in 2:13) {
  fcs[[i]][fcs[[i]] %in% c("Enrolled prior to COVID", "ESP'd before 12 wks", "Newly enrolled")] <- "Censored"
  fcs[[i]][fcs[[i]] %in% c("tb_positive_contaminated")] <- "tb_positive"
  fcs[[i]][fcs[[i]] %in% c("tb_negative_contaminated", "Lab closed", "Poor Quality", "Couldn't find participant", "Hospitalized")] <- NA
  fcs[[i]][fcs[[i]] %in% c("No sputum produced")] <- "tb_negative"
}

for(j in 1:dim(fcs)[1]){
  fcs$censor_week[j] = 0
  for (i in 2:13) {
    
    if(fcs$censor_week[j] == 0 & is.na(fcs[j,i]) ==F & (fcs[j,i]=="Censored")){
      fcs$censor_week[j] = i-1
    }
  }
}

for (i in 2:13) {
  fcs[[i]][fcs[[i]] %in% c("Censored")] <- NA
  
}

fcs$cxr_cavity_chest_radiograph_1[fcs$cxr_cavity_chest_radiograph_1 =="Unknown"] <- NA
fcs$age <- ifelse(fcs$screen_years < 20, 1, ifelse(fcs$screen_years < 30, 2, ifelse(fcs$screen_years < 40, 3, ifelse(fcs$screen_years < 50, 4, ifelse(fcs$screen_years < 60, 5,6)))))
fcs$age <- factor(fcs$age, levels = c(1,2,3,4,5,6), labels = c("<20", "20-29", "30-39", "40-49", "50-59", "60+"))


samples <- fcs %>% dplyr::select(paste("culture_conversion_sputum_specimen_", 1:12, sep = ""), paste("s_concafb_sputum_specimen_", 1:12, sep = ""), screen_sex, age, bl_hiv, bl_prevtb,
                                           fstrom1_baseline, bmi, smoked_substance_use,  cxr_cavity_chest_radiograph_1, problem_alcohol, censor_week, pid)



#set everything to factors
samples <- samples %>% mutate_at(c(paste0("culture_conversion_sputum_specimen_", 1:12),paste0("s_concafb_sputum_specimen_", 1:12),  "screen_sex", "bl_hiv", "bl_prevtb", "cxr_cavity_chest_radiograph_1", "problem_alcohol", "bmi", "smoked_substance_use", "fstrom1_baseline"), factor)


init <-  mice(samples,  maxit = 0)

#methods vector made in excel

meth <-init$method
meth[14:24] <- "polyreg"
meth[27] <- "logreg"
meth[32] <- "logreg"

#predictor matrix made in excel 
pred = read.csv("C:\\Users\\smala\\Documents\\RA\\TRUST\\AIM 1\\culture_conversion\\data/pred.csv", row.names = 1)
pred = as.matrix(pred)

#set seed
set.seed(123)

#run full imputation with updated methods and predictor matrix, creates 5 data sets
#predictors: age, sex, hiv, tobacco use, bmi, previous TB, smoked substance use, problem alcohol
imputed <-  mice(samples, method=meth, predictorMatrix = pred, m=20, printFlag=F, mincor=0, remove.collinear = F)

#transform imputed data object into one long data set
long1 <- complete(imputed, action='long', include = T)

long1$culture_conversion_sputum_specimen_1 <- as.character(long1$culture_conversion_sputum_specimen_1)
long1$culture_conversion_sputum_specimen_2 <- as.character(long1$culture_conversion_sputum_specimen_2)
long1$culture_conversion_sputum_specimen_3 <- as.character(long1$culture_conversion_sputum_specimen_3)
long1$culture_conversion_sputum_specimen_4 <- as.character(long1$culture_conversion_sputum_specimen_4)
long1$culture_conversion_sputum_specimen_5 <- as.character(long1$culture_conversion_sputum_specimen_5)
long1$culture_conversion_sputum_specimen_6 <- as.character(long1$culture_conversion_sputum_specimen_6)
long1$culture_conversion_sputum_specimen_7 <- as.character(long1$culture_conversion_sputum_specimen_7)
long1$culture_conversion_sputum_specimen_8 <- as.character(long1$culture_conversion_sputum_specimen_8)
long1$culture_conversion_sputum_specimen_9 <- as.character(long1$culture_conversion_sputum_specimen_9)
long1$culture_conversion_sputum_specimen_10 <- as.character(long1$culture_conversion_sputum_specimen_10)
long1$culture_conversion_sputum_specimen_11 <- as.character(long1$culture_conversion_sputum_specimen_11)
long1$culture_conversion_sputum_specimen_12 <- as.character(long1$culture_conversion_sputum_specimen_12)

#update each imputed data set

#censor weeks for participants newly enrolled, enrolled before COVID, ESP'd
for(j in 1:dim(long1)[1]){
  for(i in 3:14){
    if(long1$censor_week[j] != 0 & (long1$censor_week[j] + 2 <= i)) {
      long1[j,i] = "Censored"
    }
  }
}

#calculate time to conversion

long1$status <- 0
long1$conversion_week <- apply(long1[,3:14], MARGIN = 1, FUN = function(x) {
  count = 0
  flag = FALSE
  for(i in 1:length(x)) { 
    if(is.na(x[i]) == F & x[i] == "tb_negative") {
      count = count + 1
      if(count == 2) {
        return(i-1)
      }
    }else{
      count = 0 #reset count
    }
  }
  return(-1) #No conversion
})



long1$status <- ifelse(long1$conversion_week == -1, 0,1)

long1$time <- 0

long1$time <- ifelse(long1$status == 1, long1$conversion_week, 0)

long1$weeks <- 0
for(j in 1:dim(long1)[1]){
  for (i in 3:14) {
    if(is.na(long1[j,i]) == T | (is.na(long1[j,i]) == F & long1[j,i] != "Censored" )){
      long1$weeks[j] <- long1$weeks[j] +1
    }
  }
}

long1$time <- ifelse(long1$status == 1,long1$conversion_week, long1$weeks )

fcs_final <- long1 %>% dplyr::select(.imp, .id, paste("culture_conversion_sputum_specimen_", 1:12, sep = ""),pid, status, time )

write.csv(fcs_final, "C:\\Users\\smala\\Documents\\RA\\TRUST\\AIM 1\\culture_conversion\\data/fcs.csv", row.names = F)














