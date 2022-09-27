library(dplyr)
library(data.table)
library(kableExtra)

data <- read.csv("C:\\Users\\smala\\OneDrive\\Documents\\TRUST\\imputeTBculture\\data\\culture_complete.csv", stringsAsFactors = F)

culture <- data %>% dplyr::select(pid, paste0("culture_conversion_sputum_specimen_", 1:12))

long <- melt(data.table(culture), id.vars = "pid")

t1 <- long %>% dplyr::group_by(variable, value) %>% dplyr::summarise(n = n()) %>% ungroup()

t2 <- t1 %>% spread(value, n)


#missing summary stats
t1 %>% dplyr::group_by(value) %>% summarize(n = sum(n))

#set up culture data to be used for analysis
culture2 <- data %>% dplyr::select(pid, paste0("culture_conversion_sputum_specimen_", 1:12)) %>% dplyr::filter( culture_conversion_sputum_specimen_1 %in% c("tb_positive", "tb_positive_contaminated"))
culture2[culture2 == "tb_positive_contaminated"] <- "tb_positive"
culture2[culture2 == "No sputum produced"] <- "tb_negative"

culture2[culture2 =="tb_negative_contaminated" | culture2 =="Hospitalized" | culture2 =="Lab closed" |  culture2 =="ESP'd before 12 wks"  | culture2 =="Poor Quality" | culture2 =="Enrolled prior to COVID" | culture2 =="Couldn't find participant"] <- NA

long2 <- melt(data.table(culture2), id.vars = "pid")

t3 <- long2 %>% dplyr::group_by(variable, value) %>% dplyr::summarise(n = n()) %>% ungroup()

t4 <- t3 %>% tidyr::spread(value, n)

t5 <- long2 %>% dplyr::group_by(pid) %>% dplyr::summarise(n = sum(is.na(value))) %>% ungroup()

#missing data summary stats

culture2$conversion_week <- apply(culture2[,2:13], MARGIN = 1, FUN = function(x) {
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



#data frame contains 12 weeks of culture results
#divide 12 weeks into 3 time windows: 1-4, 5-8, 9-12.For each window assign each participant as "C" for complete data, "I" for intermittent (having at least one sample missing), or "M" for all samples missing. 
miss_group <-culture2

#create variable miss1 for 1-4 week window. "M" not possible since we require all participants to be tb positive at week 1
miss_group$miss1 <- NA
miss_group <- miss_group %>% dplyr::mutate(miss1 = case_when(is.na(culture_conversion_sputum_specimen_2) | is.na(culture_conversion_sputum_specimen_3) | is.na(culture_conversion_sputum_specimen_4) ~ "I", is.na(culture_conversion_sputum_specimen_2) == F & is.na(culture_conversion_sputum_specimen_3) == F & is.na(culture_conversion_sputum_specimen_4) == F ~ "C"))

#create variable miss2 for 5-8 week window
miss_group$miss2 <- NA
miss_group <- miss_group %>% dplyr::mutate(miss2 = case_when(is.na(culture_conversion_sputum_specimen_5) & is.na(culture_conversion_sputum_specimen_6) & is.na(culture_conversion_sputum_specimen_7) &
                                                               is.na(culture_conversion_sputum_specimen_8) ~ "M", is.na(culture_conversion_sputum_specimen_5) == F & is.na(culture_conversion_sputum_specimen_6) == F & is.na(culture_conversion_sputum_specimen_7) == F &
                                                               is.na(culture_conversion_sputum_specimen_8) == F ~ "C", is.na(miss2) ~ "I"))

#create variable miss3 for 9-12 week window
miss_group$miss3 <- NA
miss_group <- miss_group %>% dplyr::mutate(miss3 = case_when(is.na(culture_conversion_sputum_specimen_9) & is.na(culture_conversion_sputum_specimen_10) & is.na(culture_conversion_sputum_specimen_11) &
                                                               is.na(culture_conversion_sputum_specimen_12) ~ "M", is.na(culture_conversion_sputum_specimen_9) == F & is.na(culture_conversion_sputum_specimen_10) == F & is.na(culture_conversion_sputum_specimen_11) == F &
                                                               is.na(culture_conversion_sputum_specimen_12) == F ~ "C", is.na(miss3) ~ "I"))


#create indicator for conversion
miss_group$status <- ifelse(miss_group$conversion_week == -1, 0, 1)

miss_group_table <- miss_group %>% dplyr:::group_by(miss1, miss2, miss3) %>% summarise(n = n(), perc = round((n/238)*100, digits=1), convert = sum(status, na.rm=T)) %>% dplyr::mutate(total = paste0(n, "(", perc, ")"), conv_perc = round((convert/141)*100, digits=1), total_conv = paste0(convert, "(", conv_perc, ")") ) %>% dplyr::arrange(miss1, miss2, miss3) %>% dplyr::select(miss1, miss2, miss3, total, total_conv) 
colnames(miss_group_table) <- c("Weeks 1-4", "Weeks 5-8", "Weeks 9-12", "Total N(%)", "Conversion N(%)")

print(miss_group_table)  %>%  kable(caption = "Missing data patterns (C = no samples missing, I = at least 1 sample missing, M = all samples missing)") %>% kable_styling()
write.csv(data.frame(miss_group_table), "C:\\Users\\smala\\OneDrive\\Documents\\TRUST\\imputeTBculture\\output/missing_patterns.csv", row.names = F)

