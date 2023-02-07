library(parallel)
library(data.table)
library(tidyr)
library(dplyr)
library(splines)
library(fst)


# scenarios
scenarios <- expand.grid(dual = c("both", "high", "low"),
                         race = c("all", "white", "black"))
# scenarios <- expand.grid(dual = c("high", "low"),
#                          race = c("white", "black", "asian", "hispanic", "other"))
scenarios$dual <- as.character(scenarios$dual)
scenarios$race <- as.character(scenarios$race)
a.vals <- seq(2, 31, length.out = 146)

### Fit Outcome Models

# Load/Save models
dir_data = '/n/dominici_nsaph_l3/Lab/projects/analytic/erc_strata/qd/'
dir_mod = '/n/dominici_nsaph_l3/projects/kjosey-erc-strata/Output/Strata_Data/'

options(digits = 5)
for (i in c(1,4,7,6,5,9,8)) {
  
  scenario <- scenarios[i,]
  load(paste0(dir_data, scenario$dual, "_", scenario$race, ".RData"))
  zip_data <- new_data$x
  mean(zip_data$mean_bmi,na.rm=T); sd(zip_data$mean_bmi,na.rm=T)
  mean(zip_data$smoke_rate,na.rm=T)*100;  sd(zip_data$smoke_rate*100,na.rm=T)
  mean(zip_data$poverty,na.rm=T)*100; sd(zip_data$poverty*100,na.rm=T)
  mean(zip_data$education,na.rm=T)*100; sd(zip_data$education*100,na.rm=T)
  mean(zip_data$pct_owner_occ,na.rm=T)*100; sd(zip_data$pct_owner_occ*100,na.rm=T)
  mean(zip_data$pct_blk,na.rm=T)*100; sd(zip_data$pct_blk*100,na.rm=T)
  mean(zip_data$hispanic,na.rm=T)*100; sd(zip_data$hispanic*100,na.rm=T)
  mean(zip_data$medhouseholdincome,na.rm=T)/1000; sd(zip_data$medhouseholdincome/1000,na.rm=T)
  mean(zip_data$medianhousevalue,na.rm=T)/1000; sd(zip_data$medianhousevalue/1000,na.rm=T)
  mean(zip_data$popdensity,na.rm=T); sd(zip_data$popdensity,na.rm=T)
  mean(zip_data$summer_tmmx,na.rm=T)-273.15; sd(zip_data$summer_tmmx,na.rm=T)
  mean(zip_data$winter_tmmx,na.rm=T)-273.15; sd(zip_data$winter_tmmx,na.rm=T)
  mean(zip_data$summer_rmax,na.rm=T); sd(zip_data$summer_rmax,na.rm=T)
  mean(zip_data$winter_rmax,na.rm=T); sd(zip_data$winter_rmax,na.rm=T)
  mean(zip_data$pm25,na.rm=T); sd(zip_data$pm25,na.rm=T)
  mean(zip_data$regionNORTHEAST)*100
  mean(zip_data$regionSOUTH)*100
  (1-mean(zip_data$regionNORTHEAST)-mean(zip_data$regionSOUTH)-mean(zip_data$regionWEST))*100
  mean(zip_data$regionWEST)*100
  
}


f <- list.files("/n/dominici_nsaph_l3/Lab/data/ci3_health_data/medicare/mortality/1999_2016/wu/cache_data/merged_by_year_v2",
                pattern = "\\.fst", full.names = TRUE)
# 
myvars <- c("qid","year","zip","sex","race","age","dual","entry_age_break","statecode","followup_year",
            "followup_year_plus_one","dead","pm25_ensemble","mean_bmi","smoke_rate","hispanic","pct_blk",
            "medhouseholdincome","medianhousevalue","poverty","education","popdensity", "pct_owner_occ",
            "summer_tmmx","winter_tmmx","summer_rmax","winter_rmax")
# 
national_merged2016 <- rbindlist(lapply(f, read_fst, columns = myvars, as.data.table = TRUE))
national_merged2016$zip <- sprintf("%05d", national_merged2016$zip)
# 
NORTHEAST = c("NY","MA","PA","RI","NH","ME","VT","CT","NJ")
SOUTH = c("DC","VA","NC","WV","KY","SC","GA","FL","AL","TN","MS","AR","MD","DE","OK","TX","LA")
MIDWEST = c("OH","IN","MI","IA","MO","WI","MN","SD","ND","IL","KS","NE")
WEST = c("MT","CO","WY","ID","UT","NV","CA","OR","WA","AZ","NM")
# 
# # create region
national_merged2016$region=ifelse(national_merged2016$state %in% NORTHEAST, "NORTHEAST",
                                  ifelse(national_merged2016$state %in% SOUTH, "SOUTH",
                                         ifelse(national_merged2016$state %in% MIDWEST, "MIDWEST",
                                                ifelse(national_merged2016$state %in% WEST, "WEST", NA))))
# 
national_merged2016 <- national_merged2016[complete.cases(national_merged2016[,c(1:28)]) ,]

### Build Clean Aggregate Data Set - Outcomes, Exposures, Covariates, and Offsets

national_merged2016$time_count <- rep(1, nrow(national_merged2016))
national_merged2016$sex <- national_merged2016$sex - 1
colnames(national_merged2016)[3] <- c("female")

# label unknown and american/alaska natives as "other"
national_merged2016$race[national_merged2016$race == 6] <- 3
national_merged2016$race[national_merged2016$race == 0] <- 3

# collapse entry age breaks
national_merged2016$age_break[national_merged2016$entry_age_break %in% c(1,2)] <- "[65,75)" 
national_merged2016$age_break[national_merged2016$entry_age_break %in% c(3,4)] <- "[75,85)" 
national_merged2016$age_break[national_merged2016$entry_age_break %in% c(5,6)] <- "[85,95)" 
national_merged2016$age_break[national_merged2016$entry_age_break %in% c(7,8)] <- "[95,125)"

national_merged2016$entry_age_break <- national_merged2016$age_break
national_merged2016$age_break <- NULL

qid_data <- national_merged2016 %>% group_by(qid) %>% filter(row_number()==1)

national_merged2016_sub <- subset(national_merged2016, race == 1 & dual == 1)
qid_data_sub <- subset(qid_data, race == 1 & dual == 1)
nrow(qid_data_sub)
nrow(national_merged2016_sub)
nrow(qid_data_sub)/nrow(qid_data)*100
nrow(national_merged2016_sub)/nrow(national_merged2016)*100
#mean(qid_data_sub$entry_age);sd(qid_data$entry_age)
table(qid_data_sub$entry_age_break)/nrow(qid_data_sub)*100
table(qid_data_sub$sex)/nrow(qid_data_sub)*100
table(qid_data_sub$dual)/nrow(qid_data_sub)*100
sum(national_merged2016_sub$dead)
sum(national_merged2016_sub$dead)/sum(national_merged2016$dead)*100

# calculate median follow-up year
followup_data_sub <- national_merged2016_sub %>% group_by(qid) %>% summarise(max_followup = max(followup_year))
median(followup_data_sub$max_followup)
