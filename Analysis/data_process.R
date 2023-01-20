# Generate all required data from different data sources (see Table S1).
library(fst)
library(data.table)
library(tidyr)
library(dplyr)
library(parallel)
library(foreign)

### Big Data Cleaning (BEWARE!)

f <- list.files("/nfs/nsaph_ci3/ci3_health_data/medicare/mortality/1999_2016/wu/cache_data/merged_by_year_v2",
                pattern = "\\.fst", full.names = TRUE)

myvars <- c("qid","year","zip","sex","race","age","dual","entry_age_break","statecode","followup_year",
            "followup_year_plus_one","dead","pm25_ensemble","mean_bmi","smoke_rate","hispanic","pct_blk",
            "medhouseholdincome","medianhousevalue","poverty","education","popdensity", "pct_owner_occ",
            "summer_tmmx","winter_tmmx","summer_rmax","winter_rmax")

national_merged2016 <- rbindlist(lapply(f, read_fst, columns = myvars, as.data.table = TRUE))
national_merged2016$zip <- sprintf("%05d", national_merged2016$zip)

NORTHEAST = c("NY","MA","PA","RI","NH","ME","VT","CT","NJ")
SOUTH = c("DC","VA","NC","WV","KY","SC","GA","FL","AL","TN","MS","AR","MD","DE","OK","TX","LA")
MIDWEST = c("OH","IN","MI","IA","MO","WI","MN","SD","ND","IL","KS","NE")
WEST = c("MT","CO","WY","ID","UT","NV","CA","OR","WA","AZ","NM")

# create region
national_merged2016$region=ifelse(national_merged2016$state %in% NORTHEAST, "NORTHEAST",
                                     ifelse(national_merged2016$state %in% SOUTH, "SOUTH",
                                            ifelse(national_merged2016$state %in% MIDWEST, "MIDWEST",
                                                   ifelse(national_merged2016$state %in% WEST, "WEST", NA))))

national_merged2016 <- national_merged2016[complete.cases(national_merged2016[,c(1:27)]) ,]
save(national_merged2016, file = "/n/dominici_nsaph_l3/Lab/projects/analytic/erc_strata/national_merged2016.RData")

### Build Clean Aggregate Data Set - Outcomes, Exposures, Covariates, and Offsets

load("/n/dominici_nsaph_l3/Lab/projects/analytic/erc_strata/national_merged2016.RData")
national_merged2016$time_count <- rep(1, nrow(national_merged2016))
national_merged2016$age_break <- cut(national_merged2016$age, c(65,75,85,95,125), right = FALSE)
national_merged2016$female <- national_merged2016$sex - 1
colnames(national_merged2016)[12] <- "pm25"
national_merged2016$race[national_merged2016$race == 6] <- 3

dead_personyear <- aggregate(data.frame(dead = national_merged2016$dead,
                                        time_count = national_merged2016$time_count),
                             by=list(zip = national_merged2016$zip,
                                     year = national_merged2016$year,
                                     female = national_merged2016$female,
                                     race = national_merged2016$race,
                                     dual = national_merged2016$dual,
                                     age_break = national_merged2016$age_break),
                             FUN=sum)

new_data <- national_merged2016 %>% distinct(zip, year, female, race, dual, age_break, .keep_all = TRUE)
confounders <- new_data[,c(1,2,4,6,12:27,29,30)]

rm(national_merged2016, new_data); gc()

aggregate_data <- merge(dead_personyear, confounders, by = c("zip","year","female","race","dual","age_break"))
aggregate_data <- aggregate_data[complete.cases(aggregate_data),]

save(aggregate_data, file = "/n/dominici_nsaph_l3/Lab/projects/analytic/erc_strata/aggregate_data.RData")

rm(dead_personyear, confounders, aggregate_data); gc()

# Generate all required data from different data sources (see Table S1).
library(fst)
library(data.table)
library(tidyr)
library(dplyr)
library(parallel)
library(foreign)

### Big Data Cleaning (BEWARE! fasse_bigmem required!)

f <- list.files("/nfs/nsaph_ci3/ci3_health_data/medicare/mortality/1999_2016/wu/cache_data/merged_by_year_v2",
                pattern = "\\.fst", full.names = TRUE)

myvars <- c("qid","year","zip","sex","race","age","dual","entry_age_break","statecode","followup_year",
            "followup_year_plus_one","dead","pm25_ensemble","mean_bmi","smoke_rate","hispanic","pct_blk",
            "medhouseholdincome","medianhousevalue","poverty","education","popdensity", "pct_owner_occ",
            "summer_tmmx","winter_tmmx","summer_rmax","winter_rmax")

national_merged2016 <- rbindlist(lapply(f, read_fst, columns = myvars, as.data.table = TRUE))
national_merged2016$zip <- sprintf("%05d", national_merged2016$zip)

NORTHEAST = c("NY","MA","PA","RI","NH","ME","VT","CT","NJ")
SOUTH = c("DC","VA","NC","WV","KY","SC","GA","FL","AL","TN","MS","AR","MD","DE","OK","TX","LA")
MIDWEST = c("OH","IN","MI","IA","MO","WI","MN","SD","ND","IL","KS","NE")
WEST = c("MT","CO","WY","ID","UT","NV","CA","OR","WA","AZ","NM")

# create region
national_merged2016$region=ifelse(national_merged2016$state %in% NORTHEAST, "NORTHEAST",
                                  ifelse(national_merged2016$state %in% SOUTH, "SOUTH",
                                         ifelse(national_merged2016$state %in% MIDWEST, "MIDWEST",
                                                ifelse(national_merged2016$state %in% WEST, "WEST", NA))))

national_merged2016 <- national_merged2016[complete.cases(national_merged2016[,c(1:27)]) ,]
save(national_merged2016, file = "/n/dominici_nsaph_l3/Lab/projects/analytic/erc_strata/national_merged2016.RData")

### Build Clean Aggregate Data Set - Outcomes, Exposures, Covariates, and Offsets

load("/n/dominici_nsaph_l3/Lab/projects/analytic/erc_strata/national_merged2016.RData")
national_merged2016$time_count <- rep(1, nrow(national_merged2016))
national_merged2016$age_break <- cut(national_merged2016$age, c(65,75,85,95,125), right = FALSE)
national_merged2016$female <- national_merged2016$sex - 1
colnames(national_merged2016)[12] <- "pm25"
national_merged2016$race[national_merged2016$race == 6] <- 3

dead_personyear <- aggregate(data.frame(dead = national_merged2016$dead,
                                        time_count = national_merged2016$time_count),
                             by=list(zip = national_merged2016$zip,
                                     year = national_merged2016$year,
                                     female = national_merged2016$female,
                                     race = national_merged2016$race,
                                     dual = national_merged2016$dual,
                                     age_break = national_merged2016$age_break),
                             FUN=sum)

new_data <- national_merged2016 %>% distinct(zip, year, female, race, dual, age_break, .keep_all = TRUE)
confounders <- new_data[,c(1,2,4,6,12:27,29,30)]

rm(national_merged2016, new_data); gc()

aggregate_data <- merge(dead_personyear, confounders, by = c("zip","year","female","race","dual","age_break"))
aggregate_data <- aggregate_data[complete.cases(aggregate_data),]

save(aggregate_data, file = "/n/dominici_nsaph_l3/Lab/projects/analytic/erc_strata/aggregate_data.RData")

rm(dead_personyear, confounders, aggregate_data); gc()
