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

# creates region
national_merged2016$region=ifelse(national_merged2016$state %in% NORTHEAST, "NORTHEAST",
                                     ifelse(national_merged2016$state %in% SOUTH, "SOUTH",
                                            ifelse(national_merged2016$state %in% MIDWEST, "MIDWEST",
                                                   ifelse(national_merged2016$state %in% WEST, "WEST", NA))))

national_merged2016 <- national_merged2016[complete.cases(national_merged2016[,c(1:27)]) ,]
save(national_merged2016, file = "~/shared_space/ci3_analysis/josey_erc_strata/Data/national_merged2016.RData")

### Build Clean Aggregate Data Set

## Outcomes, Exposures, Covariates, and Offsets

load("~/shared_space/ci3_analysis/josey_erc_strata/Data/national_merged2016.RData")
national_merged2016$time_count <- rep(1, nrow(national_merged2016))
national_merged2016$age_break <- cut(national_merged2016$age, c(65,70,75,80,85,90,95,125), right = FALSE)
national_merged2016$female <- national_merged2016$sex - 1
colnames(national_merged2016)[12] <- "pm25"

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

save(aggregate_data, file = "~/shared_space/ci3_analysis/josey_erc_strata/Data/aggregate_data.RData")

rm(dead_personyear, confounders, aggregate_data); gc()

### Create Strata Data

create_strata <- function(data, dual = c(0,1,2), race = c("all", "white", "black", "asian", "hispanic")) {
  
  zip_cov <- c("pm25", "mean_bmi", "smoke_rate", "hispanic", "pct_blk", "medhouseholdincome", "medianhousevalue", "poverty", "education",
               "popdensity", "pct_owner_occ", "summer_tmmx", "winter_tmmx", "summer_rmax", "winter_rmax", "region")
  
  if (dual == 0) {
    dual0 <- 0
  } else if (dual == 1) {
    dual0 <- 1
  } else {
    dual0 <- c(0,1)
  }
  
  if (race == "white") {
    race0 <- 1
  } else if (race == "black") {
    race0 <- 2
  } else if (race == "asian") {
    race0 <- 4
  } else if (race == "hispanic") {
    race0 <- 5
  } else {
    race0 <- c(1,2,3,4,5)
  }
  
  sub_data <- subset(data, race %in% race0 & dual %in% dual0)
  
  # Covariates and Outcomes
  w <- data.table(zip = sub_data$zip, year = sub_data$year, race = sub_data$race,
                  female = sub_data$female, dual = sub_data$dual, age_break = sub_data$age_break,
                  dead = sub_data$dead, time_count = sub_data$time_count)[
                    ,lapply(.SD, sum), by = c("zip", "year", "race", "female", "dual", "age_break")]
  
  x <- data.table(zip = sub_data$zip, year = sub_data$year, model.matrix(~ ., data = sub_data[,zip_cov])[,-1])[
                    ,lapply(.SD, min), by = c("zip", "year")]
  
  return(list(w = w, x = x))
  
}

# scenarios
scenarios <- expand.grid(dual = c(0, 1, 2), race = c("white", "black", "asian", "hispanic", "all"))
scenarios$dual <- as.numeric(scenarios$dual)
scenarios$race <- as.character(scenarios$race)

# Save Location
dir_data = '/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Data/qd/'

## QD Strata

load(paste0(dir_data,"aggregate_data.RData"))
aggregate_data$zip <- factor(aggregate_data$zip)
aggregate_data$year <- factor(aggregate_data$year)
aggregate_data$female <- as.numeric(aggregate_data$female)
aggregate_data$dual <- as.numeric(aggregate_data$dual)
aggregate_data$region <- factor(aggregate_data$region)
aggregate_data$age_break <- factor(aggregate_data$age_break)
aggregate_data$race <- factor(aggregate_data$race)

lapply(1:nrow(scenarios), function(i, ...) {
  
  scenario <- scenarios[i,]
  new_data <- create_strata(data = aggregate_data, dual = scenario$dual, race = scenario$race)
  save(new_data, file = paste0(dir_data, scenario$dual, "_", scenario$race, ".RData"))
  
})
