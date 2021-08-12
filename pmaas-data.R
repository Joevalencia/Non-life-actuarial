#
# Auto GLM chapter for
# Predictive Modeling Applications in Actuarial Science
# Volume II
#
# Ernesto Schirmacher
#

#
# Source functions needed
#
source("./src/fns.R")

#
# DATA IMPORT AND MANIPULATION
#

# Define the data path and filename
data.path <- "./data"
data.fn <- "sim-modeling-dataset.csv"

#
# IMPORT DATA
#

# Define columnn class for dataset
colCls <- c("integer",         # row id
            "character",       # analysis year
            "numeric",         # exposure
            "character",       # new business / renewal business
            "numeric",         # driver age (continuous)
            "character",       # driver age (categorical)
            "character",       # driver gender
            "character",       # marital status
            "numeric",         # years licensed (continuous)
            "character",       # years licensed (categorical)
            "character",       # ncd level
            "character",       # region
            "character",       # body code
            "numeric",         # vehicle age (continuous)
            "character",       # vehicle age (categorical)
            "numeric",         # vehicle value
            "character",       # seats
            rep("numeric", 6), # ccm, hp, weight, length, width, height (all continuous)
            "character",       # fuel type
            rep("numeric", 3)  # prior claims, claim count, claim incurred (all continuous)
)

# read in the data with the appropriate column classes
dta <- read.csv(paste(data.path, data.fn, sep = "/"),
                colClasses = colCls)

# Clean workspace
rm(data.path, data.fn, colCls)

# response and miscelaneous variables
res.vars <- c("clm.count", "clm.incurred")
mis.vars <- c("row.id", "exposure")

# add a random variable to every record and include
# the variable in the list of miscelaneous variables
set.seed(1092387456)
N <- dim(dta)[1]
dta$rnd <- runif(N, min = 0, max = 1)
mis.vars <- c(mis.vars, "rnd")

# clean up the workspace
rm(N)

# add indicators for training and validation datasets.
# we want about 60% in the training set and 40% in the
# set. Add train and valid variables names to the
# miscelaneous list
set.seed(564738291)
u <- runif(dim(dta)[1], min = 0, max = 1)
dta$train <- u < 0.6
dta$valid <- !(dta$train)
mis.vars <- c(mis.vars, "train", "valid")

# clean the workspace
rm(u)

#
# Roll high values of years licensed (categorical variable) into a lower value
idx <- which(dta$yrs.licensed > 7)
dta$yrs.lic[idx] <- "8+"

# clean up
rm(idx)

#
# Roll high values of veh.age (categorical variable) into a lower value
idx <- which(dta$vehicle.age > 13)
dta$veh.age[idx] <- "14+"

# clean up
rm(idx)

#
# Roll high values of seats (categorical variable) into a lower value
idx <- dta$seats %in% c("6", "7", "8", "9")
dta$seats[idx] <- "6+"

# clean up
rm(idx)

# Set the character columns to class factor for the whole dataset
dta$year <- as.factor(dta$year)
dta$nb.rb <- as.factor(dta$nb.rb)
dta$drv.age <- as.factor(dta$drv.age)
dta$driver.gender <- as.factor(dta$driver.gender)
dta$marital.status <- as.factor(dta$marital.status)
dta$yrs.lic <- as.factor(dta$yrs.lic)
dta$ncd.level <- as.factor(dta$ncd.level)
dta$region <- as.factor(dta$region)
dta$body.code <- as.factor(dta$body.code)
dta$veh.age <- as.factor(dta$veh.age)
dta$seats <- as.factor(dta$seats)
dta$fuel.type <- as.factor(dta$fuel.type)


# Set the base level of categorical variables to
# the level with largest exposure for the whole dataset
dta$nb.rb <- re.lev(dta$nb.rb, expo = dta$exposure)
dta$drv.age <- re.lev(dta$drv.age, expo = dta$exposure)
dta$driver.gender <- re.lev(dta$driver.gender, expo = dta$exposure)
dta$marital.status <- re.lev(dta$marital.status, expo = dta$exposure)
dta$yrs.lic <- re.lev(dta$yrs.lic, expo = dta$exposure)
dta$ncd.level <- re.lev(dta$ncd.level, expo = dta$exposure)
dta$region <- re.lev(dta$region, expo = dta$exposure)
dta$body.code <- re.lev(dta$body.code, expo = dta$exposure)
dta$veh.age <- re.lev(dta$veh.age, expo = dta$exposure)
dta$seats <- re.lev(dta$seats, expo = dta$exposure)
dta$fuel.type <- re.lev(dta$fuel.type, expo = dta$exposure)


# Set the base level for analysis period to the last year
dta$year <- relevel(dta$year, ref = "2013")


# Create a new variable for length, height, width in units of decimeters
dta$len.dm <- dta$length * 10
dta$hei.dm <- dta$height * 10
dta$wid.dm <- dta$width  * 10
mis.vars <- c(mis.vars, "length", "height", "width")


# Create a new grouping for a categorical variable
dta$drv.age.gr1 <- ""

dta$drv.age.gr1[dta$driver.age %in% 18:33] <- "18-33"
dta$drv.age.gr1[dta$driver.age %in% 34:38] <- "34-38"
dta$drv.age.gr1[dta$driver.age %in% 39:43] <- "39-43"
dta$drv.age.gr1[dta$driver.age %in% 44:50] <- "44-50"
dta$drv.age.gr1[dta$driver.age %in% 51:60] <- "51-60"
dta$drv.age.gr1[dta$driver.age %in% 61:69] <- "61-69"
dta$drv.age.gr1[dta$driver.age %in% 70:93] <- "70+"

dta$drv.age.gr1 <- as.factor(dta$drv.age.gr1)
dta$drv.age.gr1 <- relevel(dta$drv.age.gr1, "44-50")

# Create another grouping for a categorical variable
dta$drv.age.gr2 <- ""

dta$drv.age.gr2[dta$driver.age %in% 18:22] <- "18-22"
dta$drv.age.gr2[dta$driver.age %in% 23:27] <- "23-27"
dta$drv.age.gr2[dta$driver.age %in% 28:32] <- "28-32"
dta$drv.age.gr2[dta$driver.age %in% 33:37] <- "33-37"
dta$drv.age.gr2[dta$driver.age %in% 38:42] <- "38-42"
dta$drv.age.gr2[dta$driver.age %in% 43:47] <- "43-47"
dta$drv.age.gr2[dta$driver.age %in% 48:52] <- "48-52"
dta$drv.age.gr2[dta$driver.age %in% 53:57] <- "53-57"
dta$drv.age.gr2[dta$driver.age %in% 58:62] <- "58-62"
dta$drv.age.gr2[dta$driver.age %in% 63:99] <- "63+"

dta$drv.age.gr2 <- as.factor(dta$drv.age.gr2)
dta$drv.age.gr2 <- relevel(dta$drv.age.gr2, "38-42")

# Define a new grouping of region
#
dta$region.g1 <- ""

rg <- c("21","26","17")
dta$region.g1[dta$region %in% rg] <- "R0"

rg <- c("27","37","34","28")
dta$region.g1[dta$region %in% rg] <- "R1"

rg <- c("6", "29", "30")
dta$region.g1[dta$region %in% rg] <- "R2"

rg <- c("12","35","4","11","5","3","13")
dta$region.g1[dta$region %in% rg] <- "R3"

rg <- c("38","16","10","9","15")
dta$region.g1[dta$region %in% rg] <- "R4"

rg <- c("22","1","31","25")
dta$region.g1[dta$region %in% rg] <- "R5"

rg <- c("14","36","32","24","19","23","33","18")
dta$region.g1[dta$region %in% rg] <- "R6"

rg <- c("7", "2", "20")
dta$region.g1[dta$region %in% rg] <- "R7"

rg <- c("8")
dta$region.g1[dta$region %in% rg] <- "R8"

dta$region.g1 <- as.factor(dta$region.g1)
dta$region.g1 <- relevel(dta$region.g1, ref = "R0")

# Clean workspace
rm(rg)


# Add severity variable
dta$sev <- NA
sv <- dta$clm.count > 0
dta[sv, "sev"] <- dta[sv, "clm.incurred"] / dta[sv, "clm.count"]
mis.vars <- c(mis.vars, "sev")

# Clean workspace
rm(sv)


# Create two linear basis elements for vehicle.value
dta$veh.val.q15 <- pmax(0, dta$vehicle.value - 15)
dta$veh.val.q35 <- pmax(0, dta$vehicle.value - 35)

# Create linear basis elements for dirver.age
dta$driver.age.q35 <- pmax(0, dta$driver.age - 35)
dta$driver.age.q49 <- pmax(0, dta$driver.age - 49)
dta$driver.age.q59 <- pmax(0, dta$driver.age - 59)


# Create a factor out of hp and relevel to max exposure
hp.cat <- cut(dta$hp, breaks = c(0,65,70,75,85,90,100,105,110,120,140,200))
hp.cat <- re.lev(hp.cat, expo = dta$exposure)
dta$hp.cat <- hp.cat

# Clean workspace
rm(hp.cat)

#
# From the Pure Premium Section define final 
# frequency and severity models, modify their
# parameters, and score the entire dataset.
#

# Final indicated frequency model
fq.m <- glm(clm.count ~ year + ncd.level + drv.age.gr2 + yrs.lic + region.g1 + prior.claims,
            family = poisson(link = "log"),
            data = dta,
            subset = train,
            offset = log(exposure))

# Final indicated severity model
sv.m <- glm(clm.incurred ~ year + marital.status + driver.gender + weight + body.code,
            family = Gamma(link = "log"),
            data = dta,
            subset = train & clm.count > 0,
            offset = log(clm.count))

# grab frequency coefficients and make modifications
fq.c <- grab.coef(fq.m)

fq.c[["ncd.level"]]["2"] <- -0.05

fq.c[["drv.age.gr2"]]["58-62"] <- 0.20
fq.c[["drv.age.gr2"]]["63+"] <- 0.25

fq.c[["yrs.lic"]]["4"]  <- -0.255
fq.c[["yrs.lic"]]["6"]  <- -0.500
fq.c[["yrs.lic"]]["7"]  <- -0.500
fq.c[["yrs.lic"]]["8+"] <- -0.500

# Score the entire dataset for frequency and determine
# the adjustment on training data only, and modify 
# scoring parameters
dta$efq <- score(fq.c, newdata = dta, offset = log(dta[,"exposure"]))$mu
adj <- -log(sum(dta[dta$train, "efq"])/ sum(dta[dta$train, "clm.count"]))
fq.c[["Base"]] <- fq.c[["Base"]] + adj

# grab severity coefficients (no modifications for severity)
sv.c <- grab.coef(sv.m)

# score the entire dataset for frequency, severity, and calculate pure premium
dta$efq <- score(fq.c, newdata = dta, offset = log(dta[,"exposure"]))$mu
dta$esv <- score(sv.c, newdata = dta)$mu
dta$epp <- dta$efq * dta$esv

#
# SAVE THE DATA in an RObject
#
save(dta, file = "./data/all-data.RData")
