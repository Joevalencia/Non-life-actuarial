#
# Auto GLM chapter for
# Predictive Modeling Applications in Actuarial Science
# Volume II
#
# Ernesto Schirmacher
#

#
# Source some utility functions needed later in the text
#
source("./src/fns.R")

# load the prepared dataset
load("./data/all-data.RData")

####-####+####-####+####-####+####-####+####-####+####-####+####-####+####-####+
# DATA CHARACTERISTICS
####-####+####-####+####-####+####-####+####-####+####-####+####-####+####-####+

# dimensions of the dataset
dim(dta)

# Create a summary table of frequency and severity
# by analysis period
yr.expo <- with(dta, tapply(exposure, year, sum))
yr.clm.count <- with(dta, tapply(clm.count, year, sum))
yr.clm.incr <- with(dta, tapply(clm.incurred, year, sum))
yr.summary <- cbind(
  exposure = round(yr.expo,1),
  clm.count = yr.clm.count,
  clm.incurred = round(yr.clm.incr,0),
  frequency = round(yr.clm.count / yr.expo, 3),
  severity = round(yr.clm.incr / yr.clm.count, 1))

yr.summary <- rbind(yr.summary,
                    total = c(
                      round(sum(yr.expo),1),
                      sum(yr.clm.count),
                      round(sum(yr.clm.incr),0),
                      round(sum(yr.clm.count)/sum(yr.expo),3),
                      round(sum(yr.clm.incr)/sum(yr.clm.count),1)))

yr.summary <- yr.summary[c(2,3,4,1,5),] #reorder the rows


# Clean the global environment
rm(yr.expo, yr.clm.count, yr.clm.incr, yr.summary)

#
# Volume of business over the years
expo <- with(dta, tapply(exposure, year, sum))
round(expo["2012"] / expo["2010"] - 1, 2)
round(expo["2013"] / expo["2012"] - 1, 2)

# Clean the workspace
rm(expo)


####-####+####-####+####-####+####-####+####-####+####-####+####-####+####-####+
# EXPLORATORY DATA ANALYSIS FOR FREQUENCY
####-####+####-####+####-####+####-####+####-####+####-####+####-####+####-####+

#
# Frequency for the training set
tmp <- subset(dta, subset = train, select = c("exposure", "clm.count"))
ex <- with(tmp, sum(exposure))
cc <- with(tmp, sum(clm.count))
fq <- cc / ex

# clean up
rm(tmp, ex, cc, fq)

#
# Empirical frequency by nb.rb  in training dataset
tmp <- subset(dta, subset = train, select = c("year", "nb.rb", "exposure", "clm.count"))
ex <- with(tmp, tapply(exposure, list(year, nb.rb), sum))[c(2,3,4,1),]
ex <- rbind(ex, total = apply(ex, 2, sum))
cc <- with(tmp, tapply(clm.count, list(year, nb.rb), sum))[c(2,3,4,1),]
cc <- rbind(cc, total = apply(cc, 2, sum))
fq <- round((cc / ex)*100,1)

#clean up the workspace
rm(tmp, ex, cc, fq)

#
# Empirical frequency across entire dataset by region
ex <- with(dta, tapply(exposure, region, sum))
cc <- with(dta, tapply(clm.count, region, sum))
f <- sort(round((cc / ex)*100, 1))

# Clean up the workspace
rm(ex, cc, f)

#
# Driver age simple characteristics over training dataset
tmp <- sort(unique(dta$driver.age[dta$train]))
diff(tmp)
length(tmp)

# clean up
rm(tmp)

# Unique ages over entire dataset
length(unique(dta$driver.age))
setdiff(unique(dta$driver.age), unique(dta$driver.age[dta$train]))


#
# Empirical frequency for driver age across the training dataset
tmp <- subset(dta, subset = train, select = c("driver.age", "exposure", "clm.count"))
ex <- with(tmp, tapply(exposure, driver.age, sum))
cc <- with(tmp, tapply(clm.count, driver.age, sum))
fr <- round((cc / ex)*100, 1)
summary(fr)

# Which ages have frequencies greater than 100%?
nm <- names(fr[which(fr > 100)])

# What is the data?
tmp[tmp$driver.age %in% nm,]

# sort freq
sort(fr, decreasing = TRUE)[1:5]

# Data with frequency greater than 50% but smaller than 90%
nm <- names(fr[which((fr > 90) & (fr < 100))])
tbl <- tmp[tmp$driver.age %in% nm, c("exposure", "clm.count")]
tbl
apply(tbl, 2, sum)

# show zig-zag pattern
fr[as.character(30:34)]
fr[as.character(50:70)]
fr[as.character(43:53)]

# Clean up the workspace
rm(tmp, ex, cc, fr, nm, tbl)


#
# Empirical frequency by size of engine for the entire dataset
bk <- unique(quantile(dta$ccm, probs = seq(0, 1, by = 0.05)))
bk[1] <- bk[1] - 1
ccm.d <- cut(dta$ccm, breaks = bk)
expo <- with(dta, tapply(exposure, ccm.d, sum))
clms <- with(dta, tapply(clm.count, ccm.d, sum))
freq <- round(clms / expo*100,2)
freq

# Clean the workspace
rm(bk, ccm.d, expo, clms, freq)

#
# Empirical frequency for driver gender
tmp <- subset(dta, subset = train, select = c("driver.gender", "year", "exposure", "clm.count"))
ex <- with(tmp, tapply(exposure, driver.gender, sum))
cc <- with(tmp, tapply(clm.count, driver.gender, sum))
fr <- round((cc / ex)*100, 1)
fr

# add year as a dimension
ex <- with(tmp, tapply(exposure, list(year, driver.gender), sum))[c(2,3,4,1),]
ex <- rbind(ex, total = apply(ex, 2, sum))
cc <- with(tmp, tapply(clm.count, list(year, driver.gender), sum))[c(2,3,4,1),]
cc <- rbind(cc, total = apply(cc, 2, sum))
fr <- round((cc / ex)*100, 1)
fr

# Clean the workspace
rm(tmp, ex, cc, fr)

#
# Empirical frequency for marital status
tmp <- subset(dta, subset = train, select = c("marital.status", "year", "exposure", "clm.count"))
ex <- with(tmp, tapply(exposure, list(year, marital.status), sum))[c(2,3,4,1),]
ex <- rbind(ex, total = apply(ex, 2, sum))
cc <- with(tmp, tapply(clm.count, list(year, marital.status), sum))[c(2,3,4,1),]
cc <- rbind(cc, total = apply(cc, 2, sum))
fr <- round((cc / ex)*100, 2)
fr[,c("Single", "Married", "Divorced", "Widow")]

# Clean the workspace
rm(tmp, ex, cc, fr)

#
# Is difference between married and widow significant?
N <- 10000
tmp <- subset(dta,
              subset = train &
                marital.status %in% c("Married", "Widow"),
              select = c("marital.status",
                         "exposure", "clm.count"))
f <- tmp$marital.status == "Married"
d <- numeric(N)
for(i in 1:N) {
  g <- sample(f, length(f))
  married.fq <- sum(tmp$clm.count[g]) / sum(tmp$exposure[g])
  widow.fq <- sum(tmp$clm.count[!g]) / sum(tmp$exposure[!g])
  d[i] <- married.fq - widow.fq
}
quantile(d, c(0.025, 0.05, 0.1, 0.25, 0.5,
              0.75, 0.9, 0.95, 0.975))
rm(N, tmp, f, d, i, g, married.fq, widow.fq) # clean-up

#
# Is difference between married and single significant?
N <- 10000
tmp <- subset(dta,
              subset = train &
                marital.status %in% c("Married", "Single"),
              select = c("marital.status",
                         "exposure", "clm.count"))
f <- tmp$marital.status == "Married"
d <- numeric(N)
for(i in 1:N) {
  g <- sample(f, length(f))
  married.fq <- sum(tmp$clm.count[g]) / sum(tmp$exposure[g])
  single.fq <- sum(tmp$clm.count[!g]) / sum(tmp$exposure[!g])
  d[i] <- married.fq - single.fq
}
quantile(d, c(0.005, 0.025, 0.05, 0.1, 0.25, 0.5,
              0.75, 0.9, 0.95, 0.975, 0.995))
rm(N, tmp, f, d, i, g, married.fq, single.fq) # clean-up



#
# Horsepower over the training dataset
length(unique(dta$hp[dta$train]))
range(dta$hp[dta$train])

#
# Empirical frequency for horsepower
tmp <- subset(dta, subset = train, select = c("hp", "year", "exposure", "clm.count"))
ex <- with(tmp, tapply(exposure, hp, sum))
cc <- with(tmp, tapply(clm.count, hp, sum))
fr <- round((cc / ex), 3)
sort(fr, decreasing = TRUE)[6:1]
sort(round(ex,0), decreasing = TRUE)[5:1]
round(100*sum(sort(round(ex,0), decreasing = TRUE)[5:1])/sum(ex),0)

# Clean the workspace
rm(tmp, ex, cc, fr)

#
# Empirical frequency for length of vehicle
length(unique(dta$len.dm[dta$train]))
range(dta$len.dm[dta$train])
round(mean(dta$len.dm[dta$train]), 2)
quantile(dta$len.dm[dta$train], prob = c(0, 0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975, 1))
round(quantile(dta$len.dm[dta$train], prob = 0.8)/10,1)


#
# using length in units of m or dm does not make a difference
m.dm <- glm(clm.count ~ len.dm, data = dta, family = poisson(link = "log"), offset = log(exposure))
m.m  <- glm(clm.count ~ length, data = dta, family = poisson(link = "log"), offset = log(exposure))

coef(m.dm)["len.dm"] * 10 - coef(m.m)["length"]

# clean up
rm(m.dm, m.m)


####-####+####-####+####-####+####-####+####-####+####-####+####-####+####-####+
# EXPLORATORY DATA ANALYSIS FOR SEVERITY
####-####+####-####+####-####+####-####+####-####+####-####+####-####+####-####+

#
# How many records have a claim? And how many claims are there?
sum(dta$clm.count > 0)
sum(dta$clm.count)
table(dta$clm.count)

#
# Average severity across training dataset?
tmp <- subset(dta, subset = train & clm.count > 0, select = c("year", "clm.count", "clm.incurred"))
cc <- with(tmp, tapply(clm.count, year, sum))[c(2,3,4,1)]
ci <- with(tmp, tapply(clm.incurred, year, sum))[c(2,3,4,1)]
sv <- ci / cc
c(round(sv, 2), total = sum(ci)/sum(cc))

# Clean the worspace
rm(tmp,cc, ci, sv)

#
# Average severity by new/renewal indicator
tmp <- subset(dta, subset = train & clm.count > 0, select = c("nb.rb", "clm.count", "clm.incurred"))
cc <- with(tmp, tapply(clm.count, nb.rb, sum))
ci <- with(tmp, tapply(clm.incurred, nb.rb, sum))
sv <- ci / cc
c(round(sv, 2), total = sum(ci)/sum(cc))

# Clean the worspace
rm(tmp,cc, ci, sv)

#
# Number of claims in training dataset and amount of exposure for fuel type
sum(dta$clm.count[dta$train])
with(dta[dta$train,], tapply(clm.count, fuel.type, sum))
round(with(dta[dta$train,], tapply(clm.count, fuel.type, sum)/sum(clm.count))*100,1)

round(with(dta[dta$train,], tapply(exposure, fuel.type, sum)),0)
round(with(dta[dta$train,], tapply(exposure, fuel.type, sum)/sum(exposure))*100,1)


#
# Average severity by fuel type indicator
tmp <- subset(dta, subset = train & clm.count > 0, select = c("fuel.type", "clm.count", "clm.incurred"))
cc <- with(tmp, tapply(clm.count, fuel.type, sum))
ci <- with(tmp, tapply(clm.incurred, fuel.type, sum))
sv <- ci / cc
c(round(sv, 2), total = sum(ci)/sum(cc))

# Clean the worspace
rm(tmp,cc, ci, sv)


####-####+####-####+####-####+####-####+####-####+####-####+####-####+####-####+
# FREQUENCY MODELING
####-####+####-####+####-####+####-####+####-####+####-####+####-####+####-####+

#
# Frequency for training dataset
round(sum(dta$clm.count[dta$train])/sum(dta$exposure[dta$train])*100,1)

#
# Null model for frequency on training dataset
m <- glm(clm.count ~ 1,
          data = dta, subset = train,
          family = poisson(link = "log"),
          offset = log(exposure))

summary(m)

exp(coef(m))

# Clean the workspace
rm(m)

#
# Single variable frequency model for driver.age
m <- glm(clm.count ~ driver.age,
         data = dta, subset = train,
         family = poisson(link = "log"),
         offset = log(exposure))

round(c(aic = m$aic, dev = m$deviance),0)

# Clean the workspace
rm(m)

#
# Single variable frequency model for drv.age
m <- glm(clm.count ~ drv.age,
         data = dta, subset = train,
         family = poisson(link = "log"),
         offset = log(exposure))

round(c(aic = m$aic, dev = m$deviance),0)

# Clean the workspace
rm(m)

#
# Number of unique ages in training data
length(unique(dta$drv.age[dta$train]))

#
# Base level for drv.age
levels(dta$drv.age)[1]

#
# Frequency for drv.age variable
tmp <- subset(dta, subset = train, select = c("drv.age", "exposure", "clm.count"))
ex <- with(tmp, tapply(exposure, drv.age, sum))
cc <- with(tmp, tapply(clm.count, drv.age, sum))
fr <- cc / ex
c(ex = ex["29"], cc = cc["29"], fr = fr["29"])

# Clean the workspace
rm(tmp, ex, cc, fr)

#
# Single variable frequency model for drv.age
m <- glm(clm.count ~ drv.age,
         data = dta, subset = train,
         family = poisson(link = "log"),
         offset = log(exposure))

cf <- coef(m)
(age29 <- c(cf["(Intercept)"], cf["drv.age29"]))

exp(sum(age29))

#
# Drivers age 38
round(exp(cf["(Intercept)"]),4)

# Clean the workspace
rm(m, cf, age29)

#
# drv.age with zero claims
i <- (with(dta[dta$train,], tapply(clm.count, drv.age, sum)))
dta$drv.age[dta$train][i]

rm(i) #clean up

#
# Frequency for yrs.lic on the training dataset
m <- glm(clm.count ~ yrs.lic - 1,
         data = dta, subset = train,
         family = poisson(link = "log"),
         offset = log(exposure))

cf <- round(summary(m)$coefficients[,-c(3,4)],5)
ex <- with(dta[dta$train,], tapply(exposure, yrs.lic, sum))
tbl <- cbind(cf, Exposure = round(ex,1))

# Clean the workspace
rm(m, cf, ex, tbl)

#
exp(-2.5)

#
# Year included as part of the null model
#
# Variables that are NOT predictors
res.vars <- c("clm.count", "clm.incurred")
mis.vars <- c("row.id", "exposure", "rnd", "train", "valid", "sev",
              "length", "height", "width", "year")

# Null frequency model
fq.null <- glm(clm.count ~ year,
               data = dta,
               subset = train,
               family = poisson(link = "log"),
               offset = log(exposure))
aic.null <- fq.null$aic
dev.null <- fq.null$deviance
rms.null <- RMSEP(dta$clm.count[dta$train],
                  predict(fq.null, newdata = dta[dta$train,],
                          type = "response"))

# Select predictor variables
vars <- setdiff(colnames(dta),
                c(res.vars, mis.vars,
                  "region.g1", "drv.age.gr1", "hp.cat",
                  "veh.val.pl", "driver.age.q35", "driver.age.q49",
                  "driver.age.q59", "veh.val.q15", "veh.val.q35"))


aic <- numeric(0)
dev <- numeric(0)
rms <- numeric(0)

for(v in vars) {
  fmla <- as.formula(paste("clm.count ~ year", v, sep = " + "))
  f <- glm(fmla,
           data = dta[dta$train,],
           family = poisson(link = "log"),
           offset = log(exposure))
  aic[v] <- f$aic
  dev[v] <- f$deviance
  rms[v] <- RMSEP(dta$clm.count[dta$train],
                  predict(f, newdata = dta[dta$train,],
                          type = "response"))
}

# sort variables by the distance from null* model
# using aic and dev
sort(sqrt((aic.null - aic)^2 + (dev.null - dev)^2),
     decreasing = TRUE)

# sort variables by the distance from null* model
# using aic and rmsep
sort(sqrt((aic.null - aic)^2 + (rms.null - rms)^2),
     decreasing = TRUE)

# clean up
rm(res.vars,mis.vars,fq.null,aic.null,
   dev.null,vars,aic,dev,v,fmla, f, rms, rms.null)

#
# Number of regions over traininig and entire dataset
length(unique(dta$region))
length(unique(dta$region[dta$train]))

#
# Frequency for yrs.lic on the training dataset
m1 <- glm(clm.count ~ year + region.g1 + ncd.level + driver.age + yrs.lic + prior.claims,
         data = dta, subset = train,
         family = poisson(link = "log"),
         offset = log(exposure))
summary(m1)

#
# Analysis of deviance table
anova(m1, test = "Chisq")

# Coordinates aic, dev
c(aic = m1$aic, dev = m1$deviance)

#
# Frequency. Add additional variables
m2 <- update(m1, . ~ . + marital.status + body.code + nb.rb)
summary(m2)
anova(m2, test = "Chisq")

#
# Analysis of deviance table between the two models
anova(m1, m2, test = "Chisq")

#
# Analysis of deviance for m2
anova(m2, test = "Chisq")

#
# Add seats to m2 and do analysis of deviance
m3 <- update(m2, . ~ . + seats)
anova(m2,m3, test = "Chisq")


# Clean the workspace
rm(m1,m2,m3)

#
# CROSS VALIDATION FOR A FREQUENCY MODEL
#
vars <- c("year", "region.g1", "ncd.level", "drv.age.gr1", "prior.claims",
          "yrs.lic", "ccm", "hei.dm", "wid.dm", "len.dm", "nb.rb", "driver.gender",
          "marital.status", "hp.cat", "fuel.type", "veh.age", "vehicle.value", "seats")
K <- length(vars)
t.rmsep <- numeric(K)
v.rmsep <- numeric(K)

for(j in 1:K) {
  predictors <- paste(vars[1:j], sep = "", collapse = " + ")
  fmla <- as.formula(paste("clm.count ~ ", predictors, sep = ""))
  m <- glm(fmla,
           data = dta,
           subset = train,
           family = poisson(link = "log"),
           offset = log(exposure))
  L <- length(coef(m))
  M <- 20
  bk <- seq(0,1,length = M+1)
  rmsep <- numeric(length(bk) - 1)
  mt <- matrix(0, nrow = L, ncol = M)
  for(i in 2:length(bk)) {
    m <- glm(fmla,
             data = dta,
             subset = train & ((rnd <= bk[i-1]) | (rnd > bk[i])),
             family = poisson(link = "log"),
             offset = log(exposure))
    mt[,i-1] <- coef(m)
    pred <- predict(m, newdata = dta[dta$train & ((dta$rnd > bk[i-1]) & (dta$rnd <= bk[i])), ], type = "response")
    actu <- dta[dta$train & ((dta$rnd > bk[i-1]) & (dta$rnd <= bk[i])), "clm.count"]
    rmsep[i-1] <- RMSEP(actu, pred)
  }
  dimnames(mt) <- list(names(coef(m)), paste("v",1:(length(bk)-1), sep = ""))
  avg.coef <- apply(mt, 1, mean)
  m$coefficients <- avg.coef

  t.rmsep[j] <- mean(rmsep)
  o <- grab.coef(m)
  scr <- score(o, newdata = dta[dta$valid,], offset = log(dta[dta$valid, "exposure"]))
  pred <- scr$mu

  actu <- dta[dta$valid, "clm.count"]
  v.rmsep[j] <- RMSEP(actu, pred)
}
op <- par(mar = c(2,4,1,1))
yl <- range(c(t.rmsep, v.rmsep))
plot(x = 1:K, y = t.rmsep, ylim = yl, axes = FALSE, frame.plot = TRUE,
     xlab = "", ylab = "RMSEP", pch = 16)
lines(x = 1:K, y = t.rmsep)
axis(2)
points(x = 1:K, y = v.rmsep, col = "red", pch = 16)
lines(x = 1:K, y = v.rmsep, col = "red")

text(x = 1:K, y = t.rmsep + 0.001, labels = vars, adj = c(0,0.5), srt = 90)
mtext("Sequentially Adding Variables", side = 1, line = 1)



N <- 20
s <- seq(-0.1,0.2,length = N)
vrmsep <- numeric(N)
for(i in 1:N){
  o[["year"]] <- c("2013" = s[i])
  pr <- score(o, newdata = dta[dta$valid, ], offset = log(dta[dta$valid, "exposure"]))

  pred <- pr$mu
  actu <- dta[dta$valid, "clm.count"]
  vrmsep[i] <- RMSEP(actu, pred)
}

# clean up the workspace
rm(bk, rmsep, m, pred, actu, i, cof)

m <- glm(clm.count ~ year + prior.claims + driver.age,
         data = dta,
         subset = train,
         family = poisson(link = "log"),
         offset = log(exposure))
p <- predict(m, newdata = dta[dta$train,], type = "response")
a <- dta$clm.count[dta$train]
RMSEP(a,p)

a<-dta$clm.count[dta$valid]
s <- c(seq(0.218,0.228,len = 50),exp(coef(m)))
r <- numeric(length(s))
for(i in 1:length(s)) {
  p <- s[i] * dta$exposure[dta$valid]
  r[i] <- RMSEP(a,p)
}
plot(x = s, y = r, type = "b")
k <- which(r == min(r))
c(s[k],r[k])



####-####+####-####+####-####+####-####+####-####+####-####+####-####+####-####+
# SEVERITY MODELING
####-####+####-####+####-####+####-####+####-####+####-####+####-####+####-####+

#
# How many records have severity above 6,000?
sum(dta$sev > 6000, na.rm = TRUE)
sum(dta$clm.incurred > 0, na.rm = TRUE)

# Average severity across entire dataset
ci <- with(dta[dta$clm.count > 0,], sum(clm.incurred))
cc <- with(dta[dta$clm.count > 0,], sum(clm.count))
sv <- round(ci / cc,2)

# Average severity across training dataset
ci <- with(dta[dta$train & dta$clm.count > 0,], sum(clm.incurred))
cc <- with(dta[dta$train & dta$clm.count > 0,], sum(clm.count))
sv <- round(ci / cc,2)

# Clean the workspace
rm(ci, cc, sv)

#
# Six point summary for severity across training dataset
with(dta[dta$train & dta$clm.count > 0, ], summary(sev))

#
# Null severity model
m <- glm(clm.incurred ~ 1, data = dta, subset = train & clm.count > 0,
         family = Gamma(link = "log"), weights = clm.count, offset = log(clm.count))
summary(m)
exp(coef(m))

# Clean the workspace
rm(m)

# Number of records with claims in the training dataset?
sum(dta$clm.count[dta$train] > 0)

#
# Number of regions in training dataset
length(unique(dta$region[dta$train]))

#
# Severity model for year and region
m <- glm(clm.incurred ~ year + region, data = dta, subset = train & clm.count > 0,
         family = Gamma(link = "log"), weights = clm.count, offset = log(clm.count))
summary(m)
anova(m, test = "Chisq")

# Clean the workspace
rm(m)

#
# Severity model for year and region.g1
m <- glm(clm.incurred ~ year + region.g1, data = dta, subset = train & clm.count > 0,
         family = Gamma(link = "log"), weights = clm.count, offset = log(clm.count))
c(aic = m$aic, dev = m$deviance)

#
# Severity model for year and drv.age.gr1
m <- glm(clm.incurred ~ year + drv.age.gr1, data = dta, subset = train & clm.count > 0,
         family = Gamma(link = "log"), weights = clm.count, offset = log(clm.count))
c(aic = m$aic, dev = m$deviance)

#
# Base level for region variable
levels(dta$region)[1]

#
# NCD.LEVEL variable
#
# Severity model for year
m <- glm(clm.incurred ~ year, data = dta, subset = train & clm.count > 0,
         family = Gamma(link = "log"), weights = clm.count, offset = log(clm.count))
summary(m)
anova(m, test = "Chisq")
with(dta[dta$train & dta$clm.count > 0,], RMSEP(clm.incurred, predict(m, type = "response")))

# adding ncd.level to model m
m1 <- update(m, . ~ . + ncd.level)
summary(m1)
anova(m1, test = "Chisq")
with(dta[dta$train & dta$clm.count > 0,], RMSEP(clm.incurred, predict(m1, type = "response")))

# Clean the workspace
rm(m,m1)

#
# DRIVER.GENDER variable
#
# Severity model for year
m <- glm(clm.incurred ~ year, data = dta, subset = train & clm.count > 0,
         family = Gamma(link = "log"), weights = clm.count, offset = log(clm.count))
with(dta[dta$train & dta$clm.count > 0,], RMSEP(clm.incurred, predict(m, type = "response")))

# adding driver.gender to model m
m1 <- update(m, . ~ . + driver.gender)
summary(m1)
anova(m1, test = "Chisq")
with(dta[dta$train & dta$clm.count > 0,], RMSEP(clm.incurred, predict(m1, type = "response")))

# Clean the workspace
rm(m,m1)

#
# MARITAL.STATUS variable
#
# Severity model for year
m <- glm(clm.incurred ~ year, data = dta, subset = train & clm.count > 0,
         family = Gamma(link = "log"), weights = clm.count, offset = log(clm.count))
m.rms <- with(dta[dta$train & dta$clm.count > 0,], RMSEP(clm.incurred, predict(m, type = "response")))

# adding marital.status to model m
m1 <- update(m, . ~ . + marital.status)
summary(m1)
anova(m1, test = "Chisq")
m1.rms <- with(dta[dta$train & dta$clm.count > 0,], RMSEP(clm.incurred, predict(m1, type = "response")))
m.rms - m1.rms

# Clean the workspace
rm(m,m1,m.rms,m1.rms)

#
# VEHICLE.VALUE variable
#
# Severity model for year
m <- glm(clm.incurred ~ year, data = dta, subset = train & clm.count > 0,
         family = Gamma(link = "log"), weights = clm.count, offset = log(clm.count))
m.rms <- with(dta[dta$train & dta$clm.count > 0,], RMSEP(clm.incurred, predict(m, type = "response")))

# adding vehicle.value to model m
m1 <- update(m, . ~ . + vehicle.value)
summary(m1)
anova(m1, test = "Chisq")
m1.rms <- with(dta[dta$train & dta$clm.count > 0,], RMSEP(clm.incurred, predict(m1, type = "response")))
m.rms - m1.rms

# Clean the workspace
rm(m,m1,m.rms,m1.rms)


#
# Severity model for year and vehicle.value
m <- glm(clm.incurred ~ year + vehicle.value, data = dta, subset = train & clm.count > 0,
         family = Gamma(link = "log"), weights = clm.count, offset = log(clm.count))
summary(m)
anova(m, test = "Chisq")

# Clean the workspace
rm(m)

#
# What proportion of severities are below 1500?
tmp <- subset(dta, subset = train & clm.count > 0, select = c("sev"))
with(tmp, round(sum(sev < 1500)/ length(sev)*100,1))

# clean the workspace
rm(tmp)

#
# Severity model for year and piecewise linear vehicle.value
m0 <- glm(clm.incurred ~ year,
          data = dta, subset = train & clm.count > 0,
          family = Gamma(link = "log"), weights = clm.count, offset = log(clm.count))

m1 <- update(m0, . ~ . + vehicle.value + veh.val.q15 + veh.val.q35)
disp <- sum(m1$weights * m1$residuals^2)/m1$df.residual
summary(m1)
anova(m1, test = "Chisq")

c('Df' = m0$df.residual - m1$df.residual, 'Resid. Dev' = m0$deviance - m1$deviance,
  'Pr(>Chi)' = pchisq((m0$deviance - m1$deviance)/disp, m0$df.residual - m1$df.residual, lower.tail = FALSE))

# Clean the workspace
rm(m0, m1, disp)

#
# REGION
#
with(dta[dta$train & dta$clm.count >0,], length(sort(unique(as.character(region)))))

m <- glm(clm.incurred ~ year + region,
         data = dta, subset = train & clm.count > 0,
         family = Gamma(link = "log"),
         weights = clm.count, offset = log(clm.count))
summary(m)

levels(dta$region)[1]

#
# Which driver ages have no reported claims?
#
cc <- with(dta[dta$train,], tapply(clm.count, drv.age, sum))
which(cc == 0)
length(which(cc == 0))

# clean up
rm(cc)


#
# Severity model for year, marital.status, and drv.age
m <- glm(clm.incurred ~ year + marital.status + drv.age, data = dta, subset = train & clm.count > 0,
         family = Gamma(link = "log"), weights = clm.count, offset = log(clm.count))
summary(m)
anova(m, test = "Chisq")

# Clean the workspace
rm(m)

#
# Amount of exposure on ages 22 and 46
ex <- with(dta[dta$train,], tapply(exposure, drv.age, sum))[-71]
pr <- ex / sum(ex)
sum(pr)
ex[c("22","46")]
round(100*pr[c("22","46")],2)

#clean up
rm(ex, pr)

#
# Severity model for year, marital.status, and drv.age.gr2
m <- glm(clm.incurred ~ year + marital.status + drv.age.gr2, data = dta, subset = train & clm.count > 0,
         family = Gamma(link = "log"), weights = clm.count, offset = log(clm.count))
summary(m)
anova(m, test = "Chisq")

# Clean the workspace
rm(m)

#
# Severity model for year, driver.gender
m <- glm(clm.incurred ~ year + driver.gender, data = dta, subset = train & clm.count > 0,
         family = Gamma(link = "log"), weights = clm.count, offset = log(clm.count))
summary(m)
anova(m, test = "Chisq")

# Clean the workspace
rm(m)


#
# Severity model for year, ncd.level, and driver.age
m <- glm(clm.incurred ~ year + ncd.level + driver.age,
         data = dta, subset = train & clm.count > 0,
         family = Gamma(link = "log"), weights = clm.count, offset = log(clm.count))
summary(m)

# plot quantile residuals against driver.age and add a smmothing spline
#     Note how the smoothing spline is not horizontal.  There is
#     information that a linear function of driver.age is not picking up
qresid <- qresiduals(m)
plot(x = dta$driver.age[dta$train & dta$clm.count > 0], y = qresid)
lines(smooth.spline(dta$driver.age[dta$train & dta$clm.count > 0], qresid), col = "red")
abline(h = 0, col = "gray")

# zoom into the residual plot
plot(x = dta$driver.age[dta$train & dta$clm.count > 0], y = qresid, ylim = c(-0.5, 0.5))
lines(smooth.spline(dta$driver.age[dta$train & dta$clm.count > 0], qresid), col = "red")
abline(h = 0, col = "gray")
abline(v = c(35, 49, 59, 67), col = "blue")

# Severity model for year, ncd.level, and piecewise linear function of driver.age
m <- glm(clm.incurred ~ year + ncd.level + driver.age +
           driver.age.q35 + driver.age.q49 + driver.age.q59,
         data = dta, subset = train & clm.count > 0,
         family = Gamma(link = "log"), weights = clm.count, offset = log(clm.count))
summary(m)

#
# Plot the new residuals against the piecewise linear function of driver.age
qresid <- qresiduals(m)
plot(x = dta$driver.age[dta$train & dta$clm.count > 0], y = qresid, ylim = c(-0.5, 0.5))
lines(smooth.spline(dta$driver.age[dta$train & dta$clm.count > 0], qresid), col = "red")
abline(h = 0, col = "gray")

# Analysis of deviance table
anova(m, test = "Chisq")

# Clean the workspace
rm(m, qresid)






####-####+####-####+####-####+####-####+####-####+####-####+####-####+####-####+
# PURE PREMIUM MODELING
####-####+####-####+####-####+####-####+####-####+####-####+####-####+####-####+

# Frequency model
fq.m <- glm(clm.count ~ year + ncd.level + drv.age.gr2 + yrs.lic + region.g1 + prior.claims,
            family = poisson(link = "log"),
            data = dta,
            subset = train,
            offset = log(exposure))
summary(fq.m)

# Severity model
sv.m <- glm(clm.incurred ~ year + marital.status + driver.gender + weight + body.code,
            family = Gamma(link = "log"),
            data = dta,
            subset = train & clm.count > 0,
            offset = log(clm.count))
summary(sv.m)

# grab the traning dataset
tr <- subset(dta, subset = train)

# Score all training records for frequency and severity
tr$fq <- score(grab.coef(fq.m), newdata = tr, offset = log(tr[,"exposure"]))$mu
tr$sv <- score(grab.coef(sv.m), newdata = tr)$mu
tr$pp <- tr$fq * tr$sv

# Ratio of pure premium to clm.incurred
sum(tr$pp) / sum(tr$clm.incurred)

# Summary statistics on freq, sev, and pure premium
summary(round(tr$fq*100,1))
summary(round(tr$sv,1))
summary(round(tr$pp,1))


# find records with very low and very high pure premiums
lo <- which(tr$pp < quantile(tr$pp, 0.01))
hi <- which(tr$pp > quantile(tr$pp, 0.99))

lo.expo <- with(tr[lo,], tapply(exposure, ncd.level, sum)) / sum(tr[lo,"exposure"])
hi.expo <- with(tr[hi,], tapply(exposure, ncd.level, sum)) / sum(tr[hi,"exposure"])
expo <- with(tr, tapply(exposure, ncd.level, sum)) / sum(tr[,"exposure"])
round(rbind(expo, lo.expo, hi.expo)*100,2)

# grab frequency coefficients and make modifications
fq.c <- grab.coef(fq.m)

fq.c[["ncd.level"]]["2"] <- -0.05

fq.c[["drv.age.gr2"]]["58-62"] <- 0.20
fq.c[["drv.age.gr2"]]["63+"] <- 0.25

fq.c[["yrs.lic"]]["4"]  <- -0.255
fq.c[["yrs.lic"]]["6"]  <- -0.500
fq.c[["yrs.lic"]]["7"]  <- -0.500
fq.c[["yrs.lic"]]["8+"] <- -0.500

#score the entire dataset for frequency and determine the adjustment on training data only
dta$efq <- score(fq.c, newdata = dta, offset = log(dta[,"exposure"]))$mu
adj <- -log(sum(dta[dta$train, "efq"])/ sum(dta[dta$train, "clm.count"]))
fq.c[["Base"]] <- fq.c[["Base"]] + adj

# grab severity coefficients
sv.c <- grab.coef(sv.m)

# score the entire dataset for frequency, severity, and calculate pure premium
dta$efq <- score(fq.c, newdata = dta, offset = log(dta[,"exposure"]))$mu
dta$esv <- score(sv.c, newdata = dta)$mu
dta$epp <- dta$efq * dta$esv

#summary statistics for the entire dataset
with(dta,
     round(rbind(summary(efq * 100),
                 summary(esv),
                 summary(epp)),1))

#summary statistics for the training dataset
with(dta[dta$train,],
     round(rbind(summary(efq * 100),
                 summary(esv),
                 summary(epp)),1))

# sum of pure premium and claim incurred for training dataset
c(pp = sum(dta[dta$train,"epp"]), ci = sum(dta[dta$train, "clm.incurred"]))



p <- order(dta[dta$train, "epp"])
plot(x = cumsum(dta[dta$train, "epp"][p]), y = cumsum(dta[dta$train, "clm.incurred"][p]), pch = ".")
abline(a = 0, b = 1, col = "red")
plot(x = 1:24495 , y = cumsum(dta[dta$train, "epp"][p])-cumsum(dta[dta$train, "clm.incurred"][p]), pch = ".")
abline(h= 0, col = "gray")


# Frequency model with year, ncd.level, and driver.gender
fq.m <- glm(clm.count ~ year + ncd.level + driver.gender,
            family = poisson(link = "log"),
            data = dta,
            subset = train,
            offset = log(exposure))
summary(fq.m)

# create a restricted ncd.level variable
disc <- c("1" = 0, "2" = 0.10, "3" = 0.20, "4" = 0.25, "5" = 0.30, "6" = 0.35)
ncd.level.r <- log(1 - disc)

# create the offset
dta$ofst <- ncd.level.r[dta$ncd.level] + log(dta$exposure)

# fit a restricted model
fq.mr <- glm(clm.count ~ year + driver.gender,
             family = poisson(link = "log"),
             data = dta,
             subset = train,
             offset = ofst)
summary(fq.mr)

####-####+####-####+####-####+####-####+####-####+####-####+####-####+####-####+
# VALIDATION
####-####+####-####+####-####+####-####+####-####+####-####+####-####+####-####+

# Check total actual claim counts vs. predicted claims on validation dataset
a.cc <- c("train" = sum(dta[dta$train, "clm.count"]),
          "valid" = sum(dta[dta$valid, "clm.count"]))
p.fq <- c("train" = sum(dta[dta$train, "efq"]),
          "valid" = sum(dta[dta$valid, "efq"]))
a.ci <- c("train" = sum(dta[dta$train, "clm.incurred"]),
          "valid" = sum(dta[dta$valid, "clm.incurred"]))
p.pp <- c("train" = sum(dta[dta$train, "epp"]),
          "valid" = sum(dta[dta$valid, "epp"]))

#clean your workspace
rm(a.cc, p.fq, a.ci, p.pp)
#
# Check across driver.gender
#
a.cc <- c("train" = with(dta[dta$train,], tapply(clm.count, driver.gender, sum)),
          "valid" = with(dta[dta$valid,], tapply(clm.count, driver.gender, sum)))
p.fq <- c("train" = with(dta[dta$train,], tapply(efq, driver.gender, sum)),
          "valid" = with(dta[dta$valid,], tapply(efq, driver.gender, sum)))
a.ci <- c("train" = with(dta[dta$train,], tapply(clm.incurred, driver.gender, sum)),
          "valid" = with(dta[dta$valid,], tapply(clm.incurred, driver.gender, sum)))
p.pp <- c("train" = with(dta[dta$train,], tapply(epp, driver.gender, sum)),
          "valid" = with(dta[dta$valid,], tapply(epp, driver.gender, sum)))

ex <- with(dta, tapply(exposure, list(train,driver.gender), sum))

ty.pt = c(t.all = sum(dta$epp[dta$train])/sum(dta$clm.incurred[dta$train]),
  t.m = sum(dta$epp[dta$train & dta$driver.gender == "Male"])/sum(dta$clm.incurred[dta$train & dta$driver.gender == "Male"]),
  t.f = sum(dta$epp[dta$train & dta$driver.gender == "Female"])/sum(dta$clm.incurred[dta$train & dta$driver.gender == "Female"]))

vy.pt = c(t.all = sum(dta$epp[dta$valid])/sum(dta$clm.incurred[dta$valid]),
          t.m = sum(dta$epp[dta$valid & dta$driver.gender == "Male"])/sum(dta$clm.incurred[dta$valid & dta$driver.gender == "Male"]),
          t.f = sum(dta$epp[dta$valid & dta$driver.gender == "Female"])/sum(dta$clm.incurred[dta$valid & dta$driver.gender == "Female"]))

plot(x = c(1,2,2,3,4,4), y = c(ty.pt,vy.pt))


a.ci <- c("train" = with(dta[dta$train,], tapply(clm.incurred, ncd.level, sum)),
          "valid" = with(dta[dta$valid,], tapply(clm.incurred, ncd.level, sum)))
p.pp <- c("train" = with(dta[dta$train,], tapply(epp, ncd.level, sum)),
          "valid" = with(dta[dta$valid,], tapply(epp, ncd.level, sum)))


bks <- quantile(dta$epp, probs = seq(0,1,length= 11))
#bks[1] <- 0
tmp <- cut(dta$epp, bks, include.lowest = TRUE)


#
# Gini index calculations
#

# on training data
o <- with(dta[dta$train,], order(epp))
x <- with(dta[dta$train,], cumsum(exposure[o])/sum(exposure))
y <- with(dta[dta$train,], cumsum(clm.incurred[o])/sum(clm.incurred))

dx <- x[-1] - x[-length(x)]
h <- (y[-1] + y[-length(y)])/2
gini <- 2 * (0.5 - sum(h * dx))

op <- par(mar = c(5,4,0,0)+0.2)
plot(x = x, y = y, type = "n")
lines(x = c(0,1), y = c(0,1), col = "grey")
lines(x = x, y = y, col = "red")

# on validation data
o <- with(dta[dta$valid,], order(epp))
x <- with(dta[dta$valid,], cumsum(exposure[o])/sum(exposure))
y <- with(dta[dta$valid,], cumsum(clm.incurred[o])/sum(clm.incurred))

dx <- x[-1] - x[-length(x)]
h <- (y[-1] + y[-length(y)])/2
gini2 <- 2 * (0.5 - sum(h * dx))
lines(x = x, y = y, col = "red")

par(op)

#
#
#
op <- par(mar = c(5,4,0,0)+0.2)
plot(x = 1, y = 1, type = "n", xlim = c(0,1), ylim = c(0,1),
     xlab = "Cumulative Relative Exposure", ylab = "Cumulative Relative Losses")
lines(x = c(0,1), y = c(0,1))
x <- seq(0,1,length=100)
y <- x^2
lines(x=x,y=y,col = "red")

par(op)

#
# Predicted vs. observed plots
o <- order(dta$epp[dta$valid])
plot(x = (dta$epp[dta$valid])[o],
     y = (dta$clm.incurred[dta$valid])[o])



####-####+####-####+####-####+####-####+####-####+####-####+####-####+####-####+
# CONCLUSIONS
####-####+####-####+####-####+####-####+####-####+####-####+####-####+####-####+


####-####+####-####+####-####+####-####+####-####+####-####+####-####+####-####+
#
####-####+####-####+####-####+####-####+####-####+####-####+####-####+####-####+
