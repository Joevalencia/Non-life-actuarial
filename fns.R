#
# Auto GLM chapter for
# Predictive Modeling Applications in Actuarial Science
# Volume II
#
# Ernesto Schirmacher
#

re.lev <- function(x, ref, expo = NULL, ...) {
  if(is.null(expo))
    return(relevel(x, ref, ...))
  
  p <- order(tapply(expo, x, sum), decreasing = TRUE)
  lv <- levels(x)[p[1]]
  rlv <- relevel(x, ref = lv)
  return(rlv)
}

fq.tbl <- function(var) {
  yr.idx <- dta$year == "2010"
  expo <- tapply(dta[yr.idx, "exposure"], dta[yr.idx, var], sum)
  clms <- tapply(dta[yr.idx, "clm.count"], dta[yr.idx, var], sum)
  tbl <- c(expo, clms, clms/expo)
  
  yr.idx <- dta$year == "2011"
  expo <- tapply(dta[yr.idx, "exposure"], dta[yr.idx, var], sum)
  clms <- tapply(dta[yr.idx, "clm.count"], dta[yr.idx, var], sum)
  tbl <- rbind(tbl, c(expo, clms, clms/expo))
  
  yr.idx <- dta$year == "2012"
  expo <- tapply(dta[yr.idx, "exposure"], dta[yr.idx, var], sum)
  clms <- tapply(dta[yr.idx, "clm.count"], dta[yr.idx, var], sum)
  tbl <- rbind(tbl, c(expo, clms, clms/expo))
  
  tbl <- rbind(tbl, apply(tbl, 2, sum))
  tbl[4,5:6] <- tbl[4,3:4]/tbl[4,1:2]
  
  dimnames(tbl)[[1]] <- c("2010", "2011", "2012", "Total")
  
  return(tbl)
}

RMSEP <- function(actual, predicted) {
  sq.d <- (actual - predicted)^2
  ans <- sqrt(mean(sq.d))
  return(ans)
}

score <- function(o, newdata, offset = 0) {

  if("Base" %in% names(o)) {
    base <- o[["Base"]]
    o <- o[-which(names(o) == "Base")]
  }
  else
    base <- 0
  
  N <- dim(newdata)[1]
  idx <- match(names(o), colnames(newdata))
  f <- data.frame(rep(base, N), offset)

  for(i in 1:length(idx)) {
    if(class(newdata[,idx[i]])[1] == "factor")
      f <- cbind(f, o[[i]][as.character(newdata[,idx[i]])])
    else
      f <- cbind(f, o[[i]] * newdata[,idx[i]])
  }
  names(f) <- c("f.base", "f.offset", paste("f", names(o), sep = "."))

  eta <- apply(f, 1, sum)
  mu <- exp(eta)
  ans <- list(dt =  newdata[,idx], fr = f, eta = eta, mu = mu)

  return(ans)
}

grab.coef <- function(m) {
  ans <- list()
  cf <- coef(m)
  
  vars <- attr(m$terms, "term.labels")
  dtClasses <- attr(m$terms, "dataClasses")
  is.fac <- dtClasses[vars] == "factor"
  xlev <- m$xlevels
  
  fac.idx <- which(is.fac)
  for(i in fac.idx){
    v <- vars[i]
    val <- c(0, cf[grep(v, names(cf))])
    names(val) <- xlev[[v]]
    ans[[v]] <- val
  }
  
  num.idx <- which(!is.fac)
  for(i in num.idx) {
    v <- vars[i]
    val <- cf[grep(v, names(cf))]
    ans[[v]] <- val
  }
  
  if(length(grep("Intercept", names(cf))) > 0) {
    ans[["Base"]] <- cf["(Intercept)"]
  }
  return(ans)
}

num.summary <- function(x)
  c("Mean"=mean(x), "Median"=median(x), "Std.Dev"=sd(x),
    "Min"=min(x), "Max"=max(x),"NA"=sum(is.na(x)))

cat.summary <- function(x) {
  n <- length(levels(x))
  if(n < 8)
    N <- n
  else
    N <- 8
  n.lev <- as.character(length(levels(x)))
  b.lev <- levels(x)[1]
  m.lev <- names(sort(table(x), decreasing = TRUE))[1]
  s.lev <- paste(sort(levels(x))[1:N], collapse = ", ")
  ans <- c("No. Lev"=n.lev, "Base"=b.lev, "Most"=m.lev, "Sample"=s.lev)
  return(ans)
}
