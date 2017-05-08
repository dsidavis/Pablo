PrevValues = list(beta = NULL, stocks = NULL)
model <- function(time, stocks, parms) {
    st = split(stocks, rep(1:3, each = length(stocks)/3))
    names(st) = c("S", "I", "R")

if(is.null(PrevValues$beta)) {
 PrevValues[1:2] <<- list(parms$beta, stocks)
} else {
  stopifnot(all.equal(PrevValues[[2]], stocks))
  stopifnot(all.equal(PrevValues[[1]], parms$beta))
}


    lambda = parms$beta %*% st$I
    IR <- lambda * st$S
    RR <- st$I/parms$delays
    MR <- st$I * parms$mort
    SR <- st$R/parms$returns
    dS_dt <- SR - IR
    dI_dt <- IR - RR - MR
    dR_dt <- RR - SR

    list(c(dS_dt, dI_dt, dR_dt))
}
