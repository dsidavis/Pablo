model <- function(time, stocks, parms) {
    n = length(stocks)/3
    S = stocks[1:n]
    I = stocks[(n+1):(2*n)]
    R = stocks[(2*n+1):(3*n)]
#    st = split(stocks, rep(1:3, each = length(stocks)/3))
#    names(st) = c("S", "I", "R")

    lambda = parms$beta %*% I
    IR <- lambda * S
    RR <- I/parms$delays
    MR <- I * parms$mort
    SR <- R/parms$returns
    dS_dt <- SR - IR
    dI_dt <- IR - RR - MR
    dR_dt <- RR - SR

    list(c(dS_dt, dI_dt, dR_dt))
}
