model <- function(time, stocks, parms) {
    st = split(stocks, rep(1:3, each = length(stocks)/3))
    names(st) = c("S", "I", "R")
# But really we know they are in positions
#   1:n1, n1+1:n2, n2+1:n3 and so we don't need substring.
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
