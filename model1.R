model <- function(time, stocks, parms) {
    S <- stocks[grep("S", names(stocks))]
    I <- stocks[grep("I", names(stocks))]
    R <- stocks[grep("R", names(stocks))]
    with(c(list(S = S, I = I, R = R), parms), {
        lambda = beta %*% I
        IR <- lambda * S
        RR <- I/delays
        MR <- I * mort
        SR <- R/returns
        dS_dt <- SR - IR
        dI_dt <- IR - RR - MR
        dR_dt <- RR - SR
        return(list(c(dS_dt, dI_dt, dR_dt)))
    })
}

