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

get_D_imm = function(k, D_imm_r, Farms) {
    if (k == 1) {
        D_imm = rep(D_imm_r[1], nrow(Farms))
    }
    else {
        if (k == 2) {
            D_imm = rep(D_imm_r[2], nrow(Farms))
        }
        else {
            D_imm = runif(nrow(Farms), min = D_imm_r[1], max = D_imm_r[2])
        }
    }
    return(D_imm)
}

get_D = function(k, D_s_r, D_p_r, Farms) {
    if (k == 1) {
        D = ifelse(Farms$sow == 1, D_s_r[1], D_p_r[1])
    }
    else {
        if (k == 2) {
            D = ifelse(Farms$sow == 1, D_s_r[2], D_p_r[2])
        }
        else {
            D = ifelse(Farms$sow == 1, runif(1, min = D_s_r[1], 
                max = D_s_r[2]), runif(1, min = D_p_r[1], max = D_p_r[2]))
        }
    }
    return(D)
}

get_v = function(k, v_s_r, v_p_r, Farms) {
    if (k == 1) {
        v = ifelse(Farms$sow == 1, v_s_r[1], v_p_r[1])
    }
    else {
        if (k == 2) {
            v = ifelse(Farms$sow == 1, v_s_r[2], v_p_r[2])
        }
        else {
            v = ifelse(Farms$sow == 1, runif(1, min = v_s_r[1], 
                max = v_s_r[2]), runif(1, min = v_p_r[1], max = v_p_r[2]))
        }
    }
    return(v)
}

get_D25 = function(k, D_s_r, D_p_r, Dv_s_r, Farms) {
    if (k == 1) {
        D25 = ifelse(Farms$sowV25 == 1, Dv_s_r[1], ifelse(Farms$sow == 
            1, D_s_r[1], D_p_r[1]))
    }
    else {
        if (k == 2) {
            D25 = ifelse(Farms$sowV25 == 1, Dv_s_r[2], ifelse(Farms$sow == 
                1, D_s_r[2], D_p_r[2]))
        }
        else {
            D25 = ifelse(Farms$sowV25 == 1, runif(1, min = Dv_s_r[1], 
                max = Dv_s_r[2]), ifelse(Farms$sow == 1, runif(1, 
                min = D_s_r[1], max = D_s_r[2]), runif(1, min = D_p_r[1], 
                max = D_p_r[2])))
        }
    }
    return(D25)
}

get_D50 = function(k, D_s_r, D_p_r, Dv_s_r, Farms) {
    if (k == 1) {
        D50 = ifelse(Farms$sowV50 == 1, Dv_s_r[1], ifelse(Farms$sow == 
            1, D_s_r[1], D_p_r[1]))
    }
    else {
        if (k == 2) {
            D50 = ifelse(Farms$sowV50 == 1, Dv_s_r[2], ifelse(Farms$sow == 
                1, D_s_r[2], D_p_r[2]))
        }
        else {
            D50 = ifelse(Farms$sowV50 == 1, runif(1, min = Dv_s_r[1], 
                max = Dv_s_r[2]), ifelse(Farms$sow == 1, runif(1, 
                min = D_s_r[1], max = D_s_r[2]), runif(1, min = D_p_r[1], 
                max = D_p_r[2])))
        }
    }
    return(D50)
}

get_D75 = function(k, D_s_r, D_p_r, Dv_s_r, Farms) {
    if (k == 1) {
        D75 = ifelse(Farms$sowV75 == 1, Dv_s_r[1], ifelse(Farms$sow == 
            1, D_s_r[1], D_p_r[1]))
    }
    else {
        if (k == 2) {
            D75 = ifelse(Farms$sowV75 == 1, Dv_s_r[2], ifelse(Farms$sow == 
                1, D_s_r[2], D_p_r[2]))
        }
        else {
            D75 = ifelse(Farms$sowV75 == 1, runif(1, min = Dv_s_r[1], 
                max = Dv_s_r[2]), ifelse(Farms$sow == 1, runif(1, 
                min = D_s_r[1], max = D_s_r[2]), runif(1, min = D_p_r[1], 
                max = D_p_r[2])))
        }
    }
    return(D75)
}

get_D100 = function(k, D_s_r, D_p_r, Dv_s_r, Farms) {
    if (k == 1) {
        D100 = ifelse(Farms$sowV100 == 1, Dv_s_r[1], ifelse(Farms$sow == 
            1, D_s_r[1], D_p_r[1]))
    }
    else {
        if (k == 2) {
            D100 = ifelse(Farms$sowV100 == 1, Dv_s_r[2], ifelse(Farms$sow == 
                1, D_s_r[2], D_p_r[2]))
        }
        else {
            D100 = ifelse(Farms$sowV100 == 1, runif(1, min = Dv_s_r[1], 
                max = Dv_s_r[2]), ifelse(Farms$sow == 1, runif(1, 
                min = D_s_r[1], max = D_s_r[2]), runif(1, min = D_p_r[1], 
                max = D_p_r[2])))
        }
    }
    return(D100)
}
