 #Add colClasses for these to specify the types.
matrix <- read.csv("matrix_to_D.csv" )
Farms <- read.csv("Farms_to_D.csv" )


#library(deSolve)
library(ggplot2)
library(dplyr)
library(plyr)
library(reshape2)
library(base)
library(MASS)
library(FME)
library(lhs)
library(zoo)

library(igraph)
library(EpiContactTrace)
library(RColorBrewer)
library(maps)
library(doBy)

library(grid)
library(gridExtra)

library(nlme)
library(multcomp)


# Rho as matrices ------------------------------------------------

rho1 <- acast(matrix, destination ~ origin, value.var = 'rho1') 
rho2 <- acast(matrix, destination ~ origin, value.var = 'rho2')
rho3 <- acast(matrix, destination ~ origin, value.var = 'rho3')
rho4 <- acast(matrix, destination ~ origin, value.var = 'rho4')
rho5 <- acast(matrix, destination ~ origin, value.var = 'rho5')
rho6 <- acast(matrix, destination ~ origin, value.var = 'rho6')
rho7 <- acast(matrix, destination ~ origin, value.var = 'rho7')
rho8 <- acast(matrix, destination ~ origin, value.var = 'rho8')
rho9 <- acast(matrix, destination ~ origin, value.var = 'rho9')

diag(rho1) <- 1;diag(rho2) <- 1;diag(rho3) <- 1;diag(rho4) <- 1;diag(rho5) <- 1;diag(rho6) <- 1;diag(rho7) <- 1;diag(rho8) <- 1;diag(rho9) <- 1 # Fill diagonal with 1, due to rho_ii=1

# PARAMETERS  -------------------------------------------------------------
# Besic reproduction number (R0)
# Jeong (2014)
R0_l_s <- 0.14    # (R0 of low virulence in sows. min)
R0_m_s <- 3       # (R0 of medium virulence in sows. mean)
R0_h_s <- 3.22    # (R0 of high virulence in sows. max)
R0_l_p <- 7.26     # (R0 of low virulence in piglets. min)
R0_m_p <- 9.26     # (R0 of medium virulence in piglets. mean)
R0_h_p <- 13.13    # (R0 of high virulence in piglets. max)

# Charpin (2012)
R0_l_p2 <- 1.8     # (R0 of low virulence in piglets. min)
R0_m_p2 <- 2.6     # (R0 of medium virulence in piglets. mean)
R0_h_p2 <- 3.3    # (R0 of high virulence in piglets. max)

# Days of infection (D)
# Jeong (2014)
D_us <- 4 # Days unvaccinated sows
D_vs <- 2.8 # Days vaccinated sows
D_up <- 8 # Days unvaccinated pigs

# Linhares (2012)
D_us2 <- 10 # Days unvaccinated sows
D_vs2 <- 6.4 # Days vaccinated sows

# Weeks of active immunity
D_imm <- 36 # (min=26, mean=36, max=52)

# Beta (without using vaccine) in sow farms
beta_l_s <- R0_l_s/D_us
beta_m_s <- R0_m_s/D_us
beta_h_s <- R0_h_s/D_us

# Beta (using vaccine) in sow farms
betaV_l_s <- R0_l_s/D_vs
betaV_m_s <- R0_m_s/D_vs
betaV_h_s <- R0_h_s/D_vs

# Beta (without using vaccine) in pigs farms
beta_l_p <- R0_l_p2/D_up
beta_m_p <- R0_m_p2/D_up
beta_h_p <- R0_h_p2/D_up

# Mortality in sow farms
v_s = 0.001 #(min=0, max=0.002)

# Mortality in pigs farms
v_p = 0.07 #(min=0.02, max=0.18)

# INFECT ONE FARM AT TIME=0 -----------------------------------------------
# Farms without vaccination or other action to control PRRS (low virulance)
# Assign posistive animals (~1% of the population infeceted) in one farm
set.seed(0)
Farms$X <- NA; Farms$Y <- NA; Farms$Z <- 0 # Create new columns for number of susceptible (X), infected (Y) and recovered (Z) within each farm

which(Farms$within3k>5 & Farms$links>10) # Criteria to seed the first infected (10% of inventory)
#which(Farms$type=="Fa" & Farms$county=="Stevens") 
Farms$Y[which(Farms$within3k>5 & Farms$links>10)] <- sample(c(rep(0,nrow(Farms[which(Farms$within3k>5 & 
                                                                                       Farms$links>10),])-1),1),
                                                            nrow(Farms[which(Farms$within3k>5 & Farms$links>10),])) # Fill with 0 and 1 farms under criteria of above

#Farms$Y[which(Farms$type=="Fa" & Farms$county=="Stevens")] <- sample(c(rep(0,nrow(Farms[which(Farms$type=="Fa" & Farms$county=="Stevens"),])-1),1),nrow(Farms[which(Farms$type=="Fa" & Farms$county=="Stevens"),])) # Fill with 0 and 1 farms under criteria of above
Farms$Y[is.na(Farms$Y)] <- 0 # Fill with 0 other farms under infected vector

Farms$Y <- ifelse(Farms$Y==1 , Farms$inventory*0.01, 0) # selected farm will have 1% of inventory infected
Farms$X <- ifelse(Farms$Y == 0, Farms$inventory, Farms$inventory-Farms$Y) 

Farms$S <- NA; Farms$I <- NA; Farms$R <- NA # Create new columns with proportion of susceptible (S), infected (I) and recovered (R) within each farm
Farms$S <- Farms$X/Farms$inventory
Farms$I <- Farms$Y/Farms$inventory
Farms$R <- Farms$Z/Farms$inventory

summary(Farms)
# END Farms 

# BETA parameter ----------------------------------------------------------

Farms$beta_l <- ifelse(Farms$sow==1, beta_l_s,beta_l_p)
Farms$beta_m <- ifelse(Farms$sow==1, beta_m_s,beta_m_p)
Farms$beta_h <- ifelse(Farms$sow==1, beta_h_s,beta_h_p)

Farms$beta_l_N <- Farms$beta_l/Farms$inventory #Beta divided by N from each farm
Farms$beta_m_N <- Farms$beta_m/Farms$inventory
Farms$beta_h_N <- Farms$beta_h/Farms$inventory

# List all possibilities of BETA and RHO ----------------------------------

beta_list = list(Farms$beta_l_N, Farms$beta_m_N, Farms$beta_h_N)
rho_list = list(rho1, rho2, rho3, rho4, rho5, rho6, rho7, rho8, rho9)

# Function SIR dissagregated model ----------------------------------------

start <- 0; finish <- 26; step <- .05 # When I increase the step time, increase the speed of my analyses. Unfortunately, if I increase step beyond than 0.05, the number of Susceptibles, Infected and Recovered goes to + or - Inf. Or I get NA values.
time <- seq(start, finish, step) #frame time period (1/2 year)

model <- function(time, stocks, parms){
  
  #states <- Farms[,c("X", "Y", "Z")]
  
  S <- stocks[grep("S",names(stocks))]
  I <- stocks[grep("I",names(stocks))]
  R <- stocks[grep("R",names(stocks))]
  
  with(c(list("S"=S,"I"=I,"R"=R), parms), {
    
    lambda = beta %*% I
    
    IR <- lambda*S
    RR <- I/delays
    MR <- I*mort
    SR <- R/returns
    
    dS_dt <- SR - IR
    dI_dt <- IR - RR - MR
    dR_dt <- RR - SR
    
    return(list(c(dS_dt, dI_dt, dR_dt)))
    
  })
}

# CONTROL STRATEGIES 1 ----------------------------------------------------
###########################################################################
# Set number of simulation  (to save time I have reduced to 10 instead 20)
nsims = 10 # (27 * 10 = 270 simulations total)

# Range of parameters (days of immunity [D_imm], increase in mortality [v] and days of infection [D])
D_imm_r = c(26, 52)
v_s_r = c(0, 0.002)
v_p_r = c(0.02, 0.18)
D_s_r = c(1, 6)
D_p_r = c(4, 12)

get_D_imm = function(k, D_imm_r, Farms){
  if(k == 1){
    D_imm = rep(D_imm_r[1], nrow(Farms))
  } else {
    if(k == 2) {
      D_imm = rep(D_imm_r[2], nrow(Farms))
    } else {
      D_imm = runif(nrow(Farms), min = D_imm_r[1], max = D_imm_r[2])
    }
  }
  return(D_imm)
}
D_imm = lapply(1:nsims, get_D_imm, D_imm_r, Farms)

get_D = function(k, D_s_r, D_p_r, Farms){
  if(k == 1){
    D = ifelse(Farms$sow == 1, D_s_r[1], D_p_r[1])
  } else {
    if(k == 2) {
      D = ifelse(Farms$sow == 1, D_s_r[2], D_p_r[2])
    } else {
      D = ifelse(Farms$sow == 1, 
                 runif(1, min = D_s_r[1], max = D_s_r[2]),
                 runif(1, min = D_p_r[1], max = D_p_r[2]))
    }
  }
  return(D)
}
D = lapply(1:nsims, get_D, D_s_r, D_p_r, Farms)

get_v = function(k, v_s_r, v_p_r, Farms){
  if(k == 1){
    v = ifelse(Farms$sow == 1, v_s_r[1], v_p_r[1])
  } else {
    if(k == 2) {
      v = ifelse(Farms$sow == 1, v_s_r[2], v_p_r[2])
    } else {
      v = ifelse(Farms$sow == 1, 
                 runif(1, min = v_s_r[1], max = v_s_r[2]),
                 runif(1, min = v_p_r[1], max = v_p_r[2]))
    }
  }
  return(v)
}
v = lapply(1:nsims, get_v, v_s_r, v_p_r, Farms)

set.seed(1) # I only can impose any strategy of control in some farms (sow farms), so I am teasting the disease dyamics by implementing 2 control strategies in those farms with different level of regional coverage (25%, 50%, 75% and 100%)
Farms3 <- subset(Farms, sow == 1)[sample(nrow(subset(Farms, sow == 1)), round(sum(Farms$sow)*0.25)), ]
Farms3$sowV25 <- NA 
Farms3$sowV25 <- 1
Farms <- merge(Farms, subset(Farms3[ , c("id", "sowV25")]), by.x = "id", by.y = "id", all.x = T)
Farms$sowV25[is.na(Farms$sowV25)] <- 0
set.seed(1)
Farms3 <- subset(Farms, sow == 1)[sample(nrow(subset(Farms, sow == 1)), round(sum(Farms$sow)*0.5)), ]
Farms3$sowV50 <- NA 
Farms3$sowV50 <- 1
Farms <- merge(Farms, subset(Farms3[ , c("id", "sowV50")]), by.x = "id", by.y = "id", all.x = T)
Farms$sowV50[is.na(Farms$sowV50)] <- 0
set.seed(1)
Farms3 <- subset(Farms, sow == 1)[sample(nrow(subset(Farms, sow == 1)), round(sum(Farms$sow)*0.75)), ]
Farms3$sowV75 <- NA 
Farms3$sowV75 <- 1
Farms <- merge(Farms, subset(Farms3[ , c("id", "sowV75")]), by.x = "id", by.y = "id", all.x = T)
Farms$sowV75[is.na(Farms$sowV75)] <- 0
set.seed(1)
Farms3 <- subset(Farms, sow == 1)[sample(nrow(subset(Farms, sow == 1)), round(sum(Farms$sow)*1)), ]
Farms3$sowV100<- NA 
Farms3$sowV100 <- 1
Farms <- merge(Farms, subset(Farms3[ , c("id", "sowV100")]), by.x = "id", by.y = "id", all.x = T)
Farms$sowV100[is.na(Farms$sowV100)] <- 0

# Vaccine efficacy / Beta*(1-E)
E = c(0.1, 0.5)

# 25% coverage
Farms$beta_l_NV25 <- ifelse(Farms$sowV25 ==1, Farms$beta_l_N*(1-E[2]), Farms$beta_l_N)
Farms$beta_m_NV25 <- ifelse(Farms$sowV25 ==1, Farms$beta_m_N*(1-E[2]), Farms$beta_m_N)
Farms$beta_h_NV25 <- ifelse(Farms$sowV25 ==1, Farms$beta_h_N*(1-E[2]), Farms$beta_h_N)
beta_listV25 = list(Farms$beta_l_NV25, Farms$beta_m_NV25, Farms$beta_h_NV25)

# 50% coverage
Farms$beta_l_NV50 <- ifelse(Farms$sowV50 ==1, Farms$beta_l_N*(1-E[2]), Farms$beta_l_N)
Farms$beta_m_NV50 <- ifelse(Farms$sowV50 ==1, Farms$beta_m_N*(1-E[2]), Farms$beta_m_N)
Farms$beta_h_NV50 <- ifelse(Farms$sowV50 ==1, Farms$beta_h_N*(1-E[2]), Farms$beta_h_N)
beta_listV50 = list(Farms$beta_l_NV50, Farms$beta_m_NV50, Farms$beta_h_NV50)

# 75%
Farms$beta_l_NV75 <- ifelse(Farms$sowV75 ==1, Farms$beta_l_N*(1-E[2]), Farms$beta_l_N)
Farms$beta_m_NV75 <- ifelse(Farms$sowV75 ==1, Farms$beta_m_N*(1-E[2]), Farms$beta_m_N)
Farms$beta_h_NV75 <- ifelse(Farms$sowV75 ==1, Farms$beta_h_N*(1-E[2]), Farms$beta_h_N)
beta_listV75 = list(Farms$beta_l_NV75, Farms$beta_m_NV75, Farms$beta_h_NV75)

# 100% coverage
Farms$beta_l_NV100 <- ifelse(Farms$sowV100 ==1, Farms$beta_l_N*(1-E[2]), Farms$beta_l_N)
Farms$beta_m_NV100 <- ifelse(Farms$sowV100 ==1, Farms$beta_m_N*(1-E[2]), Farms$beta_m_N)
Farms$beta_h_NV100 <- ifelse(Farms$sowV100 ==1, Farms$beta_h_N*(1-E[2]), Farms$beta_h_N)
beta_listV100 = list(Farms$beta_l_NV100, Farms$beta_m_NV100, Farms$beta_h_NV100)

# Set Rhos for filters
F_s = c(.6,.9)
B_s = c(.6)
matrix <- merge(matrix, Farms[,c("id","sowV25","sowV50","sowV75","sowV100")], 
                by.x = "destination", by.y = "id", all.x = T )
matrix <- merge(matrix, Farms[,c("id","sowV25","sowV50","sowV75","sowV100")], 
                by.x = "origin", by.y = "id", all.x = T )
matrix$sowV25 = matrix$sowV25.x+matrix$sowV25.y 
matrix$sowV50 = matrix$sowV50.x+matrix$sowV50.y 
matrix$sowV75 = matrix$sowV75.x+matrix$sowV75.y 
matrix$sowV100 = matrix$sowV100.x+matrix$sowV100.y 
matrix[,c("sowV25.x","sowV25.y","sowV50.x","sowV50.y","sowV75.x","sowV75.y","sowV100.x","sowV100.y")] <- list(NULL) 

# Reduce K_ij by X% ramndlnly using filtering at a given protection
# 80% protection
matrix$k_ij_F80_25 <- ifelse(matrix$sowV25 >0, matrix$k_ij*(1-F_s[2]),matrix$k_ij)
matrix$k_ij_F80_50 <- ifelse(matrix$sowV50 >0, matrix$k_ij*(1-F_s[2]),matrix$k_ij)
matrix$k_ij_F80_75 <- ifelse(matrix$sowV75 >0, matrix$k_ij*(1-F_s[2]),matrix$k_ij)
matrix$k_ij_F80_100 <- ifelse(matrix$sowV100 >0, matrix$k_ij*(1-F_s[2]),matrix$k_ij)

matrix$k_ij_up_F80_25 <- ifelse(matrix$sowV25 >0, matrix$k_ij_up*(1-F_s[2]),matrix$k_ij_up)
matrix$k_ij_up_F80_50 <- ifelse(matrix$sowV50 >0, matrix$k_ij_up*(1-F_s[2]),matrix$k_ij_up)
matrix$k_ij_up_F80_75 <- ifelse(matrix$sowV75 >0, matrix$k_ij_up*(1-F_s[2]),matrix$k_ij_up)
matrix$k_ij_up_F80_100 <- ifelse(matrix$sowV100 >0, matrix$k_ij_up*(1-F_s[2]),matrix$k_ij_up)

matrix$k_ij_low_F80_25 <- ifelse(matrix$sowV25 >0, matrix$k_ij_low*(1-F_s[2]),matrix$k_ij_low)
matrix$k_ij_low_F80_50 <- ifelse(matrix$sowV50 >0, matrix$k_ij_low*(1-F_s[2]),matrix$k_ij_low)
matrix$k_ij_low_F80_75 <- ifelse(matrix$sowV75 >0, matrix$k_ij_low*(1-F_s[2]),matrix$k_ij_low)
matrix$k_ij_low_F80_100 <- ifelse(matrix$sowV100 >0, matrix$k_ij_low*(1-F_s[2]),matrix$k_ij_low)

# 40% protection
matrix$k_ij_F40_25 <- ifelse(matrix$sowV25 >0, matrix$k_ij*(1-F_s[1]),matrix$k_ij)
matrix$k_ij_F40_50 <- ifelse(matrix$sowV50 >0, matrix$k_ij*(1-F_s[1]),matrix$k_ij)
matrix$k_ij_F40_75 <- ifelse(matrix$sowV75 >0, matrix$k_ij*(1-F_s[1]),matrix$k_ij)
matrix$k_ij_F40_100 <- ifelse(matrix$sowV100 >0, matrix$k_ij*(1-F_s[1]),matrix$k_ij)

matrix$k_ij_up_F40_25 <- ifelse(matrix$sowV25 >0, matrix$k_ij_up*(1-F_s[1]),matrix$k_ij_up)
matrix$k_ij_up_F40_50 <- ifelse(matrix$sowV50 >0, matrix$k_ij_up*(1-F_s[1]),matrix$k_ij_up)
matrix$k_ij_up_F40_75 <- ifelse(matrix$sowV75 >0, matrix$k_ij_up*(1-F_s[1]),matrix$k_ij_up)
matrix$k_ij_up_F40_100 <- ifelse(matrix$sowV100 >0, matrix$k_ij_up*(1-F_s[1]),matrix$k_ij_up)

matrix$k_ij_low_F40_25 <- ifelse(matrix$sowV25 >0, matrix$k_ij_low*(1-F_s[1]),matrix$k_ij_low)
matrix$k_ij_low_F40_50 <- ifelse(matrix$sowV50 >0, matrix$k_ij_low*(1-F_s[1]),matrix$k_ij_low)
matrix$k_ij_low_F40_75 <- ifelse(matrix$sowV75 >0, matrix$k_ij_low*(1-F_s[1]),matrix$k_ij_low)
matrix$k_ij_low_F40_100 <- ifelse(matrix$sowV100 >0, matrix$k_ij_low*(1-F_s[1]),matrix$k_ij_low)

# Increase in biosecurity 
matrix$p_ij1_25 <- ifelse(matrix$sowV25 >0, matrix$p_ij1*(1-B_s[1]),matrix$p_ij1)
matrix$p_ij1_50 <- ifelse(matrix$sowV50 >0, matrix$p_ij1*(1-B_s[1]),matrix$p_ij1)
matrix$p_ij1_75 <- ifelse(matrix$sowV75 >0, matrix$p_ij1*(1-B_s[1]),matrix$p_ij1)
matrix$p_ij1_100 <- ifelse(matrix$sowV100 >0, matrix$p_ij1*(1-B_s[1]),matrix$p_ij1)

matrix$p_ij.6_25 <- ifelse(matrix$sowV25 >0, matrix$p_ij.6*(1-B_s[1]),matrix$p_ij.6)
matrix$p_ij.6_50 <- ifelse(matrix$sowV50 >0, matrix$p_ij.6*(1-B_s[1]),matrix$p_ij.6)
matrix$p_ij.6_75 <- ifelse(matrix$sowV75 >0, matrix$p_ij.6*(1-B_s[1]),matrix$p_ij.6)
matrix$p_ij.6_100 <- ifelse(matrix$sowV100 >0, matrix$p_ij.6*(1-B_s[1]),matrix$p_ij.6)

matrix$p_ij.3_25 <- ifelse(matrix$sowV25 >0, matrix$p_ij.3*(1-B_s[1]),matrix$p_ij.3)
matrix$p_ij.3_50 <- ifelse(matrix$sowV50 >0, matrix$p_ij.3*(1-B_s[1]),matrix$p_ij.3)
matrix$p_ij.3_75 <- ifelse(matrix$sowV75 >0, matrix$p_ij.3*(1-B_s[1]),matrix$p_ij.3)
matrix$p_ij.3_100 <- ifelse(matrix$sowV100 >0, matrix$p_ij.3*(1-B_s[1]),matrix$p_ij.3)

# List Rho -80- 25%
matrix$rho1_F80_25 <- matrix$p_ij1_25 + matrix$k_ij_F80_25 - matrix$p_ij1_25*matrix$k_ij_F80_25 
matrix$rho2_F80_25 <- matrix$p_ij1_25 + matrix$k_ij_up_F80_25 - matrix$p_ij1_25*matrix$k_ij_up_F80_25 
matrix$rho3_F80_25 <- matrix$p_ij1_25 + matrix$k_ij_low_F80_25 - matrix$p_ij1_25*matrix$k_ij_low_F80_25 
matrix$rho4_F80_25 <- matrix$p_ij.3_25 + matrix$k_ij_F80_25 - matrix$p_ij.3_25*matrix$k_ij_F80_25
matrix$rho5_F80_25 <- matrix$p_ij.3_25 + matrix$k_ij_up_F80_25 - matrix$p_ij.3_25*matrix$k_ij_up_F80_25
matrix$rho6_F80_25 <- matrix$p_ij.3_25 + matrix$k_ij_low_F80_25 - matrix$p_ij.3_25*matrix$k_ij_low_F80_25 
matrix$rho7_F80_25 <- matrix$p_ij.6_25 + matrix$k_ij_F80_25 - matrix$p_ij.6_25*matrix$k_ij_F80_25
matrix$rho8_F80_25 <- matrix$p_ij.6_25 + matrix$k_ij_up_F80_25 - matrix$p_ij.6_25*matrix$k_ij_up_F80_25
matrix$rho9_F80_25 <- matrix$p_ij.6_25 + matrix$k_ij_low_F80_25 - matrix$p_ij.6_25*matrix$k_ij_low_F80_25 

rho1_F80_25 <- acast(matrix, destination ~ origin, value.var = 'rho1_F80_25') 
rho2_F80_25 <- acast(matrix, destination ~ origin, value.var = 'rho2_F80_25')
rho3_F80_25 <- acast(matrix, destination ~ origin, value.var = 'rho3_F80_25')
rho4_F80_25 <- acast(matrix, destination ~ origin, value.var = 'rho4_F80_25')
rho5_F80_25 <- acast(matrix, destination ~ origin, value.var = 'rho5_F80_25')
rho6_F80_25 <- acast(matrix, destination ~ origin, value.var = 'rho6_F80_25')
rho7_F80_25 <- acast(matrix, destination ~ origin, value.var = 'rho7_F80_25')
rho8_F80_25 <- acast(matrix, destination ~ origin, value.var = 'rho8_F80_25')
rho9_F80_25 <- acast(matrix, destination ~ origin, value.var = 'rho9_F80_25')

diag(rho1_F80_25) <- 1;diag(rho2_F80_25) <- 1;diag(rho3_F80_25) <- 1;diag(rho4_F80_25) <- 1;diag(rho5_F80_25) <- 1;diag(rho6_F80_25) <- 1;diag(rho7_F80_25) <- 1;diag(rho8_F80_25) <- 1;diag(rho9_F80_25) <- 1 

rho_list_F80_25 = list(rho1_F80_25, rho2_F80_25, rho3_F80_25, rho4_F80_25, rho5_F80_25,
                       rho6_F80_25, rho7_F80_25, rho8_F80_25, rho9_F80_25)

# List Rho -80- 50%
matrix$rho1_F80_50 <- matrix$p_ij1_50 + matrix$k_ij_F80_50 - matrix$p_ij1_50*matrix$k_ij_F80_50 
matrix$rho2_F80_50 <- matrix$p_ij1_50 + matrix$k_ij_up_F80_50 - matrix$p_ij1_50*matrix$k_ij_up_F80_50 # rho max
matrix$rho3_F80_50 <- matrix$p_ij1_50 + matrix$k_ij_low_F80_50 - matrix$p_ij1_50*matrix$k_ij_low_F80_50 
matrix$rho4_F80_50 <- matrix$p_ij.3_50 + matrix$k_ij_F80_50 - matrix$p_ij.3_50*matrix$k_ij_F80_50
matrix$rho5_F80_50 <- matrix$p_ij.3_50 + matrix$k_ij_up_F80_50 - matrix$p_ij.3_50*matrix$k_ij_up_F80_50
matrix$rho6_F80_50 <- matrix$p_ij.3_50 + matrix$k_ij_low_F80_50 - matrix$p_ij.3_50*matrix$k_ij_low_F80_50 # rho min
matrix$rho7_F80_50 <- matrix$p_ij.6_50 + matrix$k_ij_F80_50 - matrix$p_ij.6_50*matrix$k_ij_F80_50
matrix$rho8_F80_50 <- matrix$p_ij.6_50 + matrix$k_ij_up_F80_50 - matrix$p_ij.6_50*matrix$k_ij_up_F80_50
matrix$rho9_F80_50 <- matrix$p_ij.6_50 + matrix$k_ij_low_F80_50 - matrix$p_ij.6_50*matrix$k_ij_low_F80_50 

rho1_F80_50 <- acast(matrix, destination ~ origin, value.var = 'rho1_F80_50') 
rho2_F80_50 <- acast(matrix, destination ~ origin, value.var = 'rho2_F80_50')
rho3_F80_50 <- acast(matrix, destination ~ origin, value.var = 'rho3_F80_50')
rho4_F80_50 <- acast(matrix, destination ~ origin, value.var = 'rho4_F80_50')
rho5_F80_50 <- acast(matrix, destination ~ origin, value.var = 'rho5_F80_50')
rho6_F80_50 <- acast(matrix, destination ~ origin, value.var = 'rho6_F80_50')
rho7_F80_50 <- acast(matrix, destination ~ origin, value.var = 'rho7_F80_50')
rho8_F80_50 <- acast(matrix, destination ~ origin, value.var = 'rho8_F80_50')
rho9_F80_50 <- acast(matrix, destination ~ origin, value.var = 'rho9_F80_50')

diag(rho1_F80_50) <- 1;diag(rho2_F80_50) <- 1;diag(rho3_F80_50) <- 1;diag(rho4_F80_50) <- 1;
diag(rho5_F80_50) <- 1;diag(rho6_F80_50) <- 1;diag(rho7_F80_50) <- 1;diag(rho8_F80_50) <- 1;
diag(rho9_F80_50) <- 1 # Fill diagonal with 1, due to rho_ii=1

rho_list_F80_50 = list(rho1_F80_50, rho2_F80_50, rho3_F80_50, rho4_F80_50, rho5_F80_50,
                       rho6_F80_50, rho7_F80_50, rho8_F80_50, rho9_F80_50)

# List Rho -80- 75%
matrix$rho1_F80_75 <- matrix$p_ij1_75 + matrix$k_ij_F80_75 - matrix$p_ij1_75*matrix$k_ij_F80_75 
matrix$rho2_F80_75 <- matrix$p_ij1_75 + matrix$k_ij_up_F80_75 - matrix$p_ij1_75*matrix$k_ij_up_F80_75 # rho max
matrix$rho3_F80_75 <- matrix$p_ij1_75 + matrix$k_ij_low_F80_75 - matrix$p_ij1_75*matrix$k_ij_low_F80_75 
matrix$rho4_F80_75 <- matrix$p_ij.3_75 + matrix$k_ij_F80_75 - matrix$p_ij.3_75*matrix$k_ij_F80_75
matrix$rho5_F80_75 <- matrix$p_ij.3_75 + matrix$k_ij_up_F80_75 - matrix$p_ij.3_75*matrix$k_ij_up_F80_75
matrix$rho6_F80_75 <- matrix$p_ij.3_75 + matrix$k_ij_low_F80_75 - matrix$p_ij.3_75*matrix$k_ij_low_F80_75 # rho min
matrix$rho7_F80_75 <- matrix$p_ij.6_75 + matrix$k_ij_F80_75 - matrix$p_ij.6_75*matrix$k_ij_F80_75
matrix$rho8_F80_75 <- matrix$p_ij.6_75 + matrix$k_ij_up_F80_75 - matrix$p_ij.6_75*matrix$k_ij_up_F80_75
matrix$rho9_F80_75 <- matrix$p_ij.6_75 + matrix$k_ij_low_F80_75 - matrix$p_ij.6_75*matrix$k_ij_low_F80_75 

rho1_F80_75 <- acast(matrix, destination ~ origin, value.var = 'rho1_F80_75') 
rho2_F80_75 <- acast(matrix, destination ~ origin, value.var = 'rho2_F80_75')
rho3_F80_75 <- acast(matrix, destination ~ origin, value.var = 'rho3_F80_75')
rho4_F80_75 <- acast(matrix, destination ~ origin, value.var = 'rho4_F80_75')
rho5_F80_75 <- acast(matrix, destination ~ origin, value.var = 'rho5_F80_75')
rho6_F80_75 <- acast(matrix, destination ~ origin, value.var = 'rho6_F80_75')
rho7_F80_75 <- acast(matrix, destination ~ origin, value.var = 'rho7_F80_75')
rho8_F80_75 <- acast(matrix, destination ~ origin, value.var = 'rho8_F80_75')
rho9_F80_75 <- acast(matrix, destination ~ origin, value.var = 'rho9_F80_75')

diag(rho1_F80_75) <- 1;diag(rho2_F80_75) <- 1;diag(rho3_F80_75) <- 1;diag(rho4_F80_75) <- 1;
diag(rho5_F80_75) <- 1;diag(rho6_F80_75) <- 1;diag(rho7_F80_75) <- 1;diag(rho8_F80_75) <- 1;
diag(rho9_F80_75) <- 1 # Fill diagonal with 1, due to rho_ii=1

rho_list_F80_75 = list(rho1_F80_75, rho2_F80_75, rho3_F80_75, rho4_F80_75, rho5_F80_75,
                       rho6_F80_75, rho7_F80_75, rho8_F80_75, rho9_F80_75)

# List Rho -80- 100%
matrix$rho1_F80_100 <- matrix$p_ij1_100 + matrix$k_ij_F80_100 - matrix$p_ij1_100*matrix$k_ij_F80_100 
matrix$rho2_F80_100 <- matrix$p_ij1_100 + matrix$k_ij_up_F80_100 - matrix$p_ij1_100*matrix$k_ij_up_F80_100 # rho max
matrix$rho3_F80_100 <- matrix$p_ij1_100 + matrix$k_ij_low_F80_100 - matrix$p_ij1_100*matrix$k_ij_low_F80_100 
matrix$rho4_F80_100 <- matrix$p_ij.3_100 + matrix$k_ij_F80_100 - matrix$p_ij.3_100*matrix$k_ij_F80_100
matrix$rho5_F80_100 <- matrix$p_ij.3_100 + matrix$k_ij_up_F80_100 - matrix$p_ij.3_100*matrix$k_ij_up_F80_100
matrix$rho6_F80_100 <- matrix$p_ij.3_100 + matrix$k_ij_low_F80_100 - matrix$p_ij.3_100*matrix$k_ij_low_F80_100 # rho min
matrix$rho7_F80_100 <- matrix$p_ij.6_100 + matrix$k_ij_F80_100 - matrix$p_ij.6_100*matrix$k_ij_F80_100
matrix$rho8_F80_100 <- matrix$p_ij.6_100 + matrix$k_ij_up_F80_100 - matrix$p_ij.6_100*matrix$k_ij_up_F80_100
matrix$rho9_F80_100 <- matrix$p_ij.6_100 + matrix$k_ij_low_F80_100 - matrix$p_ij.6_100*matrix$k_ij_low_F80_100 

rho1_F80_100 <- acast(matrix, destination ~ origin, value.var = 'rho1_F80_100') 
rho2_F80_100 <- acast(matrix, destination ~ origin, value.var = 'rho2_F80_100')
rho3_F80_100 <- acast(matrix, destination ~ origin, value.var = 'rho3_F80_100')
rho4_F80_100 <- acast(matrix, destination ~ origin, value.var = 'rho4_F80_100')
rho5_F80_100 <- acast(matrix, destination ~ origin, value.var = 'rho5_F80_100')
rho6_F80_100 <- acast(matrix, destination ~ origin, value.var = 'rho6_F80_100')
rho7_F80_100 <- acast(matrix, destination ~ origin, value.var = 'rho7_F80_100')
rho8_F80_100 <- acast(matrix, destination ~ origin, value.var = 'rho8_F80_100')
rho9_F80_100 <- acast(matrix, destination ~ origin, value.var = 'rho9_F80_100')

diag(rho1_F80_100) <- 1;diag(rho2_F80_100) <- 1;diag(rho3_F80_100) <- 1;diag(rho4_F80_100) <- 1;
diag(rho5_F80_100) <- 1;diag(rho6_F80_100) <- 1;diag(rho7_F80_100) <- 1;diag(rho8_F80_100) <- 1;
diag(rho9_F80_100) <- 1 # Fill diagonal with 1, due to rho_ii=1

rho_list_F80_100 = list(rho1_F80_100, rho2_F80_100, rho3_F80_100, rho4_F80_100, rho5_F80_100,
                        rho6_F80_100, rho7_F80_100, rho8_F80_100, rho9_F80_100)

# List Rho -40- 25%
matrix$rho1_F40_25 <- matrix$p_ij1_25 + matrix$k_ij_F40_25 - matrix$p_ij1_25*matrix$k_ij_F40_25 
matrix$rho2_F40_25 <- matrix$p_ij1_25 + matrix$k_ij_up_F40_25 - matrix$p_ij1_25*matrix$k_ij_up_F40_25 
matrix$rho3_F40_25 <- matrix$p_ij1_25 + matrix$k_ij_low_F40_25 - matrix$p_ij1_25*matrix$k_ij_low_F40_25 
matrix$rho4_F40_25 <- matrix$p_ij.3_25 + matrix$k_ij_F40_25 - matrix$p_ij.3_25*matrix$k_ij_F40_25
matrix$rho5_F40_25 <- matrix$p_ij.3_25 + matrix$k_ij_up_F40_25 - matrix$p_ij.3_25*matrix$k_ij_up_F40_25
matrix$rho6_F40_25 <- matrix$p_ij.3_25 + matrix$k_ij_low_F40_25 - matrix$p_ij.3_25*matrix$k_ij_low_F40_25 
matrix$rho7_F40_25 <- matrix$p_ij.6_25 + matrix$k_ij_F40_25 - matrix$p_ij.6_25*matrix$k_ij_F40_25
matrix$rho8_F40_25 <- matrix$p_ij.6_25 + matrix$k_ij_up_F40_25 - matrix$p_ij.6_25*matrix$k_ij_up_F40_25
matrix$rho9_F40_25 <- matrix$p_ij.6_25 + matrix$k_ij_low_F40_25 - matrix$p_ij.6_25*matrix$k_ij_low_F40_25 

rho1_F40_25 <- acast(matrix, destination ~ origin, value.var = 'rho1_F40_25') 
rho2_F40_25 <- acast(matrix, destination ~ origin, value.var = 'rho2_F40_25')
rho3_F40_25 <- acast(matrix, destination ~ origin, value.var = 'rho3_F40_25')
rho4_F40_25 <- acast(matrix, destination ~ origin, value.var = 'rho4_F40_25')
rho5_F40_25 <- acast(matrix, destination ~ origin, value.var = 'rho5_F40_25')
rho6_F40_25 <- acast(matrix, destination ~ origin, value.var = 'rho6_F40_25')
rho7_F40_25 <- acast(matrix, destination ~ origin, value.var = 'rho7_F40_25')
rho8_F40_25 <- acast(matrix, destination ~ origin, value.var = 'rho8_F40_25')
rho9_F40_25 <- acast(matrix, destination ~ origin, value.var = 'rho9_F40_25')

diag(rho1_F40_25) <- 1;diag(rho2_F40_25) <- 1;diag(rho3_F40_25) <- 1;diag(rho4_F40_25) <- 1;diag(rho5_F40_25) <- 1;diag(rho6_F40_25) <- 1;diag(rho7_F40_25) <- 1;diag(rho8_F40_25) <- 1;diag(rho9_F40_25) <- 1 

rho_list_F40_25 = list(rho1_F40_25, rho2_F40_25, rho3_F40_25, rho4_F40_25, rho5_F40_25,
                       rho6_F40_25, rho7_F40_25, rho8_F40_25, rho9_F40_25)

# List Rho -40- 50%
matrix$rho1_F40_50 <- matrix$p_ij1_50 + matrix$k_ij_F40_50 - matrix$p_ij1_50*matrix$k_ij_F40_50 
matrix$rho2_F40_50 <- matrix$p_ij1_50 + matrix$k_ij_up_F40_50 - matrix$p_ij1_50*matrix$k_ij_up_F40_50 # rho max
matrix$rho3_F40_50 <- matrix$p_ij1_50 + matrix$k_ij_low_F40_50 - matrix$p_ij1_50*matrix$k_ij_low_F40_50 
matrix$rho4_F40_50 <- matrix$p_ij.3_50 + matrix$k_ij_F40_50 - matrix$p_ij.3_50*matrix$k_ij_F40_50
matrix$rho5_F40_50 <- matrix$p_ij.3_50 + matrix$k_ij_up_F40_50 - matrix$p_ij.3_50*matrix$k_ij_up_F40_50
matrix$rho6_F40_50 <- matrix$p_ij.3_50 + matrix$k_ij_low_F40_50 - matrix$p_ij.3_50*matrix$k_ij_low_F40_50 # rho min
matrix$rho7_F40_50 <- matrix$p_ij.6_50 + matrix$k_ij_F40_50 - matrix$p_ij.6_50*matrix$k_ij_F40_50
matrix$rho8_F40_50 <- matrix$p_ij.6_50 + matrix$k_ij_up_F40_50 - matrix$p_ij.6_50*matrix$k_ij_up_F40_50
matrix$rho9_F40_50 <- matrix$p_ij.6_50 + matrix$k_ij_low_F40_50 - matrix$p_ij.6_50*matrix$k_ij_low_F40_50 

rho1_F40_50 <- acast(matrix, destination ~ origin, value.var = 'rho1_F40_50') 
rho2_F40_50 <- acast(matrix, destination ~ origin, value.var = 'rho2_F40_50')
rho3_F40_50 <- acast(matrix, destination ~ origin, value.var = 'rho3_F40_50')
rho4_F40_50 <- acast(matrix, destination ~ origin, value.var = 'rho4_F40_50')
rho5_F40_50 <- acast(matrix, destination ~ origin, value.var = 'rho5_F40_50')
rho6_F40_50 <- acast(matrix, destination ~ origin, value.var = 'rho6_F40_50')
rho7_F40_50 <- acast(matrix, destination ~ origin, value.var = 'rho7_F40_50')
rho8_F40_50 <- acast(matrix, destination ~ origin, value.var = 'rho8_F40_50')
rho9_F40_50 <- acast(matrix, destination ~ origin, value.var = 'rho9_F40_50')

diag(rho1_F40_50) <- 1;diag(rho2_F40_50) <- 1;diag(rho3_F40_50) <- 1;diag(rho4_F40_50) <- 1;
diag(rho5_F40_50) <- 1;diag(rho6_F40_50) <- 1;diag(rho7_F40_50) <- 1;diag(rho8_F40_50) <- 1;
diag(rho9_F40_50) <- 1 

rho_list_F40_50 = list(rho1_F40_50, rho2_F40_50, rho3_F40_50, rho4_F40_50, rho5_F40_50,
                       rho6_F40_50, rho7_F40_50, rho8_F40_50, rho9_F40_50)

# List Rho -40- 75%
matrix$rho1_F40_75 <- matrix$p_ij1_75 + matrix$k_ij_F40_75 - matrix$p_ij1_75*matrix$k_ij_F40_75 
matrix$rho2_F40_75 <- matrix$p_ij1_75 + matrix$k_ij_up_F40_75 - matrix$p_ij1_75*matrix$k_ij_up_F40_75 # rho max
matrix$rho3_F40_75 <- matrix$p_ij1_75 + matrix$k_ij_low_F40_75 - matrix$p_ij1_75*matrix$k_ij_low_F40_75 
matrix$rho4_F40_75 <- matrix$p_ij.3_75 + matrix$k_ij_F40_75 - matrix$p_ij.3_75*matrix$k_ij_F40_75
matrix$rho5_F40_75 <- matrix$p_ij.3_75 + matrix$k_ij_up_F40_75 - matrix$p_ij.3_75*matrix$k_ij_up_F40_75
matrix$rho6_F40_75 <- matrix$p_ij.3_75 + matrix$k_ij_low_F40_75 - matrix$p_ij.3_75*matrix$k_ij_low_F40_75 # rho min
matrix$rho7_F40_75 <- matrix$p_ij.6_75 + matrix$k_ij_F40_75 - matrix$p_ij.6_75*matrix$k_ij_F40_75
matrix$rho8_F40_75 <- matrix$p_ij.6_75 + matrix$k_ij_up_F40_75 - matrix$p_ij.6_75*matrix$k_ij_up_F40_75
matrix$rho9_F40_75 <- matrix$p_ij.6_75 + matrix$k_ij_low_F40_75 - matrix$p_ij.6_75*matrix$k_ij_low_F40_75 

rho1_F40_75 <- acast(matrix, destination ~ origin, value.var = 'rho1_F40_75') 
rho2_F40_75 <- acast(matrix, destination ~ origin, value.var = 'rho2_F40_75')
rho3_F40_75 <- acast(matrix, destination ~ origin, value.var = 'rho3_F40_75')
rho4_F40_75 <- acast(matrix, destination ~ origin, value.var = 'rho4_F40_75')
rho5_F40_75 <- acast(matrix, destination ~ origin, value.var = 'rho5_F40_75')
rho6_F40_75 <- acast(matrix, destination ~ origin, value.var = 'rho6_F40_75')
rho7_F40_75 <- acast(matrix, destination ~ origin, value.var = 'rho7_F40_75')
rho8_F40_75 <- acast(matrix, destination ~ origin, value.var = 'rho8_F40_75')
rho9_F40_75 <- acast(matrix, destination ~ origin, value.var = 'rho9_F40_75')

diag(rho1_F40_75) <- 1;diag(rho2_F40_75) <- 1;diag(rho3_F40_75) <- 1;diag(rho4_F40_75) <- 1;
diag(rho5_F40_75) <- 1;diag(rho6_F40_75) <- 1;diag(rho7_F40_75) <- 1;diag(rho8_F40_75) <- 1;
diag(rho9_F40_75) <- 1 # Fill diagonal with 1, due to rho_ii=1

rho_list_F40_75 = list(rho1_F40_75, rho2_F40_75, rho3_F40_75, rho4_F40_75, rho5_F40_75,
                       rho6_F40_75, rho7_F40_75, rho8_F40_75, rho9_F40_75)

# List Rho -40- 100%
matrix$rho1_F40_100 <- matrix$p_ij1_100 + matrix$k_ij_F40_100 - matrix$p_ij1_100*matrix$k_ij_F40_100 
matrix$rho2_F40_100 <- matrix$p_ij1_100 + matrix$k_ij_up_F40_100 - matrix$p_ij1_100*matrix$k_ij_up_F40_100 # rho max
matrix$rho3_F40_100 <- matrix$p_ij1_100 + matrix$k_ij_low_F40_100 - matrix$p_ij1_100*matrix$k_ij_low_F40_100 
matrix$rho4_F40_100 <- matrix$p_ij.3_100 + matrix$k_ij_F40_100 - matrix$p_ij.3_100*matrix$k_ij_F40_100
matrix$rho5_F40_100 <- matrix$p_ij.3_100 + matrix$k_ij_up_F40_100 - matrix$p_ij.3_100*matrix$k_ij_up_F40_100
matrix$rho6_F40_100 <- matrix$p_ij.3_100 + matrix$k_ij_low_F40_100 - matrix$p_ij.3_100*matrix$k_ij_low_F40_100 # rho min
matrix$rho7_F40_100 <- matrix$p_ij.6_100 + matrix$k_ij_F40_100 - matrix$p_ij.6_100*matrix$k_ij_F40_100
matrix$rho8_F40_100 <- matrix$p_ij.6_100 + matrix$k_ij_up_F40_100 - matrix$p_ij.6_100*matrix$k_ij_up_F40_100
matrix$rho9_F40_100 <- matrix$p_ij.6_100 + matrix$k_ij_low_F40_100 - matrix$p_ij.6_100*matrix$k_ij_low_F40_100 

rho1_F40_100 <- acast(matrix, destination ~ origin, value.var = 'rho1_F40_100') 
rho2_F40_100 <- acast(matrix, destination ~ origin, value.var = 'rho2_F40_100')
rho3_F40_100 <- acast(matrix, destination ~ origin, value.var = 'rho3_F40_100')
rho4_F40_100 <- acast(matrix, destination ~ origin, value.var = 'rho4_F40_100')
rho5_F40_100 <- acast(matrix, destination ~ origin, value.var = 'rho5_F40_100')
rho6_F40_100 <- acast(matrix, destination ~ origin, value.var = 'rho6_F40_100')
rho7_F40_100 <- acast(matrix, destination ~ origin, value.var = 'rho7_F40_100')
rho8_F40_100 <- acast(matrix, destination ~ origin, value.var = 'rho8_F40_100')
rho9_F40_100 <- acast(matrix, destination ~ origin, value.var = 'rho9_F40_100')

diag(rho1_F40_100) <- 1;diag(rho2_F40_100) <- 1;diag(rho3_F40_100) <- 1;diag(rho4_F40_100) <- 1;
diag(rho5_F40_100) <- 1;diag(rho6_F40_100) <- 1;diag(rho7_F40_100) <- 1;diag(rho8_F40_100) <- 1;
diag(rho9_F40_100) <- 1 # Fill diagonal with 1, due to rho_ii=1

rho_list_F40_100 = list(rho1_F40_100, rho2_F40_100, rho3_F40_100, rho4_F40_100, rho5_F40_100,
                        rho6_F40_100, rho7_F40_100, rho8_F40_100, rho9_F40_100)

# Set new parameters for Dv (Days of infection) in vaccinated animals
Dv_s_r = c(0.7, 2.8)

# D- 25% Coverage 
get_D25 = function(k, D_s_r, D_p_r, Dv_s_r, Farms){ # Add Dv_s_r
  if(k == 1){
    D25 = ifelse(Farms$sowV25 == 1, Dv_s_r[1], ifelse(Farms$sow == 1, D_s_r[1], D_p_r[1]))
  } else {
    if(k == 2) {
      D25 = ifelse(Farms$sowV25 == 1, Dv_s_r[2], ifelse(Farms$sow == 1, D_s_r[2], D_p_r[2]))
    } else {
      D25 = ifelse(Farms$sowV25 == 1,
                   runif(1, min = Dv_s_r[1], max = Dv_s_r[2]),
                   ifelse(Farms$sow == 1, 
                          runif(1, min = D_s_r[1], max = D_s_r[2]),
                          runif(1, min = D_p_r[1], max = D_p_r[2])))
    }
  }
  return(D25)
} 
D25 = lapply(1:nsims, get_D25, D_s_r, D_p_r, Dv_s_r, Farms)

# D- 50% Coverage 
get_D50 = function(k, D_s_r, D_p_r, Dv_s_r, Farms){ # Add Dv_s_r
  if(k == 1){
    D50 = ifelse(Farms$sowV50 == 1, Dv_s_r[1], ifelse(Farms$sow == 1, D_s_r[1], D_p_r[1]))
  } else {
    if(k == 2) {
      D50 = ifelse(Farms$sowV50 == 1, Dv_s_r[2], ifelse(Farms$sow == 1, D_s_r[2], D_p_r[2]))
    } else {
      D50 = ifelse(Farms$sowV50 == 1,
                   runif(1, min = Dv_s_r[1], max = Dv_s_r[2]),
                   ifelse(Farms$sow == 1, 
                          runif(1, min = D_s_r[1], max = D_s_r[2]),
                          runif(1, min = D_p_r[1], max = D_p_r[2])))
    }
  }
  return(D50)
} 
D50 = lapply(1:nsims, get_D50, D_s_r, D_p_r, Dv_s_r, Farms)

# D- 75% Coverage
get_D75 = function(k, D_s_r, D_p_r, Dv_s_r, Farms){ # Add Dv_s_r
  if(k == 1){
    D75 = ifelse(Farms$sowV75 == 1, Dv_s_r[1], ifelse(Farms$sow == 1, D_s_r[1], D_p_r[1]))
  } else {
    if(k == 2) {
      D75 = ifelse(Farms$sowV75 == 1, Dv_s_r[2], ifelse(Farms$sow == 1, D_s_r[2], D_p_r[2]))
    } else {
      D75 = ifelse(Farms$sowV75 == 1,
                   runif(1, min = Dv_s_r[1], max = Dv_s_r[2]),
                   ifelse(Farms$sow == 1, 
                          runif(1, min = D_s_r[1], max = D_s_r[2]),
                          runif(1, min = D_p_r[1], max = D_p_r[2])))
    }
  }
  return(D75)
} 
D75 = lapply(1:nsims, get_D75, D_s_r, D_p_r, Dv_s_r, Farms)

# D- 100% Coverage
get_D100 = function(k, D_s_r, D_p_r, Dv_s_r, Farms){ # Add Dv_s_r
  if(k == 1){
    D100 = ifelse(Farms$sowV100 == 1, Dv_s_r[1], ifelse(Farms$sow == 1, D_s_r[1], D_p_r[1]))
  } else {
    if(k == 2) {
      D100 = ifelse(Farms$sowV100 == 1, Dv_s_r[2], ifelse(Farms$sow == 1, D_s_r[2], D_p_r[2]))
    } else {
      D100 = ifelse(Farms$sowV100 == 1,
                    runif(1, min = Dv_s_r[1], max = Dv_s_r[2]),
                    ifelse(Farms$sow == 1, 
                           runif(1, min = D_s_r[1], max = D_s_r[2]),
                           runif(1, min = D_p_r[1], max = D_p_r[2])))
    }
  }
  return(D100)
} 
D100 = lapply(1:nsims, get_D100, D_s_r, D_p_r, Dv_s_r, Farms)

# I have been trying to create a single function, where I can vary only 3 componentes (beta_list, rho_list and D) of these loops to obtain different scenatios. Unfortunatelly I have had error messages again and again. 

# Only 25% vaccine ------------------------------------------------------
sim_res1 = vector("list", nsims)
sim_res2 = vector("list", nsims)
names(sim_res1) = 1:nsims
names(sim_res2) = 1:nsims
result_1 = vector("list",27)
result_2 = vector("list",27)
names(result_1) = c(paste0("beta_l_",1:9), paste0("beta_m_",1:9),paste0("beta_h_",1:9) )
names(result_2) = c(paste0("beta_l_",1:9), paste0("beta_m_",1:9),paste0("beta_h_",1:9) )

for( k in 1:(nsims)) {
  loop=0
  
  for(i in 1:3){
    beta = beta_listV25[[i]]
    
    for(j in 1:9){
      print(paste0("k = ", k, "; i = ", i, "; j = ", j))
      loop = loop+1
      
      rho = rho_list[[j]]
      beta_matrix = rho * beta 
      
      stocks <- c(S=Farms$X, I=Farms$Y, R=Farms$Z)
      parameters <- list("beta"=beta_matrix, "delays"=c(D25[[k]]), 
                         "mort"=c(v[[k]]), "returns"=c(D_imm[[k]])) #remove farms.  
      
      out <- data.frame(ode(y=stocks, times=time, func=model, parms=parameters, 
                            method = "euler")) 
      out <- out[out$time %in% seq(0,26,1), ]
      
      o <-data.frame(time=out$time)
      o$S <- apply(out[,c(grep('S', names(out), value=TRUE))], 1, sum)
      o$I <- apply(out[,c(grep('I', names(out), value=TRUE))], 1, sum)
      o$R <- apply(out[,c(grep('R', names(out), value=TRUE))], 1, sum)
      o$animals <- o$S+o$R+o$I
      o$prop <- o$I/o$animals
      
      result_1[[loop]] = out
      result_2[[loop]] = o
    }
  }
  sim_res1[[k]] = result_1
  sim_res2[[k]] = result_2
}
#saveRDS(sim_res1, file="/Volumes/PVD2/Davis/PhD/DISSERTATION/Chapter4/R4/PRRSTM/newset3/sim_res1V25.rds")
#saveRDS(sim_res2, file="/Volumes/PVD2/Davis/PhD/DISSERTATION/Chapter4/R4/PRRSTM/newset3/sim_res2V25.rds")


stop("Run first simulation")

# 25% vaccine 25 filter (80) -----------------------------------------------
sim_res1 = vector("list", nsims)
sim_res2 = vector("list", nsims)
names(sim_res1) = 1:nsims
names(sim_res2) = 1:nsims
result_1 = vector("list",27)
result_2 = vector("list",27)
names(result_1) = c(paste0("beta_l_",1:9), paste0("beta_m_",1:9),paste0("beta_h_",1:9) )
names(result_2) = c(paste0("beta_l_",1:9), paste0("beta_m_",1:9),paste0("beta_h_",1:9) )

for( k in 1:(nsims)) {
  loop=0
  
  for(i in 1:3){
    beta = beta_listV25[[i]]
    
    for(j in 1:9){
      print(paste0("k = ", k, "; i = ", i, "; j = ", j))
      loop = loop+1
      
      rho = rho_list_F80_25[[j]]
      beta_matrix = rho * beta 
      
      stocks <- c(S=Farms$X, I=Farms$Y, R=Farms$Z)
      parameters <- list("beta"=beta_matrix, "delays"=c(D25[[k]]), 
                         "mort"=c(v[[k]]), "returns"=c(D_imm[[k]])) #remove farms.  
      
      out <- data.frame(ode(y=stocks, times=time, func=model, parms=parameters, 
                            method = "euler")) 
      out <- out[out$time %in% seq(0,26,1), ]
      
      o <-data.frame(time=out$time)
      o$S <- apply(out[,c(grep('S', names(out), value=TRUE))], 1, sum)
      o$I <- apply(out[,c(grep('I', names(out), value=TRUE))], 1, sum)
      o$R <- apply(out[,c(grep('R', names(out), value=TRUE))], 1, sum)
      o$animals <- o$S+o$R+o$I
      o$prop <- o$I/o$animals
      
      result_1[[loop]] = out
      result_2[[loop]] = o
    }
  }
  sim_res1[[k]] = result_1
  sim_res2[[k]] = result_2
}
#saveRDS(sim_res1, file="/Volumes/PVD2/Davis/PhD/DISSERTATION/Chapter4/R4/PRRSTM/newset3/sim_res1V25_F80_25.rds")
#saveRDS(sim_res2, file="/Volumes/PVD2/Davis/PhD/DISSERTATION/Chapter4/R4/PRRSTM/newset3/sim_res2V25_F80_25.rds")

# 25% vaccine 25 filter (40) -----------------------------------------------
sim_res1 = vector("list", nsims)
sim_res2 = vector("list", nsims)
names(sim_res1) = 1:nsims
names(sim_res2) = 1:nsims
result_1 = vector("list",27)
result_2 = vector("list",27)
names(result_1) = c(paste0("beta_l_",1:9), paste0("beta_m_",1:9),paste0("beta_h_",1:9) )
names(result_2) = c(paste0("beta_l_",1:9), paste0("beta_m_",1:9),paste0("beta_h_",1:9) )

for( k in 1:(nsims)) {
  loop=0
  
  for(i in 1:3){
    beta = beta_listV25[[i]]
    
    for(j in 1:9){
      print(paste0("k = ", k, "; i = ", i, "; j = ", j))
      loop = loop+1
      
      rho = rho_list_F40_25[[j]]
      beta_matrix = rho * beta 
      
      stocks <- c(S=Farms$X, I=Farms$Y, R=Farms$Z)
      parameters <- list("beta"=beta_matrix, "delays"=c(D25[[k]]), 
                         "mort"=c(v[[k]]), "returns"=c(D_imm[[k]])) #remove farms.  
      
      out <- data.frame(ode(y=stocks, times=time, func=model, parms=parameters, 
                            method = "euler")) 
      out <- out[out$time %in% seq(0,26,1), ]
      
      o <-data.frame(time=out$time)
      o$S <- apply(out[,c(grep('S', names(out), value=TRUE))], 1, sum)
      o$I <- apply(out[,c(grep('I', names(out), value=TRUE))], 1, sum)
      o$R <- apply(out[,c(grep('R', names(out), value=TRUE))], 1, sum)
      o$animals <- o$S+o$R+o$I
      o$prop <- o$I/o$animals
      
      result_1[[loop]] = out
      result_2[[loop]] = o
    }
  }
  sim_res1[[k]] = result_1
  sim_res2[[k]] = result_2
}
#saveRDS(sim_res1, file="/Volumes/PVD2/Davis/PhD/DISSERTATION/Chapter4/R4/PRRSTM/newset3/sim_res1V25_F40_25.rds")
#saveRDS(sim_res2, file="/Volumes/PVD2/Davis/PhD/DISSERTATION/Chapter4/R4/PRRSTM/newset3/sim_res2V25_F40_25.rds")

# Only 50% vaccine --------------------------------------------------------
sim_res1 = vector("list", nsims)
sim_res2 = vector("list", nsims)
names(sim_res1) = 1:nsims
names(sim_res2) = 1:nsims
result_1 = vector("list",27)
result_2 = vector("list",27)
names(result_1) = c(paste0("beta_l_",1:9), paste0("beta_m_",1:9),paste0("beta_h_",1:9) )
names(result_2) = c(paste0("beta_l_",1:9), paste0("beta_m_",1:9),paste0("beta_h_",1:9) )

for( k in 1:(nsims)) {
  loop=0
  
  for(i in 1:3){
    beta = beta_listV50[[i]]
    
    for(j in 1:9){
      print(paste0("k = ", k, "; i = ", i, "; j = ", j))
      loop = loop+1
      
      rho = rho_list[[j]]
      beta_matrix = rho * beta 
      
      stocks <- c(S=Farms$X, I=Farms$Y, R=Farms$Z)
      parameters <- list("beta"=beta_matrix, "delays"=c(D50[[k]]), 
                         "mort"=c(v[[k]]), "returns"=c(D_imm[[k]])) #remove farms.  
      
      out <- data.frame(ode(y=stocks, times=time, func=model, parms=parameters, 
                            method = "euler")) 
      out <- out[out$time %in% seq(0,26,1), ]
      
      o <-data.frame(time=out$time)
      o$S <- apply(out[,c(grep('S', names(out), value=TRUE))], 1, sum)
      o$I <- apply(out[,c(grep('I', names(out), value=TRUE))], 1, sum)
      o$R <- apply(out[,c(grep('R', names(out), value=TRUE))], 1, sum)
      o$animals <- o$S+o$R+o$I
      o$prop <- o$I/o$animals
      
      result_1[[loop]] = out
      result_2[[loop]] = o
    }
  }
  sim_res1[[k]] = result_1
  sim_res2[[k]] = result_2
}
#saveRDS(sim_res1, file="/Volumes/PVD2/Davis/PhD/DISSERTATION/Chapter4/R4/PRRSTM/newset3/sim_res1V50.rds")
#saveRDS(sim_res2, file="/Volumes/PVD2/Davis/PhD/DISSERTATION/Chapter4/R4/PRRSTM/newset3/sim_res2V50.rds")


# 50% vaccine 50% filter (80) ---------------------------------------------
sim_res1 = vector("list", nsims)
sim_res2 = vector("list", nsims)
names(sim_res1) = 1:nsims
names(sim_res2) = 1:nsims
result_1 = vector("list",27)
result_2 = vector("list",27)
names(result_1) = c(paste0("beta_l_",1:9), paste0("beta_m_",1:9),paste0("beta_h_",1:9) )
names(result_2) = c(paste0("beta_l_",1:9), paste0("beta_m_",1:9),paste0("beta_h_",1:9) )

for( k in 1:(nsims)) {
  loop=0
  
  for(i in 1:3){
    beta = beta_listV50[[i]]
    
    for(j in 1:9){
      print(paste0("k = ", k, "; i = ", i, "; j = ", j))
      loop = loop+1
      
      rho = rho_list_F80_50[[j]]
      beta_matrix = rho * beta 
      
      stocks <- c(S=Farms$X, I=Farms$Y, R=Farms$Z)
      parameters <- list("beta"=beta_matrix, "delays"=c(D50[[k]]), 
                         "mort"=c(v[[k]]), "returns"=c(D_imm[[k]])) #remove farms.  
      
      out <- data.frame(ode(y=stocks, times=time, func=model, parms=parameters, 
                            method = "euler")) 
      out <- out[out$time %in% seq(0,26,1), ]
      
      o <-data.frame(time=out$time)
      o$S <- apply(out[,c(grep('S', names(out), value=TRUE))], 1, sum)
      o$I <- apply(out[,c(grep('I', names(out), value=TRUE))], 1, sum)
      o$R <- apply(out[,c(grep('R', names(out), value=TRUE))], 1, sum)
      o$animals <- o$S+o$R+o$I
      o$prop <- o$I/o$animals
      
      result_1[[loop]] = out
      result_2[[loop]] = o
    }
  }
  sim_res1[[k]] = result_1
  sim_res2[[k]] = result_2
}
#saveRDS(sim_res1, file="/Volumes/PVD2/Davis/PhD/DISSERTATION/Chapter4/R4/PRRSTM/newset3/sim_res1V50_F80_50.rds")
#saveRDS(sim_res2, file="/Volumes/PVD2/Davis/PhD/DISSERTATION/Chapter4/R4/PRRSTM/newset3/sim_res2V50_F80_50.rds")

# 50% vaccine 50% filter (40) ---------------------------------------------
sim_res1 = vector("list", nsims)
sim_res2 = vector("list", nsims)
names(sim_res1) = 1:nsims
names(sim_res2) = 1:nsims
result_1 = vector("list",27)
result_2 = vector("list",27)
names(result_1) = c(paste0("beta_l_",1:9), paste0("beta_m_",1:9),paste0("beta_h_",1:9) )
names(result_2) = c(paste0("beta_l_",1:9), paste0("beta_m_",1:9),paste0("beta_h_",1:9) )

for( k in 1:(nsims)) {
  loop=0
  
  for(i in 1:3){
    beta = beta_listV50[[i]]
    
    for(j in 1:9){
      print(paste0("k = ", k, "; i = ", i, "; j = ", j))
      loop = loop+1
      
      rho = rho_list_F40_50[[j]]
      beta_matrix = rho * beta 
      
      stocks <- c(S=Farms$X, I=Farms$Y, R=Farms$Z)
      parameters <- list("beta"=beta_matrix, "delays"=c(D50[[k]]), 
                         "mort"=c(v[[k]]), "returns"=c(D_imm[[k]])) #remove farms.  
      
      out <- data.frame(ode(y=stocks, times=time, func=model, parms=parameters, 
                            method = "euler")) 
      out <- out[out$time %in% seq(0,26,1), ]
      
      o <-data.frame(time=out$time)
      o$S <- apply(out[,c(grep('S', names(out), value=TRUE))], 1, sum)
      o$I <- apply(out[,c(grep('I', names(out), value=TRUE))], 1, sum)
      o$R <- apply(out[,c(grep('R', names(out), value=TRUE))], 1, sum)
      o$animals <- o$S+o$R+o$I
      o$prop <- o$I/o$animals
      
      result_1[[loop]] = out
      result_2[[loop]] = o
    }
  }
  sim_res1[[k]] = result_1
  sim_res2[[k]] = result_2
}
#saveRDS(sim_res1, file="/Volumes/PVD2/Davis/PhD/DISSERTATION/Chapter4/R4/PRRSTM/newset3/sim_res1V50_F40_50.rds")
#saveRDS(sim_res2, file="/Volumes/PVD2/Davis/PhD/DISSERTATION/Chapter4/R4/PRRSTM/newset3/sim_res2V50_F40_50.rds")

# Only 75% vaccine --------------------------------------------------------
sim_res1 = vector("list", nsims)
sim_res2 = vector("list", nsims)
names(sim_res1) = 1:nsims
names(sim_res2) = 1:nsims
result_1 = vector("list",27)
result_2 = vector("list",27)
names(result_1) = c(paste0("beta_l_",1:9), paste0("beta_m_",1:9),paste0("beta_h_",1:9) )
names(result_2) = c(paste0("beta_l_",1:9), paste0("beta_m_",1:9),paste0("beta_h_",1:9) )

for( k in 1:(nsims)) {
  loop=0
  
  for(i in 1:3){
    beta = beta_listV75[[i]]
    
    for(j in 1:9){
      print(paste0("k = ", k, "; i = ", i, "; j = ", j))
      loop = loop+1
      
      rho = rho_list[[j]]
      beta_matrix = rho * beta 
      
      stocks <- c(S=Farms$X, I=Farms$Y, R=Farms$Z)
      parameters <- list("beta"=beta_matrix, "delays"=c(D75[[k]]), 
                         "mort"=c(v[[k]]), "returns"=c(D_imm[[k]])) #remove farms.  
      
      out <- data.frame(ode(y=stocks, times=time, func=model, parms=parameters, 
                            method = "euler")) 
      out <- out[out$time %in% seq(0,26,1), ]
      
      o <-data.frame(time=out$time)
      o$S <- apply(out[,c(grep('S', names(out), value=TRUE))], 1, sum)
      o$I <- apply(out[,c(grep('I', names(out), value=TRUE))], 1, sum)
      o$R <- apply(out[,c(grep('R', names(out), value=TRUE))], 1, sum)
      o$animals <- o$S+o$R+o$I
      o$prop <- o$I/o$animals
      
      result_1[[loop]] = out
      result_2[[loop]] = o
    }
  }
  sim_res1[[k]] = result_1
  sim_res2[[k]] = result_2
}
#saveRDS(sim_res1, file="/Volumes/PVD2/Davis/PhD/DISSERTATION/Chapter4/R4/PRRSTM/newset3/sim_res1V75.rds")
#saveRDS(sim_res2, file="/Volumes/PVD2/Davis/PhD/DISSERTATION/Chapter4/R4/PRRSTM/newset3/sim_res2V75.rds")

# 75% vaccine 75% filter (80) ---------------------------------------------
sim_res1 = vector("list", nsims)
sim_res2 = vector("list", nsims)
names(sim_res1) = 1:nsims
names(sim_res2) = 1:nsims
result_1 = vector("list",27)
result_2 = vector("list",27)
names(result_1) = c(paste0("beta_l_",1:9), paste0("beta_m_",1:9),paste0("beta_h_",1:9) )
names(result_2) = c(paste0("beta_l_",1:9), paste0("beta_m_",1:9),paste0("beta_h_",1:9) )

for( k in 1:(nsims)) {
  loop=0
  
  for(i in 1:3){
    beta = beta_listV75[[i]]
    
    for(j in 1:9){
      print(paste0("k = ", k, "; i = ", i, "; j = ", j))
      loop = loop+1
      
      rho = rho_list_F80_75[[j]]
      beta_matrix = rho * beta 
      
      stocks <- c(S=Farms$X, I=Farms$Y, R=Farms$Z)
      parameters <- list("beta"=beta_matrix, "delays"=c(D75[[k]]), 
                         "mort"=c(v[[k]]), "returns"=c(D_imm[[k]])) #remove farms.  
      
      out <- data.frame(ode(y=stocks, times=time, func=model, parms=parameters, 
                            method = "euler")) 
      out <- out[out$time %in% seq(0,26,1), ]
      
      o <-data.frame(time=out$time)
      o$S <- apply(out[,c(grep('S', names(out), value=TRUE))], 1, sum)
      o$I <- apply(out[,c(grep('I', names(out), value=TRUE))], 1, sum)
      o$R <- apply(out[,c(grep('R', names(out), value=TRUE))], 1, sum)
      o$animals <- o$S+o$R+o$I
      o$prop <- o$I/o$animals
      
      result_1[[loop]] = out
      result_2[[loop]] = o
    }
  }
  sim_res1[[k]] = result_1
  sim_res2[[k]] = result_2
}
#saveRDS(sim_res1, file="/Volumes/PVD2/Davis/PhD/DISSERTATION/Chapter4/R4/PRRSTM/newset3/sim_res1V75_F80_75.rds")
#saveRDS(sim_res2, file="/Volumes/PVD2/Davis/PhD/DISSERTATION/Chapter4/R4/PRRSTM/newset3/sim_res2V75_F80_75.rds")

# 75% vaccine 75% filter (40) ---------------------------------------------
sim_res1 = vector("list", nsims)
sim_res2 = vector("list", nsims)
names(sim_res1) = 1:nsims
names(sim_res2) = 1:nsims
result_1 = vector("list",27)
result_2 = vector("list",27)
names(result_1) = c(paste0("beta_l_",1:9), paste0("beta_m_",1:9),paste0("beta_h_",1:9) )
names(result_2) = c(paste0("beta_l_",1:9), paste0("beta_m_",1:9),paste0("beta_h_",1:9) )

for( k in 1:(nsims)) {
  loop=0
  
  for(i in 1:3){
    beta = beta_listV75[[i]]
    
    for(j in 1:9){
      print(paste0("k = ", k, "; i = ", i, "; j = ", j))
      loop = loop+1
      
      rho = rho_list_F40_75[[j]]
      beta_matrix = rho * beta 
      
      stocks <- c(S=Farms$X, I=Farms$Y, R=Farms$Z)
      parameters <- list("beta"=beta_matrix, "delays"=c(D75[[k]]), 
                         "mort"=c(v[[k]]), "returns"=c(D_imm[[k]])) #remove farms.  
      
      out <- data.frame(ode(y=stocks, times=time, func=model, parms=parameters, 
                            method = "euler")) 
      out <- out[out$time %in% seq(0,26,1), ]
      
      o <-data.frame(time=out$time)
      o$S <- apply(out[,c(grep('S', names(out), value=TRUE))], 1, sum)
      o$I <- apply(out[,c(grep('I', names(out), value=TRUE))], 1, sum)
      o$R <- apply(out[,c(grep('R', names(out), value=TRUE))], 1, sum)
      o$animals <- o$S+o$R+o$I
      o$prop <- o$I/o$animals
      
      result_1[[loop]] = out
      result_2[[loop]] = o
    }
  }
  sim_res1[[k]] = result_1
  sim_res2[[k]] = result_2
}
#saveRDS(sim_res1, file="/Volumes/PVD2/Davis/PhD/DISSERTATION/Chapter4/R4/PRRSTM/newset3/sim_res1V75_F40_75.rds")
#saveRDS(sim_res2, file="/Volumes/PVD2/Davis/PhD/DISSERTATION/Chapter4/R4/PRRSTM/newset3/sim_res2V75_F40_75.rds")

# Only 100% vaccine --------------------------------------------------------
sim_res1 = vector("list", nsims)
sim_res2 = vector("list", nsims)
names(sim_res1) = 1:nsims
names(sim_res2) = 1:nsims
result_1 = vector("list",27)
result_2 = vector("list",27)
names(result_1) = c(paste0("beta_l_",1:9), paste0("beta_m_",1:9),paste0("beta_h_",1:9) )
names(result_2) = c(paste0("beta_l_",1:9), paste0("beta_m_",1:9),paste0("beta_h_",1:9) )

for( k in 1:(nsims)) {
  loop=0
  
  for(i in 1:3){
    beta = beta_listV100[[i]]
    
    for(j in 1:9){
      print(paste0("k = ", k, "; i = ", i, "; j = ", j))
      loop = loop+1
      
      rho = rho_list[[j]]
      beta_matrix = rho * beta 
      
      stocks <- c(S=Farms$X, I=Farms$Y, R=Farms$Z)
      parameters <- list("beta"=beta_matrix, "delays"=c(D100[[k]]), 
                         "mort"=c(v[[k]]), "returns"=c(D_imm[[k]])) #remove farms.  
      
      out <- data.frame(ode(y=stocks, times=time, func=model, parms=parameters, 
                            method = "euler")) 
      out <- out[out$time %in% seq(0,26,1), ]
      
      o <-data.frame(time=out$time)
      o$S <- apply(out[,c(grep('S', names(out), value=TRUE))], 1, sum)
      o$I <- apply(out[,c(grep('I', names(out), value=TRUE))], 1, sum)
      o$R <- apply(out[,c(grep('R', names(out), value=TRUE))], 1, sum)
      o$animals <- o$S+o$R+o$I
      o$prop <- o$I/o$animals
      
      result_1[[loop]] = out
      result_2[[loop]] = o
    }
  }
  sim_res1[[k]] = result_1
  sim_res2[[k]] = result_2
}
#saveRDS(sim_res1, file="/Volumes/PVD2/Davis/PhD/DISSERTATION/Chapter4/R4/PRRSTM/newset3/sim_res1V100.rds")
#saveRDS(sim_res2, file="/Volumes/PVD2/Davis/PhD/DISSERTATION/Chapter4/R4/PRRSTM/newset3/sim_res2V100.rds")

# 100% vaccine 100% filter (80) ---------------------------------------------
sim_res1 = vector("list", nsims)
sim_res2 = vector("list", nsims)
names(sim_res1) = 1:nsims
names(sim_res2) = 1:nsims
result_1 = vector("list",27)
result_2 = vector("list",27)
names(result_1) = c(paste0("beta_l_",1:9), paste0("beta_m_",1:9),paste0("beta_h_",1:9) )
names(result_2) = c(paste0("beta_l_",1:9), paste0("beta_m_",1:9),paste0("beta_h_",1:9) )

for( k in 1:(nsims)) {
  loop=0
  
  for(i in 1:3){
    beta = beta_listV100[[i]]
    
    for(j in 1:9){
      print(paste0("k = ", k, "; i = ", i, "; j = ", j))
      loop = loop+1
      
      rho = rho_list_F80_100[[j]]
      beta_matrix = rho * beta 
      
      stocks <- c(S=Farms$X, I=Farms$Y, R=Farms$Z)
      parameters <- list("beta"=beta_matrix, "delays"=c(D100[[k]]), 
                         "mort"=c(v[[k]]), "returns"=c(D_imm[[k]])) #remove farms.  
      
      out <- data.frame(ode(y=stocks, times=time, func=model, parms=parameters, 
                            method = "euler")) 
      out <- out[out$time %in% seq(0,26,1), ]
      
      o <-data.frame(time=out$time)
      o$S <- apply(out[,c(grep('S', names(out), value=TRUE))], 1, sum)
      o$I <- apply(out[,c(grep('I', names(out), value=TRUE))], 1, sum)
      o$R <- apply(out[,c(grep('R', names(out), value=TRUE))], 1, sum)
      o$animals <- o$S+o$R+o$I
      o$prop <- o$I/o$animals
      
      result_1[[loop]] = out
      result_2[[loop]] = o
    }
  }
  sim_res1[[k]] = result_1
  sim_res2[[k]] = result_2
}
#saveRDS(sim_res1, file="/Volumes/PVD2/Davis/PhD/DISSERTATION/Chapter4/R4/PRRSTM/newset3/sim_res1V100_F80_100.rds")
#saveRDS(sim_res2, file="/Volumes/PVD2/Davis/PhD/DISSERTATION/Chapter4/R4/PRRSTM/newset3/sim_res2V100_F80_100.rds")

# 100% vaccine 100% filter (40) ---------------------------------------------
sim_res1 = vector("list", nsims)
sim_res2 = vector("list", nsims)
names(sim_res1) = 1:nsims
names(sim_res2) = 1:nsims
result_1 = vector("list",27)
result_2 = vector("list",27)
names(result_1) = c(paste0("beta_l_",1:9), paste0("beta_m_",1:9),paste0("beta_h_",1:9) )
names(result_2) = c(paste0("beta_l_",1:9), paste0("beta_m_",1:9),paste0("beta_h_",1:9) )

for( k in 1:(nsims)) {
  loop=0
  
  for(i in 1:3){
    beta = beta_listV100[[i]]
    
    for(j in 1:9){
      print(paste0("k = ", k, "; i = ", i, "; j = ", j))
      loop = loop+1
      
      rho = rho_list_F40_100[[j]]
      beta_matrix = rho * beta 
      
      stocks <- c(S=Farms$X, I=Farms$Y, R=Farms$Z)
      parameters <- list("beta"=beta_matrix, "delays"=c(D100[[k]]), 
                         "mort"=c(v[[k]]), "returns"=c(D_imm[[k]])) #remove farms.  
      
      out <- data.frame(ode(y=stocks, times=time, func=model, parms=parameters, 
                            method = "euler")) 
      out <- out[out$time %in% seq(0,26,1), ]
      
      o <-data.frame(time=out$time)
      o$S <- apply(out[,c(grep('S', names(out), value=TRUE))], 1, sum)
      o$I <- apply(out[,c(grep('I', names(out), value=TRUE))], 1, sum)
      o$R <- apply(out[,c(grep('R', names(out), value=TRUE))], 1, sum)
      o$animals <- o$S+o$R+o$I
      o$prop <- o$I/o$animals
      
      result_1[[loop]] = out
      result_2[[loop]] = o
    }
  }
  sim_res1[[k]] = result_1
  sim_res2[[k]] = result_2
}
#saveRDS(sim_res1, file="/Volumes/PVD2/Davis/PhD/DISSERTATION/Chapter4/R4/PRRSTM/newset3/sim_res1V100_F40_100.rds")
#saveRDS(sim_res2, file="/Volumes/PVD2/Davis/PhD/DISSERTATION/Chapter4/R4/PRRSTM/newset3/sim_res2V100_F40_100.rds")

# 25% filter (80) ---------------------------------------------
sim_res1 = vector("list", nsims)
sim_res2 = vector("list", nsims)
names(sim_res1) = 1:nsims
names(sim_res2) = 1:nsims
result_1 = vector("list",27)
result_2 = vector("list",27)
names(result_1) = c(paste0("beta_l_",1:9), paste0("beta_m_",1:9),paste0("beta_h_",1:9) )
names(result_2) = c(paste0("beta_l_",1:9), paste0("beta_m_",1:9),paste0("beta_h_",1:9) )

for( k in 1:(nsims)) {
  loop=0
  
  for(i in 1:3){
    beta = beta_list[[i]]
    
    for(j in 1:9){
      print(paste0("k = ", k, "; i = ", i, "; j = ", j))
      loop = loop+1
      
      rho = rho_list_F80_25[[j]]
      beta_matrix = rho * beta 
      
      stocks <- c(S=Farms$X, I=Farms$Y, R=Farms$Z)
      parameters <- list("beta"=beta_matrix, "delays"=c(D[[k]]), 
                         "mort"=c(v[[k]]), "returns"=c(D_imm[[k]])) #remove farms.  
      
      out <- data.frame(ode(y=stocks, times=time, func=model, parms=parameters, 
                            method = "euler")) 
      out <- out[out$time %in% seq(0,26,1), ]
      
      o <-data.frame(time=out$time)
      o$S <- apply(out[,c(grep('S', names(out), value=TRUE))], 1, sum)
      o$I <- apply(out[,c(grep('I', names(out), value=TRUE))], 1, sum)
      o$R <- apply(out[,c(grep('R', names(out), value=TRUE))], 1, sum)
      o$animals <- o$S+o$R+o$I
      o$prop <- o$I/o$animals
      
      result_1[[loop]] = out
      result_2[[loop]] = o
    }
  }
  sim_res1[[k]] = result_1
  sim_res2[[k]] = result_2
}
#saveRDS(sim_res1, file="/Volumes/PVD2/Davis/PhD/DISSERTATION/Chapter4/R4/PRRSTM/newset3/sim_res1_F80_25.rds")
#saveRDS(sim_res2, file="/Volumes/PVD2/Davis/PhD/DISSERTATION/Chapter4/R4/PRRSTM/newset3/sim_res2_F80_25.rds")

# 50% filter (80) ---------------------------------------------
sim_res1 = vector("list", nsims)
sim_res2 = vector("list", nsims)
names(sim_res1) = 1:nsims
names(sim_res2) = 1:nsims
result_1 = vector("list",27)
result_2 = vector("list",27)
names(result_1) = c(paste0("beta_l_",1:9), paste0("beta_m_",1:9),paste0("beta_h_",1:9) )
names(result_2) = c(paste0("beta_l_",1:9), paste0("beta_m_",1:9),paste0("beta_h_",1:9) )

for( k in 1:(nsims)) {
  loop=0
  
  for(i in 1:3){
    beta = beta_list[[i]]
    
    for(j in 1:9){
      print(paste0("k = ", k, "; i = ", i, "; j = ", j))
      loop = loop+1
      
      rho = rho_list_F80_50[[j]]
      beta_matrix = rho * beta 
      
      stocks <- c(S=Farms$X, I=Farms$Y, R=Farms$Z)
      parameters <- list("beta"=beta_matrix, "delays"=c(D[[k]]), 
                         "mort"=c(v[[k]]), "returns"=c(D_imm[[k]])) #remove farms.  
      
      out <- data.frame(ode(y=stocks, times=time, func=model, parms=parameters, 
                            method = "euler")) 
      out <- out[out$time %in% seq(0,26,1), ]
      
      o <-data.frame(time=out$time)
      o$S <- apply(out[,c(grep('S', names(out), value=TRUE))], 1, sum)
      o$I <- apply(out[,c(grep('I', names(out), value=TRUE))], 1, sum)
      o$R <- apply(out[,c(grep('R', names(out), value=TRUE))], 1, sum)
      o$animals <- o$S+o$R+o$I
      o$prop <- o$I/o$animals
      
      result_1[[loop]] = out
      result_2[[loop]] = o
    }
  }
  sim_res1[[k]] = result_1
  sim_res2[[k]] = result_2
}
#saveRDS(sim_res1, file="/Volumes/PVD2/Davis/PhD/DISSERTATION/Chapter4/R4/PRRSTM/newset3/sim_res1_F80_50.rds")
#saveRDS(sim_res2, file="/Volumes/PVD2/Davis/PhD/DISSERTATION/Chapter4/R4/PRRSTM/newset3/sim_res2_F80_50.rds")

# 75% filter (80) ---------------------------------------------
sim_res1 = vector("list", nsims)
sim_res2 = vector("list", nsims)
names(sim_res1) = 1:nsims
names(sim_res2) = 1:nsims
result_1 = vector("list",27)
result_2 = vector("list",27)
names(result_1) = c(paste0("beta_l_",1:9), paste0("beta_m_",1:9),paste0("beta_h_",1:9) )
names(result_2) = c(paste0("beta_l_",1:9), paste0("beta_m_",1:9),paste0("beta_h_",1:9) )

for( k in 1:(nsims)) {
  loop=0
  
  for(i in 1:3){
    beta = beta_list[[i]]
    
    for(j in 1:9){
      print(paste0("k = ", k, "; i = ", i, "; j = ", j))
      loop = loop+1
      
      rho = rho_list_F80_75[[j]]
      beta_matrix = rho * beta 
      
      stocks <- c(S=Farms$X, I=Farms$Y, R=Farms$Z)
      parameters <- list("beta"=beta_matrix, "delays"=c(D[[k]]), 
                         "mort"=c(v[[k]]), "returns"=c(D_imm[[k]])) #remove farms.  
      
      out <- data.frame(ode(y=stocks, times=time, func=model, parms=parameters, 
                            method = "euler")) 
      out <- out[out$time %in% seq(0,26,1), ]
      
      o <-data.frame(time=out$time)
      o$S <- apply(out[,c(grep('S', names(out), value=TRUE))], 1, sum)
      o$I <- apply(out[,c(grep('I', names(out), value=TRUE))], 1, sum)
      o$R <- apply(out[,c(grep('R', names(out), value=TRUE))], 1, sum)
      o$animals <- o$S+o$R+o$I
      o$prop <- o$I/o$animals
      
      result_1[[loop]] = out
      result_2[[loop]] = o
    }
  }
  sim_res1[[k]] = result_1
  sim_res2[[k]] = result_2
}
#saveRDS(sim_res1, file="/Volumes/PVD2/Davis/PhD/DISSERTATION/Chapter4/R4/PRRSTM/newset3/sim_res1_F80_75.rds")
#saveRDS(sim_res2, file="/Volumes/PVD2/Davis/PhD/DISSERTATION/Chapter4/R4/PRRSTM/newset3/sim_res2_F80_75.rds")

# 100% filter (80) ---------------------------------------------
sim_res1 = vector("list", nsims)
sim_res2 = vector("list", nsims)
names(sim_res1) = 1:nsims
names(sim_res2) = 1:nsims
result_1 = vector("list",27)
result_2 = vector("list",27)
names(result_1) = c(paste0("beta_l_",1:9), paste0("beta_m_",1:9),paste0("beta_h_",1:9) )
names(result_2) = c(paste0("beta_l_",1:9), paste0("beta_m_",1:9),paste0("beta_h_",1:9) )

for( k in 1:(nsims)) {
  loop=0
  
  for(i in 1:3){
    beta = beta_list[[i]]
    
    for(j in 1:9){
      print(paste0("k = ", k, "; i = ", i, "; j = ", j))
      loop = loop+1
      
      rho = rho_list_F80_100[[j]]
      beta_matrix = rho * beta 
      
      stocks <- c(S=Farms$X, I=Farms$Y, R=Farms$Z)
      parameters <- list("beta"=beta_matrix, "delays"=c(D[[k]]), 
                         "mort"=c(v[[k]]), "returns"=c(D_imm[[k]])) #remove farms.  
      
      out <- data.frame(ode(y=stocks, times=time, func=model, parms=parameters, 
                            method = "euler")) 
      out <- out[out$time %in% seq(0,26,1), ]
      
      o <-data.frame(time=out$time)
      o$S <- apply(out[,c(grep('S', names(out), value=TRUE))], 1, sum)
      o$I <- apply(out[,c(grep('I', names(out), value=TRUE))], 1, sum)
      o$R <- apply(out[,c(grep('R', names(out), value=TRUE))], 1, sum)
      o$animals <- o$S+o$R+o$I
      o$prop <- o$I/o$animals
      
      result_1[[loop]] = out
      result_2[[loop]] = o
    }
  }
  sim_res1[[k]] = result_1
  sim_res2[[k]] = result_2
}
#saveRDS(sim_res1, file="/Volumes/PVD2/Davis/PhD/DISSERTATION/Chapter4/R4/PRRSTM/newset3/sim_res1_F80_100.rds")
#saveRDS(sim_res2, file="/Volumes/PVD2/Davis/PhD/DISSERTATION/Chapter4/R4/PRRSTM/newset3/sim_res2_F80_100.rds")

# 25% filter (40) ---------------------------------------------
sim_res1 = vector("list", nsims)
sim_res2 = vector("list", nsims)
names(sim_res1) = 1:nsims
names(sim_res2) = 1:nsims
result_1 = vector("list",27)
result_2 = vector("list",27)
names(result_1) = c(paste0("beta_l_",1:9), paste0("beta_m_",1:9),paste0("beta_h_",1:9) )
names(result_2) = c(paste0("beta_l_",1:9), paste0("beta_m_",1:9),paste0("beta_h_",1:9) )

for( k in 1:(nsims)) {
  loop=0
  
  for(i in 1:3){
    beta = beta_list[[i]]
    
    for(j in 1:9){
      print(paste0("k = ", k, "; i = ", i, "; j = ", j))
      loop = loop+1
      
      rho = rho_list_F40_25[[j]]
      beta_matrix = rho * beta 
      
      stocks <- c(S=Farms$X, I=Farms$Y, R=Farms$Z)
      parameters <- list("beta"=beta_matrix, "delays"=c(D[[k]]), 
                         "mort"=c(v[[k]]), "returns"=c(D_imm[[k]])) #remove farms.  
      
      out <- data.frame(ode(y=stocks, times=time, func=model, parms=parameters, 
                            method = "euler")) 
      out <- out[out$time %in% seq(0,26,1), ]
      
      o <-data.frame(time=out$time)
      o$S <- apply(out[,c(grep('S', names(out), value=TRUE))], 1, sum)
      o$I <- apply(out[,c(grep('I', names(out), value=TRUE))], 1, sum)
      o$R <- apply(out[,c(grep('R', names(out), value=TRUE))], 1, sum)
      o$animals <- o$S+o$R+o$I
      o$prop <- o$I/o$animals
      
      result_1[[loop]] = out
      result_2[[loop]] = o
    }
  }
  sim_res1[[k]] = result_1
  sim_res2[[k]] = result_2
}
#saveRDS(sim_res1, file="/Volumes/PVD2/Davis/PhD/DISSERTATION/Chapter4/R4/PRRSTM/newset3/sim_res1_F40_25.rds")
#saveRDS(sim_res2, file="/Volumes/PVD2/Davis/PhD/DISSERTATION/Chapter4/R4/PRRSTM/newset3/sim_res2_F40_25.rds")

# 50% filter (80) ---------------------------------------------
sim_res1 = vector("list", nsims)
sim_res2 = vector("list", nsims)
names(sim_res1) = 1:nsims
names(sim_res2) = 1:nsims
result_1 = vector("list",27)
result_2 = vector("list",27)
names(result_1) = c(paste0("beta_l_",1:9), paste0("beta_m_",1:9),paste0("beta_h_",1:9) )
names(result_2) = c(paste0("beta_l_",1:9), paste0("beta_m_",1:9),paste0("beta_h_",1:9) )

for( k in 1:(nsims)) {
  loop=0
  
  for(i in 1:3){
    beta = beta_list[[i]]
    
    for(j in 1:9){
      print(paste0("k = ", k, "; i = ", i, "; j = ", j))
      loop = loop+1
      
      rho = rho_list_F40_50[[j]]
      beta_matrix = rho * beta 
      
      stocks <- c(S=Farms$X, I=Farms$Y, R=Farms$Z)
      parameters <- list("beta"=beta_matrix, "delays"=c(D[[k]]), 
                         "mort"=c(v[[k]]), "returns"=c(D_imm[[k]])) #remove farms.  
      
      out <- data.frame(ode(y=stocks, times=time, func=model, parms=parameters, 
                            method = "euler")) 
      out <- out[out$time %in% seq(0,26,1), ]
      
      o <-data.frame(time=out$time)
      o$S <- apply(out[,c(grep('S', names(out), value=TRUE))], 1, sum)
      o$I <- apply(out[,c(grep('I', names(out), value=TRUE))], 1, sum)
      o$R <- apply(out[,c(grep('R', names(out), value=TRUE))], 1, sum)
      o$animals <- o$S+o$R+o$I
      o$prop <- o$I/o$animals
      
      result_1[[loop]] = out
      result_2[[loop]] = o
    }
  }
  sim_res1[[k]] = result_1
  sim_res2[[k]] = result_2
}
#saveRDS(sim_res1, file="/Volumes/PVD2/Davis/PhD/DISSERTATION/Chapter4/R4/PRRSTM/newset3/sim_res1_F40_50.rds")
#saveRDS(sim_res2, file="/Volumes/PVD2/Davis/PhD/DISSERTATION/Chapter4/R4/PRRSTM/newset3/sim_res2_F40_50.rds")

# 75% filter (80) ---------------------------------------------
sim_res1 = vector("list", nsims)
sim_res2 = vector("list", nsims)
names(sim_res1) = 1:nsims
names(sim_res2) = 1:nsims
result_1 = vector("list",27)
result_2 = vector("list",27)
names(result_1) = c(paste0("beta_l_",1:9), paste0("beta_m_",1:9),paste0("beta_h_",1:9) )
names(result_2) = c(paste0("beta_l_",1:9), paste0("beta_m_",1:9),paste0("beta_h_",1:9) )

for( k in 1:(nsims)) {
  loop=0
  
  for(i in 1:3){
    beta = beta_list[[i]]
    
    for(j in 1:9){
      print(paste0("k = ", k, "; i = ", i, "; j = ", j))
      loop = loop+1
      
      rho = rho_list_F40_75[[j]]
      beta_matrix = rho * beta 
      
      stocks <- c(S=Farms$X, I=Farms$Y, R=Farms$Z)
      parameters <- list("beta"=beta_matrix, "delays"=c(D[[k]]), 
                         "mort"=c(v[[k]]), "returns"=c(D_imm[[k]])) #remove farms.  
      
      out <- data.frame(ode(y=stocks, times=time, func=model, parms=parameters, 
                            method = "euler")) 
      out <- out[out$time %in% seq(0,26,1), ]
      
      o <-data.frame(time=out$time)
      o$S <- apply(out[,c(grep('S', names(out), value=TRUE))], 1, sum)
      o$I <- apply(out[,c(grep('I', names(out), value=TRUE))], 1, sum)
      o$R <- apply(out[,c(grep('R', names(out), value=TRUE))], 1, sum)
      o$animals <- o$S+o$R+o$I
      o$prop <- o$I/o$animals
      
      result_1[[loop]] = out
      result_2[[loop]] = o
    }
  }
  sim_res1[[k]] = result_1
  sim_res2[[k]] = result_2
}
#saveRDS(sim_res1, file="/Volumes/PVD2/Davis/PhD/DISSERTATION/Chapter4/R4/PRRSTM/newset3/sim_res1_F40_75.rds")
#saveRDS(sim_res2, file="/Volumes/PVD2/Davis/PhD/DISSERTATION/Chapter4/R4/PRRSTM/newset3/sim_res2_F40_75.rds")

# 100% filter (80) ---------------------------------------------
sim_res1 = vector("list", nsims)
sim_res2 = vector("list", nsims)
names(sim_res1) = 1:nsims
names(sim_res2) = 1:nsims
result_1 = vector("list",27)
result_2 = vector("list",27)
names(result_1) = c(paste0("beta_l_",1:9), paste0("beta_m_",1:9),paste0("beta_h_",1:9) )
names(result_2) = c(paste0("beta_l_",1:9), paste0("beta_m_",1:9),paste0("beta_h_",1:9) )

for( k in 1:(nsims)) {
  loop=0
  
  for(i in 1:3){
    beta = beta_list[[i]]
    
    for(j in 1:9){
      print(paste0("k = ", k, "; i = ", i, "; j = ", j))
      loop = loop+1
      
      rho = rho_list_F40_100[[j]]
      beta_matrix = rho * beta 
      
      stocks <- c(S=Farms$X, I=Farms$Y, R=Farms$Z)
      parameters <- list("beta"=beta_matrix, "delays"=c(D[[k]]), 
                         "mort"=c(v[[k]]), "returns"=c(D_imm[[k]])) #remove farms.  
      
      out <- data.frame(ode(y=stocks, times=time, func=model, parms=parameters, 
                            method = "euler")) 
      out <- out[out$time %in% seq(0,26,1), ]
      
      o <-data.frame(time=out$time)
      o$S <- apply(out[,c(grep('S', names(out), value=TRUE))], 1, sum)
      o$I <- apply(out[,c(grep('I', names(out), value=TRUE))], 1, sum)
      o$R <- apply(out[,c(grep('R', names(out), value=TRUE))], 1, sum)
      o$animals <- o$S+o$R+o$I
      o$prop <- o$I/o$animals
      
      result_1[[loop]] = out
      result_2[[loop]] = o
    }
  }
  sim_res1[[k]] = result_1
  sim_res2[[k]] = result_2
}
#saveRDS(sim_res1, file="/Volumes/PVD2/Davis/PhD/DISSERTATION/Chapter4/R4/PRRSTM/newset3/sim_res1_F40_100.rds")
#saveRDS(sim_res2, file="/Volumes/PVD2/Davis/PhD/DISSERTATION/Chapter4/R4/PRRSTM/newset3/sim_res2_F40_100.rds")

# Baseline ---------------------------------------------
sim_res1 = vector("list", nsims)
sim_res2 = vector("list", nsims)
names(sim_res1) = 1:nsims
names(sim_res2) = 1:nsims
result_1 = vector("list",27)
result_2 = vector("list",27)
names(result_1) = c(paste0("beta_l_",1:9), paste0("beta_m_",1:9),paste0("beta_h_",1:9) )
names(result_2) = c(paste0("beta_l_",1:9), paste0("beta_m_",1:9),paste0("beta_h_",1:9) )

for( k in 1:(nsims)) {
  loop=0
  
  for(i in 1:3){
    beta = beta_list[[i]]
    
    for(j in 1:9){
      print(paste0("k = ", k, "; i = ", i, "; j = ", j))
      loop = loop+1
      
      rho = rho_list[[j]]
      beta_matrix = rho * beta 
      
      stocks <- c(S=Farms$X, I=Farms$Y, R=Farms$Z)
      parameters <- list("beta"=beta_matrix, "delays"=c(D[[k]]), 
                         "mort"=c(v[[k]]), "returns"=c(D_imm[[k]])) #remove farms.  
      
      out <- data.frame(ode(y=stocks, times=time, func=model, parms=parameters, 
                            method = "euler")) 
      out <- out[out$time %in% seq(0,26,1), ]
      
      o <-data.frame(time=out$time)
      o$S <- apply(out[,c(grep('S', names(out), value=TRUE))], 1, sum)
      o$I <- apply(out[,c(grep('I', names(out), value=TRUE))], 1, sum)
      o$R <- apply(out[,c(grep('R', names(out), value=TRUE))], 1, sum)
      o$animals <- o$S+o$R+o$I
      o$prop <- o$I/o$animals
      
      result_1[[loop]] = out
      result_2[[loop]] = o
    }
  }
  sim_res1[[k]] = result_1
  sim_res2[[k]] = result_2
}
#saveRDS(sim_res1, file="/Volumes/PVD2/Davis/PhD/DISSERTATION/Chapter4/R4/PRRSTM/newset3/sim_res1.rds")
#saveRDS(sim_res2, file="/Volumes/PVD2/Davis/PhD/DISSERTATION/Chapter4/R4/PRRSTM/newset3/sim_res2.rds")







