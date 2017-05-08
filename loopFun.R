#
# Took the first loop. Compared it to all of the others to identify what variables changed.
# Then took that first for loop code and put it in 25.R
# Copied it to here.
# Then modified this version.
# Start with no parameters.
# After creating the function, use findGlobals() to find missing variables which need to be parameters
# findGlobals(doSim, FALSE)
# [1] "beta_listV25" "D_imm"        "D25"          "Farms"        "model"       
# [6] "nsims"        "paste0"       "rho_list"     "sum"          "time"        
#[11] "v"   
# paste0 is not a function, or sum, or model o

doSim = 
function(beta_list, rho_list, D, time, v, D_imm, Farms, nsims = 2, verbose = TRUE)
{

   sim_res1 = sim_res2 = structure(vector("list", nsims), names = 1:nsims)

   result_1 = result_2 = structure(vector("list",27) , 
                                   names = apply(expand.grid(sprintf("beta_%s_", c("l", "m", "h")), 1:9), 1, paste0, collapse = ""))

   stocks <- c(S = Farms$X, I = Farms$Y, R = Farms$Z)
   IDX = matrix(2:(3*nrow(Farms) + 1),, 3)

for(k in 1:nsims) {
  loop=0
 
    # doesn't depend on i or j
  parameters <- list("beta" = beta_matrix, "delays"=D[[k]], 
                     "mort"= v[[k]], "returns"=D_imm[[k]]) #remove farms.  
  
  for(i in 1:3){
    beta = beta_list[[i]]
    
    for(j in 1:9){
      if(verbose)
         cat(paste0("k = ", k, "; i = ", i, "; j = ", j), "\n")

      loop = loop + 1
      
      rho = rho_list[[j]]
      beta_matrix = rho * beta 
     
      tmp = ode(y = stocks, times = time, func = model, parms = parameters, method = "euler")
        # Can this be just a simple tmp[,1] < 26
      tmp = tmp[ tmp[, "time"] %in% seq(0,26,1), ]
      out <- data.frame(tmp) 
      
      o <- data.frame(time = tmp[, 1]) # "time"
         # It appears (verify) that the S, I and R rows can be determined directly
         # from the nrow(Farms). There are three sets of equal length for each
         # and the S starts at 2 since time is first.
      o$S <- rowSums(tmp[, IDX[,1]])
      o$I <- rowSums(tmp[, IDX[,2]])
      o$R <- rowSums(tmp[, IDX[,3]])

      o$animals <- o$S + o$R + o$I
      o$prop <- o$I/o$animals
      
      result_1[[loop]] = out  # what about leaving as a matrix?
      result_2[[loop]] = o
    }
  }
  sim_res1[[k]] = result_1
  sim_res2[[k]] = result_2

 }
 return(list(sim_res1 = sim_res1, sim_res2 = sim_res2))
}
#saveRDS(sim_res1, file="/Volumes/PVD2/Davis/PhD/DISSERTATION/Chapter4/R4/PRRSTM/newset3/sim_res1V25.rds")
#saveRDS(sim_res2, file="/Volumes/PVD2/Davis/PhD/DISSERTATION/Chapter4/R4/PRRSTM/newset3/sim_res2V25.rds")