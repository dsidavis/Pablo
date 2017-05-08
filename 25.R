nsims = 1
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