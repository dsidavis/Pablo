source("setup.R")
source("functions.R")
source("model5.R")

Rprof("prof")
system.time(doSim(beta_list, rho_list, D, time, v, D_imm, Farms, nsims = 1, verbose = TRUE))
Rprof(NULL)
summaryRprof('prof')$by.self


system.time(doSim(beta_listV25, rho_list_F80_25, D25, time, v, D_imm, Farms, nsims = 1, verbose = TRUE))


