

y <- as.vector(read.table("simulation_experiments_scripts/section_6/SW/USdata_updated.txt")$x)

T <- length(y)

#cc <- stanc_builder(file="sw_drhmc_tauc.stan",allow_undefined=TRUE, isystem=CIP_header_path())
#drhmc.mdlt <- stan_model(stanc_ret=cc,allow_undefined=TRUE,includes=CIP_include())


standta <- list(T=T,y=y,alpha=5.0,beta=0.5,zmean=0.0*y[1:(T-1)],xmean=0*y,lambdamean=0.0,lambdaprec=5+T-1.5)


#drhmc <- sampling(drhmc.mdlt,data=standta,chains=10,seed=2,init=inits,control=list(stepsize_jitter=0.5))

stan.fit <- rstan::stan("simulation_experiments_scripts/section_6/SW/sw_drhmc_tauc.stan",data=standta,iter=10000)



