require(pdphmc)
require(rstan)
l <- read.table("simulation_experiments_scripts/german/set2.dat")

y <- as.vector(l[,1])
X <- as.matrix(l[,3:26])
glmfit<- glm.fit(x=l[,2:26],y=y,family = binomial())  # for reference only
glmfit_new <- glm(V1 ~ . - 1, data = l, family = "binomial")
n <- dim(X)[1]
p <- dim(X)[2]

Tmax <- 5000.0

fit.stan <- stan("logistic.stan",data=list(n=n,p=p,y=as.integer(y),X=X),chains=10)

mdl <- build("logistic.cpp",lambda="lambda_constant",massMatrix = "diagMassVARI")
fit.c <- run(mdl,data=list(y=as.integer(y),X=X),Tmax=Tmax,chains=10,
             control=list(absTol=1.0e-3,relTol=1.0e-3,lambdaAdapt=c(1.0,1.0,5.0)))
fit.cl <- run(mdl,data=list(y=as.integer(y),X=X),Tmax=Tmax,chains=10,
              control=list(absTol=1.0e-3,relTol=1.0e-3,lambdaAdapt=c(1.0,1.0,20.0)))

mdl.a <- build("logistic.cpp",lambda="lambda_arclength",massMatrix = "diagMassVARI")
fit.a <- run(mdl.a,data=list(y=as.integer(y),X=X),Tmax=Tmax,chains=10,
             control=list(absTol=1.0e-3,relTol=1.0e-3,lambdaAdapt=c(1.0,1.0,5.0)))
fit.al <- run(mdl.a,data=list(y=as.integer(y),X=X),Tmax=Tmax,chains=10,
              control=list(absTol=1.0e-3,relTol=1.0e-3,lambdaAdapt=c(1.0,1.0,20.0)))

save.image("Computations")
