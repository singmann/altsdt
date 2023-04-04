
library("tidyverse")
library("MPTinR")

load("EmpiricalData/MWWData/MWW2007_preprocesseddata.rds")
# mww_wide <- MWW2007 %>%
#   group_by(exp, id, oldnew, response) %>%
#   count() %>%
#   ungroup() %>%
#   arrange(oldnew, n) %>%
#   pivot_wider(id_cols = c(exp, id), names_from = c(oldnew, response), values_from = n, values_fill = 0) %>%
#   select(exp, id, Old_1, Old_2, Old_3, Old_4, Old_5, Old_6,  New_1, New_2, New_3, New_4, New_5, New_6)


uvsd <- "
pnorm(cr1, mu, sigma)
pnorm(cr1+cr2, mu, sigma) - pnorm(cr1, mu, sigma)
pnorm(cr3+cr2+cr1, mu, sigma) - pnorm(cr2+cr1, mu, sigma)
pnorm(cr4+cr3+cr2+cr1, mu, sigma) - pnorm(cr3+cr2+cr1, mu, sigma)
pnorm(cr5+cr4+cr3+cr2+cr1, mu, sigma) - pnorm(cr4+cr3+cr2+cr1, mu, sigma)
1 - pnorm(cr5+cr4+cr3+cr2+cr1, mu, sigma)

pnorm(cr1)
pnorm(cr2+cr1) - pnorm(cr1)
pnorm(cr3+cr2+cr1) - pnorm(cr2+cr1)
pnorm(cr4+cr3+cr2+cr1) - pnorm(cr3+cr2+cr1)
pnorm(cr5+cr4+cr3+cr2+cr1) - pnorm(cr4+cr3+cr2+cr1)
1 - pnorm(cr5+cr4+cr3+cr2+cr1)
"

check.mpt(textConnection(uvsd))

uvsd_fit <- fit.model(mww2007_6point[,-(1:2)], textConnection(uvsd),
            lower.bound=c(-Inf, rep(0, 5), 0.001), upper.bound=Inf, n.optim = 20, use.gradient = FALSE)
uvsd_fit$goodness.of.fit

uvsd_fit$goodness.of.fit$aggregated

uvsd_fit_agg <- fit.model(colSums(mww2007_6point[,-(1:2)]), textConnection(uvsd),
            lower.bound=c(-Inf, rep(0, 5), 0.001), upper.bound=Inf, n.optim = 20, use.gradient = FALSE)
uvsd_fit_agg$goodness.of.fit

uvsd_fit_rest <- fit.model(mww2007_6point[,-(1:2)], textConnection(uvsd),
            lower.bound=c(-Inf, rep(0, 2), 0.001), upper.bound=Inf,
            restrictions.filename = list("cr2 = cr3 = cr4 = cr5"),
            n.optim = 20)
uvsd_fit_rest$goodness.of.fit

#plot(uvsd_fit$goodness.of.fit$individual$G.Squared, uvsd_fit_rest$goodness.of.fit$individual$G.Squared)

hist(uvsd_fit$goodness.of.fit$individual$G.Squared - uvsd_fit_rest$goodness.of.fit$individual$G.Squared)

mean(uvsd_fit$goodness.of.fit$individual$G.Squared - uvsd_fit_rest$goodness.of.fit$individual$G.Squared < -qchisq(0.95,4))

gaussian_uvsdt_opt <- function(par, data, param.names, n.params, tmp.env, lower.bound, upper.bound,
                               nconf = 6){
  d     <- as.numeric(par[1])
  sigo  <- as.numeric(par[2])
  # start_crit <- as.numeric(par[3])
  # crit_dist <- as.numeric(par[4])
  c     <- cumsum(as.numeric(c(par[3:n.params])))
  sign  <- 1
  I <- c(-Inf,c,Inf) # put criteria into larger array
  #I <- c(-Inf,rep(start_crit, nconf-1) + seq(from = 0, to = nconf-1, by = 1)*crit_dist,Inf) # put criteria into larger array

  # Likelihood of every trial
  # New items
  pNlikJ <- vector("numeric", nconf)
  for (i in 1:nconf){
    pNlikJ[i] <- pnorm(I[i+1],mean=0,sd=sign)-pnorm(I[i],mean=0,sd=sign)
  }

  # Old items
  pOlikJ <- vector("numeric", nconf)
  for (i in 1:nconf){

    pOlikJ[i] <- pnorm(I[i+1],mean=d,sd=sigo)-pnorm(I[i],mean=d,sd=sigo)
  }
  e <- c(pOlikJ, pNlikJ)

  llk <- -sum(data[data!=0]*log(e[data!=0]))

  if (is.na(llk)) llk <- 1e10
  if (llk == Inf) llk <- 1e10

  return(llk)

}

gaussian_uvsdt_opt_rest <- function(par, data, param.names, n.params, tmp.env, lower.bound, upper.bound,
                               nconf = 6){
  d     <- as.numeric(par[1])
  sigo  <- as.numeric(par[2])
  start_crit <- as.numeric(par[3])
  crit_dist <- as.numeric(par[4])
  #c     <- cumsum(as.numeric(c(par[3:n.params])))
  sign  <- 1
  I <- c(-Inf,rep(start_crit, nconf-1) + seq(from = 0, to = nconf-2, by = 1)*crit_dist,Inf) # put criteria into larger array

  # Likelihood of every trial
  # New items
  pNlikJ <- vector("numeric", nconf)
  for (i in 1:nconf){
    pNlikJ[i] <- pnorm(I[i+1],mean=0,sd=sign)-pnorm(I[i],mean=0,sd=sign)
  }

  # Old items
  pOlikJ <- vector("numeric", nconf)
  for (i in 1:nconf){

    pOlikJ[i] <- pnorm(I[i+1],mean=d,sd=sigo)-pnorm(I[i],mean=d,sd=sigo)
  }
  e <- c(pOlikJ, pNlikJ)

  llk <- -sum( data[data!=0]*log(e[data!=0]))

  if (is.na(llk)) llk <- 1e10
  if (llk == Inf) llk <- 1e10

  return(llk)

}

gaussian_uvsdt_pred<- function(par, data, param.names, n.params, tmp.env, lower.bound, upper.bound,
                               nconf = 6){
  d     <- as.numeric(par[1])
  sigo  <- as.numeric(par[2])
  # start_crit <- as.numeric(par[3])
  # crit_dist <- as.numeric(par[4])
  c     <- cumsum(as.numeric(c(par[3:n.params])))
  sign  <- 1
  I <- c(-Inf,c,Inf) # put criteria into larger array
  #I <- c(-Inf,rep(start_crit, nconf-1) + seq(from = 0, to = nconf-1, by = 1)*crit_dist,Inf) # put criteria into larger array

  # Likelihood of every trial
  # New items
  pNlikJ <- vector("numeric", nconf)
  for (i in 1:(length(c)+1)){
    pNlikJ[i] <- pnorm(I[i+1],mean=0,sd=sign)-pnorm(I[i],mean=0,sd=sign)
  }

  # Old items
  pOlikJ <- vector("numeric", nconf)
  for (i in 1:(length(c)+1)){

    pOlikJ[i] <- pnorm(I[i+1],mean=d,sd=sigo)-pnorm(I[i],mean=d,sd=sigo)
  }

  return(c(pOlikJ, pNlikJ))
}


uvsd_fit2 <- fit.mptinr(mww2007_6point[,-(1:2)], objective = gaussian_uvsdt_opt, prediction = gaussian_uvsdt_pred,
                        param.names = c("mu", "sigma", paste0("c", 1:5)), categories.per.type = c(6, 6),
            lower.bound=c(0, 0.0001, -Inf, rep(0, 4)), upper.bound=Inf, n.optim = 20)

uvsd_fit2$goodness.of.fit$individual
all.equal(uvsd_fit$goodness.of.fit, uvsd_fit2$goodness.of.fit)

uvsd_fit2_rest <- fit.mptinr(mww2007_6point[,-(1:2)], objective = gaussian_uvsdt_opt_rest,
                        param.names = c("mu", "sigma", "cstart", "cincr"), categories.per.type = c(6, 6),
            lower.bound=c(0, 0.0001, -Inf, 0), upper.bound=Inf, n.optim = 20)

all.equal(uvsd_fit_rest$goodness.of.fit, uvsd_fit2_rest$goodness.of.fit)

uvsd_fit_rest$goodness.of.fit$individual
uvsd_fit2_rest$goodness.of.fit$individual

uvsd_fit_rest$parameters$individual[,,1]

uvsd_fit2_rest$parameters$individual[,1,]

uvsd_20scale <- fit.mptinr(mww2007_20point[,-(1:2)], objective = gaussian_uvsdt_opt,
                        param.names = c("mu", "sigma", paste0("c", 1:19)), categories.per.type = c(20, 20),
            lower.bound=c(0, 0.0001, -Inf, rep(0, 18)), upper.bound=Inf, n.optim = 20,
            starting.values = list(
              c(0.3, 1, -1, rep(0, 18)),
              c(1, 2, 0, rep(0.1, 18))
            ), nconf = 20)

uvsd_20scale$goodness.of.fit$individual

uvsd_20scale_rest <- fit.mptinr(mww2007_20point[,-(1:2)], objective = gaussian_uvsdt_opt_rest,
                        param.names = c("mu", "sigma", "cstart", "cincr"), categories.per.type = c(20, 20),
            lower.bound=c(0, 0.0001, -Inf, 0), upper.bound=Inf, n.optim = 20,
            starting.values = list(
              c(0.3, 1, -1, 0),
              c(1, 2, 0, 1)
            ), nconf = 20)

uvsd_20scale_rest$goodness.of.fit$individual$G.Squared

uvsd_20scale_rest$goodness.of.fit$individual$G.Squared - uvsd_20scale$goodness.of.fit$individual$G.Squared

qchisq(0.95, 17)

######### 99 scale #####
uvsd_99scale <- fit.mptinr(mww2007_99point[,-(1:2)], objective = gaussian_uvsdt_opt,
                        param.names = c("mu", "sigma", paste0("c", 1:98)), categories.per.type = c(99, 99),
            lower.bound=c(0, 0.0001, -Inf, rep(0, 97)), upper.bound=Inf, n.optim = 20,
            starting.values = list(
              c(0.3, 1, -1, rep(0, 97)),
              c(1, 2, 0, rep(0.1, 97))
            ), nconf = 99)

uvsd_99scale$goodness.of.fit

uvsd_99scale_rest <- fit.mptinr(mww2007_99point[,-(1:2)], objective = gaussian_uvsdt_opt_rest,
                        param.names = c("mu", "sigma", "cstart", "cincr"), categories.per.type = c(99, 99),
            lower.bound=c(0, 0.0001, -Inf, 0), upper.bound=Inf, n.optim = 20,
            starting.values = list(
              c(0.3, 1, -1, 0),
              c(1, 2, 0, 0.1)
            ), nconf = 99)

diff_99 <- uvsd_99scale_rest$goodness.of.fit$individual$G.Squared - uvsd_99scale$goodness.of.fit$individual$G.Squared

diff_99
min(diff_99)
qchisq(0.95, 96)
## all rejected with p < .0001
pchisq(min(diff_99), 0.05, 96, lower.tail = FALSE)

#####