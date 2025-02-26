### AR1" with INLA LDS
source("My working/ThesisFunctions.R") # <--All packages and custom functions stored here


# Model type: Autoregressive model with lag 1
# The INLA model is fit on the AR(1) model where we estimate 2 hyperparameters
# HP 1: Precision of the obeservations
# HP 2: Rho

### MODEL PARAMETERS
set.seed(1242)
n = 100
rho = 0.65
rho_z = (log(1+rho)/(1-rho)) # Fishers z transformation
s = 0.1
x = arima.sim(n, model=list(ar=rho))
x=scale(x)
y = x + rnorm(n, sd = s)

data = data.frame(y=y, x=x, idx=1:n)

formula = y ~ -1 + f(idx, model = "ar1",
                     hyper = list(
                       rho = list(initial = rho_z,
                                  fixed = FALSE,
                                  pior = "betacorrelation",
                                  param = c(5,1)),
                       prec = list(initial = 1,
                                   fixed = TRUE,
                                   prior = "loggamma",
                                   param=c(1,1)
                       )
                     ))

result = inla(formula, data = data, family="gaussian",
              control.family = list(hyper = list(
                prec = list(
                  fixed = FALSE,
                  initial = log(1 / s^2),
                  piror = "loggamma",
                  param = c(100,1)
                )
              )),
              control.inla = list(int.strategy = "grid"))

result.hyperpar = inla.hyperpar(result, diff.logdens=20)
