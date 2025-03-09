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
                                  prior = "betacorrelation",
                                  param = c(5,1)),
                       prec = list(initial = 1,
                                   fixed = FALSE,
                                   prior = "loggamma",
                                   param=c(1,1)
                       )
                     ))

result = inla(formula, data = data, family="gaussian",
              control.family = list(hyper = list(
                prec = list(
                  fixed = TRUE,
                  initial = log(1 / s^2),
                  piror = "loggamma",
                  param = c(100,1)
                )
              )),
              control.inla = list(int.strategy = "grid"))

result.hyperpar = inla.hyperpar(result, diff.logdens=20)

summary(result)
summary(result.hyperpar)

names(result$marginals.hyperpar)

# Support bounds
limit = 0.001

# HP1 Precision
theta.inla.hyperpar1.list = inla.smarginal(result.hyperpar$internal.marginals.hyperpar[[1]])
theta.inla.hyperpar1.mat = cbind(theta.inla.hyperpar1.list$x, theta.inla.hyperpar1.list$y)
theta.inla.hyperpar1 = subset(theta.inla.hyperpar1.mat, theta.inla.hyperpar1.mat[,2] >= limit)
inla_hyperpar1 = result.hyperpar$marginals.hyperpar[[1]]

# HP2 Rho
theta.inla.hyperpar2.list = inla.smarginal(result.hyperpar$internal.marginals.hyperpar[[2]])
theta.inla.hyperpar2.mat = cbind(theta.inla.hyperpar2.list$x, theta.inla.hyperpar2.list$y)
theta.inla.hyperpar2 = subset(theta.inla.hyperpar2.mat, theta.inla.hyperpar2.mat[,2] >= limit)
inla_hyperpar2 = result.hyperpar$marginals.hyperpar[[2]]

# HP Details
theta.names = names(result$marginals.hyperpar)
len.theta = length(theta.names)
est.theta = result$mode$theta # Mode estimates from inla
cov.intern = result$misc$cov.intern # Covariance matrix estimate from inla


# compute the reference
marglik.estimator = 2
marglik.ref = inla.evaluate(est.theta, result)$mlik[marglik.estimator] + 
  (dgamma(exp(est.theta[1]), shape = 1, rate = 1, log=TRUE)+ est.theta[1]) +
  (scaledBeta(((exp(est.theta[2])-1)/(exp(est.theta[2])+1)),5,1,log=TRUE)+
     log(abs((2*exp(est.theta[2]))/((exp(est.theta[2])+1)^2))))

#### LDS
npoints = 64

# Generating constant from LatticeBuilder (a = 27)
korobov.non = shift.koro(npoints, len.theta, 27, FALSE)
x11(width=30,height=30); pairs(korobov.non, pch=16)
# vs 19, which seems like the best choice so far
korobov.non = shift.koro(npoints, len.theta, 19, FALSE)
x11(width=30,height=30); pairs(korobov.non, pch=16)

# Low and high limits
low = c(min(theta.inla.hyperpar1[,1]), min(theta.inla.hyperpar2[,1]))
high = c(max(theta.inla.hyperpar1[,1]), max(theta.inla.hyperpar2[,1]))

# Transform pointset from unit hypercube to support 
korobov = trans.pset(korobov.non, low, high)

#Sobol sequence
sobol.seq = sobol(npoints, len.theta)
x11(30,30); pairs(sobol.seq, pch=16)
sobol.pset = trans.pset(sobol.seq, low, high)
#Halton
halton.seq = halton(npoints, len.theta)
x11(30,30); pairs(halton.seq, pch=16)
halton.pset = trans.pset(halton.seq, low, high)

# Matrix that will store the function evaluations and the points 
#eval.results = matrix(NA, npoints, len.theta+1)

######CONVERT THIS TO PROCESS IN PARALLEL
# Compute function evaluations and store in eval.results
# system.time({
#   for (i in 1:npoints){
#     theta.new = korobov[i,]
#     result.new = inla.evaluate(theta.new, result)
#     logpost.rel = result.new$mlik[marglik.estimator] - marglik.ref +
#       (dgamma(exp(theta.new[1]), shape = 100, rate = 1, log=TRUE)+ theta.new[1]) +
#       (scaledBeta(((exp(theta.new[2])-1)/(exp(theta.new[2])+1)),5,1,log=TRUE)+
#          log(abs((2*exp(theta.new[2]))/((exp(theta.new[2])+1)^2))))
#     
#     cat(i, "Evaluate theta.new", theta.new,"with relative posterior",logpost.rel,"\n")
#     eval.results[i, ] = c(theta.new, logpost.rel)
#   }
# })


comp_function_eval <- function(x_i,i){
  #return(sum(x_i))
  theta.new = x_i
  result.new = inla.evaluate(theta.new, result)
  logpost.rel = result.new$mlik[marglik.estimator] - marglik.ref +
    (dgamma(exp(theta.new[1]), shape = 1, rate = 1, log=TRUE)+ theta.new[1]) +
    (scaledBeta(((exp(theta.new[2])-1)/(exp(theta.new[2])+1)),5,1,log=TRUE)+
       log(abs((2*exp(theta.new[2]))/((exp(theta.new[2])+1)^2))))
  return(c(theta.new, logpost.rel))
  # cat(i, "Evaluate theta.new", theta.new,"with relative posterior",logpost.rel,"\n")
  # return(c(theta.new, logpost.rel))
}

task_function_eval <- function(x){
  result.new <- matrix(NA, nrow(x), ncol(x)+1)
  for(i in 1:nrow(x)){
    result.new[i,]<- comp_function_eval(x[i,],i)
  }
  return(result.new)
}
n_dfs <- 3
rows_per_split <- ceiling(nrow(korobov) / n_dfs)
split_df_list <- vector("list", n_dfs)
for (i in 1:n_dfs) {
  start_row <- (i - 1) * rows_per_split + 1
  end_row <- min(i * rows_per_split, nrow(korobov))
  split_df_list[[i]] <- korobov[start_row:end_row, ]
}

ptset_list <- split_df_list

cl <- makeCluster(16)
clusterEvalQ(cl, {
  source("My working/ThesisFunctions.R")
})
clusterExport(cl
              , list("task_function_eval"
                     ,"ptset_list"
                     ,"comp_function_eval"
                     ,"inla.evaluate"
                     ,"result"
                     ,"rho_z","marglik.estimator","marglik.ref","scaledBeta"))
system.time({
  result_list <- parLapply(cl, ptset_list, task_function_eval)
})


stopCluster(cl)

final_Eval_res <- do.call(rbind, result_list)
final_Eval_res
######

# Quadratic Approximation
# partitions
parts = 5
position = 2
poly.degree = 2

# Finding marginals when function evaluations are in log-scale
psi1 = final_Eval_res[,c(1, len.theta+1)]
psi1 = psi1[order(psi1[,1]),]
x11();plot(psi1)
psi2 = final_Eval_res[,c(2, len.theta+1)]
psi2 = psi2[order(psi2[,1]),]
x11();plot(psi2)

# Partition and find pointwise means
psi1.parts.exp = partition(cbind(psi1[,1],exp(psi1[,2])),c(low[1], high[1]),parts,position)
psi2.parts.exp = partition(cbind(psi2[,1],exp(psi2[,2])),c(low[2], high[2]),parts,position)

# Co-ordinates for pointwise means
psi1.parts = cbind(psi1.parts.exp[,1], log(psi1.parts.exp[,2]))
psi2.parts = cbind(psi2.parts.exp[,1], log(psi2.parts.exp[,2]))

# Data matrix
van.partM1 = vandermonde.matrix(psi1.parts[,1],poly.degree)
van.partM2 = vandermonde.matrix(psi2.parts[,1],poly.degree)

# X'X
gramian.part1 = t(van.partM1)%*%van.partM1
gramian.part2 = t(van.partM2)%*%van.partM2

# Least squares approximation to quadratic coefficients
marg1.beta.part = qr.solve(gramian.part1, tol = 1e-32)%*%t(van.partM1)%*%psi1.parts[,2]
marg2.beta.part = qr.solve(gramian.part2, tol = 1e-32)%*%t(van.partM2)%*%psi2.parts[,2]

# Create a new set of values over the range for polynomial
l.o = 1000
pt.part1 = seq(low[1], high[1], length.out = l.o)
pt.part2 = seq(low[2], high[2], length.out = l.o)

# This is the quadratic approximation to the log(true) marginal
marg1.part.pred = vandermonde.matrix(pt.part1, poly.degree)%*%marg1.beta.part
marg2.part.pred = vandermonde.matrix(pt.part2, poly.degree)%*%marg2.beta.part

theta.poly.part1 = cbind(pt.part1, marg1.part.pred)
theta.poly.part1 = theta.poly.part1[order(theta.poly.part1[,1]),]
theta.poly.part2 = cbind(pt.part2, marg2.part.pred)
theta.poly.part2 = theta.poly.part2[order(theta.poly.part2[,1]),]

# Partition lines for graphs
lines.part1 = seq(low[1], high[1], length.out = parts + 1)
lines.part2 = seq(low[2], high[2], length.out = parts + 1)

### Plots for function evaluations, partitions, pointwise means and quadratic approximation
x11(30,10); par(mfrow=c(1,2))
plot(psi1, pch=16, main="function evaluations", col="light blue")
points(psi1.parts, pch=16, col="red")
lines(theta.poly.part1, col="black",lwd=1.5)
abline(v=lines.part1, col="green")

plot(psi2, pch=16, main="function evaluations", col="light blue")
points(psi2.parts, pch=16, col="red")
lines(theta.poly.part2, col="black",lwd=1.5)
abline(v=lines.part2, col="green")

### Plots above, without function evaluations - Can see difference between polynomial and pointwise means
### For next step - treat these as "errors/residuals", to help fit a correction..
x11(30,10); par(mfrow=c(1,2))
plot(psi1.parts, pch=16, col="red")
lines(theta.poly.part1, col="black", lwd=1.5)
abline(v = lines.part1, col = "green")

plot(psi2.parts, pch=16, col="red")
lines(theta.poly.part2, col="black", lwd=1.5)
abline(v = lines.part2, col = "green")

### Proper theta approximations --> HP approximations
theta.lds.part1 = cbind(pt.part1, exp(marg1.part.pred))
theta.ppn.const1 = diff(range(theta.lds.part1[,1]))*mean(theta.lds.part1[,2])
theta.lds.part1 = cbind(theta.lds.part1[,1], theta.lds.part1[,2]/theta.ppn.const1)

lds.part1 = cbind(exp(theta.lds.part1[,1]), theta.lds.part1[,2])
hp.ppn.const1 = diff(range(lds.part1[,1]))*mean(lds.part1[,2])
lds.part1 = cbind(lds.part1[,1], lds.part1[,2]/hp.ppn.const1)

theta.lds.part2 = cbind(pt.part2, exp(marg2.part.pred))
theta.ppn.const2 = diff(range(theta.lds.part2[,1]))*mean(theta.lds.part2[,2])
theta.lds.part2 = cbind(theta.lds.part2[,1], theta.lds.part2[,2]/theta.ppn.const2)

lds.part2 = cbind(exp(theta.lds.part2[,1]), theta.lds.part2[,2])
hp.ppn.const2 = diff(range(lds.part2[,1]))*mean(lds.part2[,2])
lds.part2 = cbind(lds.part2[,1], lds.part2[,2]/hp.ppn.const2)

#### Some plotting stuff ###
max.theta.ylim1 = max(c(theta.inla.hyperpar1[,2], theta.lds.part1[,2]))
max.theta.ylim2 = max(c(theta.inla.hyperpar2[,2], theta.lds.part2[,2]))

max.hp.ylim1 = max(c(inla_hyperpar1[,2], lds.part1[,2]))
max.hp.ylim2 = max(c(inla_hyperpar2[,2], lds.part2[,2]))

### Approximations - Our quadratic approximation vs INLA's Dense Grid
x11(20,10); par(mfrow = c(1,3))
plot(theta.inla.hyperpar1, type="l", lwd = 2, ylim = c(0,max.theta.ylim1*1.05), cex.lab = 1.25,
     xlab = "theta_1_prec")
lines(theta.lds.part1, type="l", col = "red", lwd = 2)
abline(v = 1, col="green")

plot(theta.inla.hyperpar2, type="l", lwd = 2, ylim = c(0,max.theta.ylim2*1.05), cex.lab = 1.25,
     xlab = "theta_2_rho")
lines(theta.lds.part2, type="l", col = "red", lwd = 2)
#abline(v = log((1+0.65)/(1-0.65)), col="green")
abline(v = rho_z, col="green")

plot(0,0,col="white")
text(0,0, cex = 1.5,
     "THETAS \n Black is inla.hyperpar \n Red is LDS-polynomial \n Green is the true value")

########################
### Cubic Correction ###
########################

### Polynomial fit at abscissa points
quad.poly1 = van.partM1%*%marg1.beta.part
quad.poly2 = van.partM2%*%marg2.beta.part

### residuals = poly - pointwise means
resids1 = quad.poly1 - psi1.parts[,2]
resids2 = quad.poly2 - psi2.parts[,2]

ppdeg = 3

### get cubic coefficients using least squares through residuals
### Vandermonde matrix here is cubic..
vR1 = vandermonde.matrix(psi1.parts[,1], ppdeg)
vR2 = vandermonde.matrix(psi2.parts[,1], ppdeg)

another.beta1 = qr.solve(t(vR1)%*%vR1, tol = 1e-32)%*%t(vR1)%*%resids1
another.beta2 = qr.solve(t(vR2)%*%vR2, tol = 1e-32)%*%t(vR2)%*%resids2

corrected.beta1 = rbind(marg1.beta.part,0) - another.beta1
corrected.beta2 = rbind(marg2.beta.part,0) - another.beta2

## Cubic correction
marg1.corr.pred = vandermonde.matrix(pt.part1, ppdeg)%*%corrected.beta1
marg2.corr.pred = vandermonde.matrix(pt.part2, ppdeg)%*%corrected.beta2

l.o = 1000

pt.part1 = seq(low[1], high[1], length.out = l.o)
pt.part2 = seq(low[2], high[2], length.out = l.o)

theta.poly.corr1 = cbind(pt.part1, marg1.corr.pred)
theta.poly.corr1 = theta.poly.corr1[order(theta.poly.corr1[,1]),]
theta.poly.corr2 = cbind(pt.part2, marg2.corr.pred)
theta.poly.corr2 = theta.poly.corr2[order(theta.poly.corr2[,1]),]

lines.part1 = seq(low[1], high[1], length.out = parts + 1)
lines.part2 = seq(low[2], high[2], length.out = parts + 1)

### Plots of function evals, partitions, pointwise means and cubic correction
x11(20,10); par(mfrow=c(1,2))
plot(psi1, pch=16, main = "function evals", col="light blue")
points(psi1.parts, pch=16, col="red")
lines(theta.poly.corr1, col="black", lwd=1.5)
abline(v = lines.part1, col = "green")

plot(psi2, pch=16, main = "function evals", col="light blue")
points(psi2.parts, pch=16, col="red")
lines(theta.poly.corr2, col="black", lwd=1.5)
abline(v = lines.part2, col = "green")

### Cubic correction vs Initial (quadratic) approximation
x11(20,10); par(mfrow=c(1,2))
plot(psi1.parts, pch=16, col="red", cex = 1.5, xlab="", ylab="", main = "Theta(Tau_y))", cex.main=2)
lines(theta.poly.part1, col="black", lwd=2)
lines(theta.poly.corr1, col='blue', lty=2, lwd=2)
legend("topleft",
       c("Initial", "Cubic-Cor"), col = c("black","blue"), lwd=c(2,2), lty=c(1,2), cex=1.5)

plot(psi2.parts, pch=16, col="red", cex = 1.5, xlab="", ylab="", main = "Theta(Rho))", cex.main=2)
lines(theta.poly.part2, col="black", lwd=2)
lines(theta.poly.corr2, col='blue', lty=2, lwd=2)
legend("topleft",
       c("Initial", "Cubic-Cor"), col = c("black","blue"), lwd=c(2,2), lty=c(1,2), cex=1.5)

## Proper correction plots 
theta.lds.corr1 = cbind(pt.part1, exp(marg1.corr.pred))
theta.corr.const1 = diff(range(theta.lds.corr1[,1]))*mean(theta.lds.corr1[,2])
theta.lds.corr1 = cbind(theta.lds.corr1[,1], theta.lds.corr1[,2]/theta.corr.const1)

lds.corr1 = cbind(exp(theta.lds.corr1[,1]), theta.lds.corr1[,2])
hp.corr.const1 = diff(range(lds.corr1[,1]))*mean(lds.corr1[,2])
lds.corr1 = cbind(lds.corr1[,1], lds.corr1[,2]/hp.corr.const1)

theta.lds.corr2 = cbind(pt.part2, exp(marg2.corr.pred))
theta.corr.const2 = diff(range(theta.lds.corr2[,1]))*mean(theta.lds.corr2[,2])
theta.lds.corr2 = cbind(theta.lds.corr2[,1], theta.lds.corr2[,2]/theta.corr.const2)

lds.corr2 = cbind(exp(theta.lds.corr2[,1]), theta.lds.corr2[,2])
hp.corr.const2 = diff(range(lds.corr2[,1]))*mean(lds.corr2[,2])
lds.corr2 = cbind(lds.corr2[,1], lds.corr2[,2]/hp.corr.const2)

#### Some plotting stuff ###
max.theta.ylim1 = max(c(theta.inla.hyperpar1[,2], theta.lds.part1[,2]))
max.theta.ylim2 = max(c(theta.inla.hyperpar2[,2], theta.lds.part2[,2]))

max.hp.ylim1 = max(c(inla_hyperpar1[,2], lds.part1[,2]))
max.hp.ylim2 = max(c(inla_hyperpar2[,2], lds.part2[,2]))

### Plots - Thetas and hyperpars
## INLA dense grid vs Initial quadratic approx vs cubic correction
x11(35,20); par(mfrow = c(2,2))
plot(theta.inla.hyperpar1, type="l", lwd = 2, ylim = c(0,max.theta.ylim1*1.05), cex.lab = 1.25,
     xlab = "", ylab="", main = "Theta(Tau_y)", cex.main = 2)
lines(theta.lds.part1, type="l", col = "red", lwd = 2)
lines(theta.lds.corr1, type="l", col = "blue", lty =2, lwd = 2)
abline(v = log(1 / s^2), col="green")
legend("topleft",
       c("INLA-grid", "Initial", "C-Corr"), lty = c(1,1,2), lwd = c(2,2,2), 
       col = c("black","red", "blue"), cex = 1.5)

plot(theta.inla.hyperpar2, type="l", lwd = 2, ylim = c(0,max.theta.ylim2*1.05), cex.lab = 1.25,
     xlab = "", ylab="", main = "Theta(Rho)", cex.main = 2)
lines(theta.lds.part2, type="l", col = "red", lwd = 2)
lines(theta.lds.corr2, type="l", col = "blue", lty =2, lwd = 2)
abline(v = rho_z, col="green")
legend("topleft",
       c("INLA-grid", "Initial", "C-Corr"), lty = c(1,1,2), lwd = c(2,2,2), 
       col = c("black","red", "blue"), cex = 1.5)

plot(inla.smarginal(result.hyperpar$marginals.hyperpar[[1]]), 
     type="l", lwd = 2, cex.lab = 1.25, xlim = c(min(lds.part1[,1]), max(lds.part1[,1])),
     xlab = "", ylab="", main = "Tau_y", cex.main = 2)
lines(lds.part1, type="l", col = "red", lwd = 2)
lines(lds.corr1, type="l", col = "blue", lty =2, lwd = 2)
abline(v = 100, col="green")
legend("topleft",
       c("INLA-grid", "Initial", "C-Corr"), lty = c(1,1,2), lwd = c(2,2,2), 
       col = c("black","red", "blue"), cex = 1.5)

plot(inla.smarginal(result.hyperpar$marginals.hyperpar[[2]]), 
     type="l", lwd = 2, cex.lab = 1.25, xlim = c(min(lds.part2[,1]), max(lds.part2[,1])),
     ylim = c(min(lds.part2[1,]), max(lds.part2[1,])),
     xlab = "", ylab="", main = "Rho", cex.main = 2)
lines(lds.part2, type="l", col = "red", lwd = 2)
lines(lds.corr2, type="l", col = "blue", lty =2, lwd = 2)
abline(v = rho, col="green")
legend("topright",
       c("INLA-grid", "Initial", "C-Corr"), lty = c(1,1,2), lwd = c(2,2,2), 
       col = c("black","red", "blue"), cex = 1.5)


#########################
#### LATENT ESTIMATION
#########################

result_lat <- inla(formula, data = data, family="gaussian",
                   control.family = list(hyper = list(
                     prec = list(
                       fixed = FALSE,
                       initial = log(1 / s^2),
                       prior = "loggamma",
                       param=c(100,1)
                     ))),
                   control.compute = list(config=TRUE),
                   control.inla = list(int.strategy = "grid"))
grid.points <- result_lat$joint.hyper
dim(grid.points)

a.min1 = min(grid.points[,1])
a.min2 = min(grid.points[,2])

b.max1 = max(grid.points[,1])
b.max2 = max(grid.points[,2])

low = c(a.min1, a.min2)
high = c(b.max1, b.max2)

#Halton
halt <- halton(npoints, len.theta)
trans.halt <- trans.pset(halt, low, high)
#Sobol
sobol.seq <- sobol(npoints, len.theta)
trans.sobol <- trans.pset(sobol.seq, low, high)
#Korobov
korobov.non = shift.koro(npoints, len.theta, 23, FALSE)
korobov <- trans.pset(korobov.non, low, high)

result.korobov <- inla(formula, data = data, family="gaussian",
                       control.family = list(hyper = list(
                         prec = list(
                           fixed = FALSE,
                           initial = log(1 / s^2),
                           prior = "loggamma",
                           param=c(100,1)
                         ))),
                       control.compute = list(config=TRUE),
                       control.inla = list(int.strategy = "user",
                                           int.design = cbind(korobov,1),
                                           reordering = "metis"))
result.halt <- inla(formula, data = data, family="gaussian",
                    control.family = list(hyper = list(
                      prec = list(
                        fixed = FALSE,
                        initial = log(1 / s^2),
                        prior = "loggamma",
                        param=c(100,1)
                      ))),
                    control.compute = list(config=TRUE),
                    control.inla = list(int.strategy = "user",
                                        int.design = cbind(trans.halt,1),
                                        reordering = "metis"))
result.sobol <- inla(formula, data = data, family="gaussian",
                     control.family = list(hyper = list(
                       prec = list(
                         fixed = FALSE,
                         initial = log(1 / s^2),
                         prior = "loggamma",
                         param=c(100,1)
                       ))),
                     control.compute = list(config=TRUE),
                     control.inla = list(int.strategy = "user",
                                         int.design = cbind(trans.sobol,1),
                                         reordering = "metis"))
summary(result_lat)
summary(result.korobov)
summary(result.halt)
summary(result.sobol)

pairs(result_lat$joint.hyper[,c(1:2)], pch=16)
pairs(result.korobov$joint.hyper[,c(1:2)], pch=16)
pairs(result.halt$joint.hyper[,c(1:2)], pch=16)
pairs(result.sobol$joint.hyper[,c(1:2)], pch=16)

# Latent posteriors


# Observed
x11(width=20,height=20)
plot(data$idx,data$y, main=paste("LDS npoints = ",npoints))
lines(data$idx,data$y)
# Korobov
post_samples <- inla.posterior.sample(1,result.korobov) 
kor_latent_post <- post_samples[[1]]$latent[,1]
lines(data$idx,head(kor_latent_post,100))
# Halton
post_samples <- inla.posterior.sample(1,result.halt) 
halt_latent_post <- post_samples[[1]]$latent[,1]
lines(data$idx,head(halt_latent_post,100),col="red")
# Sobol
post_samples <- inla.posterior.sample(1,result.sobol) 
sob_latent_post <- post_samples[[1]]$latent[,1]
lines(data$idx,head(sob_latent_post,100),col="purple")

