#### AR1 Example
### Approximate hyperparameter posterior marginals using:
### initial approximation (LDS-QA), then
### update with a Cubic Correction

source("My working/ThesisFunctions.R") # <--All packages and custom functions stored here


### AR1 model
### 3 hyperparameters, tau, kappa, rho
n = 100
s = 0.1
set.seed(781984); x = arima.sim(n, model = list(ar = 0.65))
x = scale(x)
y = x + rnorm(n, sd = s)
data = data.frame(y, x, idx = 1:n)
r = inla(y ~ -1 + f(idx, model="ar1",
                    ## this is 0.65 in the internal scale
                    hyper = list(rho = list(initial = 1.550597, fixed=FALSE,
                                            prior = "betacorrelation",
                                            param = c(5, 1)),
                                 prec = list(prior = "loggamma",
                                             param = c(1, 1)))), 
         data = data,
         family = "gaussian",
         control.family = list(hyper = list(
           prec = list(
             fixed = FALSE,
             initial = log(1/s^2),
             prior = "loggamma",
             param = c(100, 1)))), 
         control.inla = list(int.strategy = "eb", reordering = "metis"))

r2 = inla.hyperpar(r, diff.logdens = 20)

### SUPPORT BOUNDS ###
### Using INLA approximations

### 1) Fit inla.spline to inla theta approximation
### 2) Take a marginal support where a_i and b_i contain the bulk of the density
### 3) construct box out of these "bounds"

### INLAs (inla.hyperpar) approximations to theta and HP for support bounds!
limit = 0.001

### INLA Theta/HP 1
theta.inla.hyperpar1.list = inla.smarginal(r2$internal.marginals.hyperpar[[1]])
theta.inla.hyperpar1.mat = cbind(theta.inla.hyperpar1.list$x, theta.inla.hyperpar1.list$y)
theta.inla.hyperpar1 = subset(theta.inla.hyperpar1.mat, theta.inla.hyperpar1.mat[,2] >= limit)
inla_hyperpar1 = r2$marginals.hyperpar[[1]]

### INLA Theta/HP2
theta.inla.hyperpar2.list = inla.smarginal(r2$internal.marginals.hyperpar[[2]])
theta.inla.hyperpar2.mat = cbind(theta.inla.hyperpar2.list$x, theta.inla.hyperpar2.list$y)
theta.inla.hyperpar2 = subset(theta.inla.hyperpar2.mat, theta.inla.hyperpar2.mat[,2] >= limit)
inla_hyperpar2 = r2$marginals.hyperpar[[2]]

### INLA Theta/HP3
theta.inla.hyperpar3.list = inla.smarginal(r2$internal.marginals.hyperpar[[3]])
theta.inla.hyperpar3.mat = cbind(theta.inla.hyperpar3.list$x, theta.inla.hyperpar3.list$y)
theta.inla.hyperpar3 = subset(theta.inla.hyperpar3.mat, theta.inla.hyperpar3.mat[,2] >= limit)
inla_hyperpar3 = r2$marginals.hyperpar[[3]]

### Details
theta.names = names(r$marginals.hyperpar)
len.theta = length(theta.names)
est.theta = r$mode$theta   # Mode value(s) estimate (from inla function)
cov.intern = r$misc$cov.intern # Covariance matrix estimate from inla.

cat("We have ", len.theta, "hyperparmaeters\n")
for(nm in theta.names) {
  cat("\t- ", nm, "\n")
}
cat("Estimates are ")
cat(est.theta, "\n")
cat("With Covariance matrix\n")
print(cov.intern)

### compute the reference (without the priors on the hyperparameters)
mlik.estimator = 2
mlik.ref = inla.evaluate(r$mode$theta, r)$mlik[mlik.estimator] + 
  (dgamma(exp(r$mode$theta[1]), shape = 100, rate = 1, log=TRUE) + r$mode$theta[1]) +
  (dgamma(exp(r$mode$theta[2]), shape = 1, rate = 1, log=TRUE) + r$mode$theta[2]) +
  (scaledBeta(((exp(r$mode$theta[3])-1)/(exp(r$mode$theta[3])+1)), 5, 1, log=TRUE) + 
     log(abs((2*exp(r$mode$theta[3]))/((exp(r$mode$theta[3])+1)^2))))

### Generate Korobov point set and calculate function evaluations
### Need - shift.koro(...), inla.evaluate(...)
### Number of points to evaluate function
npoints = 512

### Generating constant from LatticeBuilder (a = 19)
korobov.non = shift.koro(npoints, len.theta, 19, FALSE)
x11(50,50); pairs(korobov.non, pch=16)

#Low and High based on INLA "where the bulk lies" limits
low = c(min(theta.inla.hyperpar1[,1]), min(theta.inla.hyperpar2[,1]), min(theta.inla.hyperpar3[,1]))
high = c(max(theta.inla.hyperpar1[,1]), max(theta.inla.hyperpar2[,1]), max(theta.inla.hyperpar3[,1]))

### Transform pointset from unit hypercube to support
korobov = trans.pset(korobov.non, low, high)
res1 = matrix(NA, npoints, len.theta+1)

### Compute function evaluations
for(i in 1:npoints) {
  theta.new = korobov[i,]
  ## compute the relative posterior based on the Laplace approximation. ## Since we are rerunning
  ## with fixed=TRUE (above), the priors for the hyperparameters are 
  ## not included, so we need to add those.
  r.new = inla.evaluate(theta.new, r)
  lpost.rel = r.new$mlik[mlik.estimator] - mlik.ref +
    (dgamma(exp(theta.new[1]), shape = 100, rate = 1, log=TRUE) + theta.new[1]) +
    (dgamma(exp(theta.new[2]), shape = 1, rate = 1, log=TRUE) + theta.new[2]) +
    (scaledBeta(((exp(theta.new[3])-1)/(exp(theta.new[3])+1)), 5, 1, log=TRUE) + 
       log(abs((2*exp(theta.new[3]))/((exp(theta.new[3])+1)^2))))
  
  ## lpost.rel are the function evaluations (log) for the posterior 
  ## adding the log hyperpriors and subtracting the log marginal likelihood ## computed before.
  cat(i, "Evaluate theta.new", theta.new, "with relative posterior",
      lpost.rel, "\n")
  res1[i, ] = c(theta.new, lpost.rel)
  ## res1 contains augmented matrix (points and function evaluations)
}

comp_function_eval <- function(x_i,i){
  #return(sum(x_i))
  theta.new = x_i
  result.new = inla.evaluate(theta.new, r)
  logpost.rel = result.new$mlik[marglik.estimator] - marglik.ref +
    (dgamma(exp(theta.new[1]), shape = 100, rate = 1, log=TRUE)+ theta.new[1]) +
    (dgamma(exp(theta.new[2]), shape = 1, rate = 1, log=TRUE) + theta.new[2]) +
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
n_dfs <- 4
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
                     ,"r"
                     ,"marglik.estimator","marglik.ref","scaledBeta"))
system.time({
  result_list <- parLapply(cl, ptset_list, task_function_eval)
})


stopCluster(cl)

final_Eval_res <- do.call(rbind, result_list)
final_Eval_res

#########################################
### Initial Approximation (Quadratic) ###
#########################################
### Number of partitions, partition position and degree (2)
parts = 15
position = 2
poly.degree = 2

### Finding marginals when function evaluations are in log-scale...
psi1 = final_Eval_res[,c(1,len.theta+1)]
psi1 = psi1[order(psi1[,1]),]
psi2 = final_Eval_res[,c(2,len.theta+1)]
psi2 = psi2[order(psi2[,1]),]
psi3 = final_Eval_res[,c(3,len.theta+1)]
psi3 = psi3[order(psi3[,1]),]

### Partition and find pointwise means
psi1.parts.exp = partition(cbind(psi1[,1], exp(psi1[,2])), c(low[1], high[1]), parts, position)
psi2.parts.exp = partition(cbind(psi2[,1], exp(psi2[,2])), c(low[2], high[2]), parts, position)
psi3.parts.exp = partition(cbind(psi3[,1], exp(psi3[,2])), c(low[3], high[3]), parts, position)

### Co-ordinates for pointwise means
psi1.parts = cbind(psi1.parts.exp[,1], log(psi1.parts.exp[,2]))
psi2.parts = cbind(psi2.parts.exp[,1], log(psi2.parts.exp[,2]))
psi3.parts = cbind(psi3.parts.exp[,1], log(psi3.parts.exp[,2]))

## Data matrix
van.partM1 = vandermonde.matrix(psi1.parts[,1], poly.degree)
van.partM2 = vandermonde.matrix(psi2.parts[,1], poly.degree)
van.partM3 = vandermonde.matrix(psi3.parts[,1], poly.degree)

## X'X
gramian.part1 = t(van.partM1)%*%van.partM1
gramian.part2 = t(van.partM2)%*%van.partM2
gramian.part3 = t(van.partM3)%*%van.partM3

## Least squares approximation to quadratic coefficients
marg1.beta.part = qr.solve(gramian.part1, tol = 1e-32)%*%t(van.partM1)%*%psi1.parts[,2]
marg2.beta.part = qr.solve(gramian.part2, tol = 1e-32)%*%t(van.partM2)%*%psi2.parts[,2]
marg3.beta.part = qr.solve(gramian.part3, tol = 1e-32)%*%t(van.partM3)%*%psi3.parts[,2]

## Create a new set of values over the range for polynomial
l.o = 1000
pt.part1 = seq(low[1], high[1], length.out = l.o)
pt.part2 = seq(low[2], high[2], length.out = l.o)
pt.part3 = seq(low[3], high[3], length.out = l.o)

## This is the quadratic approximation to the log(true) marginal
marg1.part.pred = vandermonde.matrix(pt.part1, poly.degree)%*%marg1.beta.part
marg2.part.pred = vandermonde.matrix(pt.part2, poly.degree)%*%marg2.beta.part
marg3.part.pred = vandermonde.matrix(pt.part3, poly.degree)%*%marg3.beta.part

theta.poly.part1 = cbind(pt.part1, marg1.part.pred)
theta.poly.part1 = theta.poly.part1[order(theta.poly.part1[,1]),]
theta.poly.part2 = cbind(pt.part2, marg2.part.pred)
theta.poly.part2 = theta.poly.part2[order(theta.poly.part2[,1]),]
theta.poly.part3 = cbind(pt.part3, marg3.part.pred)
theta.poly.part3 = theta.poly.part3[order(theta.poly.part3[,1]),]

## Partition lines for graphs
lines.part1 = seq(low[1], high[1], length.out = parts + 1)
lines.part2 = seq(low[2], high[2], length.out = parts + 1)
lines.part3 = seq(low[3], high[3], length.out = parts + 1)

### Plots for function evaluations, partitions, pointwise means and quadratic approximation
x11(30,10); par(mfrow=c(1,3))
plot(psi1, pch=16, main = "function evals", col="light blue")
points(psi1.parts, pch=16, col="red")
lines(theta.poly.part1, col="black", lwd=1.5)
abline(v = lines.part1, col = "green")

plot(psi2, pch=16, main = "function evals", col="light blue")
points(psi2.parts, pch=16, col="red")
lines(theta.poly.part2, col="black", lwd=1.5)
abline(v = lines.part2, col = "green")

plot(psi3, pch=16, main = "function evals", col="light blue")
points(psi3.parts, pch=16, col="red")
lines(theta.poly.part3, col="black", lwd=1.5)
abline(v = lines.part3, col = "green")

### Plots above, without function evaluations - Can see difference between polynomial and pointwise means
### For next step - treat these as "errors/residuals", to help fit a correction..
x11(30,10); par(mfrow=c(1,3))
plot(psi1.parts, pch=16, col="red")
lines(theta.poly.part1, col="black", lwd=1.5)
abline(v = lines.part1, col = "green")

plot(psi2.parts, pch=16, col="red")
lines(theta.poly.part2, col="black", lwd=1.5)
abline(v = lines.part2, col = "green")

plot(psi3.parts, pch=16, col="red")
lines(theta.poly.part3, col="black", lwd=1.5)
abline(v = lines.part3, col = "green")

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

theta.lds.part3 = cbind(pt.part3, exp(marg3.part.pred))
theta.ppn.const3 = diff(range(theta.lds.part3[,1]))*mean(theta.lds.part3[,2])
theta.lds.part3 = cbind(theta.lds.part3[,1], theta.lds.part3[,2]/theta.ppn.const3)

lds.part3 = cbind(exp(theta.lds.part3[,1]), theta.lds.part3[,2])
hp.ppn.const3 = diff(range(lds.part3[,1]))*mean(lds.part3[,2])
lds.part3 = cbind(lds.part3[,1], lds.part3[,2]/hp.ppn.const3)

#### Some plotting stuff ###
max.theta.ylim1 = max(c(theta.inla.hyperpar1[,2], theta.lds.part1[,2]))
max.theta.ylim2 = max(c(theta.inla.hyperpar2[,2], theta.lds.part2[,2]))
max.theta.ylim3 = max(c(theta.inla.hyperpar3[,2], theta.lds.part3[,2]))

max.hp.ylim1 = max(c(inla_hyperpar1[,2], lds.part1[,2]))
max.hp.ylim2 = max(c(inla_hyperpar2[,2], lds.part2[,2]))
max.hp.ylim3 = max(c(inla_hyperpar3[,2], lds.part3[,2]))

### Approximations - Our quadratic pproximation vs INLA's Dense Grid
x11(30,30); par(mfrow = c(2,2))
plot(theta.inla.hyperpar1, type="l", lwd = 2, ylim = c(0,max.theta.ylim1*1.05), cex.lab = 1.25,
     xlab = "theta_1 kappa")
lines(theta.lds.part1, type="l", col = "red", lwd = 2)
abline(v = log(100), col="green")

plot(theta.inla.hyperpar2, type="l", lwd = 2, ylim = c(0,max.theta.ylim2*1.05), cex.lab = 1.25,
     xlab = "theta_2 tau")
lines(theta.lds.part2, type="l", col = "red", lwd = 2)
abline(v = log(1), col="green")

plot(theta.inla.hyperpar3, type="l", lwd = 2, ylim = c(0,max.theta.ylim3*1.05), cex.lab = 1.25,
     xlab = "theta_3 rho")
lines(theta.lds.part3, type="l", col = "red", lwd = 2)
abline(v = log((1+0.65)/(1-0.65)), col="green")

plot(0,0,col="white")
text(0,0, cex = 1.5,
     "THETAS \n Black is inla.hyperpar \n Red is LDS-polynomial \n Green is the true value")

########################
### Cubic Correction ###
########################

### Polynomial fit at abscissa points
quad.poly1 = van.partM1%*%marg1.beta.part
quad.poly2 = van.partM2%*%marg2.beta.part
quad.poly3 = van.partM3%*%marg3.beta.part

### residuals = poly - pointwise means
resids1 = quad.poly1 - psi1.parts[,2]
resids2 = quad.poly2 - psi2.parts[,2]
resids3 = quad.poly3 - psi3.parts[,2]

ppdeg = 3

### get cubic coefficients using least squares through residuals
### Vandermonde matrix here is cubic..
vR1 = vandermonde.matrix(psi1.parts[,1], ppdeg)
vR2 = vandermonde.matrix(psi2.parts[,1], ppdeg)
vR3 = vandermonde.matrix(psi3.parts[,1], ppdeg)

another.beta1 = qr.solve(t(vR1)%*%vR1, tol = 1e-32)%*%t(vR1)%*%resids1
another.beta2 = qr.solve(t(vR2)%*%vR2, tol = 1e-32)%*%t(vR2)%*%resids2
another.beta3 = qr.solve(t(vR3)%*%vR3, tol = 1e-32)%*%t(vR3)%*%resids3

corrected.beta1 = rbind(marg1.beta.part,0) - another.beta1
corrected.beta2 = rbind(marg2.beta.part,0) - another.beta2
corrected.beta3 = rbind(marg3.beta.part,0) - another.beta3

## Cubic correction
marg1.corr.pred = vandermonde.matrix(pt.part1, ppdeg)%*%corrected.beta1
marg2.corr.pred = vandermonde.matrix(pt.part2, ppdeg)%*%corrected.beta2
marg3.corr.pred = vandermonde.matrix(pt.part3, ppdeg)%*%corrected.beta3

l.o = 1000

pt.part1 = seq(low[1], high[1], length.out = l.o)
pt.part2 = seq(low[2], high[2], length.out = l.o)
pt.part3 = seq(low[3], high[3], length.out = l.o)

theta.poly.corr1 = cbind(pt.part1, marg1.corr.pred)
theta.poly.corr1 = theta.poly.corr1[order(theta.poly.corr1[,1]),]
theta.poly.corr2 = cbind(pt.part2, marg2.corr.pred)
theta.poly.corr2 = theta.poly.corr2[order(theta.poly.corr2[,1]),]
theta.poly.corr3 = cbind(pt.part3, marg3.corr.pred)
theta.poly.corr3 = theta.poly.corr3[order(theta.poly.corr3[,1]),]

lines.part1 = seq(low[1], high[1], length.out = parts + 1)
lines.part2 = seq(low[2], high[2], length.out = parts + 1)
lines.part3 = seq(low[3], high[3], length.out = parts + 1)

### Plots of function evals, partitions, pointwise means and cubic correction
x11(30,10); par(mfrow=c(1,3))
plot(psi1, pch=16, main = "function evals", col="light blue")
points(psi1.parts, pch=16, col="red")
lines(theta.poly.corr1, col="black", lwd=1.5)
abline(v = lines.part1, col = "green")

plot(psi2, pch=16, main = "function evals", col="light blue")
points(psi2.parts, pch=16, col="red")
lines(theta.poly.corr2, col="black", lwd=1.5)
abline(v = lines.part2, col = "green")

plot(psi3, pch=16, main = "function evals", col="light blue")
points(psi3.parts, pch=16, col="red")
lines(theta.poly.corr3, col="black", lwd=1.5)
abline(v = lines.part3, col = "green")

### Cubic correction vs Initial (quadratic) approximation
x11(30,10); par(mfrow=c(1,3))
plot(psi1.parts, pch=16, col="red", cex = 1.5, xlab="", ylab="", main = "Theta(Tau_y))", cex.main=2)
lines(theta.poly.part1, col="black", lwd=2)
lines(theta.poly.corr1, col='blue', lty=2, lwd=2)
legend("topleft",
       c("Initial", "Cubic-Cor"), col = c("black","blue"), lwd=c(2,2), lty=c(1,2), cex=1.5)

plot(psi2.parts, pch=16, col="red", cex = 1.5, xlab="", ylab="", main = "Theta(Kappa))", cex.main=2)
lines(theta.poly.part2, col="black", lwd=2)
lines(theta.poly.corr2, col='blue', lty=2, lwd=2)
legend("topleft",
       c("Initial", "Cubic-Cor"), col = c("black","blue"), lwd=c(2,2), lty=c(1,2), cex=1.5)


plot(psi3.parts, pch=16, col="red", cex = 1.5, xlab="", ylab="", main = "Theta(Rho))", cex.main=2)
lines(theta.poly.part3, col="black", lwd=2)
lines(theta.poly.corr3, col='blue', lty=2, lwd=2)
legend("topleft",
       c("Initial", "Cubic-Cor"), col = c("black","blue"), lwd=c(2,2), lty=c(1,2), cex=1.5)

### Proper correction plots 
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

theta.lds.corr3 = cbind(pt.part3, exp(marg3.corr.pred))
theta.corr.const3 = diff(range(theta.lds.corr3[,1]))*mean(theta.lds.corr3[,2])
theta.lds.corr3 = cbind(theta.lds.corr3[,1], theta.lds.corr3[,2]/theta.corr.const3)

lds.corr3 = cbind((exp(theta.lds.corr3[,1])-1)/(exp(theta.lds.corr3[,1])+1), theta.lds.corr3[,2])
hp.corr.const3 = diff(range(lds.corr3[,1]))*mean(lds.corr3[,2])
lds.corr3 = cbind(lds.corr3[,1], lds.corr3[,2]/hp.corr.const3)

#### Some plotting stuff ###
max.theta.ylim1 = max(c(theta.inla.hyperpar1[,2], theta.lds.part1[,2]))
max.theta.ylim2 = max(c(theta.inla.hyperpar2[,2], theta.lds.part2[,2]))
max.theta.ylim3 = max(c(theta.inla.hyperpar3[,2], theta.lds.part3[,2]))

max.hp.ylim1 = max(c(inla_hyperpar1[,2], lds.part1[,2]))
max.hp.ylim2 = max(c(inla_hyperpar2[,2], lds.part2[,2]))
max.hp.ylim3 = max(c(inla_hyperpar3[,2], lds.part3[,2]))

### Plots - Thetas and hyperpars
## INLA dense grid vs Initial quadratic approx vs cubic correction
x11(35,20); par(mfrow = c(2,3))
plot(theta.inla.hyperpar1, type="l", lwd = 2, ylim = c(0,max.theta.ylim1*1.05), cex.lab = 1.25,
     xlab = "", ylab="", main = "Theta(Tau_y)", cex.main = 2)
lines(theta.lds.part1, type="l", col = "red", lwd = 2)
lines(theta.lds.corr1, type="l", col = "blue", lty =2, lwd = 2)
abline(v = log(100), col="green")
legend("topleft",
       c("INLA-grid", "Initial", "C-Corr"), lty = c(1,1,2), lwd = c(2,2,2), 
       col = c("black","red", "blue"), cex = 1.5)

plot(theta.inla.hyperpar2, type="l", lwd = 2, ylim = c(0,max.theta.ylim2*1.05), cex.lab = 1.25,
     xlab = "", ylab="", main = "Theta(Kappa)", cex.main = 2)
lines(theta.lds.part2, type="l", col = "red", lwd = 2)
lines(theta.lds.corr2, type="l", col = "blue", lty =2, lwd = 2)
abline(v = log(1), col="green")
legend("topleft",
       c("INLA-grid", "Initial", "C-Corr"), lty = c(1,1,2), lwd = c(2,2,2), 
       col = c("black","red", "blue"), cex = 1.5)

plot(theta.inla.hyperpar3, type="l", lwd = 2, ylim = c(0,max.theta.ylim3*1.05), cex.lab = 1.25,
     xlab = "", ylab="", main = "Theta(Rho)", cex.main = 2)
lines(theta.lds.part3, type="l", col = "red", lwd = 2)
lines(theta.lds.corr3, type="l", col = "blue", lty =2, lwd = 2)
abline(v = log((1+0.65)/(1-0.65)), col="green")
legend("topleft",
       c("INLA-grid", "Initial", "C-Corr"), lty = c(1,1,2), lwd = c(2,2,2), 
       col = c("black","red", "blue"), cex = 1.5)

plot(inla.smarginal(r2$marginals.hyperpar[[1]]), 
     type="l", lwd = 2, cex.lab = 1.25, xlim = c(min(lds.part1[,1]), max(lds.part1[,1])),
     xlab = "", ylab="", main = "Tau_y", cex.main = 2)
lines(lds.part1, type="l", col = "red", lwd = 2)
lines(lds.corr1, type="l", col = "blue", lty =2, lwd = 2)
abline(v = 100, col="green")
legend("topleft",
       c("INLA-grid", "Initial", "C-Corr"), lty = c(1,1,2), lwd = c(2,2,2), 
       col = c("black","red", "blue"), cex = 1.5)

plot(inla.smarginal(r2$marginals.hyperpar[[2]]), 
     type="l", lwd = 2, cex.lab = 1.25, xlim = c(min(lds.part2[,1]), max(lds.part2[,1])),
     xlab = "", ylab="", main = "Kappa", cex.main = 2)
lines(lds.part2, type="l", col = "red", lwd = 2)
lines(lds.corr2, type="l", col = "blue", lty =2, lwd = 2)
abline(v = 1, col="green")
legend("topright",
       c("INLA-grid", "Initial", "C-Corr"), lty = c(1,1,2), lwd = c(2,2,2), 
       col = c("black","red", "blue"), cex = 1.5)

plot(inla.smarginal(r2$marginals.hyperpar[[3]]), 
     type="l", lwd = 2, cex.lab = 1.25, xlim = c(min(lds.part3[,1]), max(lds.part3[,1])),
     xlab = "", ylab="", main = "Rho", cex.main = 2)
lines(lds.part3, type="l", col = "red", lwd = 2)
lines(lds.corr3, type="l", col = "blue", lty =2, lwd = 2)
abline(v = 0.65, col="green")
legend("topleft",
       c("INLA-grid", "Initial", "C-Corr"), lty = c(1,1,2), lwd = c(2,2,2), 
       col = c("black","red", "blue"), cex = 1.5)