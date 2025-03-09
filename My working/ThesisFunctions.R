## Functions storage
library(randtoolbox)
library(INLA)
library(parallel)
### Korobov Lattice
shift.koro = function(n.pts, dimension, gen_constant, shift){
  n = n.pts; s = dimension; a = gen_constant; v = runif(s,0,1)
  z.korob = vector(length = s)
  z.korob[1] = 1
  for (j in 1:(s-1)){
    z.korob[j+1] = (z.korob[j] * a) %% n
  }
  korob = array(dim = c(n,s))
  if(shift == TRUE){
    for (k in 0:(n-1)){
      korob[k+1,] = (((k/n)*z.korob) + v) %% 1
    }
  } else {
    for (k in 0:(n-1)){
      korob[k+1,] = ((k/n)*z.korob) %% 1
    }
  }
  return(korob)
}

#Sobol
gen_sobol <- function(n,s){
  return(sobol(n,s,seed=1242))
}

#Halton
gen_halton <- function(n,s){
  return(halton(n,s))
}


### Transform pointset to square integration region
trans.pset = function(points, low, high){
  vec = vector(length=length(low))
  for (i in 1:length(vec)){
    vec[i] = high[i]-low[i]
  }
  t.pset = array(dim = dim(points))
  for(j in 1:nrow(points)){
    for(i in 1:ncol(points)){
      t.pset[j,i] = vec[i]*points[j,i] + low[i]
    }}
  return(t.pset)
}

### Kullback-Liebler divergence - for comparisons of distributions
kl.div = function(true.dist,approx.dist)
{	if(length(true.dist) != length(approx.dist))
{message("Try again mate!")}
  else{
    n = length(true.dist)
    try = cbind(true.dist,approx.dist)
    h = apply(try,2,sum)
    try.again = cbind(try[,1]/h[1], try[,2]/h[2])	#Normalising
    kld = vector(length=n)
    for(i in 1:n)
    {
      kld[i] = log(try.again[i,1]/try.again[i,2])*try.again[i,1]
    }
    return(sum(kld))
  }} 

# For function evaluations using INLA
inla.evaluate = function(theta, r){
  ## re-evaluate the model at fixed configuration 'theta'. return the new result.
  r$.args$control.mode = list(theta = theta, x = r$mode$x, fixed = TRUE)
  return (do.call("inla", args = r$.args))
}

### Scaled-Beta [Beta with support (-1,1)] - for prior of rho
scaledBeta = function(x, a, b, log=FALSE)
{
  ## the input 'x' is in the interval (-1, 1)
  xx = (x + 1)/2
  f = dbeta(xx, shape1=a, shape=b, log=TRUE) - log(2)
  return (if (log) f else exp(f))
}

### Partitioning functions
meanBetween = function(vec){
  meanVec = vector(length = length(vec) - 1)
  for(i in 1:length(meanVec)){
    meanVec[i] = 0.5*(vec[i] + vec[i+1])
  }
  return(meanVec)
}

partition = function(mat, support, nPart, position){
  parts = seq(support[1], support[2], length.out = nPart+1)
  m.vec = vector(length = nPart)
  m.pos = vector(length = nPart) 
  for(i in 1:nPart-1){		###FINDING AVERAGES
    m.vec[i] = mean(subset(mat[,2], mat[,1] >= parts[i] & mat[,1] < parts[i+1]))
  }
  m.vec[nPart] = mean(subset(mat[,2], mat[,1] >= parts[nPart] & mat[,1] <= parts[nPart+1]))
  ###POSITIONING
  if(position == 0){
    m.pos = parts[1:nPart]
  } else if(position == 1){
    m.pos = parts[2:(nPart+1)]
  } else{
    m.pos = meanBetween(parts)
  }
  return(cbind(m.pos, m.vec))
}

### Vandermonde Matrix
vandermonde.matrix = function(x.pts, degree){
  x = x.pts; y = 0:degree
  v = outer(x,y,"^")
  return(v)
}
