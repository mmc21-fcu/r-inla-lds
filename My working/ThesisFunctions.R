## Functions storage
library(randtoolbox)

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

