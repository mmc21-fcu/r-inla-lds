
# -------------------------------
# Simulate example climate dataset
# -------------------------------
set.seed(123)
n <- 200
climate_data <- data.frame(
  lon = runif(n, 0, 10),             # Longitude coordinate
  lat = runif(n, 0, 10),             # Latitude coordinate
  time = sample(2000:2010, n, replace = TRUE),  # Years as a proxy for time
  covariate = rnorm(n)               # Some environmental covariate (e.g., solar radiation)
)
# Simulate a response variable (e.g., temperature anomaly)
# with an intercept and effect from the covariate.
climate_data$y <- 0.5 + 0.8 * climate_data$covariate + rnorm(n)

# -------------------------------
# 1. Create Spatial Mesh for the SPDE
# -------------------------------
mesh <- inla.mesh.2d(
  loc = cbind(climate_data$lon, climate_data$lat),
  max.edge = c(1, 3),   # Adjust depending on the density of your locations
  cutoff = 0.1          # Minimum distance between points
)

# -------------------------------
# 2. Define the SPDE Model for the Spatial Field
# -------------------------------
# Here we use the popular PC-Matérn formulation.
# The pc-priors are set so that P(range < 2) = 0.01 and P(sigma > 1) = 0.01.
spde <- inla.spde2.pcmatern(
  mesh = mesh,
  prior.range = c(2, 0.01),   # Prior for spatial range
  prior.sigma = c(1, 0.01)    # Prior for spatial standard deviation
)

# Create a projector matrix 'A' for mapping observed locations onto the mesh nodes.
A <- inla.spde.make.A(mesh = mesh, loc = cbind(climate_data$lon, climate_data$lat))

# Generate an index for the spatial random field.
s.index <- inla.spde.make.index("spatial.field", n.spde = spde$n.spde)

# -------------------------------
# 3. Define the Temporal Effect
# -------------------------------
# Create a time index (ensuring it is numeric and sequential)
climate_data$time.index <- as.numeric(as.factor(climate_data$time))
n_time <- length(unique(climate_data$time.index))

# -------------------------------
# 4. Prepare the Data Stack for INLA
# -------------------------------
stack <- inla.stack(
  data = list(y = climate_data$y),
  A = list(A, 1),  # The first matrix 'A' links the spatial field; the second is the identity for fixed effects.
  effects = list(
    spatial.field = s.index$spatial.field,
    data = data.frame(
      Intercept = 1, 
      covariate = climate_data$covariate,
      time.index = climate_data$time.index
    )
  ),
  tag = "est"
)

# -------------------------------
# 5. Build the Model Formula
# -------------------------------
# Notice the inclusion of:
#  - Fixed effects: intercept and covariate.
#  - f(spatial.field): the spatial effect via the SPDE model.
#  - f(time.index): an AR(1) temporal effect with custom PC priors.
formula <- y ~ 0 + Intercept + covariate +
  f(spatial.field, model = spde) + 
  f(time.index, model = "ar1",
    hyper = list(
      # AR1: precision hyperparameter; here pc.prec prior is used.
      prec = list(prior = "pc.prec", param = c(1, 0.01)),
      # AR1: autocorrelation hyperparameter, here the PC prior for AR1 is set.
      rho  = list(prior = "pc.ar", param = 0.9)
    )
  )

# -------------------------------
# 6. Fit the Model Using INLA
# -------------------------------
result <- inla(
  formula,
  data = inla.stack.data(stack),
  family = "gaussian",  # Gaussian likelihood (noise precision associated to this likelihood will be estimated)
  control.predictor = list(A = inla.stack.A(stack), compute = TRUE),
  control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE)
)

# -------------------------------
# 7. Explore the Results
# -------------------------------
summary(result)