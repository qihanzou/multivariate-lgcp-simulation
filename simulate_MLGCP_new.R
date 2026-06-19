library(spatstat.geom)     
library(spatstat.explore)  
library(spatstat.random)   
library(spatstat.model)
library(ggplot2)


simulate_grf = function(nx, ny, L, range, var) {
  mu0 = as.im(0, owin(xrange = c(0, L), yrange = c(0, L)), dimyx = c(ny, nx))
  sim_field = rLGCP(model = "exponential", mu = mu0, var = var, scale = range, win = as.owin(mu0), saveLambda = TRUE)
  field_im  = attr(sim_field, "Lambda")
  field_im  = eval.im(log(field_im))   
  as.matrix(field_im)
}

simulate_LGCP_point_pattern <- function(params) {
  Y = simulate_grf(params$nx, params$ny, params$L, range = params$Y_scale, var = 1)
  Z = simulate_grf(params$nx, params$ny, params$L, range = params$Z_scale, var = 1)
  V = simulate_grf(params$nx, params$ny, params$L, range = params$V_scale, var = 1)
  
  
  U_list <- vector("list", params$p)
  for (j in seq_len(params$p)) {
    U_list[[j]] <- simulate_grf(params$nx, params$ny, params$L, range = params$U_scale[j], var = 1)
  }
  
  lambda0 <- exp(params$V_var*V - (params$V_var^2)/2)
  Lambda_list <- vector("list", params$p)
  X_list <- vector("list", params$p)
  for (i in seq_len(params$p)) {
    Lambda_list[[i]] <- lambda0*exp(params$beta0[i] + params$beta1[i]*Z)*exp(params$Y_var[i]*Y + params$U_var[i]*U_list[[i]] - (params$Y_var[i]^2)/2 - (params$U_var[i]^2)/2)
    Lambda_im = im(mat = Lambda_list[[i]], xrange = c(0, params$L), yrange = c(0, params$L))
    X_list[[i]] = rpoispp(Lambda_im)
  }
  
  
  list(X_list = X_list, 
       Lambda_list = Lambda_list, 
       lambda_base = lambda0,
       fields = list(Y = Y, Z = Z, V = V, U_list = U_list),
       params = params)
}



params = list(p = 4, L = 1, nx = 256, ny = 256,  
              beta0 = c(5.17, 5.44, 5.88, 6.13),
              beta1 = c(0.0, 0.3, -0.6, 0.6),
              
              Y_var = c(0.5, -0.4, 0.6, -0.3),
              V_var = 0.5,
              U_var = c(sqrt(0.5), sqrt(0.5), sqrt(0.5), sqrt(0.5)),
              
              Y_scale = 0.1,
              V_scale = 0.05,  
              Z_scale = 0.05,
              U_scale = c(0.05, 0.05, 0.05, 0.05)
)

sim <- simulate_LGCP_point_pattern(params)



# plots the plots
df <- do.call(
  rbind,
  lapply(seq_along(sim$X_list), function(i) {
    Xi <- sim$X_list[[i]]
    if (Xi$n == 0) return(NULL)
    data.frame(
      x = Xi$x,
      y = Xi$y,
      type = factor(i)
    )
  })
)

ggplot(df, aes(x, y, colour = type)) +
  geom_point(size = 1, alpha = 1) +
  coord_equal(xlim = c(0, params$L), ylim = c(0, params$L)) +
  theme_minimal() +
  labs(title = "Simulated multitype LGCP points", colour = "Type")


# Plot U fields
op1 <- par(mfrow = c(2, 2), mar = c(3, 3, 3, 3))
for (i in 1:4) {
  plot(spatstat.geom::as.im(sim$fields$U_list[[i]]),
       main = paste0("U_", i))
}
par(op1)
# Plot Y, Z, V fields
op2 <- par(mfrow = c(1, 3), mar = c(3, 3, 3, 3))
plot(spatstat.geom::as.im(sim$fields$Y), main = "Y")
plot(spatstat.geom::as.im(sim$fields$Z), main = "Z")
plot(spatstat.geom::as.im(sim$fields$V), main = "V")
par(op2)






X <- superimpose(
  Type1 = sim$X_list[[1]],
  Type2 = sim$X_list[[2]],
  Type3 = sim$X_list[[3]],
  Type4 = sim$X_list[[4]],
  W = Window(sim$X_list[[1]])
)

marks(X) <- factor(marks(X))  
levels(marks(X))





types <- levels(factor(marks(X)))
Kcross_list <- list()
for (i in types) {
  for (j in types) {
    if (i != j) {
      name <- paste0(i, "_", j)
      Kcross_list[[name]] <- Kcross.inhom(X, i = i, j = j, correction = "isotropic")
    }
  }
}
names(Kcross_list)



Kinhom_list <- list()
for (i in types) {
  name <- paste0(i)
  Xi <- X[marks(X) == i]
  Kinhom_list[[i]] <- Kinhom(Xi, correction = "isotropic")
}
names(Kinhom_list)


# print number of points
point_counts <- sapply(sim$X_list, function(x) x$n)
names(point_counts) <- paste0("Type", seq_along(point_counts))
print(point_counts)

