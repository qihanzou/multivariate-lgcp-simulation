
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



