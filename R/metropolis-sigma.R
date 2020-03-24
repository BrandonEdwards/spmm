#' Metropolis sigma
#'
#' Does the thing for sigma
#'

metropolis_sigma <- function(Z = NULL,
                             Y = NULL,
                             X = NULL,
                             D = NULL,
                             W = NULL,
                             at_risk = NULL,
                             disease = NULL,
                             region = NULL,
                             Sigma = NULL,
                             sigma_scale = NULL,
                             v = NULL,
                             k = NULL,
                             B = NULL,
                             Beta = NULL,
                             phi = NULL,
                             u_matrix = NULL,
                             pi = NULL)
{
  C <- dim(Sigma)[1]
  A <- chol(Sigma)

  for (i in 1:C)
  {
    for (j in i:C)
    {
      if (runif(1) > 0.5)
      {
        A[i,j] <- A[i,j] + (sigma_scale * rnorm(1))
      }
    }
  }

  Sigma_proposed <- A %*% t(A)#MCMCpack::riwish(v = v, S = Sigma)

  phi_u_proposed <- mcar(B = B,
                         Sigma = Sigma_proposed,
                         D = D,
                         W = W,
                         disease = disease,
                         region = region)

  u <- runif(1)

  loglik_1 <- loglik(Z = Z,
                     Y = Y,
                     X = X,
                     disease = disease,
                     region = region,
                     Beta = Beta,
                     at_risk = at_risk,
                     phi = phi_u_proposed[["phi"]],
                     pi = pi)

  prior_1 <- target_sigma(Sigma = Sigma_proposed,
                          v = v,
                          C = C,
                          nu = eigen(Sigma_proposed)$values,
                          k = k)

  loglik_2 <- loglik(Z = Z,
                     Y = Y,
                     X = X,
                     disease = disease,
                     region = region,
                     Beta = Beta,
                     at_risk = at_risk,
                     phi = phi,
                     pi = pi)

  prior_2 <- target_sigma(Sigma = Sigma,
                          v = v,
                          C = C,
                          nu = eigen(Sigma)$values,
                          k = k)

  log_ratio_sigma <- (loglik_1 + prior_1) - (loglik_2 + prior_2)

  if (u > 0)
  {
    if (log_ratio_sigma >= log(u))
    {
      return(list(Sigma = Sigma_proposed,
                  phi = phi_u_proposed[["phi"]],
                  u = phi_u_proposed[["u"]]))
    }else
    {
      return(list(Sigma = Sigma,
                  phi = phi,
                  u = u_matrix))
    }
  }else
  {
    return(list(Sigma = Sigma,
                phi = phi,
                u = u_matrix))
  }

  # for (i in 1:C)
  # {
  #   for (j in i:C)
  #   {
  #     Sigma_proposed[i,j] <- Sigma[i,j] + (sigma_scale + rnorm(1))
  #
  #     u <- runif(1)
  #
  #     loglik_1 <- loglik(Z = Z,
  #                        Y = Y,
  #                        X = X,
  #                        disease = disease,
  #                        region = region,
  #                        Beta = Beta,
  #                        at_risk = at_risk,
  #                        phi = phi,
  #                        pi = pi)
  #   }
  # }
}
