#' Metropolis sigma
#'
#' Does the thing for sigma
#'
#' @importFrom MCMCpack rwish
#'
#' @export
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
  R <- dim(D)[1]

  Sigma_proposed <- MCMCpack::rwish(v = v,
                                    S = Sigma/v)

  A <- chol(Sigma_proposed)
  phi_proposed <- matrix(0, nrow = C, ncol = R)

  for (i in 1:C)
  {
    for (j in 1:R)
    {
      for (l in 1:C)
      {
        phi_proposed[i,j] <- phi_proposed[i,j] + A[i,l]*u_matrix[l,j]
      }
    }
  }

  u <- runif(1)

  loglik_1 <- loglik(Z = Z,
                     Y = Y,
                     X = X,
                     disease = disease,
                     region = region,
                     Beta = Beta,
                     at_risk = at_risk,
                     phi = phi_proposed,
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
                  phi = phi_proposed,
                  u = u_matrix))
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

}
