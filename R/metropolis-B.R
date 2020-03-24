#' Metropolis B
#'
#' Do thing for stuff
#'

metropolis_B <- function(Z = NULL,
                         Y = NULL,
                         X = NULL,
                         D = NULL,
                         W = NULL,
                         at_risk = NULL,
                         disease = NULL,
                         region = NULL,
                         Sigma = NULL,
                         v = NULL,
                         k = NULL,
                         B = NULL,
                         B_scale = NULL,
                         u_matrix = NULL,
                         Beta = NULL,
                         phi = NULL,
                         pi = NULL)
{
  C <- dim(B)[1]
  R <- dim(D)[1]

  G <- chol(B)

  for (i in 1:C)
  {
    for (j in i:C)
    {
      if (runif(1) > 0.5)
      {
        G[i,j] <- G[i,j] + (B_scale * rnorm(1))
      }
    }
  }

  B_proposed <- G %*% t(G)
  #B_proposed <- MCMCpack::riwish(v = v,S = B)

  phi_proposed_list <- mcar(B = B_proposed,
                            Sigma = Sigma,
                            D = D,
                            W = W,
                            disease = disease,
                            region = region)

  phi_proposed <- phi_proposed_list[["phi"]]
  u_matrix_proposed <- phi_proposed_list[["u"]]

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

  prior_1 <- target_B(D = D,
                      W = W,
                      v = v,
                      C = C,
                      R = R,
                      kappa = eigen(B_proposed)$values,
                      k = k,
                      B = B_proposed,
                      u = u_matrix_proposed)
  if (is.infinite(prior_1))
  {
    return(list(B = B,
                phi = phi,
                u = u_matrix))
  }

  loglik_2 <- loglik(Z = Z,
                     Y = Y,
                     X = X,
                     disease = disease,
                     region = region,
                     Beta = Beta,
                     at_risk = at_risk,
                     phi = phi,
                     pi = pi)

  prior_2 <- target_B(D = D,
                      W = W,
                      v = v,
                      C = C,
                      R = R,
                      kappa = eigen(B)$values,
                      k = k,
                      B = B,
                      u = u_matrix)

  log_ratio_B <- (loglik_1 + prior_1) - (loglik_2 + prior_2)

  if (u > 0)
  {
    if (log_ratio_B >= log(u))
    {
      return(list(B = B_proposed,
                  phi = phi_proposed,
                  u = u_matrix_proposed))
    }else
    {
      return(list(B = B,
                  phi = phi,
                  u = u_matrix))
    }
  }else
  {
    return(list(B = B,
                phi = phi,
                u = u_matrix))
  }
}
