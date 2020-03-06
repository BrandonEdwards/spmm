#' Metropolis Z
#'
#' Does M-H algo for Z, the allocation vector
#'

metropolis_z <- function(Y = NULL,
                         Z = NULL,
                         z = NULL,
                         X = NULL,
                         at_risk = NULL,
                         disease = NULL,
                         Beta = NULL,
                         phi = NULL,
                         pi = NULL,
                         N = NULL)
{
  C <- dim(Z)[2]
  Z_proposed <- matrix(0, nrow = N, ncol = C)

  N_u <- N - sum(z)
  repeat{
    n_u <- round(runif(C-1, min = 0, max = N_u))
    if (sum(n_u) <= N_u) break
  }
  n_u[C] <- N_u - sum(n_u)

  for (k in 1:N)
  {
    if (z[k] == 1)
    {
      Z_proposed[k, disease[k]] <- 1
    }else
    {
      r_allocation <- as.vector(t(rmultinom(1, 1, n_u / sum(n_u))))
      Z_proposed[k,] <- r_allocation
      n_u <- n_u - r_allocation
    }
  }

  pi_proposed <- colSums(Z_proposed) / sum(colSums(Z_proposed))
  # loglik_1 <- loglik_2 <- 0
  # for (i in 1:C)
  # {
  #   loglik_1 <- loglik_1 +
  #     sum((Z_proposed[,i] * y) * (X %*% beta[,i]) -
  #           (Z_proposed[,i] * at_risk) * exp(X %*% beta[,i]) +
  #           (Z_proposed[,i] * log(pi_proposed)))
  #
  #   loglik_2 <- loglik_2 +
  #     sum((Z[,i] * y) * (X %*% beta[,i]) -
  #           (Z[,i] * at_risk) * exp(X %*% beta[,i]) +
  #           (Z[,i] * log(pi)))
  # }

  u <- runif(1)
  loglik_1 <- loglik(Z = Z_proposed,
                     Y = Y,
                     X = X,
                     Beta = Beta,
                     at_risk = at_risk,
                     phi = phi,
                     pi = pi_proposed)

  loglik_2 <- loglik(Z = Z,
                     Y = Y,
                     X = X,
                     Beta = Beta,
                     at_risk = at_risk,
                     phi = phi,
                     pi = pi)
  log_ratio <- loglik_1 - loglik_2

  if (u > 0)
  {
    if (log_ratio >= log(u))
    {
      return(list(Z = Z_proposed,
                  pi = pi_proposed))
    }else
    {
      return(list(Z = Z,
                  pi = pi))
    }
  }else
  {
    return(list(Z = Z,
                pi = pi))
  }
}
