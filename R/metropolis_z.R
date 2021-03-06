#' Metropolis Z
#'
#' Does M-H algo for Z, the allocation vector
#'
#' @importFrom stats rmultinom
#'

metropolis_z <- function(Y = NULL,
                         Z = NULL,
                         alloc = NULL,
                         X = NULL,
                         at_risk = NULL,
                         disease = NULL,
                         region = NULL,
                         Beta = NULL,
                         phi = NULL,
                         pi = NULL,
                         N = NULL)
{
  C <- dim(Z)[2]
  Z_proposed <- matrix(0, nrow = N, ncol = C)

  # N_u <- N - sum(alloc)
  # repeat{
  #   n_u <- round(runif(C-1, min = 0, max = (N_u - 5)))
  #   if (sum(n_u) <= N_u) break
  # }
  # n_u[C] <- N_u - sum(n_u)

  for (k in 1:N)
  {
    if (alloc[k] == 1)
    {
      Z_proposed[k, disease[k]] <- 1
    }else
    {
      r_allocation <- as.vector(t(rmultinom(1, 1, pi)))
      Z_proposed[k,] <- r_allocation
      #n_u <- n_u - r_allocation
    }
  }

  #pi_proposed <- colSums(Z_proposed) / sum(colSums(Z_proposed))
 #print(pi_proposed)

  u <- runif(1)
  loglik_1 <- loglik(Z = Z_proposed,
                     Y = Y,
                     X = X,
                     disease = disease,
                     region = region,
                     Beta = Beta,
                     at_risk = at_risk,
                     phi = phi,
                     pi = pi)

  loglik_2 <- loglik(Z = Z,
                     Y = Y,
                     X = X,
                     disease = disease,
                     region = region,
                     Beta = Beta,
                     at_risk = at_risk,
                     phi = phi,
                     pi = pi)
  log_ratio_z <- loglik_1 - loglik_2
  #print(paste0("loglik1: ", loglik_1, " loglik2: ", loglik_2))

  if (u > 0)
  {
    if (log_ratio_z >= log(u))
    {
      return(Z_proposed)
    }else
    {
      return(Z)
    }
  }else
  {
    return(Z)
  }
}
