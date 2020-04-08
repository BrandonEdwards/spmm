#' Metropolis pi
#'
#' @export
#'

metropolis_pi <- function(Y = NULL,
                          Z = NULL,
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
  pi_proposed <- abs(pi + rnorm(n = C, mean = 0, sd = 0.5))
  pi_proposed <- pi_proposed / sum(pi_proposed)

  u <- runif(1)
  loglik_1 <- loglik(Z = Z,
                     Y = Y,
                     X = X,
                     disease = disease,
                     region = region,
                     Beta = Beta,
                     at_risk = at_risk,
                     phi = phi,
                     pi = pi_proposed)

  loglik_2 <- loglik(Z = Z,
                     Y = Y,
                     X = X,
                     disease = disease,
                     region = region,
                     Beta = Beta,
                     at_risk = at_risk,
                     phi = phi,
                     pi = pi)

  log_ratio_pi <- loglik_1 - loglik_2

  if (u > 0)
  {
    if (log_ratio_pi >= log(u))
    {
      return(pi_proposed)
    }else
    {
      return(pi)
    }
  }else
  {
    return(pi)
  }
}
