#' Log-likelihood
#'
#' \code{loglik} calculates the log-likelihood component of the SPMM
#'

loglik <- function(Z = NULL,
                   Y = NULL,
                   X = NULL,
                   Beta = NULL,
                   at_risk = NULL,
                   phi = NULL,
                   pi = NULL)
{
  ll_sum <- 0
  C <- length(pi)

  for (i in 1:C)
  {
    ll_sum <- ll_sum + sum((Z[,i] * Y) * (X %*% Beta[,i]) -
                           (Z[,i] * at_risk) * exp(X %*% Beta[,i]) +
                           (Z[,i] * log(pi)))
  }

  return(ll_sum)
}
