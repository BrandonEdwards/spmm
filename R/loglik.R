#' Log-likelihood
#'
#' \code{loglik} calculates the log-likelihood component of the SPMM
#'

loglik <- function(Z = NULL,
                   Y = NULL,
                   X = NULL,
                   disease = NULL,
                   region = NULL,
                   Beta = NULL,
                   at_risk = NULL,
                   phi = NULL,
                   pi = NULL,
                   alloc = NULL)
{
  ll_sum <- 0
  C <- length(pi)
  phi_vector <- NULL
  for (k in 1:length(disease))
  {
    phi_vector <- c(phi_vector, phi[disease[k], region[k]])
  }


  for (i in 1:C)
  {
    if (is.null(alloc))
    {
      ll_sum <- ll_sum + sum((Z[,i] * Y) * (X %*% Beta[,i] + phi_vector) -
                               (Z[,i] * at_risk) * exp(X %*% Beta[,i] + phi_vector) +
                               (Z[,i] * log(pi)))
    }else
    {
      ll_sum <- ll_sum + sum(alloc * ((Z[,i] * Y) * (X %*% Beta[,i] + phi_vector) -
                                      (Z[,i] * at_risk) * exp(X %*% Beta[,i] + phi_vector) +
                                      (Z[,i] * log(pi))))
    }
  }

  return(ll_sum)
}
