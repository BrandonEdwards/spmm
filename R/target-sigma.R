#' Target Sigma
#'
#' Prior calculation for sigma
#'
#' @export
#'

target_sigma <- function(Sigma = NULL,
                         v = NULL,
                         C = NULL,
                         nu = NULL,
                         k = NULL)
{
  i <- seq(1, C)
  A <- chol(Sigma)

  return(
    -((v + C + 1) / 2) *
      sum(log(nu)) -
      ((k*v) / 2) * sum(1 / nu) +
      sum(i * log(diag(A)))
  )
}
