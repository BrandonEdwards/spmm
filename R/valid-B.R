#' Check validity of B
#'
#' @export

valid_B <- function(zeta = NULL,
                    D = NULL,
                    W = NULL)
{
  xi <- eigen(pracma::sqrtm(D)$Binv %*% W %*% pracma::sqrtm(D)$Binv)$values

  return((1/min(xi) < zeta) && (zeta <= 1))
}
