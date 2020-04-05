#' MCAR
#'
#' Multivariate Conditionally Autoregressive
#'
#' @export
#'

mcar <- function(B = NULL,
                 Sigma = NULL,
                 D = NULL,
                 W = NULL,
                 disease = NULL,
                 region = NULL)
{
  C <- dim(Sigma)[1]
  R <- dim(D)[1]
  I_R <- diag(x = 1, nrow = R, ncol = R)
  I_C <- diag(x = 1, nrow = C, ncol = C)
  A <- chol(Sigma)

  mcar <- (diag(x = 1, nrow = C) %x% D) - (B %x% W)
  u <- as.vector(mvtnorm::rmvnorm(n = 1,
                        mean = rep(0, R*C),
                        sigma = mcar))

  u_matrix <- matrix(u, nrow = C, ncol = R, byrow = T)
  phi <- matrix(0, nrow = C, ncol = R)

  for (i in 1:C)
  {
    for (j in 1:R)
    {
      for (l in 1:C)
      {
        phi[i,j] <- phi[i,j] + A[i,l]*u_matrix[l,j]
      }
    }
  }

  # Initiate phi
  #phi <- (A %x% diag(x = 1, nrow = R) %*% u)

 # mvnorm_sigma <-

  return(list(u = u_matrix,
              phi = phi))

}
