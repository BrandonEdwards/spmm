#' Target B
#'

target_B <- function(D = NULL,
                     W = NULL,
                     v = NULL,
                     C = NULL,
                     R = NULL,
                     kappa = NULL,
                     k = NULL,
                     B = NULL,
                     u = NULL)
{
  zeta <- kappa
  xi <- eigen(pracma::sqrtm(D)$Binv %*% W %*% pracma::sqrtm(D)$Binv)$values
  #
  # eigen_check <- (1/min(xi) < zeta) && (zeta < 1)
  # if (!all(eigen_check))
  # {
  #   return(-Inf)
  # }

  if (!valid_B(zeta = kappa,
               D = D,
               W = W))
  {
    return(-Inf)
  }

  prob_u_given_B_1 <- 0

  for (i in 1:C)
  {
    for (j in 1:R)
    {
      prob_u_given_B_1 <- prob_u_given_B_1 + log(1 - (zeta[i] * xi[j]))
    }
  }

  prob_u_given_B_2 <- 0

  for (i in 1:C)
  {
    for (j in 1:R)
    {
      for (l in 1:R)
      {
        prob_u_given_B_2 <- prob_u_given_B_2 +
          B[i,i] * u[i,j] * W[j,l] * u[i,l]
      }
    }
  }

  prob_u_given_B_3 <- 0
  for (i in 1:(C - 1))
  {
    for (j in (i+1):C)
    {
      for (k in 1:R)
      {
        for (l in 1:R)
        {
          prob_u_given_B_3 <- prob_u_given_B_3 +
            B[i,i] * u[i,j] * W[j,l] * u[i,l]
        }
      }
    }
  }

  prob_u_given_B <- 0.5*prob_u_given_B_1 + 0.5*prob_u_given_B_2 + prob_u_given_B_3

  prob_B <- -((v + C + 1) / 2) *
    sum(log(kappa)) -
    ((k*v) / 2) * sum(1 / kappa)

  return(prob_u_given_B + prob_B)
}
