#' Metropolis B
#'
#' Do thing for stuff
#'
#' @importFrom MCMCpack rwish
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

  val_B <- FALSE
  while (isFALSE(val_B))
  {
    B_proposed <- MCMCpack::rwish(v = v, S = B/v)
    val_B <- valid_B(zeta = eigen(B_proposed)$values,
                     D = D,
                     W = W)
  }

  prior_1 <- target_B(D = D,
                      W = W,
                      v = v,
                      C = C,
                      R = R,
                      kappa = eigen(B_proposed)$values,
                      k = k,
                      B = B_proposed,
                      u = u_matrix)

  prior_2 <- target_B(D = D,
                      W = W,
                      v = v,
                      C = C,
                      R = R,
                      kappa = eigen(B)$values,
                      k = k,
                      B = B,
                      u = u_matrix)

  #log_ratio_B <- (loglik_1 + prior_1) - (loglik_2 + prior_2)
  log_ratio_B <- prior_1 - prior_2
  u <- runif(1)
  if (u > 0)
  {
    if (log_ratio_B >= log(u))
    {
      return(B_proposed)
    }else
    {
      return(B)
    }
  }else
  {
    return(B)
  }


#
#
#
#
#   phi_proposed_list <- mcar(B = B_proposed,
#                             Sigma = Sigma,
#                             D = D,
#                             W = W,
#                             disease = disease,
#                             region = region)
#
#   phi_proposed <- phi_proposed_list[["phi"]]
#   u_matrix_proposed <- phi_proposed_list[["u"]]
#
#   u <- runif(1)
#
#   loglik_1 <- loglik(Z = Z,
#                      Y = Y,
#                      X = X,
#                      disease = disease,
#                      region = region,
#                      Beta = Beta,
#                      at_risk = at_risk,
#                      phi = phi_proposed,
#                      pi = pi)
#
#   prior_1 <- target_B(D = D,
#                       W = W,
#                       v = v,
#                       C = C,
#                       R = R,
#                       kappa = eigen(B_proposed)$values,
#                       k = k,
#                       B = B_proposed,
#                       u = u_matrix_proposed)
#   if (is.infinite(prior_1))
#   {
#     return(list(B = B,
#                 phi = phi,
#                 u = u_matrix))
#   }
#
#   loglik_2 <- loglik(Z = Z,
#                      Y = Y,
#                      X = X,
#                      disease = disease,
#                      region = region,
#                      Beta = Beta,
#                      at_risk = at_risk,
#                      phi = phi,
#                      pi = pi)
#
#   prior_2 <- target_B(D = D,
#                       W = W,
#                       v = v,
#                       C = C,
#                       R = R,
#                       kappa = eigen(B)$values,
#                       k = k,
#                       B = B,
#                       u = u_matrix)
#
#   log_ratio_B <- (loglik_1 + prior_1) - (loglik_2 + prior_2)
#   #log_ratio_B <- prior_1 - prior_2
#
#   if (u > 0)
#   {
#     if (log_ratio_B >= log(u))
#     {
#       return(list(B = B_proposed,
#                   phi = phi_proposed,
#                   u = u_matrix_proposed))
#     }else
#     {
#       return(list(B = B,
#                   phi = phi,
#                   u = u_matrix))
#     }
#   }else
#   {
#     return(list(B = B,
#                 phi = phi,
#                 u = u_matrix))
#   }
}
