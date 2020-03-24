#' Metropolis Beta
#'
#' \code{metropolis_beta} runs the Metropolis-Hastings algorithm to
#'   update each Beta term one at a time
#'

metropolis_beta <- function(Z = NULL,
                            Y = NULL,
                            X = NULL,
                            at_risk = NULL,
                            disease = NULL,
                            region = NULL,
                            Beta = NULL,
                            beta_scale = NULL,
                            phi = NULL,
                            pi = NULL)
{
  beta_dim <- dim(Beta)
  beta_list <- as.vector(Beta)

  for (m in 1:length(beta_list))
  {
    beta_list_proposed <- beta_list

    beta_list_proposed[m] <- beta_list[m] + (beta_scale * rnorm(1))

    u <- runif(1)
    loglik_1 <- loglik(Z = Z,
                       Y = Y,
                       X = X,
                       disease = disease,
                       region = region,
                       Beta = matrix(beta_list_proposed,
                                     nrow = beta_dim[1],
                                     ncol = beta_dim[2]),
                       at_risk = at_risk,
                       phi = phi,
                       pi = pi)

    loglik_2 <- loglik(Z = Z,
                       Y = Y,
                       X = X,
                       disease = disease,
                       region = region,
                       Beta = matrix(beta_list,
                                     nrow = beta_dim[1],
                                     ncol = beta_dim[2]),
                       at_risk = at_risk,
                       phi = phi,
                       pi = pi)
    log_ratio_beta <- loglik_1 - loglik_2

    if (u > 0)
    {
      if (log_ratio_beta >= log(u))
      {
        beta_list <- beta_list_proposed
      }
    }
  }

  return(matrix(beta_list, nrow = beta_dim[1], ncol = beta_dim[2]))
}
