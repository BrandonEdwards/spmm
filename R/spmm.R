#' SPMM
#'
#' \code{spmm} is the main function to run a spatial Poisson mixture model
#'
#' @param data The actual data
#' @param C The number of diseases to model
#' @param R The number of spatial regions
#' @param D R x R neighbourhood matrix with entries D_jj, j = 1...R, i.e. the count of neighbours for region j
#' @param W R x R adjacency matrix
#' @param m Number of covariates
#' @param n_iter How many iterations to sample?
#' @param n_burnin How many iterations to burn in chain?
#' @param n_chains How many chains to run per parameter?
#'
#' @export
#'

spmm <- function(data = NULL,
                 inits = NULL,
                 n_iter = NULL,
                 n_burnin = NULL,
                 n_chains = 3,
                 constant_B = FALSE)
{
  # Do error checking for null things

  ############################################################
  # Set initial values per parameter per chain
  ############################################################
  message("Initializing model...")

  counts <- data[["counts"]]
  disease <- data[["disease"]]
  region <- data[["region"]]
  at_risk <- data[["at_risk"]]
  alloc <- data[["alloc"]]
  C <- data[["C"]]
  R <- data[["R"]]
  D <- data[["D"]]
  W <- data[["W"]]
  N <- data[["N"]]
  X <- cbind(rep(1, N), data[["X"]])
  m <- dim(X)[2]

  # Set constants
  beta_scale <- 0.05
  sigma_scale <- 10
  B_scale <- 10
  k_1 <- 0.1
  v_1 <- C + 98
  k_2 <- 0.1
  v_2 <- C + 98
  N_u <- N - sum(alloc)
  total_iterations <- round(n_iter/n_chains) + n_burnin

  deviance_sum_observed <- 0
  deviance_sum_complete <- 0

  max_observed_posterior <- 0
  MAP_observed_index <- c(0,0)
  max_obs_mar_Z <- 0
  MAP_obs_mar_Z <- c(0,0)
  max_obs_mar_beta <- 0
  MAP_obs_mar_beta <- c(0,0)
  max_obs_mar_phi <- 0
  MAP_obs_mar_phi <- c(0,0)
  max_obs_mar_pi <- 0
  MAP_obs_mar_pi <- c(0,0)

  max_complete_posterior <- 0
  MAP_complete_posterior <- c(0,0)
  max_com_mar_Z <- 0
  MAP_com_mar_Z <- c(0,0)
  max_com_mar_beta <- 0
  MAP_com_mar_beta <- c(0,0)
  max_com_mar_phi <- 0
  MAP_com_mar_phi <- c(0,0)
  max_com_mar_pi <- 0
  MAP_com_mar_pi <- c(0,0)

  # Create empty matrices
  beta <- array(dim = c(m, C, n_chains))
  B <- array(dim = c(C, C, n_chains))
  Sigma <- array(dim = c(C, C, n_chains))
  A <- array(dim = c(C, C, n_chains))
  mcar <- array(dim = c(C*R, C*R, n_chains))
  u <- array(dim = c(C, R, n_chains))
  phi <- array(dim = c(C, R, n_chains))
  Z <- array(data = 0, dim = c(N, C, n_chains))
  pi <- array(dim = c(1, C, n_chains))

  # Initiate each parameter per chain
  for (i in 1:n_chains)
  {
    # Initiate betas
    beta[,,i] <- matrix(rep(0.5, m*C), nrow = m)

    # Initiate B matrix
    if (constant_B)
    {
      # Treat B as constant, i.e. just the identity  matrix
      B[,,i] <- diag(1, nrow = C)
    }else
    {
      val_B <- FALSE
      while (isFALSE(val_B))
      {
        B[,,i] <- MCMCpack::riwish(v = v_1,
                                   S = diag(k_1 * v_1, nrow = C))
        val_B <- valid_B(zeta = eigen(B[,,i])$values,
                         D = D,
                         W = W)
      }
    }

    # Initiate Sigma and set A = Cholesky Decomposition(Sigma)
    Sigma[,,i] <- MCMCpack::riwish(v = v_2,
                                  S = diag(k_2 * v_2, nrow = C))
    A[,,i] <- chol(Sigma[,,i])

    mcar_list <- mcar(B = B[,,i],
                      Sigma = Sigma[,,i],
                      D = D,
                      W = W,
                      disease = disease,
                      region = region)

    u[,,i] <- mcar_list[["u"]]
    phi[,,i] <- mcar_list[["phi"]]

    # Initiate Z
    repeat{
      n_u <- round(runif(C-1, min = 0, max = (N_u - 5)))
      if (sum(n_u) <= N_u) break
    }
    n_u[C] <- N_u - sum(n_u)

    for (k in 1:N)
    {
      if (alloc[k] == 1)
      {
        Z[k, disease[k], i] <- 1
      }else
      {
        r_allocation <- as.vector(t(rmultinom(1, 1, n_u / sum(n_u))))
        Z[k,,i] <- r_allocation
        n_u <- n_u - r_allocation
      }
    }

    # Initiate pi
    pi[,,i] <- colSums(Z[,,i]) / sum(colSums(Z[,,i]))
  }

  # Create empty chains for each tracked parameter
  beta_chain <- array(dim = c(m, C, round(n_iter/n_chains), n_chains))
  B_chain <- array(dim = c(C, C, round(n_iter/n_chains), n_chains))
  Sigma_chain <- array(dim = c(C, C, round(n_iter/n_chains), n_chains))
  A_chain <- array(dim = c(C, C, round(n_iter/n_chains), n_chains))
  mcar_chain <- array(dim = c(C*R, C*R, round(n_iter/n_chains), n_chains))
  u_chain <- array(dim = c(C, R, round(n_iter/n_chains), n_chains))
  phi_chain <- array(dim = c(C, R, round(n_iter/n_chains), n_chains))
  Z_chain <- array(dim = c(N, C, round(n_iter/n_chains), n_chains))
  pi_chain <- array(dim = c(1, C, round(n_iter/n_chains), n_chains))

  ####################################################
  # The main sampling algorithm
  ####################################################
  message(paste0("Done! Beginning burn-in period for ", n_burnin, " iterations."))

  burnin_pb <- progress::progress_bar$new(format = "\r[:bar] :percent eta: :eta",
                                          clear = FALSE,
                                          total = n_burnin,
                                          width = 100)
  sample_pb <- NULL
  burnin_pb$tick(0)

  for (i in 1:total_iterations)
  {
    for (j in 1:n_chains)
    {
      Z[,,j] <- metropolis_z(Y = counts,
                             Z = Z[,,j],
                             alloc = alloc,
                             X = X,
                             at_risk = at_risk,
                             disease = disease,
                             region = region,
                             Beta = beta[,,j],
                             phi = phi[,,j],
                             pi = pi[,,j],
                             N = N)

      pi[,,j] <- metropolis_pi(Y = counts,
                               Z = Z[,,j],
                               X = X,
                               at_risk = at_risk,
                               disease = disease,
                               region = region,
                               Beta = beta[,,j],
                               phi = phi[,,j],
                               pi = pi[,,j],
                               N = N)

      beta[,,j] <- metropolis_beta(Y = counts,
                                   Z = Z[,,j],
                                   X = X,
                                   at_risk = at_risk,
                                   disease = disease,
                                   region = region,
                                   Beta = beta[,,j],
                                   beta_scale = beta_scale,
                                   phi = phi[,,j],
                                   pi = pi[,,j])

      Sigma_and_phi <- metropolis_sigma(Z = Z[,,j],
                                        Y = counts,
                                        X = X,
                                        D = D,
                                        W = W,
                                        at_risk = at_risk,
                                        disease = disease,
                                        region = region,
                                        Sigma = Sigma[,,j],
                                        sigma_scale = sigma_scale,
                                        v = v_2,
                                        k = k_2,
                                        B = B[,,j],
                                        Beta = beta[,,j],
                                        phi = phi[,,j],
                                        u_matrix = u[,,j],
                                        pi = pi[,,j])
      Sigma[,,j] <- Sigma_and_phi[["Sigma"]]
      phi[,,j] <- Sigma_and_phi[["phi"]]
      u[,,j] <- Sigma_and_phi[["u"]]

      if (isFALSE(constant_B))
      {
        B[,,j] <- metropolis_B(Z = Z[,,j],
                               Y = counts,
                               X = X,
                               D = D,
                               W = W,
                               at_risk = at_risk,
                               disease = disease,
                               region = region,
                               Sigma = Sigma[,,j],
                               v = v_1,
                               k = k_1,
                               B = B[,,j],
                               B_scale = B_scale,
                               u_matrix = u[,,j],
                               Beta = beta[,,j],
                               phi = phi[,,j],
                               pi = pi[,,j])
      }else
      {
        B[,,j] <- diag(1, nrow = C)
      }

      # B[,,j] <- B_and_phi[["B"]]
      # phi[,,j] <- B_and_phi[["phi"]]
      # u[,,j] <- B_and_phi[["u"]]

      phi_and_u <- metropolis_u(Z = Z[,,j],
                                Y = counts,
                                X = X,
                                D = D,
                                W = W,
                                at_risk = at_risk,
                                disease = disease,
                                region = region,
                                Sigma = Sigma[,,j],
                                v = v_1,
                                k = k_1,
                                B = B[,,j],
                                B_scale = B_scale,
                                u_matrix = u[,,j],
                                Beta = beta[,,j],
                                phi = phi[,,j],
                                pi = pi[,,j],
                                constant_B = constant_B)
      phi[,,j] <- phi_and_u[["phi"]]
      u[,,j] <- phi_and_u[["u"]]

      # Add samples to chains if we are past the burnin period
      if (i > n_burnin)
      {
        Z_chain[,, i - n_burnin, j] <- Z[,,j]
        pi_chain[,, i - n_burnin, j] <- pi[,,j]
        beta_chain[,, i - n_burnin, j] <- beta[,,j]
        Sigma_chain[,, i - n_burnin, j] <- Sigma[,,j]
        B_chain[,, i - n_burnin, j] <- B[,,j]
        phi_chain[,, i - n_burnin, j] <- phi[,,j]

        ######## Observed likelihood/posterior stuff #############
        observed_likelihood <- loglik(Z = Z[,,j],
                                      Y = counts,
                                      X = X,
                                      disease = disease,
                                      region = region,
                                      Beta = beta[,,j],
                                      at_risk = at_risk,
                                      phi = phi[,,j],
                                      pi = pi[,,j],
                                      alloc = alloc)
        deviance_sum_observed <- deviance_sum_observed +
          (-2 * observed_likelihood)

        observed_posterior <- observed_likelihood +
          target_B(D = D,
                   W = W,
                   v = v_1,
                   C = C,
                   R = R,
                   kappa = eigen(B[,,j])$values,
                   k = k_1,
                   B = B[,,j],
                   u = u[,,j]) +
          target_sigma(Sigma = Sigma[,,j],
                       v = v_2,
                       C = C,
                       nu = eigen(Sigma[,,j])$values,
                       k = k_2)

        if (observed_posterior > max_observed_posterior)
        {
          max_observed_posterior <- observed_posterior
          MAP_observed_index <- c(i - n_burnin, j)
        }

        ######## Complete likelihood/posterior stuff #############
        complete_likelihood <- loglik(Z = Z[,,j],
                                      Y = counts,
                                      X = X,
                                      disease = disease,
                                      region = region,
                                      Beta = beta[,,j],
                                      at_risk = at_risk,
                                      phi = phi[,,j],
                                      pi = pi[,,j])
        deviance_sum_complete <- deviance_sum_complete +
          ((-2) * complete_likelihood)

        complete_posterior <- complete_likelihood +
          target_B(D = D,
                   W = W,
                   v = v_1,
                   C = C,
                   R = R,
                   kappa = eigen(B[,,j])$values,
                   k = k_1,
                   B = B[,,j],
                   u = u[,,j]) +
          target_sigma(Sigma = Sigma[,,j],
                       v = v_2,
                       C = C,
                       nu = eigen(Sigma[,,j])$values,
                       k = k_2)

        if (complete_posterior > max_complete_posterior)
        {
          max_complete_posterior <- complete_posterior
          MAP_complete_index <- c(i - n_burnin, j)
        }

      }
    }

    if (i <= n_burnin)
    {
      burnin_pb$tick()

      if (i == n_burnin)
      {
        message(paste0("\nDone! Beginning MCMC sampling for ", n_iter, " iterations."))
        sample_pb <- progress::progress_bar$new(format = "\r[:bar] :percent eta: :eta",
                                                clear = FALSE,
                                                total = round(n_iter / n_chains),
                                                width = 100)
      }
    }
    else
    {
      sample_pb$tick()
    }

  }

  ############### Calculating DICs #######################
  DIC <- c(0,0,0,0,0,0,0,0)

  mean_observed_deviance <- deviance_sum_observed / n_iter
  mean_complete_deviance <- deviance_sum_complete / n_iter
  #mean_conditional_deviance <- deviance_sum_conditional / n_iter

  Z_post <- R.utils::wrap(Z_chain, map = list(1,2,NA))
  pi_post <- R.utils::wrap(pi_chain, map = list(1,2,NA))
  beta_post <- R.utils::wrap(beta_chain, map = list(1,2,NA))
  B_post <- R.utils::wrap(B_chain, map = list(1,2,NA))
  Sigma_post <- R.utils::wrap(Sigma_chain, map = list(1,2,NA))
  phi_post <- R.utils::wrap(phi_chain, map = list(1,2,NA))

  Z_mean_raw <- apply(Z_post, c(1,2), mean)
  Z_mean <- matrix(0, nrow = N, ncol = C)
  for (i in 1:N)
  {
    index <- which(Z_mean_raw[i,] == max(Z_mean_raw[i,]))
    Z_mean[i,index] <- 1
  }
  Z_obs_mode <- Z_chain[,,MAP_observed_index[1], MAP_observed_index[2]]
  Z_com_mode <- Z_chain[,,MAP_complete_index[1], MAP_complete_index[2]]
  pi_mean <- apply(pi_post, 2, mean)
  pi_obs_mode <- pi_chain[,,MAP_observed_index[1], MAP_observed_index[2]]
  pi_com_mode <- pi_chain[,,MAP_complete_index[1], MAP_complete_index[2]]
  beta_mean <- apply(beta_post, c(1,2), mean)
  beta_obs_mode <- beta_chain[,,MAP_observed_index[1], MAP_observed_index[2]]
  beta_com_mode <- beta_chain[,,MAP_complete_index[1], MAP_complete_index[2]]
  B_mean <- apply(B_post, c(1,2), mean)
  B_obs_mode <- B_chain[,,MAP_observed_index[1], MAP_observed_index[2]]
  B_com_mode <- B_chain[,,MAP_complete_index[1], MAP_complete_index[2]]
  Sigma_mean <- apply(Sigma_post, c(1,2), mean)
  Sigma_obs_mode <- Sigma_chain[,,MAP_observed_index[1], MAP_observed_index[2]]
  Sigma_com_mode <- Sigma_chain[,,MAP_complete_index[1], MAP_complete_index[2]]
  phi_mean <- apply(phi_post, c(1,2), mean)
  phi_obs_mode <- phi_chain[,,MAP_observed_index[1], MAP_observed_index[2]]
  phi_com_mode <- phi_chain[,,MAP_complete_index[1], MAP_complete_index[2]]

  DIC[1] <- (-2 * mean_observed_deviance) -
    2 * loglik(Z = Z_mean,
               Y = counts,
               X = X,
               disease = disease,
               region = region,
               Beta = beta_mean,
               at_risk = at_risk,
               phi = phi_mean,
               pi = pi_mean,
               alloc = alloc)

  DIC[2] <- (-2 * mean_observed_deviance) -
    2 * loglik(Z = Z_obs_mode,
               Y = counts,
               X = X,
               disease = disease,
               region = region,
               Beta = beta_obs_mode,
               at_risk = at_risk,
               phi = phi_obs_mode,
               pi = pi_obs_mode,
               alloc = alloc)

  #DIC[3] <- (2 * mean_observed_deviance) +

  DIC[4] <- (-2 * mean_complete_deviance) -
    2 * loglik(Z = Z_mean,
               Y = counts,
               X = X,
               disease = disease,
               region = region,
               Beta = beta_mean,
               at_risk = at_risk,
               phi = phi_mean,
               pi = pi_mean)

  DIC[5] <- (-2 * mean_complete_deviance) -
    2 * loglik(Z = Z_com_mode,
               Y = counts,
               X = X,
               disease = disease,
               region = region,
               Beta = beta_com_mode,
               at_risk = at_risk,
               phi = phi_com_mode,
               pi = pi_com_mode)

  to_return <- list(Z_samples = Z_chain,
                    pi_samples = pi_chain,
                    beta_samples = beta_chain,
                    B_samples = B_chain,
                    Sigma_samples = Sigma_chain,
                    A_samples = A_chain,
                    mcar_samples = mcar_chain,
                    u_samples = u_chain,
                    phi_samples = phi_chain,
                    mean_observed_deviance = deviance_sum_observed / n_iter,
                    mean_complete_deviance = deviance_sum_complete / n_iter,
                    MAP_observed_index = MAP_observed_index,
                    MAP_complete_index = MAP_complete_index,
                    DIC = DIC)

  return(to_return)
}
