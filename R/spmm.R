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
                 n_chains = 3)
{
  # Do error checking for null things

  message("Initializing model...")
  ###### Extract the data for ease of access #########
  counts <- data[["counts"]]
  disease <- data[["disease"]]
  region <- data[["region"]]
  at_risk <- data[["at_risk"]]
  z <- data[["z"]]
  C <- data[["C"]]
  R <- data[["R"]]
  D <- data[["D"]]
  W <- data[["W"]]
  N <- data[["N"]]
  X <- cbind(rep(1, N), data[["X"]])
  m <- dim(X)[2]

  ############# Extract or set initial values per chain #############

  # Set constants
  beta_scale <- 0.05
  k_1 <- 0.1
  v_1 <- C + 2
  k_2 <- 0.1
  v_2 <- C + 2
  N_u <- N - sum(z)
  total_iterations <- round(n_iter/n_chains) + n_burnin

  # Create empty matrices
  beta <- array(dim = c(C, m, n_chains))
  B <- array(dim = c(C, C, n_chains))
  Sigma <- array(dim = c(C, C, n_chains))
  A <- array(dim = c(C, C, n_chains))
  mcar <- array(dim = c(C*R, C*R, n_chains))
  u <- array(dim = c(C*R, 1, n_chains))
  phi <- array(dim = c(C*R, 1, n_chains))
  Z <- array(data = 0, dim = c(N, C, n_chains))
  pi <- array(dim = c(1, C, n_chains))

  for (i in 1:n_chains)
  {
    # Initiate betas
    beta[,,i] <- matrix(rep(0.5, m*C), nrow = C)

    # Initiate B matrix
    B[,,i] <- cov2cor(MCMCpack::riwish(v = v_1,
                                       S = diag(k_1 * v_1, nrow = C)))

    # Initiate Sigma and set A = Cholesky Decomposition(Sigma)
    Sigma[,,i] <- MCMCpack::riwish(v = v_2,
                                  S = diag(k_2 * v_2, nrow = C))
    A[,,i] <- chol(Sigma[,,i])

    # Initiate mcar and u
    mcar[,,i] <- (diag(x = 1, nrow = C) %x% D) - (B[,,i] %x% W)
    u[,,i] <- mvtnorm::rmvnorm(n = 1,
                             mean = rep(0, R*C),
                             sigma = mcar[,,i])

    # Initiate phi
    phi[,,i] <- (A[,,i] %x% diag(x = 1, nrow = R) %*% u[,,i])

    # Initiate Z
    repeat{
      n_u <- round(runif(C-1, min = 0, max = N_u))
      if (sum(n_u) <= N_u) break
    }
    n_u[C] <- N_u - sum(n_u)

    for (k in 1:N)
    {
      if (z[k] == 1)
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

  beta_chain <- array(dim = c(C, m, round(n_iter/n_chains), n_chains))
  B_chain <- array(dim = c(C, C, round(n_iter/n_chains), n_chains))
  Sigma_chain <- array(dim = c(C, C, round(n_iter/n_chains), n_chains))
  A_chain <- array(dim = c(C, C, round(n_iter/n_chains), n_chains))
  mcar_chain <- array(dim = c(C*R, C*R, round(n_iter/n_chains), n_chains))
  u_chain <- array(dim = c(C*R, 1, round(n_iter/n_chains), n_chains))
  phi_chain <- array(dim = c(C*R, 1, round(n_iter/n_chains), n_chains))
  Z_chain <- array(dim = c(N, C, round(n_iter/n_chains), n_chains))
  pi_chain <- array(dim = c(1, C, round(n_iter/n_chains), n_chains))

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
      Z_and_pi <- metropolis_z(Y = counts,
                               Z = Z[,,j],
                               z = z,
                               X = X,
                               at_risk = at_risk,
                               disease = disease,
                               region = region,
                               Beta = beta[,,j],
                               phi = phi[,,j],
                               pi = pi[,,j],
                               N = N)
      Z[,,j] <- Z_and_pi[["Z"]]
      pi[,,j] <- Z_and_pi[["pi"]]

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

      # Add samples to chains if we are past the burnin period
      if (i < n_burnin)
      {
        Z_chain[,, i - n_burnin, j] <- Z[,,j]
        pi_chain[,, i - n_burnin, j] <- pi[,,j]
        beta_chain[,, i - n_burnin, j] <- beta[,,j]
      }
    }

    if (i <= n_burnin){
      burnin_pb$tick()
    }else{
      sample_pb$tick()
    }

    if (i == n_burnin)
    {
      message(paste0("\nDone! Beginning MCMC sampling for ", n_iter, " iterations."))

      sample_pb <- progress::progress_bar$new(format = "\r[:bar] :percent eta: :eta",
                                              clear = FALSE,
                                              total = round(n_iter / n_chains),
                                              width = 100)
    }

  }

  to_return <- list(Z_samples = Z_chain,
                    pi_samples = pi_chain,
                    beta_samples = beta_chain,
                    B_samples = B_chain,
                    Sigma_samples = Sigma_chain,
                    A_samples = A_chain,
                    mcar_samples = mcar_chain,
                    u_samples = u_chain,
                    phi_samples = phi_chain)

  return(to_return)
}
