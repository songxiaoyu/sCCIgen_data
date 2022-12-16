#' Fit the marginal distributions for each row of a count matrix
#'
#' @keywords internal
#'
#' @param x            A matrix of shape p by n that contains count values.
#' @param marginal     Specification of the types of marginal distribution.
#'                     Default value is 'auto_choose' which chooses between ZINB, NB, ZIP
#'                     and Poisson by a likelihood ratio test (lrt) and whether there is
#'                     underdispersion.
#'                     'zinb' will fit the ZINB model. If there is underdispersion, it
#'                     will choose between ZIP and Poisson by a lrt. Otherwise, it will try to
#'                     fit the ZINB model. If in this case, there is no zero at all or an error
#'                     occurs, it will fit an NB model instead.
#'                     'nb' fits the NB model that chooses between NB and Poisson depending
#'                     on whether there is underdispersion.
#'                     'poisson' simply fits the Poisson model.
#' @param pval_cutoff  Cutoff of p-value of the lrt that determines whether
#'                     there is zero inflation.
#' @param epsilon      Threshold value for preventing the transformed quantile
#'                     to collapse to 0 or 1.
#' @param jitter       Logical, whether a random projection should be performed in the
#'                     distributional transform.
#' @param DT           Logical, whether distributional transformed should be performed.
#'                     If set to FALSE, the returned object u will be NULL.
#'
#' @return             a list with the following components:
#'\describe{
#'  \item{params}{a matrix of shape p by 3. The values of each column are: the ZI proportion,
#'  the dispersion parameter (for Poisson, it's Inf), and the mean parameter.}
#'  \item{u}{NULL or a matrix of the same shape as x, which records the transformed quantiles,
#'  by DT.}
#'}
fit_marginals <- function(x, marginal = c('auto_choose', 'zinb', 'nb', 'poisson'),
                          pval_cutoff = 0.05, epsilon = 1e-5,
                          jitter = TRUE, DT = TRUE){
  p <- nrow(x)
  n <- ncol(x)

  marginal <- match.arg(marginal)
  if(marginal == 'auto_choose'){
    params <- t(apply(x, 1, function(gene){
      m <- mean(gene)
      v <- var(gene)
      if(m >= v){
        mle_Poisson <- glm(gene ~ 1, family = poisson)
        tryCatch({
          mle_ZIP <- zeroinfl(gene ~ 1|1, dist = 'poisson')
          chisq_val <- 2 * (logLik(mle_ZIP) - logLik(mle_Poisson))
          pvalue <- as.numeric(1 - pchisq(chisq_val, 1))
          if(pvalue < pval_cutoff)
            c(plogis(mle_ZIP$coefficients$zero), Inf, exp(mle_ZIP$coefficients$count))
          else
            c(0.0, Inf, m)
        },
        error = function(cond){
          c(0.0, Inf, m)})
      }else{
        mle_NB <- MASS::glm.nb(gene ~ 1)
        if(min(gene) > 0)
          c(0.0, mle_NB$theta, exp(mle_NB$coefficients))
        else
          tryCatch({
            mle_ZINB <- zeroinfl(gene ~ 1|1, dist = 'negbin')
            chisq_val <- 2 * (logLik(mle_ZINB) - logLik(mle_NB))
            pvalue <- as.numeric(1 - pchisq(chisq_val, 1))
            if(pvalue < pval_cutoff)
              c(plogis(mle_ZINB$coefficients$zero), mle_ZINB$theta, exp(mle_ZINB$coefficients$count))
            else
              c(0.0, mle_NB$theta, exp(mle_NB$coefficients))
          },
          error = function(cond){
            c(0.0, mle_NB$theta, exp(mle_NB$coefficients))
          })
      }
    }))
  }else if(marginal == 'zinb'){
    params <- t(apply(x, 1, function(gene){
      m <- mean(gene)
      v <- var(gene)
      if(m >= v)
      {
        mle_Poisson <- glm(gene ~ 1, family = poisson)
        tryCatch({
          mle_ZIP <- zeroinfl(gene ~ 1|1, dist = 'poisson')
          chisq_val <- 2 * (logLik(mle_ZIP) - logLik(mle_Poisson))
          pvalue <- as.numeric(1 - pchisq(chisq_val, 1))
          if(pvalue < pval_cutoff)
            c(plogis(mle_ZIP$coefficients$zero), Inf, exp(mle_ZIP$coefficients$count))
          else
            c(0.0, Inf, m)
        },
        error = function(cond){
          c(0.0, Inf, m)})
      }
      else
      {
        if(min(gene) > 0)
        {
          mle_NB <- MASS::glm.nb(gene ~ 1)
          c(0.0, mle_NB$theta, exp(mle_NB$coefficients))
        }
        else
          tryCatch({
            mle_ZINB <- zeroinfl(gene ~ 1|1, dist = 'negbin')
            c(plogis(mle_ZINB$coefficients$zero), mle_ZINB$theta, exp(mle_ZINB$coefficients$count))
          },
          error = function(cond){
            mle_NB <- MASS::glm.nb(gene ~ 1)
            c(0.0, mle_NB$theta, exp(mle_NB$coefficients))
          })
      }
    }))
  }else if(marginal == 'nb'){
    params <- t(apply(x, 1, function(gene){
      m <- mean(gene)
      v <- var(gene)
      if(m >= v){
        c(0.0, Inf, m)
      }else{
        mle_NB <- MASS::glm.nb(gene ~ 1)
        c(0.0, mle_NB$theta, exp(mle_NB$coefficients))
      }
    }))
  }else if(marginal == 'poisson'){
    params <- t(apply(x, 1, function(gene){
      c(0.0, Inf, mean(gene))
    }))
  }

  if(DT){
    u <- t(sapply(1:p, function(iter){
      param <- params[iter, ]
      gene <- x[iter, ]
      prob0 <- param[1]
      u1 <- prob0 + (1 - prob0) * pnbinom(gene, size = param[2], mu = param[3])
      u2 <- (prob0 + (1 - prob0) * pnbinom(gene - 1, size = param[2], mu = param[3])) *
        as.integer(gene > 0)
      if(jitter)
        v <- runif(n)
      else
        v <- rep(0.5, n)
      r <- u1 * v + u2 * (1 - v)
      idx_adjust <- which(1-r < epsilon)
      r[idx_adjust] <- r[idx_adjust] - epsilon
      idx_adjust <- which(r < epsilon)
      r[idx_adjust] <- r[idx_adjust] + epsilon

      r
    }))
  }else{
    u <- NULL
  }

  return(list(params = params, u = u))
}



#' Fit a Gaussian copula model for a count matrix of a single cell type
#'
#' @inheritParams fit_marginals
#' @param zp_cutoff            The maximum propotion of zero allowed for a gene to be included
#'                             in the joint copula model.
#' @param min_non_zero_num     The minimum number of non-zero values required for a gene to be
#'                             fitted a marginal model.
#'
#' @return The genes of \code{x} will be partitioned into three groups. The first group contains
#' genes whose zero proportion is less than \code{zp_cutoff}. The second group contains genes
#' whose zero proportion is greater than \code{zp_cutoff} but still contains at least
#' \code{min_non_zero_num} non-zero values. The third and last group contains the rest of the
#' genes. For the first group, a joint Gaussian copula model will be fitted. For the second group,
#' only the marginal distribution of each gene will be fitted. For the last group, no model will
#' be fitted and only the index of these genes will be recorded. A list that contains the above
#' fitted model will be returned that contains the following components.
#' \describe{
#' \item{cov_mat}{The fitted covariance (or equivalently in this case, correlation) matrix of the
#' Gaussin copula model.}
#' \item{marginal_param1}{A matrix of the parameters for the marginal distributions of genes in
#' group one.}
#' \item{marginal_param2}{A matrix of the parameters for the marginal distributions of genes in
#' group two.}
#' \item{gene_sel1}{A numeric vector of the row indices of the genes in group one.}
#' \item{gene_sel2}{A numeric vector of the row indices of the genes in group two.}
#' \item{gene_sel3}{A numeric vector of the row indices of the genes in group three.}
#' \item{zp_cutoff}{Same as the input.}
#' \item{min_non_zero_num}{Same as the input.}
#' \item{sim_method}{A character string that says 'copula'. To be distinguished with the
#' (w/o copula) model.}
#' }
#'
# @export
fit_Gaussian_copula <- function(x, marginal = c('auto_choose', 'zinb', 'nb', 'poisson'),
                                jitter = TRUE, zp_cutoff = 0.8,
                                min_nonzero_num = 2){
  marginal <- match.arg(marginal)
  n <- ncol(x)
  p <- nrow(x)

  gene_zero_prop <- apply(x, 1, function(y){
    sum(y < 1e-5) / n
  })

  gene_sel1 <- which(gene_zero_prop < zp_cutoff)
  gene_sel2 <- which(gene_zero_prop < 1.0 - min_nonzero_num/n &
                       gene_zero_prop >= zp_cutoff)
  gene_sel3 <- (1:p)[-c(gene_sel1, gene_sel2)]

  if(length(gene_sel1) > 0){
    marginal_result1 <- fit_marginals(x[gene_sel1, , drop = FALSE], marginal, jitter = jitter, DT = TRUE)
    quantile_normal <- qnorm(marginal_result1$u)
    cov_mat <- cor(t(quantile_normal))
  }else{
    cov_mat = NULL
    marginal_result1 = NULL
  }

  if(length(gene_sel2) > 0){
    marginal_result2 <- fit_marginals(x[gene_sel2, , drop = FALSE], marginal, DT = FALSE)
  }else{
    marginal_result2 = NULL
  }
  return(list(cov_mat = cov_mat, marginal_param1 = marginal_result1$params,
              marginal_param2 = marginal_result2$params,
              gene_sel1 = gene_sel1, gene_sel2 = gene_sel2, gene_sel3 = gene_sel3,
              zp_cutoff = zp_cutoff, min_nonzero_num = min_nonzero_num,
              sim_method = 'copula', n_cell = n, n_read = sum(x)))
}



#' Fit a (w/o copula) model for a count matrix of a single cell type
#'
#' This function only fits the marginal distribution for each gene.
#'
#' @inheritParams fit_marginals
#' @inheritParams fit_Gaussian_copula
#' @return The genes of \code{x} will be partitioned into two groups. The first group contains
#' genes with at least \code{min_non_zero_num} non-zero values. The second group contains the
#' other genes. For the first group, the marginal distribution of each gene will be fitted.
#' For the last group, no model will be fitted and only the index of these genes will be
#' recorded. A list that contains the above fitted model will be returned that contains the
#' following components.
#' \describe{
#' \item{marginal_param1}{A matrix of the parameters for the marginal distributions of genes in
#' group one.}
#' \item{gene_sel1}{A numeric vector of the row indices of the genes in group one.}
#' \item{gene_sel2}{A numeric vector of the row indices of the genes in group two.}
#' \item{min_non_zero_num}{Same as the input.}
#' \item{sim_method}{A character string that says 'ind', short for 'independent'. To be
#' distinguished with the copula model.}
#' }
#'
# @export
fit_wo_copula <- function(x, marginal = c('auto_choose', 'zinb', 'nb', 'poisson'),
                          jitter = TRUE, min_nonzero_num = 2){
  marginal <- match.arg(marginal)
  n <- ncol(x)
  p <- nrow(x)

  gene_zero_prop <- apply(x, 1, function(y){
    sum(y < 1e-5) / n
  })

  gene_sel1 <- which(gene_zero_prop < 1.0 - min_nonzero_num/n)
  gene_sel2 <- (1:p)[-gene_sel1]

  if(length(gene_sel1) > 0){
    marginal_result1 <- fit_marginals(x[gene_sel1, ], marginal, jitter = jitter, DT = FALSE)
  }else{
    marginal_result1 = NULL
  }

  return(list(marginal_param1 = marginal_result1$params,
              gene_sel1 = gene_sel1, gene_sel2 = gene_sel2,
              min_nonzero_num = min_nonzero_num, sim_method = 'ind',
              n_cell = n, n_read = sum(x)))
}



#' Fit models for a count matrix
#'
#' @param data_mat      A matrix of shape p by n that contains count values. Each of its
#'                      column names should be the cell type names of that cell. Its column
#'                      names should also match \code{cell_type_sel}.
#' @param cell_type_sel A character vector that contains the selected cell types for which a
#'                      model will be fitted.
#' @param sim_method    Specification of the type of model.
#'                      Default value is 'copula', which selects the copula model.
#'                      'ind' will select the (w/o copula) model.
#' @inheritParams fit_Gaussian_copula
#' @param ncores        A numeric value that indicates the number of parallel cores for model
#'                      fitting. One core for each cell type.
#' @return A list with the same length as \code{cell_type_sel} that contains the fitted model
#' as each of its element.
#'
#' @export
fit_model_scDesign2 <- function(data_mat, cell_type_sel, sim_method = c('copula', 'ind'),
                                marginal = 'auto_choose',
                                jitter = TRUE, zp_cutoff = 0.8,
                                min_nonzero_num = 2, ncores = 1){
  sim_method <- match.arg(sim_method)
  marginal <- match.arg(marginal)

  if(sum(abs(data_mat - round(data_mat))) > 1e-5){
    warning('The entries in the input matrix are not integers. Rounding is performed.')
    data_mat <- round(data_mat)
  }

  if(sim_method == 'copula'){
    param <- parallel::mclapply(1:length(cell_type_sel), function(iter){
      fit_Gaussian_copula(data_mat[, colnames(data_mat) == cell_type_sel[iter]], marginal,
                          jitter = jitter, zp_cutoff = zp_cutoff,
                          min_nonzero_num = min_nonzero_num)
    }, mc.cores = ncores)
  }else if(sim_method == 'ind'){
    param <- mclapply(1:length(cell_type_sel), function(iter){
      fit_wo_copula(data_mat[, colnames(data_mat) == cell_type_sel[iter]], marginal,
                    jitter = jitter,
                    min_nonzero_num = min_nonzero_num)
    }, mc.cores = ncores)
  }

  names(param) <- cell_type_sel
  param
}
#' Simulate a count matrix for a single cell type based on a copula model
#'
#' @param copula_result A list that contains the parameters of a copula model.
#' @param n             An integer value that indicates the number of cells to generate.
#' @param marginal      A character string that indicates whether the generated values should
#'                      stay as discrete or switch to continuous. Default value is 'nb', which
#'                      should be used for generating a count marix. The alternative 'Gamma' is
#'                      only needed when this function is being called by other functions that
#'                      generate data with a user-specified sequencing depth. Normally, users
#'                      do not need to change this value.
#' @return A matrix of shape p by n that contains the simulated count values. p is derived from
#' \code{copula_result}
#'
# @export
simulate_count_copula <- function(copula_result, n = 100,
                                  marginal = c('nb', 'Gamma')){
  marginal <- match.arg(marginal)

  p1 <- length(copula_result$gene_sel1)
  if(p1 > 0){
    result1 <- MASS::mvrnorm(n = n, mu = rep(0.0, p1), Sigma = copula_result$cov_mat)
    result1 <- matrix(result1, nrow = n)
    result2 <- apply(result1, 2, pnorm)
    result2 <- matrix(result2, nrow = n)
  }


  p2 <- length(copula_result$gene_sel2)
  if(marginal == 'nb'){
    if(p1 > 0){
      result31 <- t(sapply(1:p1, function(iter){
        param <- copula_result$marginal_param1[iter, ]
        qnbinom(pmax(0.0, result2[, iter] - param[1]) / (1-param[1]),
                size = param[2], mu = param[3]) # simulate negative binomial dist
      }))
    }
    if(p2 > 0){
      result32 <- t(sapply(1:p2, function(iter){
        param <- copula_result$marginal_param2[iter, ]
        rbinom(n, 1, 1-param[1]) * rnbinom(n, size = param[2], mu = param[3])
      }))
    }
  }else if(marginal == 'Gamma'){
    if(p1 > 0){
      result31 <- t(sapply(1:p1, function(iter){
        param <- copula_result$marginal_param1[iter, ]
        qgamma(max(0.0, result2[, iter] - param[1]), shape = param[2], scale = param[3] / param[2])
      }))
    }
    if(p2 > 0){
      result32 <- t(sapply(1:p2, function(iter){
        param <- copula_result$marginal_param2[iter, ]
        rbinom(n, 1, 1-param[1]) * rgamma(n, shape = param[2], scale = param[3] / param[2])
      }))
    }
  }

  result <- matrix(0, nrow = p1 + p2 + length(copula_result$gene_sel3), ncol = n)
  if(p1 > 0){
    result[copula_result$gene_sel1, ] <- result31
  }
  if(p2 > 0){
    result[copula_result$gene_sel2, ] <- result32
  }
  result
}



#' Simulate a count matrix for a single cell type based on a (w/o copula model)
#'
#' @param model_params A list that contains the model parameters (can be either the copula
#'                     model or the (w/o copula) model).
#' @inheritParams simulate_count_copula
#' @return A matrix of shape p by n that contains the simulated count values. p is derived from
#' \code{model_params}.
#'
# @export
simulate_count_ind <- function(model_params, n = 100,
                               marginal = c('nb', 'Gamma')){
  marginal <- match.arg(marginal)

  if(model_params$sim_method == 'copula' || 'gene_sel3' %in% names(model_params)){
    p1 <- length(model_params$gene_sel1)
    p2 <- length(model_params$gene_sel2)
    if(marginal == 'nb'){
      if(p1 > 0){
        result31 <- t(sapply(1:p1, function(iter){
          param <- model_params$marginal_param1[iter, ]
          rbinom(n, 1, 1-param[1]) * rnbinom(n, size = param[2], mu = param[3])
        }))
      }
      if(p2 > 0){
        result32 <- t(sapply(1:p2, function(iter){
          param <- model_params$marginal_param2[iter, ]
          rbinom(n, 1, 1-param[1]) * rnbinom(n, size = param[2], mu = param[3])
        }))
      }
    }else if(marginal == 'Gamma'){
      if(p1 > 0){
        result31 <- t(sapply(1:p1, function(iter){
          param <- model_params$marginal_param1[iter, ]
          rbinom(n, 1, 1-param[1]) * rgamma(n, shape = param[2], scale = param[3] / param[2])
        }))
      }
      if(p2 > 0){
        result32 <- t(sapply(1:p2, function(iter){
          param <- model_params$marginal_param2[iter, ]
          rbinom(n, 1, 1-param[1]) * rgamma(n, shape = param[2], scale = param[3] / param[2])
        }))
      }
    }
    result <- matrix(0, nrow = p1 + p2 + length(model_params$gene_sel3), ncol = n)
    if(p1 > 0){
      result[model_params$gene_sel1, ] <- result31
    }
    if(p2 > 0){
      result[model_params$gene_sel2, ] <- result32
    }
  }else{
    p1 <- length(model_params$gene_sel1)
    p2 <- length(model_params$gene_sel2)
    result <- matrix(0, nrow = p1 + p2, ncol = n)
    if(p1 > 0){
      if(marginal == 'nb'){
        result31 <- t(sapply(1:p1, function(iter){
          param <- model_params$marginal_param1[iter, ]
          rbinom(n, 1, 1-param[1]) * rnbinom(n, size = param[2], mu = param[3])
        }))
      }else if(marginal == 'Gamma'){
        result31 <- t(sapply(1:p1, function(iter){
          param <- model_params$marginal_param1[iter, ]
          rbinom(n, 1, 1-param[1]) * rgamma(n, shape = param[2], scale = param[3] / param[2])
        }))
      }
      result[model_params$gene_sel1, ] <- result31
    }
  }
  result
}



#' Revise the initial scDesign2 functions to allow different sequencing depth
#'
#' @param model_params    A list with the same length as \code{cell_type_prop} that contains
#'                        the fitted model as each of its element (can be either the copula
#'                        model or the (w/o copula) model).
#' @param depth_simu_ref_ratio The (expected) sequencing depth ratio between simulated and refernece data.
#' @param n_cell_new      The total number of cells in the simulated count matrix.
#' @param cell_type_prop  The cell type proportion in the simulated count matrix.
#' @return A matrix of shape p by n that contains the simulated count values. p is derived from
#' \code{model_params}.
#'
#' @export


scDesign2.revised <- function(model_params, n_cell_new, cell_type_prop,
                              depth_simu_ref_ratio = NULL,
                              sim_method = c('copula', 'ind')){


  n_cell_vec <- sapply(model_params, function(x) x$n_cell)
  n_read_vec <- sapply(model_params, function(x) x$n_read)
  total_count_old <- sum(n_read_vec)
  n_cell_old      <- sum(n_cell_vec)

  if(length(model_params)!=length(cell_type_prop)){
    stop('Cell type proportion should have the same length as the number of models.')
  }

  n_cell_type <- length(cell_type_prop)
  cell_type_prop <- cell_type_prop / sum(cell_type_prop)
  n_cell_each <- round(cell_type_prop * n_cell_new)
  if(sum(n_cell_each) != n_cell_new){
    idx <- sample(n_cell_type, size = 1)
    n_cell_each[idx] <- n_cell_each[idx] + n_cell_new - sum(n_cell_each)
  }


  p <- length(model_params[[1]]$gene_sel1) + length(model_params[[1]]$gene_sel2) +
    length(model_params[[1]]$gene_sel3)
  new_count <- matrix(0, nrow = p, ncol = n_cell_new)


  r <- rep(depth_simu_ref_ratio, n_cell_type)

  for(iter in 1:n_cell_type)
    if(n_cell_each[iter] > 0){
      ulim <- sum(n_cell_each[1:iter])
      llim <- ulim - n_cell_each[iter] + 1
      params_new <- model_params[[iter]]
      params_new$marginal_param1[, 3] <- params_new$marginal_param1[, 3] * r[iter]
      params_new$marginal_param2[, 3] <- params_new$marginal_param2[, 3] * r[iter]
      if(sim_method == 'copula'){
        new_count[, llim:ulim] <- simulate_count_copula(model_params[[iter]],
                                                        n = n_cell_each[iter], marginal = 'Gamma')
      }else if(sim_method == 'ind'){
        new_count[, llim:ulim] <- simulate_count_ind(model_params[[iter]],
                                                     n = n_cell_each[iter], marginal = 'Gamma')
      }

    }
  if(is.null(names(model_params))){
    colnames(new_count) <- unlist(lapply(1:n_cell_type, function(x){rep(x, n_cell_each[x])}))
  }else{
    colnames(new_count) <- unlist(lapply(1:n_cell_type, function(x){
      rep(names(model_params)[x], n_cell_each[x])}))
  }
  return(new_count)

}



