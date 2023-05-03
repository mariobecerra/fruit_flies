simulate_mixture_choice_data = function(design_array, beta_vector, order, append_model_matrix = T, seed = NULL){
  
  dim_des_array = dim(design_array)
  q = dim_des_array[1]
  J = dim_des_array[2]
  S = dim_des_array[3]
  
  if(!is.null(seed)) set.seed(seed)
  
  # Only works for pairwise comparisons:
  model_tibble = map_df(1:S, function(s){
    Xs = opdesmixr::mnl_get_Xs(design_array, s = s, order = order)
    
    utilities = as.numeric(Xs %*% beta_vector)
    
    exp_utilities = exp(utilities)
    probs = exp_utilities/sum(exp_utilities)
    
    # stopifnot(sum(probs) - 1)
    if(abs(sum(probs) - 1) > 1e-8) stop("sum of probs numerically different to 1")
    
    # # Simulate Bernoulli trial with probability probs[1]
    # choose_first = rbernoulli(n = 1, p = probs[1])
    # # Vector of choices
    # choices = ifelse(choose_first, c(1, 0), c(0, 1))
    # Same as the following (but rnultinom generalizes to more classes):
    choice = as.numeric(rmultinom(n = 1, size = 1, prob = probs))
    
    out = as.data.frame(t(design_array[,,s])) %>% 
      set_names(paste0("comp_", 1:q)) %>% 
      as_tibble() %>% 
      mutate(
        choice_set = s,
        utilities = utilities,
        probs = probs, 
        choice = choice)
    
    if(append_model_matrix){
      out = out %>% 
        bind_cols(
          as.data.frame(Xs) %>% 
            set_names(paste0("model_mat_col", 1:ncol(.)))
        )
    }
    
    return(out)
  })
  
  return(model_tibble)
  
}




get_loglikelihood_repeated_respondents = function(design_array, Ny, order, beta_vec){
  
  dim_des_array = dim(design_array)
  q = dim_des_array[1]
  J = dim_des_array[2]
  S = dim_des_array[3]
  
  
  probs_vec = rep(NA_real_, J*S)
  
  for(s in 1:S){
    init = (s-1)*J + 1
    end = s*J
    probs_vec[init:end] = opdesmixr::mnl_get_Ps(design_array, beta_vec, s, order = order, transform_beta = F)
  }
  
  
  loglik = sum(Ny * log(probs_vec))
  
  return(loglik)
  
}






get_loglikelihood_gradient_repeated_respondents = function(design_array, Ny, order, beta_vec){
  
  dim_des_array = dim(design_array)
  q = dim_des_array[1]
  J = dim_des_array[2]
  S = dim_des_array[3]
  
  Ns = sapply(split(Ny, rep(1:S, each = 2)), sum)
  
  probs_vec = rep(NA_real_, J*S)
  model_matrix = matrix(NA_real_, nrow = J*S, ncol = length(beta_vec))
  
  for(s in 1:S){
    init = (s-1)*J + 1
    end = s*J
    probs_vec[init:end] = opdesmixr::mnl_get_Ps(design_array, beta_vec, s, order = order, transform_beta = F)
    model_matrix[init:end, ] = opdesmixr::mnl_get_Xs(design_array, s = s, order = order)
  }
  
  
  grad_p1 = (Ny - probs_vec*Ns)*model_matrix
  
  grad = as.numeric(apply(grad_p1, 2, sum))
  
  return(grad)
  
}




maximum_likelihood_est_repeated_respondents = function(design_array, Ny, order, starting_par, maxit = 500){
  
  minus_twice_log_likelihood_func = function(beta_param){
    return(
      -2*get_loglikelihood_repeated_respondents(
        design_array = design_array,
        Ny = Ny,
        beta_vec = beta_param,
        order = order
      )
    )
  }
  
  
  minus_twice_log_likelihood_func_grad = function(beta_param){
    return(
      -2*get_loglikelihood_gradient_repeated_respondents(
        design_array = design_array,
        Ny = Ny,
        beta_vec = beta_param,
        order = order
      )
    )
  }
  
  out = optim(par = starting_par, 
              fn = minus_twice_log_likelihood_func, 
              gr = minus_twice_log_likelihood_func_grad,
              control = list(maxit = maxit),
              method = "BFGS")
  
  return(out)
  
}









get_loglikelihood = function(design_array, y, order, beta_vec){
  # design_array: design array of dimension (q, J, S)
  # y: response (binary) vector of length J*S
  
  dim_des_array = dim(design_array)
  q = dim_des_array[1]
  J = dim_des_array[2]
  S = dim_des_array[3]
  
  # aaa = data.frame(
  #   choice_set = rep(1:S, each = J),
  #   probs = rep(NA_real_, J*S)
  # )
  
  probs_vec = rep(NA_real_, J*S)
  
  for(s in 1:S){
    init = (s-1)*J + 1
    end = s*J
    probs_vec[init:end] = opdesmixr::mnl_get_Ps(design_array, beta_vec, s, order = order, transform_beta = F)
  }
  
  
  loglik = sum(y * log(probs_vec))
  
  return(loglik)
  
}



get_loglikelihood_gradient = function(design_array, y, order, beta_vec){
  # design_array: design array of dimension (q, J, S)
  # y: response (binary) vector of length J*S
  
  dim_des_array = dim(design_array)
  q = dim_des_array[1]
  J = dim_des_array[2]
  S = dim_des_array[3]
  
  # aaa = data.frame(
  #   choice_set = rep(1:S, each = J),
  #   probs = rep(NA_real_, J*S)
  # )
  
  probs_vec = rep(NA_real_, J*S)
  model_matrix = matrix(NA_real_, nrow = J*S, ncol = length(beta_vec))
  
  for(s in 1:S){
    init = (s-1)*J + 1
    end = s*J
    probs_vec[init:end] = opdesmixr::mnl_get_Ps(design_array, beta_vec, s, order = order, transform_beta = F)
    model_matrix[init:end, ] = opdesmixr::mnl_get_Xs(design_array, s = s, order = order)
  }
  
  
  grad_p1 = (y - probs_vec)*model_matrix
  
  grad = as.numeric(apply(grad_p1, 2, sum))
  
  return(grad)
  
}




maximum_likelihood_est = function(design_array, y, order, starting_par, maxit = 500){
  minus_twice_log_likelihood_func = function(beta_param){
    return(
      -2*get_loglikelihood(
        design_array = design_array,
        y = y,
        beta_vec = beta_param,
        order = order
      )
    )
  }
  
  
  minus_twice_log_likelihood_func_grad = function(beta_param){
    return(
      -2*get_loglikelihood_gradient(
        design_array = design_array,
        y = y,
        beta_vec = beta_param,
        order = order
      )
    )
  }
  
  
  out = optim(par = starting_par, 
              fn = minus_twice_log_likelihood_func, 
              gr = minus_twice_log_likelihood_func_grad,
              control = list(maxit = maxit),
              method = "BFGS")
  
  return(out)
  
}








simulate_mixture_choice_data_no_choice = function(
  design_array, 
  beta_vector_mixture, 
  beta_no_choice,
  order, 
  append_model_matrix = T, 
  seed = NULL){
  
  dim_des_array = dim(design_array)
  q = dim_des_array[1]
  J = dim_des_array[2]
  S = dim_des_array[3]
  
  if(!is.null(seed)) set.seed(seed)
  
  # Only works for pairwise comparisons:
  model_tibble = map_df(1:S, function(s){
    
    Xs_0 = opdesmixr::mnl_get_Xs(design_array, s = s, order = order)
    
    Xs = rbind(
      cbind(Xs_0, rep(0, J)),
      c(rep(0, ncol(Xs_0)), 1)
    )
    
    
    utilities = as.numeric(Xs %*% c(beta_vector_mixture, beta_no_choice))
    
    exp_utilities = exp(utilities)
    probs = exp_utilities/sum(exp_utilities)
    
    # stopifnot(sum(probs) - 1)
    if(abs(sum(probs) - 1) > 1e-8) stop("sum of probs numerically different to 1")
    
    # # Simulate Bernoulli trial with probability probs[1]
    # choose_first = rbernoulli(n = 1, p = probs[1])
    # # Vector of choices
    # choices = ifelse(choose_first, c(1, 0), c(0, 1))
    # Same as the following (but rnultinom generalizes to more classes):
    choice = as.numeric(rmultinom(n = 1, size = 1, prob = probs))
    
    
    out = rbind(
      cbind(t(design_array[,,s]), rep(0, J)),
      c(rep(0, q), 1)
    ) %>% 
      as.data.frame() %>% 
      set_names(c(paste0("comp_", 1:q), "no_choice_ind")) %>% 
      as_tibble() %>% 
      mutate(
        choice_set = s,
        utilities = utilities,
        probs = probs, 
        choice = choice)
    
    if(append_model_matrix){
      out = out %>% 
        bind_cols(
          as.data.frame(Xs) %>% 
            set_names(paste0("model_mat_col", 1:ncol(.)))
        )
    }
    
    return(out)
  })
  
  return(model_tibble)
  
}







simulate_mixture_choice_data_no_choice_many_respondents = function(
  design_array, 
  beta_vector_mixture, 
  beta_no_choice,
  n_resp,
  order, 
  append_model_matrix = T, 
  seed = NULL){
  
  dim_des_array = dim(design_array)
  q = dim_des_array[1]
  J = dim_des_array[2]
  S = dim_des_array[3]
  
  if(!is.null(seed)) set.seed(seed)
  
  # Only works for pairwise comparisons:
  model_tibble = map_df(1:S, function(s){
    
    Xs_0 = opdesmixr::mnl_get_Xs(design_array, s = s, order = order)
    
    Xs = rbind(
      cbind(Xs_0, rep(0, J)),
      c(rep(0, ncol(Xs_0)), 1)
    )
    
    
    utilities = as.numeric(Xs %*% c(beta_vector_mixture, beta_no_choice))
    
    exp_utilities = exp(utilities)
    probs = exp_utilities/sum(exp_utilities)
    
    # stopifnot(sum(probs) - 1)
    if(abs(sum(probs) - 1) > 1e-8) stop("sum of probs numerically different to 1")
    
    # # Simulate Bernoulli trial with probability probs[1]
    # choose_first = rbernoulli(n = 1, p = probs[1])
    # # Vector of choices
    # choices = ifelse(choose_first, c(1, 0), c(0, 1))
    # Same as the following (but rmultinom generalizes to more classes):
    choices = as.numeric(rmultinom(n = 1, size = n_resp, prob = probs))
    
    
    out = rbind(
      cbind(t(design_array[,,s]), rep(0, J)),
      c(rep(0, q), 1)
    ) %>% 
      as.data.frame() %>% 
      set_names(c(paste0("comp_", 1:q), "no_choice_ind")) %>% 
      as_tibble() %>% 
      mutate(
        choice_set = s,
        utilities = utilities,
        probs = probs, 
        Ny = choices)
    
    if(append_model_matrix){
      out = out %>% 
        bind_cols(
          as.data.frame(Xs) %>% 
            set_names(paste0("model_mat_col", 1:ncol(.)))
        )
    }
    
    return(out)
  })
  
  return(model_tibble)
  
}








simulate_mixture_choice_data_no_choice_many_respondents_with_covariates = function(
  design_array, 
  beta_vector_mixture, 
  beta_no_choice,
  extra_covariates_beta,
  extra_covariates_model_matrix,
  n_resp,
  order, 
  append_model_matrix = T, 
  seed = NULL){
  
  dim_des_array = dim(design_array)
  q = dim_des_array[1]
  J = dim_des_array[2]
  S = dim_des_array[3]
  
  if(!is.null(seed)) set.seed(seed)
  
  # Only works for pairwise comparisons:
  model_tibble = map_df(1:S, function(s){
    
    Xs_0 = opdesmixr::mnl_get_Xs(design_array, s = s, order = order)
    
    Xs_1 = cbind(Xs_0, extra_covariates_model_matrix)
    
    Xs = rbind(
      cbind(Xs_1, rep(0, J)),
      c(rep(0, ncol(Xs_1)), 1)
    )
    
    
    utilities = as.numeric(Xs %*% c(beta_vector_mixture, extra_covariates_beta, beta_no_choice))
    
    exp_utilities = exp(utilities)
    probs = exp_utilities/sum(exp_utilities)
    
    # stopifnot(sum(probs) - 1)
    if(abs(sum(probs) - 1) > 1e-8) stop("sum of probs numerically different to 1")
    
    # # Simulate Bernoulli trial with probability probs[1]
    # choose_first = rbernoulli(n = 1, p = probs[1])
    # # Vector of choices
    # choices = ifelse(choose_first, c(1, 0), c(0, 1))
    # Same as the following (but rmultinom generalizes to more classes):
    choices = as.numeric(rmultinom(n = 1, size = n_resp, prob = probs))
    
    
    des_mat = cbind(t(design_array[,,s]), extra_covariates_model_matrix)
    
    out = rbind(
      cbind(des_mat, rep(0, J)),
      c(rep(0, ncol(des_mat)), 1)
    ) %>% 
      as.data.frame() %>% 
      set_names(c(paste0("comp_", 1:q), paste0("other_vars_", 1:ncol(extra_covariates_model_matrix)), "no_choice_ind")) %>% 
      as_tibble() %>% 
      mutate(
        choice_set = s,
        utilities = utilities,
        probs = probs, 
        Ny = choices)
    
    if(append_model_matrix){
      out = out %>% 
        bind_cols(
          as.data.frame(Xs) %>% 
            set_names(paste0("model_mat_col", 1:ncol(.)))
        )
    }
    
    return(out)
  })
  
  return(model_tibble)
  
}

