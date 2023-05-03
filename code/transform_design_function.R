
transform_design_flies = function(design_df, intensity_df) {
  # Inputs:
  # design_df: a data.frame (not tibble) with the design to be transformed. Its colnames should be c("R", "O", "Y", "G", "B", "P", "UV", "intensity", "cs"). The last column can be differently named, but it should identify the choice set and be numeric.
  # intensity_df: data.frame with the intensity values
  
  # Returns a list with the following objects:
  # design_orig: Original design but with one row identifying each choice set
  # designL_exp: The transformed values of the left side
  # designR_exp: The transformed values of the right side
  # design_software: A dataframe with the structure that Frank needs
  # design_with_mixture: A dataframe where each row is a choice set, and it includes the original design and the folder name where it's going to be saved. This is used for the mapping and reading of the counts later.
  
  
  q = 7
  n_pv = 1
  n_var = q + n_pv
  
  S = max(design_df[, n_var +1])
  J = nrow(design_df)/S
  
  if(J != 2) stop("J should be 2")
  if(ncol(design_df) != n_var+1) stop("Design should have 9 columns")
  
  ## Left
  
  design_origL <- design_df[seq(1, J*S, J), 1:n_var]
  colnames(design_origL) <- c('R_left', 'O_left', 'Y_left', 'G_left', 'B_left', 'P_left', 'UV1_left', 'intensityL')
  
  # calibration
  mL <- (intensity_df[n_var, 'L'] + intensity_df[n_var+1, 'L']) / 2
  deltaL <- (intensity_df[n_var, 'L'] - intensity_df[n_var+1, 'L']) / 2 
  design_origL$intensityL_actual <- design_origL$intensityL*deltaL + mL
  # experiment design left
  designL_exp <- data.frame(matrix(nrow = S, ncol = n_var-1))
  for (i in 1:S) {
    for (j in 1:(n_var-1)) {
      designL_exp[i,j] <- round(design_origL[i,j]*design_origL[i, 'intensityL_actual'] / intensity_df[j, 'L']*99, 2)
    }
  }
  colnames(designL_exp) <- c('R', 'O', 'Y', 'G', 'B', 'P', 'UV1')
  
  # original design right
  design_origR <- design_df[seq(2, J*S, J), 1:n_var]
  colnames(design_origR) <- c('R_right', 'O_right', 'Y_right', 'G_right', 'B_right', 'P_right', 'UV1_right', 'intensityR')
  # calibration right
  mR <- (intensity_df[n_var, 'R'] + intensity_df[n_var+1, 'R']) / 2
  deltaR <- (intensity_df[n_var, 'R'] - intensity_df[n_var+1, 'R']) / 2 
  design_origR$intensityR_actual <- design_origR$intensityR*deltaR + mR
  # experiment design right
  designR_exp <- data.frame(matrix(nrow = S, ncol = n_var-1))
  for (i in 1:S) {
    for (j in 1:(n_var-1)) {
      designR_exp[i,j] <- round(design_origR[i,j]*design_origR[i, 'intensityR_actual'] / intensity_df[j, 'R']*99, 2)
    }
  }
  colnames(designR_exp) <- c('R', 'O', 'Y', 'G', 'B', 'P', 'UV1')
  
  # combine 
  design_orig <- cbind(design_origL[,1:n_var], design_origR[,1:n_var])
  design_orig$chid <- seq(1,S,1)
  
  # # combine swapped
  # design_orig_swapped <- cbind(design_origR[,1:n_var], design_origL[,1:n_var])
  # colnames(design_orig_swapped) <- c('R_left', 'O_left', 'Y_left', 'G_left', 'B_left', 'P_left', 'UV1_left', 'intensityL', 'R_right', 'O_right', 'Y_right', 'G_right', 'B_right', 'P_right', 'UV1_right', 'intensityR')
  # design_orig_swapped$chid <- seq(1,S,1)
  
  
  
  
  
  
  
  
  
  expL <- data.frame(matrix(nrow = S, ncol = 2*n_var))
  colnames(expL) <- c('R', '', 'O', '', 'Y', '', 'G', '', 'B', '', 'P', '', 'UV1', '', 'UV2', '')
  # left: set values
  for (i in 1:S) {
    for (j in 1:(n_var-1)) {
      if (design_origL[i,j] == 0) {
        expL[i,2*j-1] = 'F'
      } else {
        expL[i,2*j-1] = 'T'
      }
      expL[i,2*j] = designL_exp[i,j]
    }}
  # left: drop UV2 - always F
  expL[, (2*n_var)-1] = 'F' 
  expL[, (2*n_var)] = 0
  
  # right
  expR <- data.frame(matrix(nrow = S, ncol = 2*n_var))
  colnames(expR) <- c('R', '', 'O', '', 'Y', '', 'G', '', 'B', '', 'P', '', 'UV1', '', 'UV2', '')
  for (i in 1:S) {
    for (j in 1:(n_var-1)) {
      if (design_origR[i,j] == 0) {
        expR[i,2*j-1] = 'F'
      } else {
        expR[i,2*j-1] = 'T'
      }
      expR[i,2*j] = designR_exp[i,j]
    }}
  expR[, (2*n_var)-1] = 'F' 
  expR[, (2*n_var)] = 0
  
  # convert design to fit the software: replicate rows to allow for wash-out period and set 'F 0'
  exp <- cbind(expL, expR)
  exp_rep <- exp[rep(seq_len(nrow(exp)), each = 2), ]
  for (i in seq(2, J*S, 2)) {
    for (j in seq(1, ncol(exp_rep), 2)) {
      exp_rep[i,j] = 'F'
    }
    for (k in seq(2, ncol(exp_rep), 2)) {
      exp_rep[i,k] = '0'
    }
  }
  exp_setting <- data.frame(rep(c(1200, 300), S), rep(c(10, 300), S), rep(80000, 2*S), rep(140000, 2*S))
  colnames(exp_setting) <- c('LED', 'Image', 'Exposure1', 'Exposure2')
  exp_rep <- cbind(exp_setting, exp_rep)
  
  # make connections between two designs
  # 1. round to integers
  expL_round <- data.frame(lapply(expL, function(x) if(is.numeric(x)) round(x, 0) else x))
  expR_round <- data.frame(lapply(expR, function(x) if(is.numeric(x)) round(x, 0) else x))
  colnames(expL_round) <- c('R', '', 'O', '', 'Y', '', 'G', '', 'B', '', 'P', '', 'UV1', '', 'UV2', '')
  colnames(expR_round) <- c('R', '', 'O', '', 'Y', '', 'G', '', 'B', '', 'P', '', 'UV1', '', 'UV2', '')
  exp_round <- cbind(expL_round, expR_round)
  # 2. append settings into a long string
  exp_round$name <- rep(' ', S)
  for (i in 1:S) {
    exp_round[i,4*n_var+1] <- paste(exp_round[i,], collapse = '')
  }
  exp_round$name <- substr(exp_round$name, 1, nchar(exp_round$name)-1)
  exp_round$chid <- seq(1,S,1)
  
  
  
  
  design_with_mixture <- merge(exp_round[, c('chid', 'name')], 
                               design_orig, by = 'chid')
  
  
  
  
  out <- list()
  out$design_orig <- design_orig
  # out$design_origL <- design_origL
  # out$design_origR <- design_origR
  # out$design_orig_swapped <- design_orig_swapped
  out$designL_exp <- designL_exp
  out$designR_exp <- designR_exp
  out$design_software <- exp_rep
  # out$design_exp_with_name <- exp_round
  out$design_with_mixture <- design_with_mixture
  
  
  return(out)
}

