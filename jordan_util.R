library(matrixcalc)
library(Matrix)
library(lmreg)
library(Rssa)

eps <- 5 * 1e-7
matrix_power <- function(A, N)
{
  if (N == 0) return (diag(nrow(A)))
  if (N == 1)  return (A)
  return (A %*% matrix_power(A, N - 1))
}

get_vec_basis <- function(vecs){
  eps2 <- nrow(vecs) * 2e-6
  print(paste0("precision: ", eps2))
  
  if (ncol(vecs) == 1){
    if (sum(abs(vecs)) < eps)
      return(list("basis" = numeric(0), "indices" = numeric(0)))
    return(list("basis" = vecs, "indices" = c(1)))
  }
  res <- numeric(0)
  indices <- numeric(0)
  has_vector <- FALSE
  for (i in 1:ncol(vecs))
    if (sum(abs(vecs[,i])) > eps){
      res <- as.matrix(vecs[,i], ncol = 1)
      indices <- c(indices, i)
      has_vector <- TRUE
      break
    }
  if (!has_vector)
    return(list("basis" = numeric(0), "indices" = numeric(0)))
  for (i in 1:ncol(vecs)){
    coefs <- qr.solve(res, as.matrix(vecs[,i], ncol = 1))
    if (sum(abs(vecs[,i])) > eps){
      if ( norm(abs(coefs)) < eps){
        indices <- c(indices, i)
        res <- cbind(res, vecs[,i])
      }
      else if ( norm(abs(res %*% coefs - as.matrix(vecs[,i], ncol = 1))) > eps2 ){
          print(paste0("Error: ", norm(abs(res %*% coefs - as.matrix(vecs[,i], ncol = 1)))))#
          indices <- c(indices, i)
          res <- cbind(res, vecs[,i])
        }
    }
  }
  
  return(list("basis" = res, "indices" = indices))
}

get_vec_basis2 <- function(vecs, eps2 = NA, need_vecs = NA){
  if (is.na(eps2))
    eps2 <- 2e-6 * nrow(vecs)
  if (ncol(vecs) == 1){
    if (sum(abs(vecs)) < eps)
      return(list("basis" = numeric(0), "indices" = numeric(0)))
    return(list("basis" = vecs, "indices" = c(1)))
  }
  res <- numeric(0)
  indices <- numeric(0)
  has_vector <- FALSE
  for (i in 1:ncol(vecs))
    if (sum(abs(vecs[,i])) > eps){
      res <- as.matrix(vecs[,i], ncol = 1)
      indices <- c(indices, i)
      has_vector <- TRUE
      break
    }
  if (!has_vector)
    return(list("basis" = numeric(0), "indices" = numeric(0)))
  min_error <- Inf
  for (i in 1:ncol(vecs)){
    coefs <- qr.solve(res, as.matrix(vecs[,i], ncol = 1))
    if (sum(abs(vecs[,i])) > eps){
      if ( norm(abs(coefs)) < eps){
        indices <- c(indices, i)
        res <- cbind(res, vecs[,i])
      }
      else if ( norm(abs(res %*% coefs - as.matrix(vecs[,i], ncol = 1))) > eps2 ){
        error <- norm(abs(res %*% coefs - as.matrix(vecs[,i], ncol = 1)))
        if (error < min_error)
          min_error <- error
        indices <- c(indices, i)
        res <- cbind(res, vecs[,i])
      }
    }
  }
  if (!is.na(need_vecs))
    if (length(indices) > need_vecs)
      return (get_vec_basis2(vecs, min_error * 1.01, need_vecs))
  
  return(list("basis" = res, "indices" = indices))
}

my_compbasis <- function(vecs, ndim){
  res <- numeric(0)
  if (length(vecs) == 0)
    res <- diag(1, nrow = ndim)
  else
    res <- get_coef(diag(1, nrow = nrow(vecs)), vecs, nrow(vecs) - ncol(vecs))
  return(res)
}

get_orthogonal_vec <- function(vecs, need_vecs = NA){
  return(my_compbasis(get_vec_basis2(vecs, NA, need_vecs)$basis, nrow(vecs)))
}

get_kernel_basis <- function(B, num){
  if (num == 0)
    return(numeric(0))
  if (num == 1)
    return(as.matrix(get_orthogonal_vec(t(B))[,1], ncol = 1))
  else
    return(get_orthogonal_vec(t(B))[,1:num]) 
}

#Функция ищет нужное количество векторов ортогональных базису basis2 из пространства, которое соответствует basis1 
#Используется для поиска базиса ядра оператора
#basis1 это вектора, по которым строим линейные комбинации
#basis2 это вектора которым линейная комбинация должна быть ортогональна
get_coef <- function(basis1, basis2, num_vecs){
  vecs <- basis1
  res_basis <- basis2
  
  added_basis <- numeric(0)
  S <- t(res_basis) %*% vecs
  
  types <- numeric(0)
  added_vecs <- 0
  index <- 1
  added_0_vecs <- numeric(0)
  
  while((index <= ncol(S)) & (added_vecs < num_vecs)){
    if (row_norm(S[,index]) < eps){
      added_vecs <- added_vecs + 1
      added_0_vecs <- cbind(added_0_vecs, vecs[,index])
    }
    index <- index + 1
  }
  
  if (added_vecs < num_vecs){
    while(added_vecs < num_vecs){
      S <- t(res_basis) %*% vecs
      
      rk_s <- rankMatrix(S, tol = 1e-7)
      if (rk_s == 1){
        #Берем ненулевой элемент
        n0_index <- numeric(0)
        for (i in 1:ncol(S)){
          if (abs(S[1,i]) > eps){
            n0_index <- i
            break
          }
        }
        coef <- rep(1, ncol(S))
        alpha <- -sum(S[1,-n0_index]) / S[1,n0_index]
        coef[n0_index] <- alpha
        
        
        res_basis <- cbind(res_basis, vecs %*% matrix(coef, ncol = 1))
        added_basis <- cbind(added_basis, vecs %*% matrix(coef, ncol = 1))
        added_vecs <- added_vecs + 1
        next
      }
      
      s_basis <- get_vec_basis2(S, NA, ncol(S) - (num_vecs - added_vecs))
      
      index1 <- 1
      not_basis_index <- 1
      
      while (not_basis_index <= nrow(basis1)){
        if (row_norm(S[,not_basis_index]) < eps)
          not_basis_index <- not_basis_index + 1
        else{
          if (index1 <= length(s_basis$indices)){
            if (s_basis$indices[index1] == not_basis_index){
              index1 <- index1 + 1
              not_basis_index <- not_basis_index + 1
            }
            else
              break
          }
          else
            break
        }
        
      }
      
      if (not_basis_index == nrow(basis1) + 1)
        break
      
      coef <- rep(0,ncol(S))
      coef[not_basis_index] <- 1
      
      coef_basis <- qr.solve(s_basis$basis, -S[,not_basis_index])
      for (i in 1:length(s_basis$indices)){
        coef[s_basis$indices[i]] <- coef_basis[i]
      }
      
      
      res_basis <- cbind(res_basis, vecs %*% matrix(coef, ncol = 1))
      added_basis <- cbind(added_basis, vecs %*% matrix(coef, ncol = 1))
      added_vecs <- added_vecs + 1
    }
  }
  return (cbind(added_basis,added_0_vecs))
  
}

#Эта функция берет векторы из basis1, которые линейно независимы с векторами из basis2
#Я пришел к выводу, что использовать эту функцию, а не поиск ортогональных векторов, для построения базиса более корректно
get_lin_ind <- function(basis1, basis2, num_vecs){
  added_num <- 0
  added_vecs <- numeric(0)
  for (i in 1:ncol(basis1)){
    if (added_num == num_vecs)
      break
    
    row_basis <- get_vec_basis2(t(basis2))
    basis2_corr <- t(row_basis$basis)
    b_corr <- basis1[row_basis$indices,i]
    
    #Линейная независимость
    coef <- matrix(qr.solve(basis2_corr, b_corr), ncol = 1)
    
    if (row_norm(basis2 %*% coef - basis1[,i]) > eps){
      added_vecs <- cbind(added_vecs, basis1[,i])
      basis2 <- cbind(basis2, basis1[,i])
      added_num <- added_num + 1
    }
  }
  return(added_vecs)
}


get_vecs <- function(B, res_basis, cell_sizes, ranks, k, m){
  
  B_curr <- correct_matrix(matrix_power(B, k))
  rk_curr <- rankMatrix(B_curr, tol=1e-7)
  
  B_prev <- correct_matrix(matrix_power(B, k - 1))
  rk_prev <- rankMatrix(B_prev, tol=1e-7)
  
  num_joined_curr <- numeric(0)
  if (length(res_basis) == 0){
    ranks <- c(rk_curr)
    num_joined_curr <- rk_prev - rk_curr
  }
  else{
    ranks <- c(ranks, rk_curr)
    num_joined_curr <- ranks[length(ranks)-1] - 2 * rk_curr + rk_prev
  }
  
  num_joined_prev <- numeric(0)
  if (length(res_basis) > 0){
    num_joined_prev <- ncol(res_basis)
    if (num_joined_prev == m)
      return(list("basis" = res_basis, "cell_size" = numeric(0), "ranks" = ranks))
  }
  else
    num_joined_prev <- 0 
  
  res_matrix <- numeric(0)
  if (k == 1){
    num_vecs <- min(num_joined_curr, m - num_joined_prev)
    
    if (num_vecs <= 0)
      ans <- list("basis" = res_basis, "cell_size" = numeric(0), "ranks" = ranks)
    else{
      vecs <- get_kernel_basis(B, nrow(B) - rk_curr)
      
      if (length(res_basis) == 0)
        return(list("basis" = vecs[,1:num_vecs], "cell_size" = rep(1, num_vecs), "ranks" = ranks))
      
      res_basis <- cbind(res_basis, get_lin_ind(vecs, res_basis, num_vecs))
      
      ans <- list("basis" = res_basis, "cell_size" = rep(1, num_vecs), "ranks" = ranks)
    }
    return(ans)
  }
  
  added_num <- 0
  full_block <- FALSE
  cell_size <- numeric(0)
  
  #Случай, когда ранг 0
  if ((length(res_basis) == 0) & (rk_curr == 0)){
    
    vecs <- get_kernel_basis(B_prev, nrow(B) - rk_prev)
    ort_vecs <- get_orthogonal_vec(vecs)
    
    for (i in 1:(rk_prev - rk_curr)){
      curr_block <- numeric(0)
      curr_block <- cbind(ort_vecs[,i], curr_block)
      added_num <- added_num + 1
      curr_cell_size <- 1
      
      basis_num <- numeric(0)
      if (length(res_basis) > 0)
        basis_num <- ncol(res_basis)
      else
        basis_num <- 0
      
      if (added_num + basis_num == m){
        full_block <- TRUE
        res_basis <- cbind(res_basis, curr_block)
        cell_size <- c(cell_size, curr_cell_size)
        break
      }
      
      for (j in 1:(k-1)){
        curr_block <- cbind(B %*% curr_block[,1], curr_block)
        curr_cell_size <- curr_cell_size + 1
        added_num <- added_num + 1
        if (added_num + basis_num == m){
          full_block <- TRUE
          break
        }
      }
      res_basis <- cbind(res_basis, curr_block)
      cell_size <- c(cell_size, curr_cell_size)
      
      if (full_block)
        break
    }
  }
  else if (length(res_basis) == 0){
    curr_vecs <- get_kernel_basis(B_curr, nrow(B) - rk_curr)
    prev_vecs <- get_kernel_basis(B_prev, nrow(B) - rk_prev)
    
    num_basis <- 0
    
    if (num_joined_curr > 0){
      for (i in 1:num_joined_curr){
        curr_block <- numeric(0)
        
        vec <- get_lin_ind(curr_vecs, prev_vecs, 1)
        
        curr_block <- cbind(vec, curr_block)
        prev_vecs <- cbind(prev_vecs, vec)
        added_num <- added_num + 1
        curr_cell_size <- 1
        
        if (length(res_basis) > 0)
          num_basis <- ncol(res_basis)
        else
          num_basis <- 0
        
        if (added_num + num_basis == m){
          full_block <- TRUE
          res_basis <- cbind(res_basis, curr_block)
          cell_size <- c(cell_size, curr_cell_size)
          break
        }
        
        for (j in 1:(k-1)){
          curr_block <- cbind(B %*% curr_block[,1], curr_block)
          curr_cell_size <- curr_cell_size + 1
          added_num <- added_num + 1
          if (added_num + num_basis == m){
            full_block <- TRUE
            break
          }
        }
        
        res_basis <- cbind(res_basis, curr_block)
        cell_size <- c(cell_size, curr_cell_size)
        
        if (full_block)
          break
        
      }
      
    }
  }
  else
  {
    #Здесь базис не пустой
    #Сейчас надо передать вектора с текущего уровня из базиса
    curr_vecs <- get_kernel_basis(B_curr, nrow(B) - rk_curr)#базис ядра текущей степени оператора
    prev_vecs <- get_kernel_basis(B_prev, nrow(B) - rk_prev)#базис ядра предыдущей степени оператора
    
    #Новый фрагмент
    #Достаем векторы текущего уровня
    indices <- numeric(0)
    shift <- 0
    for (i in 1:length(cell_sizes)){
      if (cell_sizes[i] > k - 1){
        indices <- c(indices, shift + k)
        shift <- shift + cell_sizes[i]
      }
      else 
        break
    }
    level_vectors <- res_basis[,indices]
    if (length(indices) == 1)
      level_vectors <- matrix(level_vectors, ncol = 1)
    prev_vecs <- cbind(prev_vecs, level_vectors)
    
    #Конец нового фрагмента
    
    if (num_joined_curr > 0)
      for (i in 1:num_joined_curr){
        curr_block <- numeric(0)
        
        vec <- get_lin_ind(curr_vecs, prev_vecs, 1)
        
        if (length(vec) == 0){
          break
        }
        curr_block <- cbind(vec, curr_block)
        prev_vecs <- cbind(prev_vecs, vec)
        added_num <- added_num + 1
        curr_cell_size <- 1
        
        if (added_num + ncol(res_basis) == m){
          full_block <- TRUE
          res_basis <- cbind(res_basis, curr_block)
          cell_size <- c(cell_size, curr_cell_size)
          break
        }
        
        for (j in 1:(k - 1)){
          curr_block <- cbind(B %*% curr_block[,1], curr_block)
          added_num <- added_num + 1
          curr_cell_size <- curr_cell_size + 1
          if (added_num + ncol(res_basis) == m){
            full_block <- TRUE
            break
          }
        }
        
        res_basis <- cbind(res_basis, curr_block)
        cell_size <- c(cell_size, curr_cell_size)
        if(full_block)
          break
      }
    
    
  }
  ans <- list("basis" = res_basis, "cell_size" = cell_size, "ranks" = ranks)
  return(ans)
}


# Возникла проблема, если оптимально вычислять здесь степени, то матрицы сильно отличаются от настоящих. Поэтому здесь неоптимально ищутся степени.
find_k <- function(B){
  k <- 1
  C <- correct_matrix(B)
  rk <- rankMatrix(C, tol=1e-7)
  prev_rk <- rk + 1
  B_curr <- B
  while(rk - prev_rk < 0){
    k <- k + 1
    B_curr <- correct_matrix(matrix_power(B, k))
    
    if (k == 50){
      print("Reached matrix power limit (50).")
      break
    }
    prev_rk <- rk
    rk <- rankMatrix(B_curr, tol = 1e-7)
  }
  return (k - 1)
}

row_norm <- function(x){
  return(norm(abs(as.matrix(x)), type="f"))
}

get_basis <- function(B, m){
  k <- find_k(B)
  res_basis <- numeric(0)
  cell_size <- numeric(0)
  ranks <- numeric(0)
  while(k > 0){
    tmp <- get_vecs(B, res_basis, cell_size, ranks, k, m)
    res_basis <- tmp$basis
    cell_size <- c(cell_size, tmp$cell_size)
    ranks <- tmp$ranks
    k <- k - 1
  }
  ans <- list("basis" = res_basis, "cell_size" = cell_size)
  return(ans)
}

correct_matrix <- function(B){
  eps_corr <- 0.001732
  bool_vec_row <- apply(B, 1, row_norm) < eps_corr
  bool_vec_col <- apply(B, 2, row_norm) < eps_corr
  C <- B
  C[bool_vec_row,] <- rep(0, ncol(C))
  C[,bool_vec_col] <- rep(0, nrow(C))
  return (C)
}

get_all_basis <- function(A, roots, m){
  basis <- numeric(0)
  cell_sizes <- list()
  for (i in 1:length(roots)){
    B <- correct_matrix(A - diag(nrow(A)) * roots[i])
    res <- get_basis(B, m[i])
    basis <- cbind(basis, res$basis)
    cell_sizes[[i]] <- res$cell_size
  }
  return(list("basis" = basis, "cell_sizes" = cell_sizes))
}

clusterize_roots <- function(roots, cl_eps = 3e-4) {
  N <- length(roots)
  are_used <- numeric(N)
  cluster_num <- 0
  next_center <- 1
  res_roots <- numeric(0)
  res_multiplicities <- numeric(0)
  indexes <- list()
  k <- 0
  
  
  while (!is.na(next_center)) {
    are_used[next_center] <- 1
    cluster_num <- cluster_num + 1
    res_roots <- c(res_roots, roots[next_center])
    k <- k + 1
    indexes[[k]] <- next_center
    res_multiplicities <- c(res_multiplicities, 1)
    center <- next_center
    next_center <- NA
    for (i in 1:N) {
      if ((! are_used[i]) & (Mod(roots[i] - roots[center]) < cl_eps)) {
        indexes[[k]] <- c(indexes[[k]], i)
        are_used[i] <- 1
        res_roots[cluster_num] <- res_roots[cluster_num] + roots[i]
        res_multiplicities[cluster_num] <- res_multiplicities[cluster_num] + 1
      }
      else if (!are_used[i])
        next_center <- i
    }
  }
  for (i in 1:cluster_num)
    res_roots[i] <- res_roots[i] / res_multiplicities[i]
  
  return (list("roots" = res_roots, "multiplicities" = res_multiplicities, "indexes" = indexes, "cluster_num" = cluster_num))
}  

get_roots_par <- function(ssa_obj, num_roots){
  roots <- parestimate(ssa_obj, groups = list(1:num_roots), method = "esprit")$roots
  res <- clusterize_roots(roots)
  return(res)
}

get_roots_eig <- function(ssa_obj, num_roots){
  w_length <- ssa_obj$window
  U <- ssa_obj$U[,1:num_roots]
  V <- ssa_obj$V[,1:num_roots]
  vals <- diag(ssa_obj$sigma[1:num_roots])
  
  P_low <- U[1:(w_length - 1),]
  P_up <- U[2:w_length,]
  M <- Rssa:::tls.solve(P_low, P_up) # Матрица сдвига
  
  dec <- eigen(M)
  roots <- dec$values
  res <- clusterize_roots(roots)
  return(res)
}