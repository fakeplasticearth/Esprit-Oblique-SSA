library(matrixcalc)
library(Matrix)
library(lmreg)

eps <- 5 * 1e-7
matrix_power <- function(A, N)
{
  if (N == 0) return (diag(nrow(A)))
  if (N == 1)  return (A)
  return (A %*% matrix_power(A, N - 1))
}

get_vec_basis <- function(vecs){
  
  if (ncol(vecs) == 1){
    if (sum(abs(vecs)) < eps)
      return(list("basis" = numeric(0), "indices" = numeric(0)))
    return(list("basis" = vecs, "indices" = c(1)))
  }
  res <- numeric(0)
  indices <- numeric(0)
  for (i in 1:ncol(vecs))
    if (sum(abs(vecs[,i])) > eps){
      res <- as.matrix(vecs[,i], ncol = 1)
      indices <- c(indices, i)
      break
    }
  for (i in 1:ncol(vecs)){
    coefs <- qr.solve(res, as.matrix(vecs[,i], ncol = 1))
    if (sum(abs(vecs[,i])) > eps){
      if ( norm(abs(coefs)) < eps){
        indices <- c(indices, i)
        res <- cbind(res, vecs[,i])
      }
      else if ( norm(abs(res %*% coefs - as.matrix(vecs[,i], ncol = 1))) > eps ){
        indices <- c(indices, i)
        N <- norm(abs(res %*% coefs - as.matrix(vecs[,i], ncol = 1)))
        res <- cbind(res, vecs[,i])
      }
    }
  }
  
  return(list("basis" = res, "indices" = indices))
}

my_compbasis <- function(vecs){
  res <- get_coef(diag(1, nrow = nrow(vecs)), vecs, nrow(vecs) - ncol(vecs))
  return(res)
}

get_orthogonal_vec <- function(vecs){
  return(my_compbasis(get_vec_basis(vecs)$basis))
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
      
      
      s_basis <- get_vec_basis(S)
      
      
      not_basis_index <- 1
      
      for (i in 1:length(s_basis$indices)){
        if (s_basis$indices[i] == not_basis_index)
          not_basis_index <- not_basis_index + 1
      }
      
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
    
    row_basis <- get_vec_basis(t(basis2))
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


get_vecs <- function(B, res_basis, ranks, k, m){
  
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
      
      #Старая версия, ищутся именно ортогональные вектора
      #res_basis <- cbind(res_basis, get_coef(vecs, res_basis, num_vecs))
      #Новая
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
    #Ситуация такая: базис еще пустой, ранг матрицы не 0
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
    curr_vecs <- get_kernel_basis(B_curr, nrow(B) - rk_curr)#базис ядра текущей степени оператора
    prev_vecs <- get_kernel_basis(B_prev, nrow(B) - rk_prev)#базис ядра предыдущей степени оператора
    
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

find_k <- function(B){
  k <- 1
  C <- correct_matrix(B)
  rk <- rankMatrix(C, tol=1e-7)
  prev_rk <- rk + 1
  while(rk - prev_rk < 0){
    k <- k + 1
    C <- correct_matrix(C %*% B)
    if (k == 50){
      print("Reached matrix power limit (50).")
      break
    }
    prev_rk <- rk
    rk <- rankMatrix(C, tol = 1e-7)
    
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
    tmp <- get_vecs(B, res_basis, ranks, k, m)
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