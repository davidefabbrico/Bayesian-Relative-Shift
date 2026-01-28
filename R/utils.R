# Functions and MCMC algorithm for Bayesian Relative Shift Regression Model

load_or_install <- function(pkgs) {
  for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      ans <- readline(
        paste0("Package '", pkg, "' is not installed. Do you want to install it now? [y/n]: ")
      )
      if (tolower(ans) == "y") {
        install.packages(pkg, dependencies = TRUE)
      } else {
        stop(paste("Package", pkg, "is required but not installed. Stopping."))
      }
    }
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
}

required_pkgs <- c(
  "BayesLogit",
  "ggplot2",
  "dplyr",
  "tidyr",
  "igraph",
  "posterior",
  "doParallel",
  "progress"
)

load_or_install(required_pkgs)

abundance_plot <- function(count_matrix, title = ""){
  if(!require(ggplot2)){
    stop("please ensure the ggplot2 package is installed")
  }
  # sort matrix by column means - make pretty
  ss = count_matrix#[, order(colMeans(count_matrix), decreasing = T)]
  # coordinates
  xyz = cbind(expand.grid(1:dim(ss)[1], 1:dim(ss)[2]), as.vector(ss), as.vector(ss) > 0)
  names(xyz) = c("Sample.ID","Phylotype","Counts","Presence")
  print(ggplot2::ggplot(xyz, ggplot2::aes(y = Sample.ID, x = Phylotype, fill = log(Counts))) +
          ggplot2::geom_raster() + ggplot2::theme_bw() + ggplot2::ggtitle(title))
}

threshold <- function(X, t){
  n <- nrow(X)
  p <- ncol(X)
  
  X <- diag(1/apply(X, 1, sum), n, n)%*%X
  if((sum(t> apply(X, 1, max)) != 0)){
    stop("Threshold too large, some rows will be truncated to a zero row.")
  }
  if((sum(t > apply(X, 1, min)) == 0)){
    stop("Threshold too small, no zero entries.")
  }
  if(min(X) < 0){
    stop("Input is not non-negative")
  }
  X_t <- X*as.numeric(X >= t)
  X_t <- diag(1/apply(X_t, 1, sum), n, n)%*%X_t
  
  prop <- mean(apply(matrix(X_t==0, n, p), 1, sum)/p)
  cat("average zero proportion in each sample", prop, "\n")
  return(X_t)
}


gendata_scenario <- function(n = 500, p = 100, t = 0.005, prop = 0.8, scenario = 1,
                             family = "gaussian") {
  temp <- matrix(rnorm(n*(p-1)), nrow = n, ncol = (p-1))
  temp <- cbind(exp(temp), rep(1, n))
  X_true <- temp/apply(temp, 1, sum)
  X_obs <- threshold(X_true, t)
  X_tree <- NULL
  A <- NULL
  TT <- NULL
  a <- 1:10
  b <- rep(1, 10)
  c <- 1:5
  d <- rep(1, 20)
  taxonomy <- cbind(c(1:100), 
                    a%x%b, 
                    c%x%d, 
                    c(rep(1, 40), rep(2, 40), rep(3, 20)))
  if (scenario == 1) {
    beta <- c(rep(2, 40), rep(-2, 20), rep(1, 20), 
              c(1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5))
    binary_true_tree <- c(
      rep(0, 80),  
      rep(1, 20),  
      rep(0, 8),   
      1, 1,
      0, 0, 1, 1, 1,
      1, 1, 1
    )
  } else if (scenario == 2) {
    beta <- c(rep(c(1,2,3,4,5), 4), rep(2, 10), 
              rep(c(1,2,3,4,5), 6), rep(-2, 10),
              rep(c(1,2,3,4,5), 6))
    binary_true_tree <- c(
      rep(1, 20),
      rep(0, 10),
      rep(1, 30),
      rep(0, 10),
      rep(1, 48)
    )
  } else if (scenario == 3) {
    beta <- c(rep(0.5, 20), rep(c(0.5,0.6,0.7,0.8,0.9), 2), 
              rep(0.2, 10), rep(0.3, 20), 
              rep(0.4, 20), rep(c(0.5,0.6,0.7,0.8,0.9), 4))
    binary_true_tree <- c(rep(0, 20), rep(1, 10), rep(0, 10), 
                          rep(0, 40), rep(1, 20),
                          0, 0, 1, 1, 0, 0, 0, 0, 1, 1,
                          rep(1, 8))
  } else {
    stop("Invalid scenario number. Choose 1, 2, or 3.")
  }
  A <- expmat(taxonomy)$A
  M2 <- matnode(taxonomy)$M2
  TT <- treemat(M2)
  d <- ncol(TT)
  X_tree <- X_true%*%A
  if (family == "gaussian") {
    signal <- X_true %*% beta
    sigma2 <- var(signal)
    err <- rnorm(n, 0, sqrt(sigma2))
    y <- X_true%*%beta + rnorm(n, 0, sqrt(sigma2))
  } else if (family == "binomial") {
    eta <- X_true %*% beta
    p_y <- 1 / (1 + exp(-eta))
    y <- rbinom(n, size = 1, prob = p_y)
  }
  num_to_select <- round(n * prop)
  random_sample <- sample(1:n, num_to_select)
  return(list(y=y[random_sample], x=X_obs[random_sample,], beta = beta, xnc=temp[random_sample,],
              xt = X_tree[random_sample,], TT = TT, A = A, M2 = M2, ytest=y[-random_sample], 
              xtest=X_obs[-random_sample,], xnctest=temp[-random_sample,], 
              xttest = X_tree[-random_sample,], true_tree = binary_true_tree))
}


matnode <- function(M){
  p <- nrow(M)
  L <- ncol(M)
  
  #check if the finest level has unique taxa
  if(length(unique(M[,1])) < p){
    stop("the first level has overlapping taxa \n")
  }
  
  #check if M is really a tree
  for(j in 1:(L-1)){
    uniqid <- unique(M[,j])
    for(k in 1:length(uniqid)){
      if(length(unique(M[M[,j]==uniqid[k], j+1])) > 1){
        stop("the input matrix does not have a tree structure!")
      }
    }
  }
  
  if(length(unique(M[, L])) != 1){
    #if the last level is not all equal, add a root node
    M <- cbind(M, rep(1, p))
    L <- L+1
  }
  
  M1 <- M #node index matrix, index goes from 1 to # of nodes
  M1[, 1] <- c(1:p)
  for(j in 2:L){
    M1[,j] <-  match(M[,j], unique(M[,j]), nomatch = 0)
    #each column contains 1,2,3,... # unique groups; 1st column has 1~p, last column is all 1
    M1[,j] <- M1[,j]+max(M1[,j-1])
    #now the numbers are unique node indices, the last column values should be the total number of nodes in the tree
  }
  total_redun <- M1[1,L]#total number of nodes in the redundant tree
  #convert to a parent node vector
  node_redun <- rep(0, total_redun)
  
  for(i in 1:(total_redun - 1)){
    #[ind_i,ind_j] = ind2sub(size(M1),find(M1==i,1));
    #[r,c] = ind2sub(size(A), ind)
    ind <- which((M1==i)==1)[1]
    m <- dim(M1)[1]
    ind_i <- ((ind-1) %% m) + 1
    ind_j <- floor((ind-1) / m) + 1
    node_redun[i] <- M1[ind_i,ind_j+1] #the parent of node i is the node index to its right
  }
  
  # A trimmed/standardized version, parent node vector without redundancy
  M2 <- M # node index matrix, index goes from 1 to # of nodes
  M2[,1] <- c(1:p)
  for(j in 2:L){
    M2[,j] <-  match(M[,j], unique(M[,j]), nomatch = 0)
    tempseq <- M2[,j]
    for(k in max(tempseq):1){
      if(length(unique(M2[(tempseq==k), j-1])) == 1){
        tempseq[(which(tempseq==k)==1)] <- 0
        tempseq[tempseq > k] <- tempseq[tempseq > k] - 1
      }
    }
    tempseq[tempseq != 0] <- tempseq[tempseq != 0] + max(M2[, j-1]) #unique new node index
    tempseq[tempseq == 0] <- M2[tempseq == 0, j-1]
    M2[,j] <- tempseq
  }
  total <- M2[1, L]# total #of nodes
  # convert to a parent node vector
  node <- rep(0, total)
  for(i in 1:(total-1)){
    ind <- tail(which((M2==i)==1), n=1)
    m <- dim(M2)[1]
    ind_i <- ((ind-1) %% m) + 1
    ind_j <- floor((ind-1) / m) + 1
    node[i] <- M2[ind_i, ind_j] #the parent of node i is the node index to its right
  }
  return(list(node = node, M2=M2))
}

#' treemat
#' Function that takes as input compact tree representation and builds adjacency matrix 
#' @param M2: compact tree representation 
#' 

treemat <- function(M2){
  # prende in input la matrice M2 di matnode
  if(length(unique(M2[, ncol(M2)])) == 1){
    #tolgo il rootnode
    M2 <- M2[, -ncol(M2)]
  }
  
  mm <- max(M2)
  pp <- ncol(M2)
  
  TT <- matrix(0, mm, mm)
  for(i in 1:mm){
    #cat(which(M2 == i, arr.ind=TRUE), "\n")
    loc <- which(M2 == i, arr.ind=TRUE)
    for(j in 1:nrow(loc)){
      if(loc[j,2]==1){
        TT[i, (M2[loc[j, 1], (loc[j, 2] + 1)])] <- 1
      } else {
        if(loc[j,2]==pp){
          TT[i, (M2[loc[j, 1], (loc[j, 2] - 1)])] <- 1
        } else {
          TT[i, (M2[loc[j, 1], (loc[j, 2] + 1)])] <- 1
          TT[i, (M2[loc[j, 1], (loc[j, 2] - 1)])] <- 1
        }
      }
    }
  }
  return(TT)
}

#' ancmat
#' Function that takes as input compact tree representation and builds ancestors matrix 
#' @param M2: compact tree representation 
#' 

ancmat <- function(M2){
  # prende in input la matrice M2 di matnode
  if(length(unique(M2[, ncol(M2)])) == 1){
    #tolgo il rootnode
    M2 <- M2[, -ncol(M2)]
  }
  
  mm <- max(M2)
  pp <- ncol(M2)
  
  AN <- matrix(0, mm, mm)
  for(i in 1:mm){
    #cat(which(M2 == i, arr.ind=TRUE), "\n")
    loc <- which(M2 == i, arr.ind=TRUE)[1,]
    idx <- M2[loc[1], -c(1:loc[2])]
    
    AN[i, idx] <- 1
  }
  return(AN)
}

expmat <- function(M){
  mn <- matnode(M)
  node <- mn$node
  M2 <- mn$M2
  p <- nrow(M)
  
  numall <- length(node) # #of nodes in T
  numint <- numall-p # # internal node (including root)
  
  A <- matrix(0, p, numall)
  for(j in 1:p){
    parent <- unique(M2[j,])
    A[j, parent] <- 1
  }
  #remove root node
  A <- A[,-numall]
  return(list(A=A))
}

pgj <- function(gamma, MRF, gval, j, mu = -4, eta = 1.35){
  fgj <- mu + eta * c(MRF[j,]%*%gamma)
  num <- exp(gval*fgj)
  den <- 1 + exp(fgj)
  out <- num/den
  return(out)
}

mrf_prob <- function(y, X, gamma, j, MRF, mu = -4, eta = 2.0){
  p1 <- pgj(gamma, MRF, 1, j, mu = mu, eta = eta)
  p0 <- pgj(gamma, MRF, 0, j, mu = mu, eta = eta)
  out <- p1/(p1+p0)
  cat("out: ", out, "\n")
  return(out)
}

postpred_mpm <- function(x, y, xh, priors = c(5.0, 0.01, 0.01), mppi){
  n <- nrow(x)
  p <- sum(mppi>.5)
  M <- nrow(xh)
  
  x <- x[,c(mppi>.5)]
  xh <- xh[,c(mppi>.5)]
  
  sigma <- priors[1]
  a0 <- priors[2]
  b0 <- priors[3]
  mu0 <- rep(0, p)
  Lambda0 <- diag(sigma, p, p)
  
  Lambdan <- t(x)%*%x + Lambda0
  mun <- ginv(Lambdan)%*%(Lambda0%*%mu0 + t(x)%*%y)
  an <- a0 + n/2
  bn <- b0 + .5*(t(y)%*%y + mu0%*%Lambda0%*%mu0 - t(mun)%*%Lambdan%*%mun)
  V <- diag(1, M, M) + xh%*%Lambdan%*%t(xh)
  Sigma <- c(bn/an)*V
  m <- xh%*%mun
  ni <- 2*an
  
  #yhat <- mvtnorm::rmvt(M, mean = m, sigma = Sigma, df = Inf)
  
  yhat <- c(mvtnorm::rmvt(1, sigma=Sigma, df=ni)) + c(m)
  #matrix(rep(m, M), ncol=M, byrow=TRUE)
  return(yhat)
}

loglikelihood <- function(y, mu, sigma){
  return(sum(dnorm(y, mu, sqrt(sigma), log = TRUE)))
}

check_children <- function(node, M) {
  # funzione che mi ritorna i figli di un nodo
  indT <- which(M==node, arr.ind = T)
  ch <- c()
  for (c in 1:dim(indT)[1]) {
    if (indT[c, 2] <= 1) {
      return(NULL)
    } else {
      indChild <- indT[c, ] - c(0, 1)
      indChildRow <- indChild[1]
      indChildCol <- indChild[2]
      if (M[indChildRow, indChildCol] == 0) {
        return(NULL)
      } else {
        ch <- c(ch, M[indChildRow, indChildCol])
      } 
    }
  }
  return(ch)
}


leaf_parents <- function(matrice) {
  # funzione che mi dice chi sono i genitori delle foglie
  lp <- apply(matrice, 1, function(riga) {
    if (all(riga == 0)) {
      return(NA)
    }
    indici_non_zero <- which(riga != 0)
    if (length(indici_non_zero) >= 2) {
      return(riga[indici_non_zero[2]])
    } else {
      return(NA)
    }
  })
  return(lp)
}

current_model_leaf <- function(matrice) {
  # funzione che mi ritorna le foglie del modello corrente
  n_col <- ncol(matrice)
  leaves <- numeric(n_col)
  for (i in 1:nrow(matrice)) {
    first_leaf <- which(matrice[i, ] == 1)[1]
    if (!is.na(first_leaf)) {
      leaves[first_leaf] <- 1
    }
  }
  return(leaves)
}

leaf_children <- function(matrice) {
  # funzione che mi dice chi sono i figli delle foglie
  # o chi erano i genitori delle foglie e adesso sono diventati foglie
  lc <- apply(matrice, 1, function(riga) {
    if (!any(riga == 0)) {
      return(NA)
    }
    indice_non_zero <- which(riga != 0)[1]
    if (!is.na(indice_non_zero)) {
      return(riga[indice_non_zero])
    } else {
      return(NA)
    }
  })
  return(lc)
}

current_model <- function(matrice) {
  # funzione per ricreare i gamma
  res <- as.numeric(colSums(matrice) > 0)
  return(res)
}

get_siblings <- function(node, M) {
  # funzione che mi ritorna i fratelli di un dato nodo
  indNodoT <- which(M==node, arr.ind = T)
  if (dim(indNodoT)[1] == 0) {
    return("Node not found!")
  }
  for (s in dim(indNodoT)[1]) {
    indNodo <- indNodoT[s,]
    indPar <- indNodo + c(0, 1)
    if (indPar[2] > dim(M)[2]) {
      # se sono già alla radice, questi non hanno genitori
      Sib <- M[,indNodo[2]]
      return(unique(Sib[Sib != node]))
    } else {
      indParRow <- indPar[1]
      indParCol <- indPar[2]
      NodePar <- M[indParRow, indParCol]
      Sib <- M[M[, indParCol] == NodePar, indNodo[2]]
      return(unique(Sib[Sib != node]))
    }
  }
}

check_siblings_children <- function(node, M) {
  # funzione che mi ritorna true se selezionato un nodo qualsiasi
  # questo nodo ha fratelli che hanno figli, no altrimenti
  sib <- get_siblings(node, M)
  # questo if si potrebbe eliminare dal momento in cui non consideriamo 
  # catene, quindi nodi che non hanno fratelli.
  if (length(sib) == 0) {
    # se il nodo non ha fratelli
    return(FALSE)
  }
  for (s in sib) {
    indSib <- which(M==s, arr.ind = T)
    indSib <- indSib[1,]
    indChild <- indSib - c(0, 1)
    if (indChild[2] < 1) {
      return(FALSE)
    } else {
      indChildRow <- indChild[1]
      indChildCol <- indChild[2]
      FindChild <- M[indChildRow, indChildCol]
      if (FindChild != 0) {
        return(TRUE)
      } else {
        return(FALSE)
      }
    }
  }
}

count_rows_with_zero <- function(mat) {
  sum(apply(mat, 1, function(row) any(row == 0)))
}

move_tree <- function(Acurr, Mcurr, M_fix, N, move = "add") {
  ###########################################################################
  # devo mettere queste cose fuori l'if altrimenti non posso ritornare numAdd
  # chi erano i figli delle foglie o chi erano i genitori delle foglie 
  # e adesso sono diventati foglie
  leafChild <- leaf_children(Mcurr) # (Acurr, M_fix) # ADD
  # chi sono gli unici elementi?
  uniq_leaf_c <- unique(na.omit(leafChild))
  numAddOld <- length(uniq_leaf_c)
  # chi sono i genitori delle foglie?
  leafPar <- leaf_parents(Mcurr) # DELETE
  # guardo quanti elementi unici ci sono
  uniq_leaf_p <- unique(na.omit(leafPar))
  for (node in uniq_leaf_p) {
    children <- check_children(node, Mcurr)
    for (ch in children) {
      if (check_siblings_children(ch, Mcurr) | length(check_children(ch, Mcurr)) > 0) {
        uniq_leaf_p <- uniq_leaf_p[uniq_leaf_p != node]
      }
    }
  }
  numDeleteOld <- length(uniq_leaf_p)
  ##########################################################################
  if (move == "add") {
    # ADD
    # se l'insieme uniq_leaf è vuoto allora non posso fare add
    # altrimenti, se l'insieme ha solo un elemento allora prendo quello
    # altrimenti campiono N elementi in modo casuale
    if (numAddOld != 0) {
      if (numAddOld == 1) {
        samp <- uniq_leaf_c
      } else {
        # qui campiono casualmente il nodo da proporre
        depths <- sapply(uniq_leaf_c, function(j) {
          idx <- which(apply(Mcurr[,ncol(Mcurr):1], 2, function(col) j %in% col))
          if (length(idx) == 0) return(NA_integer_) else return(min(idx))
        })
        depths <- na.omit(depths)
        samp <- sample(uniq_leaf_c, N, replace = FALSE) #, prob = 1-(depths/sum(depths)))
      }
      # per ciascun elemento nel campione vado a modificare
      # la matrice M e la matrice A
      for (elem in samp) {
        # ritorno gli indici associati all'elemento elem nella matrice M_fix
        newSampind <- which(M_fix==elem, arr.ind = T)
        # prendo solo le righe
        idxrow <- newSampind[,1]
        # prendo le colonne ma diminuisco di 1 perché voglio inserire i figli
        colM <- (newSampind[,2]-1)[1]
        idxcol <- M_fix[idxrow, colM]
        figli <- M_fix[which(leafChild %in% elem), colM] 
        Mcurr[idxrow, colM] <- idxcol
        for (i in 1:length(idxrow)) {
          Acurr[idxrow[i], idxcol[i]] <- 1
        }
      }
    } else {print("No leaf children!")}
  } else if (move == "delete") {
    # DELETE
    # come prima, se l'insieme uniq_leaf è vuoto allora 
    # non posso fare nessun delete, altrimenti, se l'insieme contiene
    # almeno un elemento allora se ne contiene esattamente uno, lo prendo
    # altrimenti ne campiono N.
    if (numDeleteOld != 0) { 
      if (numDeleteOld == 1) {
        samp <- uniq_leaf_p
      } else {
        # qui campiono casualmente il nodo da proporre
        depths <- sapply(uniq_leaf_p, function(j) {
          idx <- which(apply(Mcurr[,ncol(Mcurr):1], 2, function(col) j %in% col))
          if (length(idx) == 0) return(NA_integer_) else return(min(idx))
        })
        depths <- na.omit(depths)
        samp <- sample(uniq_leaf_p, N, replace = FALSE) #, prob = depths/sum(depths))
      }
      # per ciascun elemento nel campione vado a modificare
      # la matrice M e la matrice A
      for (elem in samp) {
        newSampind <- which(M_fix == elem, arr.ind = T)
        idxrow <- newSampind[,1]
        colM <- (newSampind[,2] - 1)[1]
        idxcol <- M_fix[idxrow, colM]
        colM <- newSampind[1,2] - 1
        Mcurr[idxrow, colM] <- 0
        for (i in 1:length(idxrow)) {
          Acurr[idxrow[i], idxcol[i]] <- 0
        }
      }
    }
  }
  else {stop("Error: The provided move is invalid. Allowed moves are 'add' and 'delete'.")}
  # devo adesso ricalcolare quanti add e delete posso fare per
  # il calcolo delle probabilità di transizione
  leafChild <- leaf_children(Mcurr) # (Acurr, M_fix) # ADD
  # chi sono gli unici elementi?
  uniq_leaf_c <- unique(na.omit(leafChild))
  numAdd <- length(uniq_leaf_c)
  # chi sono i genitori delle foglie?
  leafPar <- leaf_parents(Mcurr) # DELETE
  # guardo quanti elementi unici ci sono
  uniq_leaf_p <- unique(na.omit(leafPar))
  for (node in uniq_leaf_p) {
    children <- check_children(node, Mcurr)
    for (ch in children) {
      if (check_siblings_children(ch, Mcurr) | length(check_children(ch, Mcurr)) > 0) {
        uniq_leaf_p <- uniq_leaf_p[uniq_leaf_p != node]
      }
    }
  }
  numDelete <- length(uniq_leaf_p)
  aD <- count_rows_with_zero(M_fix)
  return(list(Acurr, Mcurr, numAdd - aD, numDelete, numAddOld - aD, numDeleteOld))
}

tree_plot <- function(M, mppi) {
  # da inserire la radice nella matrice M
  edges <- c()
  for (col in 1:(ncol(M)-1)) {
    for (row in 1:nrow(M)) {
      edges <- c(edges, M[row, col], M[row, col + 1])
    }
  }
  edges <- matrix(edges, nrow = 2)
  edges <- edges[, !duplicated(t(apply(edges, 2, sort)))]
  g <- graph_from_edgelist(t(edges), directed = FALSE)
  V(g)$name <- as.character(1:vcount(g))
  node_probs <- mppi
  gray_palette <- colorRampPalette(c("gray90", "black"))
  node_colors <- gray_palette(100)[as.numeric(cut(node_probs, breaks = 100))]
  get_luminance <- function(color) {
    rgb <- col2rgb(color) / 255
    0.299 * rgb[1, ] + 0.587 * rgb[2, ] + 0.114 * rgb[3, ]
  }
  label_colors <- ifelse(get_luminance(node_colors) > 0.5, "black", "white")
  
  root_node <- M[, ncol(M)][1]
  layout <- layout_as_tree(g, root = root_node, circular = FALSE)
  
  edge_targets <- as.integer(ends(g, es = E(g))[,2])
  edge_alphas <- pmax(node_probs[edge_targets], 0.1)
  edge_colors <- alpha("black", edge_alphas)
  
  plot(g,
       layout = layout,
       vertex.color = node_colors,
       vertex.label = V(g)$name,
       vertex.label.color = label_colors,
       edge.color = edge_colors,
       edge.arrow.size = 0.5,
       vertex.size = 10
  )
}

tree_plot_minimal <- function(M, mppi, threshold = 0.5, size_node = 10) {
  edges <- c()
  for (col in 1:(ncol(M) - 1)) {
    for (row in 1:nrow(M)) {
      if (!is.na(M[row, col]) && !is.na(M[row, col + 1])) {
        edges <- c(edges, M[row, col], M[row, col + 1])
      }
    }
  }
  
  if (length(edges) == 0) {
    stop("No edges!")
  }
  
  edges <- matrix(edges, nrow = 2)
  edges <- edges[, !duplicated(t(apply(edges, 2, sort))), drop = FALSE]
  
  g <- graph_from_edgelist(t(edges), directed = FALSE)
  V(g)$name <- as.character(1:vcount(g))
  
  root_node <- M[1, ncol(M)]
  layout <- calculate_intersection_layout(g, root_node, M)
  
  node_probs <- mppi[1:vcount(g)]
  gray_palette <- colorRampPalette(c("gray90", "black"))
  node_colors <- gray_palette(100)[pmax(1, pmin(100, as.numeric(cut(node_probs, breaks = 100))))]
  node_lty <- ifelse(node_probs < threshold, 3, 1)
  
  x_range <- range(layout[, 1])
  y_range <- range(layout[, 2])
  
  x_margin <- 0
  y_margin <- 0
  
  options(repr.plot.width = 10, repr.plot.height = 12)
  
  plot(NA, 
       type = "n", 
       xlim = c(x_range[1] - x_margin, x_range[2] + x_margin),
       ylim = c(y_range[1] - y_margin, y_range[2] + y_margin),
       xlab = "", ylab = "", 
       axes = FALSE)

  draw_intersection_edges(g, layout, node_probs, threshold)
  
  points(layout[, 1], layout[, 2],
         pch = 21,
         bg = node_colors,
         col = "black",
         lwd = 1.5,
         # lty = node_lty,
         cex = size_node/8)
  
  invisible(list(
    graph = g,
    layout = layout,
    node_count = vcount(g),
    edge_count = ecount(g),
    root_node = root_node
  ))
}

calculate_intersection_layout <- function(g, root_node, M) {
  n_nodes <- vcount(g)
  layout <- matrix(0, nrow = n_nodes, ncol = 2)

  levels <- get_node_levels(g, root_node)
  max_level <- max(levels)
  
  vertical_spacing <- 4 
  
  leaves <- which(levels == max_level)
  leaf_order <- unique(M[, 1])  
  leaf_order <- leaf_order[leaf_order %in% leaves]
  
  n_leaves <- length(leaf_order)
  if (n_leaves > 1) {
    x_positions <- seq(-n_leaves, n_leaves, length.out = n_leaves)
  } else {
    x_positions <- 0
  }
  
  for (i in 1:n_leaves) {
    layout[leaf_order[i], 1] <- x_positions[i]
    layout[leaf_order[i], 2] <- -max_level * vertical_spacing
  }
  
  for (level in (max_level-1):0) {
    nodes_at_level <- which(levels == level)
    
    for (node in nodes_at_level) {
      children <- neighbors(g, node)
      children <- children[levels[children] > level]
      
      if (length(children) > 0) {
        children_x <- layout[children, 1]
        layout[node, 1] <- mean(children_x)
        
        layout[node, 2] <- -level * vertical_spacing
      } else {
        layout[node, 2] <- -level * vertical_spacing
      }
    }
  }
  
  return(layout)
}


draw_intersection_edges <- function(g, layout, node_probs, threshold) {
  # edge_list <- ends(g, es = E(g))
  edge_list <- get.edgelist(g)
  # print(get.edgelist(g))
  for (i in 1:nrow(edge_list)) {
    from_node <- as.numeric(edge_list[i, 1])
    to_node <- as.numeric(edge_list[i, 2])
    
    x1 <- layout[from_node, 1]
    y1 <- layout[from_node, 2]
    x2 <- layout[to_node, 1]
    y2 <- layout[to_node, 2]
    
    if (y1 > y2) {
      parent_x <- x1
      parent_y <- y1
      child_x <- x2
      child_y <- y2
      edge_prob <- node_probs[to_node]
    } else {
      parent_x <- x2
      parent_y <- y2
      child_x <- x1
      child_y <- y1
      edge_prob <- node_probs[from_node]
    }
    
    line_style <- if (edge_prob < threshold) 3 else 1
    alpha_val <- max(edge_prob, 0.1)
    line_color <- adjustcolor("black", alpha.f = alpha_val)
    
    segments(parent_x, parent_y, child_x, parent_y, col = line_color, lty = line_style, lwd = 1.5)
    segments(child_x, parent_y, child_x, child_y, col = line_color, lty = line_style, lwd = 1.5)
  }
}

get_node_levels <- function(g, root_node) {
  n_nodes <- vcount(g)
  levels <- rep(-1, n_nodes)
  levels[root_node] <- 0
  
  queue <- root_node
  while (length(queue) > 0) {
    current <- queue[1]
    queue <- queue[-1]
    
    neighbors <- neighbors(g, current)
    for (neighbor in neighbors) {
      if (levels[neighbor] == -1) {
        levels[neighbor] <- levels[current] + 1
        queue <- c(queue, neighbor)
      }
    }
  }
  
  return(levels)
}

hierarchical_mppi_threshold <- function(mppi_leaf, M, threshold = 0.5, use_max = FALSE) {
  # vorrei attuare questa strategia... la matrice M corrisponde alla struttura
  # dell'albero. mppi_leaf sono gli mppi delle foglie. Vorrei che mi riempissi 
  # la matrice M con gli mppi_leaf, però la matrice M può avere dei duplicati:
  # per esempio se 101 è padre di 30 e 32 allora ci sarà ripetuto due volte 101.
  # Arrivati a questo punto devo calcolare i vettori di inclusione 0 e 1.
  # Ci possono essere due strategie da considerare:
  # 1) Se la somma dei figli di un dato nodo superana la threshold allora includo
  #    tutti i figli altrimenti prendo la somma dei figli e la sommo al nodo padre
  #    e cosi via propagando verso la radice.
  # 2) Se almeno un nodo figlio supera la threshold allora includo tutti i figli.
  #    Se non la supera faccio come prima andando a sommare in questo caso non la somma
  #    ma il massimo
  nodes <- unique(c(M))
  n_nodes <- length(nodes)
  
  node2idx <- setNames(1:n_nodes, nodes)
  idx2node <- nodes
  
  node_level <- rep(NA, n_nodes)
  names(node_level) <- nodes
  children <- list()
  parents <- list()
  nlev <- ncol(M)
  
  for (lev in 1:nlev) {
    unique_in_lev <- unique(M[, lev])
    node_level[as.character(unique_in_lev)] <- lev
  }
  
  for (lev in 2:nlev) {
    for (i in 1:nrow(M)) {
      parent <- M[i, lev]
      child <- M[i, lev-1]
      p_char <- as.character(parent)
      c_char <- as.character(child)
      
      if (!p_char %in% names(children)) {
        children[[p_char]] <- character(0)
      }
      if (!c_char %in% children[[p_char]]) {
        children[[p_char]] <- c(children[[p_char]], c_char)
      }
      
      if (!c_char %in% names(parents)) {
        parents[[c_char]] <- character(0)
      }
      if (!p_char %in% parents[[c_char]]) {
        parents[[c_char]] <- c(parents[[c_char]], p_char)
      }
    }
  }
  
  values <- setNames(rep(0, n_nodes), nodes)
  
  root_nodes <- unique(M[, nlev])
  non_root_nodes <- setdiff(nodes, root_nodes)
  
  if (length(mppi_leaf) != length(non_root_nodes)) {
    stop("Error: mppi_leaf must have length equal to the number of non-root nodes")
  }
  
  values[as.character(non_root_nodes)] <- mppi_leaf
  
  I <- setNames(rep(0, n_nodes), nodes)
  
  get_ancestors <- function(node) {
    anc <- character(0)
    current <- node
    while (current %in% names(parents)) {
      current <- parents[[current]][1]
      anc <- c(anc, current)
    }
    return(anc)
  }
  
  levels <- sort(unique(node_level), decreasing = FALSE)
  
  for (lev in levels) {
    nodes_at_lev <- names(node_level)[node_level == lev]
    if (lev != 1) { 
      for (node in nodes_at_lev) {
        child_nodes <- as.numeric(children[[node]])
        child_values <- values[child_nodes]
        if (use_max) {
          agg_val <- max(child_values)
        } else {
          agg_val <- mean(child_values)
        }
        if (agg_val < threshold) {
          values[node] <- values[node] + agg_val
        } else {
          values[as.numeric(children[[node]])] <- agg_val
        }
      }
    }
  }
  values[length(values)] <- 1.0
  
  for (node in 1:(n_nodes-1)) { 
    if (values[node] >= threshold) {
      idxAnc <- get_ancestors(node)
      values[idxAnc] <- 1.0
    }
  }
  I <- as.integer(values >= threshold)
  return(list(values = values, I = I))
}


# selection metrics TPR, FPR, and MCC
confusion_matrix <- function(beta_true, beta_sel) {
  
  TP <- FP <- FN <- TN <- 0
  lb <- length(beta_true)
  for (i in 1:lb) {
    # TP
    if (beta_true[i] == 1 & beta_sel[i] == 1) {
      TP <- TP + 1
    } else if (beta_true[i] == 0 & beta_sel[i] == 1) {
      FP <- FP + 1
    } else if (beta_true[i] == 1 & beta_sel[i] == 0) {
      FN <- FN + 1
    } else if (beta_true[i] == 0 & beta_sel[i] == 0) {
      TN <- TN + 1
    }
  }
  
  mat <- matrix(0, nrow = 2, ncol = 2, dimnames = list("Actual" = c("Positive", "Negative"),
                                                       "Predicted" = c("Positive", "Negative")))
  # TP FN
  # FP TN
  mat[1,1] <- TP
  mat[1,2] <- FN
  mat[2,1] <- FP
  mat[2,2] <- TN

  return(mat)
}

sel_met <- function(conf_mat) {
  TP <- conf_mat[1, 1]
  FN <- conf_mat[1, 2]
  FP <- conf_mat[2, 1]
  TN <- conf_mat[2, 2]
  
  TPR <- if ((TP + FN) > 0) TP / (TP + FN) else 0
  FPR <- if ((FP + TN) > 0) FP / (FP + TN) else 0

  denom <- sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
  MCC <- if (denom > 0) {
    (TP * TN - FP * FN) / denom
  } else {
    0 
  }
  
  return(list(tpr = TPR, fpr = FPR, mcc = MCC))
}

plot_roc_curves <- function(FPR_mat, TPR_mat, 
                            n_points_smooth = 200,
                            curve_color = "blue", curve_alpha = 0.3, curve_size = 0.5,
                            mean_color = "red", mean_size = 1.2) {
  
  n_replicas <- ncol(TPR_mat)
  curves_list <- list()
  
  for (i in 1:n_replicas) {
    x <- c(0, FPR_mat[, i], 1)
    y <- c(0, TPR_mat[, i], 1)
    
    ord <- order(x)
    x <- x[ord]
    y <- y[ord]
    
    df_unique <- data.frame(x = x, y = y) %>%
      group_by(x) %>%
      summarize(y = mean(y), .groups = "drop")
    
    spline_fit <- splinefun(df_unique$x, df_unique$y, method = "monoH.FC")
    x_smooth <- seq(0, 1, length.out = n_points_smooth)
    y_smooth <- spline_fit(x_smooth)
    
    y_smooth[1] <- 0
    y_smooth[length(y_smooth)] <- 1
    
    curves_list[[i]] <- data.frame(x = x_smooth, y = y_smooth, curve = i)
  }
  
  df_all <- bind_rows(curves_list)
  
  df_mean <- df_all %>%
    group_by(x) %>%
    summarize(y = mean(y), .groups = "drop")
  
  df_mean$y[1] <- 0
  df_mean$y[nrow(df_mean)] <- 1
  
  ggplot() +
    geom_line(data = df_all, aes(x = x, y = y, group = curve), 
              color = curve_color, alpha = curve_alpha, linewidth = curve_size) +
    geom_line(data = df_mean, aes(x = x, y = y), 
              color = mean_color, linewidth = mean_size) +
    labs(x = "FPR", y = "TPR") +
    theme_minimal()
}

plot_metric_boxplot <- function(metric_mat, 
                                method_names = NULL,
                                x_label = "Method",
                                y_label = "Metric",
                                colors = NULL) {
  
  n_methods <- nrow(metric_mat)
  n_replicas <- ncol(metric_mat)
  
  if (is.null(method_names)) {
    method_names <- paste0("Method_", 1:n_methods)
  }
  
  okabe_ito <- c(
    "#E69F00", "#56B4E9", "#009E73", "#F0E442",
    "#0072B2", "#D55E00", "#CC79A7", "#000000"
  )
  
  if (is.null(colors)) {
    if (n_methods <= length(okabe_ito)) {
      colors <- okabe_ito[1:n_methods]
    } else {
      colors <- grDevices::colorRampPalette(okabe_ito)(n_methods)
    }
  } else {
    if (length(colors) != n_methods) {
      stop("The number of colors must be equal to the number of methods.")
    }
  }
  
  df_long <- metric_mat %>%
    as.data.frame() %>%
    mutate(Method = factor(method_names, levels = method_names)) %>%
    pivot_longer(cols = -Method, names_to = "Replica", values_to = "Value")
  
  ggplot(df_long, aes(x = Method, y = Value, fill = Method)) +
    geom_boxplot(alpha = 0.7, outlier.color = "black") +
    scale_fill_manual(values = colors) +
    labs(x = x_label, y = y_label) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "none")
}


compute_auc <- function(FPR, TPR) {
  FPR <- c(0, FPR, 1)
  TPR <- c(0, TPR, 1)
  ord <- order(FPR)
  FPR <- FPR[ord]
  TPR <- TPR[ord]
  auc <- sum((FPR[-1] - FPR[-length(FPR)]) * (TPR[-1] + TPR[-length(TPR)]) / 2)
  return(auc)
}

compute_BetaPrior_nodes <- function(M, tau = 1.0, kappa = 1.0, verbose = TRUE) {
  L <- ncol(M)
  alpha0_node <- list()
  beta0_node <- list()
  
  for (ell in 1:L) {
    # norm_ell <- ell/L
    exponential_decay <- exp(-tau*(ell))
    alpha_ell <-  kappa*exponential_decay
    beta_ell <-  kappa*(1-exponential_decay)
    nodes_at_level <- unique(M[, ell])
    for (node in nodes_at_level) {
      alpha0_node[[as.character(node)]] <- alpha_ell
      beta0_node[[as.character(node)]] <- beta_ell
    }
  }
  
  alpha0_node <- unlist(alpha0_node)
  beta0_node <- unlist(beta0_node)
  
  if (verbose) {
    cat("Beta prior mean values per level (root → leafs):\n")
    for (ell in 1:L) {
      # norm_ell <- ell/L
      exponential_decay <- exp(-tau*(ell))
      alpha_ell <- kappa*exponential_decay
      beta_ell <- kappa*(1-exponential_decay)
      cat(sprintf("  Level %d: mean = %.5f\n", ell, exponential_decay))
    }
  }
  
  return(list(alpha0_node, beta0_node))
}


get_parent <- function(M, node) {
  indexNode <- which(M == node, arr.ind = TRUE)
  if (indexNode[1, 2] == ncol(M)) {
    parent <- c()
  } else {
    indexWorking <- indexNode[1,]
    parent <- M[indexWorking[1], indexWorking[2] + 1]
  }
  return(parent)
}

log_prior_tree <- function(M, model, pi) {
  d <- length(model)
  log_prior <- 0
  
  for (j in 1:d) {
    parent <- get_parent(M, j)
    if (length(parent) == 0) {
      if (model[j] == 1) {
        log_prior <- log_prior + log(pi[j])
      } else {
        log_prior <- log_prior + log(1 - pi[j])
      }
    } else {
      if (model[parent] == 0 && model[j] == 1) {
        return(-Inf)
      }
      if (model[parent] == 1) {
        if (model[j] == 1) {
          log_prior <- log_prior + log(pi[j])
        } else {
          log_prior <- log_prior + log(1 - pi[j])
        }
      }
    }
  }
  
  return(log_prior)
}

accept_standard <- function(Xold, Xnew, Z, c_beta, c_theta, y, r21, r12, a0, 
                                  b0, n, log_prior_ratio, gprior = FALSE, g = NULL) {
  
  log_marg <- function(X_comp, Z, gprior, g, c_beta, c_theta) {
    term0 <- - (n/2) * log(2 * pi)
    
    # Crea matrice design combinata
    if (is.null(X_comp) || ncol(X_comp) == 0) {
      if (is.null(Z) || ncol(Z) == 0) {
        # modello nullo
        S <- drop(crossprod(y))
        log_lik <- term0 - (a0 + (n-1)/2) * log(b0 + 0.5 * S)
        return(log_lik)
      } else {
        # Solo non-composizionali
        X_design <- Z
        prior_prec <- diag(1/c_theta, ncol(Z))
      }
    } else {
      if (is.null(Z) || ncol(Z) == 0) {
        # Solo composizionali
        X_design <- X_comp
        prior_prec <- diag(1/c_beta, ncol(X_comp))
      } else {
        # Entrambi i tipi
        X_design <- cbind(X_comp, Z)
        prior_prec <- diag(c(rep(1/c_beta, ncol(X_comp)), rep(1/c_theta, ncol(Z))))
      }
    }
    
    k_total <- ncol(X_design)
    
    if (gprior) {
      # G-PRIOR (solo per parte composizionale)
      # Per semplicità, assumiamo g-prior solo per composizionali, ridge per non-composizionali
      stop("G-prior not implemented for mixed compositional/non-compositional case")
    } else {
      # ridge prior per tutti i coefficienti
      XtX <- crossprod(X_design)
      XtY <- crossprod(X_design, y)
      
      V0_inv <- prior_prec
      Vn_inv <- XtX + V0_inv
      
      tryCatch({
        Vn_inv_chol <- chol(Vn_inv)
        log_det_Vn_inv <- 2 * sum(log(diag(Vn_inv_chol)))
        Vn <- chol2inv(Vn_inv_chol)
        
        mn <- Vn %*% XtY
        yty <- drop(crossprod(y))
        m0V0m0 <- 0  # m_0 = 0
        mnVnmn <- drop(t(mn) %*% Vn_inv %*% mn)
        bn <- b0 + 0.5 * (yty + m0V0m0 - mnVnmn)
        
        # Log-determinante della prior
        if (is.null(X_comp) || ncol(X_comp) == 0) {
          log_det_V0 <- ncol(Z) * log(c_theta)
        } else if (is.null(Z) || ncol(Z) == 0) {
          log_det_V0 <- ncol(X_comp) * log(c_beta)
        } else {
          log_det_V0 <- ncol(X_comp) * log(c_beta) + ncol(Z) * log(c_theta)
        }
        
        log_det_Vn <- -log_det_Vn_inv
        
        term1 <- 0.5 * (log_det_Vn - log_det_V0)
        term2 <- -(a0 + (n-1)/2) * log(bn)
        
        return(term0 + term1 + term2)
        
      }, error = function(e) {
        return(-Inf)
      })
    }
  }
  
  log_post_old <- log_marg(Xold, Z, gprior, g, c_beta, c_theta)
  log_post_new <- log_marg(Xnew, Z, gprior, g, c_beta, c_theta)
  
  log_accept <- (log_post_new - log_post_old) + (log(r21) - log(r12)) + log_prior_ratio
  
  return(min(log_accept, 700))
}

# accept_standard <- function(Xold, Xnew, c, y, r21, r12, a0, 
#                             b0, n, log_prior_ratio, gprior = FALSE, g = NULL) {
#   
#   log_marg <- function(X, gprior, g, c) {
#     term0 <- - (n/2) * log(2 * pi)
#     
#     if (is.null(X) || ncol(X) == 0) {
#       # modello nullo
#       S <- drop(crossprod(y))
#       log_lik <- term0 - (a0 + (n-1)/2) * log(b0 + 0.5 * S)
#       return(log_lik)
#     }
#     
#     k <- ncol(X)
#     
#     if (gprior) {
#       # G-PRIOR
#       XtX <- crossprod(X)
#       XtX_inv <- solve(XtX + diag(1e-8, k))
#       beta_ols <- XtX_inv %*% crossprod(X, y)
#       yPy <- drop(crossprod(y, X %*% beta_ols))
#       S <- drop(crossprod(y)) - (g/(g+1)) * yPy
#       log_lik <- term0 - 0.5 * k * log(1+g) - (a0 + (n-1)/2) * log(b0 + 0.5 * S)
#       return(log_lik)
#     } else {
#       # ridge prior
#       XtX <- crossprod(X)
#       XtY <- crossprod(X, y)
#       
#       V0_inv <- diag(1/c, k)
#       Vn_inv <- XtX + V0_inv
#       
#       tryCatch({
#         Vn_inv_chol <- chol(Vn_inv)
#         log_det_Vn_inv <- 2 * sum(log(diag(Vn_inv_chol)))
#         Vn <- chol2inv(Vn_inv_chol)
#         
#         mn <- Vn %*% XtY
#         yty <- drop(crossprod(y))
#         m0V0m0 <- 0  # m_0 = 0
#         mnVnmn <- drop(t(mn) %*% Vn_inv %*% mn)
#         bn <- b0 + 0.5 * (yty + m0V0m0 - mnVnmn)
#         
#         log_det_V0 <- k * log(c)
#         log_det_Vn <- -log_det_Vn_inv
#         
#         term1 <- 0.5 * (log_det_Vn - log_det_V0)
#         term2 <- -(a0 + (n-1)/2) * log(bn)
#         
#         return(term0 + term1 + term2)
#         
#       }, error = function(e) {
#         return(-Inf)
#       })
#     }
#   }
#   
#   log_post_old <- log_marg(Xold, gprior, g, c)
#   log_post_new <- log_marg(Xnew, gprior, g, c)
#   
#   log_accept <- (log_post_new - log_post_old) + (log(r21) - log(r12)) + log_prior_ratio
#   
#   return(min(log_accept, 700))
# }

fun_adaptive_moves <- function(depth, lambda) {
  if (lambda > 0) {
    out <- (exp(-lambda*depth) - exp(-lambda))/(1 - exp(-lambda))
  } else if (lambda == 0) {
    out <- 1 - depth
  }
  return(out)
}

adaptiveMoves <- function(A, M, lambda) {
  active_leaves <- which(current_model_leaf(A) == 1)
  if (length(active_leaves) > 0) {
    depths <- sapply(active_leaves, function(j) {
      idx <- which(apply(M[,ncol(M):1], 2, function(col) j %in% col))
      if (length(idx) == 0) return(NA_integer_) else return(min(idx))
    })
    depths <- na.omit(depths)
    nEntropy <- (max(depths) -  min(depths))/ncol(M) # mean(unique(depths))/ncol(M)
  } else {
    nEntropy <- 0.5
  }
  pAdd <- fun_adaptive_moves(nEntropy, lambda)
  listOne <- list(add = pAdd, metric = nEntropy)
  return(listOne)
}

relativeShiftMCMC_gaussian <- function(y, X, Z = NULL, A, M, priors = c(5.0, 10.0, 2.0, 2.0, 1.0, 1.0), 
                                       N = 1, num_iter, num_burn, num_thin, 
                                       Xnew = NULL, Znew = NULL, method = "OLS",
                                       seed = NULL, gprior = FALSE, g = NULL,
                                       adaptive_levels = TRUE, tau = 0.5, kappa = 1,
                                       adaptive_moves = TRUE, lambda = 1) {
  
  # 1. PARAMETERS INITIALISATION -----------------------------------------------
  if (!is.null(seed)) set.seed(seed)
  
  if (gprior && is.null(g)) {
    stop("If gprior=TRUE, must be specified the value of g")
  }
  
  n <- nrow(X)
  d <- ncol(X)
  
  # IF Z = NULL
  if (is.null(Z)) {
    p_theta <- 0
    Z <- matrix(0, nrow = n, ncol = 0)
  } else {
    p_theta <- ncol(Z)
  }
  
  ntest <- if (!is.null(Xnew)) nrow(Xnew) else 0
  
  c_beta <- priors[1]
  c_theta <- priors[2]
  a_sigma <- priors[3]
  b_sigma <- priors[4]
  alpha0 <- priors[5]
  beta0 <- priors[6]
  
  # INITIALISATION FROM PRIOR
  sigma <- 1 / rgamma(1, a_sigma, b_sigma)
  
  if (adaptive_levels) {
    betaPrior_vec <- compute_BetaPrior_nodes(M, tau = tau, kappa = kappa)
    alpha0_vec <- betaPrior_vec[[1]]
    beta0_vec <- betaPrior_vec[[2]]
  } else {
    alpha0_vec <- rep(alpha0, d)
    beta0_vec <- rep(beta0, d)
  }
  
  pi <- sapply(1:d, function(i) rbeta(1, alpha0_vec[i], beta0_vec[i]))
  
  # INITIALISATION THETA NON-COMPOSITIONAL COVARIATES
  if (p_theta > 0) {
    theta <- t(as.matrix(rnorm(p_theta, 0, sqrt(c_theta))))
  } else {
    theta <- numeric(0)
  }
  
  # INITIAL TREE STRUCTURE
  model <- current_model(A)
  A_fix <- A_curr <- A
  M_fix <- M_curr <- M
  
  l <- which(model == 1)
  k <- length(l)
  
  # INITIALISATION BETA COMPOSITIONAL COVARIATES
  Xsel <- X[, l, drop = FALSE]
  if (method == "OLS" && k > 0) {
    if (p_theta > 0) {
      y_res <- y - Z %*% t(theta)
    } else {
      y_res <- y
    }
    beta <- MASS::ginv(t(Xsel) %*% Xsel) %*% (t(Xsel) %*% y_res)
  } else {
    beta <- rnorm(d, mean = 0, sd = sqrt(c_beta * sigma))
  }
  
  # MOVE COUNTER
  resMove <- move_tree(A_curr, M_curr, M_fix, N, move = "delete")
  numAdd_curr <- resMove[[5]]
  numDelete_curr <- resMove[[6]]
  
  # 2. SETUP STORAGE -----------------------------------------------------------
  num_samples <- as.integer(floor((num_iter - num_burn) / num_thin))
  model_out <- model_out_leaf <- beta_out <- matrix(0, nrow = num_samples, ncol = d)
  if (p_theta > 0) {
    theta_out <- matrix(0, nrow = num_samples, ncol = p_theta)
  } else {
    theta_out <- NULL
  }
  sigma_out <- pAccept <- cMove <- meanD <- rep(0, num_samples)
  pi_out <- matrix(0, nrow = num_samples, ncol = d)
  ll_out <- rep(0, num_iter)
  Asamples <- vector("list", num_samples)
  Msamples <- vector("list", num_samples)
  
  if ((!is.null(Xnew) && p_theta == 0) || (!is.null(Xnew) && !is.null(Znew) && p_theta > 0)) {
    ypred_out <- matrix(0, nrow = ntest, ncol = num_samples)
  } else {
    ypred_out <- NULL
  }
  
  it <- 1
  betmodelac <- 0
  start_time <- proc.time()
  
  # LOADING BAR
  pb <- progress::progress_bar$new(
    format = "  MCMC Progress [:bar] :percent | ETA: :eta | Elapsed: :elapsed | Acceptance: :acc_rate",
    total = num_iter, 
    clear = FALSE,
    width = 80,
    show_after = 0
  )
  
  # 3. MCMC LOOP ---------------------------------------------------------------
  for (t in 1:num_iter) {
    
    pb$tick(tokens = list(acc_rate = sprintf("%.3f", betmodelac/max(1, t-1))))
    
    
    # if ((t %% 500) == 0) {
    #   cat("Iteration:", t, "Acceptance rate:", round(betmodelac/t, 4), "\n")
    # }
    
    # a. Residual for compositional covariates
    if (p_theta > 0) {
      y_res <- y - Z %*% t(theta)
    } else {
      y_res <- y
    }
    
    # b. New move proposal for compositional covariates
    if (adaptive_moves) {
      pEnt <- adaptiveMoves(A_curr, M_curr, lambda)
      pAdd <- pEnt$add
      metric <- pEnt$metric
    } else {
      pAdd <- 0.5
      metric <- 0
    }
    
    if (numDelete_curr == 0) {
      move <- "add"
    } else if (numAdd_curr == 0) {
      move <- "delete"
    } else {
      move <- sample(c("add", "delete"), 1, prob = c(pAdd, 1 - pAdd))
    }
    
    resMove <- move_tree(A_curr, M_curr, M_fix, N, move = move)
    A_prop <- resMove[[1]]
    M_prop <- resMove[[2]]
    numAdd_prop <- resMove[[3]]
    numDelete_prop <- resMove[[4]]
    numAdd_curr <- resMove[[5]]
    numDelete_curr <- resMove[[6]]
    
    proposemodel <- current_model(A_prop)
    
    proposal_ratio <- 1
    if (adaptive_moves) {
      pEnt_prop <- adaptiveMoves(A_prop, M_prop, lambda)
      pAdd_prop <- pEnt_prop$add
      metric_prop <- pEnt_prop$metric
      if (pAdd > 0 && pAdd < 1 && pAdd_prop > 0 && pAdd_prop < 1) {
        if (move == "add") {
          proposal_ratio <- (1 - pAdd_prop) / pAdd
        } else {
          proposal_ratio <- pAdd_prop / (1 - pAdd)
        }
      }
    }
    
    # c. ratio proposal
    r12 <- 1 / ifelse(move == "add", numAdd_curr, numDelete_curr)
    r21 <- 1 / ifelse(move == "add", numDelete_prop, numAdd_prop)
    
    # d. Evaluation of the proposed model
    newl <- which(proposemodel == 1)
    newk <- length(newl)
    Xselnew <- X[, newl, drop = FALSE] 
    
    # e. ratio prior
    log_prior_old <- log_prior_tree(M_fix, model, pi)
    log_prior_new <- log_prior_tree(M_fix, proposemodel, pi)
    log_prior_ratio <- log_prior_new - log_prior_old
    
    # f. Acceptance log-probability
    acceptMH <- accept_standard(Xsel, Xselnew, Z, c_beta, c_theta, y, 
                                r21, r12, a_sigma, b_sigma, n, log_prior_ratio,
                                gprior, g)
    # acceptMH <- accept_standard(Xsel, Xselnew, c_beta, y_res, 
    #                             r21, r12, a_sigma, b_sigma, n, log_prior_ratio,
    #                             gprior, g)
    
    if (proposal_ratio != 1) {
      acceptMH <- acceptMH + log(proposal_ratio)
    }
    
    acceptMH <- exp(min(acceptMH, 700))
    u <- runif(1)

    # g. Accept or reject the new model
    if (u <= acceptMH) {
      model <- proposemodel 
      Xsel <- Xselnew
      k <- newk
      l <- newl
      betmodelac <- betmodelac + 1
      A_curr <- A_prop
      M_curr <- M_prop
      numDelete_curr <- numDelete_prop
      numAdd_curr <- numAdd_prop
      if (adaptive_moves) {
        metric <- metric_prop
      }
    }
    
    # h. Update pi
    for (i in 1:d) {
      parent <- get_parent(M_fix, i)
      if (length(parent) == 0) {
        alpha_post <- alpha0_vec[i] + model[i]
        beta_post <- beta0_vec[i] + 1 - model[i]
        pi[i] <- min(max(rbeta(1, alpha_post, beta_post), 0.001), 0.999)
      } else {
        if (model[parent] == 0) {
          alpha_post <- alpha0_vec[i] + 0
          beta_post <- beta0_vec[i] + 1
          pi[i] <- min(max(rbeta(1, alpha_post, beta_post), 0.001), 0.999)
        } else {
          alpha_post <- alpha0_vec[i] + model[i]
          beta_post <- beta0_vec[i] + 1 - model[i]
          pi[i] <- min(max(rbeta(1, alpha_post, beta_post), 0.001), 0.999)
        }
      }
    }
    
    # i. Update regression parameters (compositional and non-compositional)
    leaf_model <- current_model_leaf(A_curr)
    l_leaf <- which(leaf_model == 1)
    k_leaf <- length(l_leaf)
    
    if (k_leaf > 0) {
      Xsel_leaf <- X[, l_leaf, drop = FALSE]
      
      # Update beta (Compositional Covariate)
      XtX <- crossprod(Xsel_leaf)
      Xty <- crossprod(Xsel_leaf, y_res)
      
      if (gprior) {
        XtX_inv <- solve(XtX + diag(1e-8, k_leaf))
        beta_ols <- XtX_inv %*% Xty
        post_mean <- (g/(g+1)) * beta_ols
        post_var <- (g/(g+1)) * sigma * XtX_inv
        beta_temp <- mvtnorm::rmvnorm(1, drop(post_mean), post_var)
        
        # yPy <- drop(crossprod(y_res, Xsel_leaf %*% beta_ols))
        # S <- drop(crossprod(y_res)) - (g/(g+1)) * yPy
        # an <- a_sigma + n/2
        # bn <- b_sigma + 0.5 * S
        # sigma <- 1 / rgamma(1, an, bn)
      } else {
        post_var <- solve(XtX + diag(1 / c_beta, k_leaf))
        post_mean <- post_var %*% Xty
        beta_temp <- mvtnorm::rmvnorm(1, drop(post_mean), post_var * sigma)
        
        # RSS <- sum((y_res - Xsel_leaf %*% post_mean)^2)
        # penalty <- (1 / c_beta) * sum(post_mean^2)
        # an <- a_sigma + (n-1) / 2
        # bn <- b_sigma + 0.5 * (RSS + penalty)
        # sigma <- 1 / rgamma(1, an, bn)
      }
      
      beta <- rep(0, d)
      beta[l_leaf] <- beta_temp
      
      # Update theta (Non-Compositional Covariate) - only if Z is not NULL
      if (p_theta > 0) {
        y_res_theta <- y - X %*% beta
        ZtZ <- crossprod(Z)
        post_var_theta <- solve(ZtZ/sigma + diag(1/c_theta, p_theta))
        post_mean_theta <- post_var_theta %*% (t(Z) %*% y_res_theta / sigma)
        theta <- mvtnorm::rmvnorm(1, post_mean_theta, post_var_theta)
      }
      
      # log-likelihood
      if (p_theta > 0) {
        ll_out[t] <- loglikelihood(y, X[, l_leaf, drop = FALSE] %*% t(beta_temp) + Z %*% t(theta), sigma)
      } else {
        ll_out[t] <- loglikelihood(y, X[, l_leaf, drop = FALSE] %*% t(beta_temp), sigma)
      }
    } else {
      # Empty model for compositional covariates
      # Update theta (non-compositional) for gaussian models
      if (p_theta > 0) {
        y_res_theta <- y
        ZtZ <- crossprod(Z)
        post_var_theta <- solve(ZtZ/sigma + diag(1/c_theta, p_theta))
        post_mean_theta <- post_var_theta %*% (t(Z) %*% y_res_theta / sigma)
        theta <- mvtnorm::rmvnorm(1, post_mean_theta, post_var_theta)
        
        # an <- a_sigma + n / 2
        # bn <- b_sigma + 0.5 * sum((y - Z %*% theta)^2)
        # sigma <- 1 / rgamma(1, an, bn)
        # 
        # ll_out[t] <- loglikelihood(y, Z %*% t(theta), sigma)
      } 
      # else {
      #   an <- a_sigma + n / 2
      #   bn <- b_sigma + 0.5 * sum(y^2)
      #   sigma <- 1 / rgamma(1, an, bn)
      #   ll_out[t] <- loglikelihood(y, rep(0, n), sigma)
      # }
    }
    
    if (k_leaf > 0) {
      if (p_theta > 0) {
        RSS <- sum((y - X[, l_leaf, drop = FALSE] %*% t(beta_temp) - Z %*% t(theta))^2)
        penalty_beta <- ifelse(gprior, 0, (1/c_beta) * sum(beta_temp^2))
        penalty_theta <- (1/c_theta) * sum(theta^2)
      } else {
        RSS <- sum((y - X[, l_leaf, drop = FALSE] %*% t(beta_temp))^2)
        penalty_beta <- ifelse(gprior, 0, (1/c_beta) * sum(beta_temp^2))
        penalty_theta <- 0
      }
    } else {
      if (p_theta > 0) {
        RSS <- sum((y - Z %*% t(theta))^2)
        penalty_beta <- 0
        penalty_theta <- (1/c_theta) * sum(theta^2)
      } else {
        RSS <- sum(y^2)
        penalty_beta <- 0
        penalty_theta <- 0
      }
    }
    
    an <- a_sigma + (n + k_leaf + p_theta) / 2
    bn <- b_sigma + 0.5 * (RSS + penalty_beta + penalty_theta)
    sigma <- 1 / rgamma(1, an, bn)
    
    if (k_leaf > 0) {
      if (p_theta > 0) {
        ll_out[t] <- loglikelihood(y, X[, l_leaf, drop = FALSE] %*% t(beta_temp) + Z %*% t(theta), sigma)
      } else {
        ll_out[t] <- loglikelihood(y, X[, l_leaf, drop = FALSE] %*% t(beta_temp), sigma)
      }
    } else {
      if (p_theta > 0) {
        ll_out[t] <- loglikelihood(y, Z %*% t(theta), sigma)
      } else {
        ll_out[t] <- loglikelihood(y, rep(0, n), sigma)
      }
    }
    
    # j. SAVE OUTPUT
    if ((t > num_burn) & (t %% num_thin) == 0) {
      model_out[it, ] <- model
      model_out_leaf[it, ] <- leaf_model
      beta_out[it, ] <- beta
      if (p_theta > 0) {
        theta_out[it, ] <- theta
      }
      sigma_out[it] <- sigma
      pi_out[it, ] <- pi
      cMove[it] <- move
      pAccept[it] <- acceptMH
      Asamples[[it]] <- A_curr
      Msamples[[it]] <- M_curr
      if (adaptive_moves) {
        meanD[it] <- metric
      }
      
      # PREDICTION
      if ((!is.null(Xnew) && p_theta == 0) || (!is.null(Xnew) && !is.null(Znew) && p_theta > 0)) {
        if (k_leaf > 0) {
          if (p_theta > 0) {
            ypred_out[, it] <- rnorm(ntest, Xnew[, l_leaf, drop = FALSE] %*% t(beta_temp) + Znew %*% t(theta), sqrt(sigma))
          } else {
            ypred_out[, it] <- rnorm(ntest, Xnew[, l_leaf, drop = FALSE] %*% t(beta_temp), sqrt(sigma))
          }
        }
      }
      it <- it + 1
    }
  }
  
  # 4. FINAL RESULTS --------------------------------------------------------
  total_time <- proc.time() - start_time
  
  pb$terminate()
  
  return(list(
    model = model_out, 
    model_leaf = model_out_leaf, 
    beta = beta_out, 
    theta = theta_out,
    sigma = sigma_out, 
    pi = pi_out,
    ypred = ypred_out, 
    move = cMove, 
    accept = pAccept, 
    ll = ll_out,
    M = Msamples, 
    A = Asamples, 
    nacc = betmodelac, 
    time = total_time,
    acceptance_rate = betmodelac / num_iter,
    metric = meanD,
    prior_type = ifelse(gprior, "g-prior", "ridge-prior"),
    g_value = ifelse(gprior, g, NA),
    tau = tau,
    family = "gaussian"
  ))
}


run_multiple_chains_gaussian <- function(y, X, Z = NULL, A, M, priors = c(5.0, 10.0, 2.0, 2.0, 1.0, 1.0),
                                         N = 1, num_iter, num_burn, num_thin, 
                                         method = "OLS", gprior = FALSE, g = NULL,
                                         adaptive_levels = TRUE, tau = 0.5, kappa = 1,
                                         adaptive_moves = TRUE, lambda = 1,
                                         n_chains, parallel = TRUE, train_prop = 0.8, seed = NULL) {
  
  # Load required libraries
  if (!require(posterior)) install.packages("posterior")
  if (!require(doParallel)) install.packages("doParallel")
  library(posterior)
  library(doParallel)
  
  # Set up parallel backend
  if (parallel) {
    n_cores <- min(n_chains, detectCores() - 1)
    cl <- makeCluster(n_cores)
    registerDoParallel(cl)
    on.exit(stopCluster(cl))
  }
  
  # Run chains in parallel or sequentially
  chains <- foreach(i = 1:n_chains, .combine = 'c', .packages = c("mvtnorm", "MASS")) %dopar% {
    
    if(file.exists("R/utils.R")) {
      source("R/utils.R")
    } else {
      warning("File 'R/utils.R' don't found in the expected path.")
    }
    
    # Set unique seed for each chain
    chain_seed <- if (!is.null(seed)) seed + i else NULL
    
    # Split data into training and test sets for this chain
    n <- nrow(X)
    n_train <- floor(train_prop * n)
    train_indices <- sample(1:n, n_train)
    
    X_train <- X[train_indices, , drop = FALSE]
    y_train <- y[train_indices]
    X_test <- X[-train_indices, , drop = FALSE]
    y_test <- y[-train_indices]
    
    if (!is.null(Z)) {
      Z_train <- Z[train_indices, , drop = FALSE]
      Z_test <- Z[-train_indices, , drop = FALSE]
    } else {
      Z_train <- NULL
      Z_test <- NULL
    }
    
    # Run MCMC on training data
    result <- relativeShiftMCMC_gaussian(
      y = y_train, X = X_train, Z = Z_train, A = A, M = M, priors = priors,
      N = N, num_iter = num_iter, num_burn = num_burn, num_thin = num_thin,
      Xnew = X_test, Znew = Z_test, method = method,
      seed = chain_seed, gprior = gprior, g = g,
      adaptive_levels = adaptive_levels, tau = tau, kappa = kappa,
      adaptive_moves = adaptive_moves, lambda = lambda
    )
    
    # Add test data and predictions to results
    result$X_test <- X_test
    result$y_test <- y_test
    result$Z_test <- Z_test
    result$train_indices <- train_indices
    
    list(result)
  }
  
  # Extract parameter samples from each chain
  extract_samples <- function(chain, param) {
    switch(param,
           "beta" = chain$beta,
           "theta" = chain$theta,
           "sigma" = chain$sigma)
  }
  
  params <- c("beta", "theta", "sigma")
  
  samples_list <- lapply(params, function(p) {
    lapply(chains, function(chain) extract_samples(chain, p))
  })
  names(samples_list) <- params
  
  compute_diagnostics <- function(samples) {
    if (is.null(samples) || length(samples) == 0 || is.null(samples[[1]])) return(NULL)
    
    chain1 <- as.matrix(samples[[1]]) 
    n_iter <- nrow(chain1)
    n_params <- ncol(chain1)
    
    samples_array <- array(NA, dim = c(n_iter, n_chains, n_params))
    
    for (chain_idx in 1:n_chains) {
      mat_data <- as.matrix(samples[[chain_idx]])
      
      if(nrow(mat_data) == n_iter && ncol(mat_data) == n_params) {
        samples_array[, chain_idx, ] <- mat_data
      }
    }
    
    param_names <- colnames(samples[[1]])
    if (is.null(param_names)) {
      if(n_params == 1) param_names <- "param" 
      else param_names <- paste0("p", 1:n_params)
    }
    dimnames(samples_array)[[3]] <- param_names
    
    draws <- posterior::as_draws_array(samples_array)
    
    tryCatch({
      res <- posterior::summarise_draws(draws, 
                                        posterior::rhat, 
                                        posterior::ess_bulk,
                                        posterior::ess_tail)
      return(res)
    }, error = function(e) {
      return(list(error = paste("Diagnostics failed:", e$message)))
    })
  }
  
  diagnostics <- lapply(samples_list, compute_diagnostics)
  
  return(list(chains = chains, diagnostics = diagnostics))
}


rpolya_gamma <- function(n, b, c) {
  if (!requireNamespace("BayesLogit", quietly = TRUE)) {
    stop("Package 'BayesLogit' is required for Polya-Gamma sampling")
  }
  BayesLogit::rpg(n, b, c)
}

accept_binary <- function(Xold, Xnew, c, kappa_pg, r21, r12, 
                          log_prior_ratio, omega_pg) {
  
  log_marg_binary <- function(X, kappa_pg, omega_pg, c) {
    n <- length(kappa_pg)
    
    # Controllo stabilità omega_pg
    if (any(omega_pg <= 1e-12) || any(!is.finite(omega_pg)) || any(omega_pg > 1e12)) {
      return(-Inf)
    }
    
    # Termini comuni
    term1 <- -0.5 * sum(kappa_pg^2 / omega_pg)
    term5 <- -0.5 * sum(log(omega_pg))
    
    # MODELLO NULLO
    if (is.null(X) || ncol(X) == 0) {
      return(term1 + term5)
    }
    
    k <- ncol(X)
    
    # MODELLO CON COVARIATE
    tryCatch({
      # Usa scaling per stabilità
      X_scaled <- X * sqrt(omega_pg)
      XtWX <- crossprod(X_scaled)
      Sigma_inv <- XtWX + diag(1/c, k)
      
      # Stabilizzazione con controllo autovalori
      eigen_vals <- eigen(Sigma_inv, symmetric = TRUE, only.values = TRUE)$values
      min_eigen <- min(eigen_vals)
      
      if (min_eigen <= 1e-8) {
        jitter_amount <- max(abs(min_eigen) + 1e-6, 1e-8)
        diag(Sigma_inv) <- diag(Sigma_inv) + jitter_amount
      }
      
      # Cholesky decomposition
      Sigma_chol <- tryCatch({
        chol(Sigma_inv)
      }, error = function(e) {
        # Se fallisce, usa decomposizione spettrale
        eig <- eigen(Sigma_inv, symmetric = TRUE)
        eig$values[eig$values < 1e-10] <- 1e-10
        Sigma_inv_fixed <- eig$vectors %*% diag(eig$values) %*% t(eig$vectors)
        chol(Sigma_inv_fixed)
      })
      
      log_det_Sigma_inv <- 2 * sum(log(diag(Sigma_chol)))
      Sigma <- chol2inv(Sigma_chol)
      
      # Calcolo di mu (più stabile)
      XtK <- crossprod(X, kappa_pg)
      mu <- Sigma %*% XtK
      
      # CORREZIONE IMPORTANTE: term2 corretto
      term2 <- 0.5 * as.numeric(crossprod(mu, XtK))
      term3 <- 0.5 * log_det_Sigma_inv          
      term4 <- -0.5 * k * log(c)                
      
      log_marg <- term1 + term2 + term3 + term4 + term5
      
      if (!is.finite(log_marg)) {
        return(-Inf)
      }
      
      return(log_marg)
      
    }, error = function(e) {
      return(-Inf)
    })
  }
  
  log_post_old <- log_marg_binary(Xold, kappa_pg, omega_pg, c)
  log_post_new <- log_marg_binary(Xnew, kappa_pg, omega_pg, c)
  
  # Gestione -Inf
  if (is.infinite(log_post_old) && is.infinite(log_post_new)) {
    return(-Inf)
  }
  
  if (is.infinite(log_post_old)) {
    return(Inf)
  }
  
  if (is.infinite(log_post_new)) {
    return(-Inf)
  }
  
  log_accept <- (log_post_new - log_post_old) + (log(r21) - log(r12)) + log_prior_ratio
  
  # Clipping per stabilità
  return(pmax(pmin(log_accept, 700), -700))
}


relativeShiftMCMC_binary <- function(y, X, Z = NULL, A, M, priors = c(5.0, 10.0, 1.0, 1.0), 
                                       N = 1, num_iter, num_burn, num_thin, 
                                       Xnew = NULL, Znew = NULL,
                                       seed = NULL, 
                                       adaptive_levels = TRUE, tau = 0.5, kappa = 1,
                                       adaptive_moves = TRUE, lambda = 1) {
  
  # 1. PARAMETERS INITIALISATION -----------------------------------------------
  if (!is.null(seed)) set.seed(seed)
  
  if (any(!unique(y) %in% c(0, 1))) {
    stop("For binomial models, y must contain only 0 and 1 values")
  }
  
  n <- nrow(X)
  d <- ncol(X)
  
  # IF Z = NULL
  if (is.null(Z)) {
    p_theta <- 0
    Z <- matrix(0, nrow = n, ncol = 0)
  } else {
    p_theta <- ncol(Z)
  }
  
  ntest <- if (!is.null(Xnew)) nrow(Xnew) else 0
  
  c_beta <- priors[1]
  c_theta <- priors[2]
  alpha0 <- priors[3] 
  beta0 <- priors[4]
  
  # INITIALISATION FROM PRIOR
  if (adaptive_levels) {
    betaPrior_vec <- compute_BetaPrior_nodes(M, tau = tau, kappa = kappa)
    alpha0_vec <- betaPrior_vec[[1]]
    beta0_vec <- betaPrior_vec[[2]]
  } else {
    alpha0_vec <- rep(alpha0, d)
    beta0_vec <- rep(beta0, d)
  }

  pi <- sapply(1:d, function(i) rbeta(1, alpha0_vec[i], beta0_vec[i]))
  
  # INITIALISATION THETA NON-COMPOSITIONAL COVARIATES
  if (p_theta > 0) {
    theta <- matrix(rnorm(p_theta, 0, sqrt(c_theta)), ncol = 1)
  } else {
    theta <- numeric(0)
  }
  
  # INITIAL TREE STRUCTURE
  model <- current_model(A)
  A_fix <- A_curr <- A
  M_fix <- M_curr <- M
  
  l <- which(model == 1)
  k <- length(l)
  
  # INITIALISATION BETA COMPOSITIONAL COVARIATES
  Xsel <- X[, l, drop = FALSE]
  beta <- matrix(0, nrow = d, ncol = 1)  
  if (length(l) > 0) {
    beta[l, 1] <- rnorm(length(l), mean = 0, sd = sqrt(c_beta)) 
  }
  
  # INITIALISATION POLYA-GAMMA
  if (p_theta > 0) {
    eta_init <- X %*% beta + Z %*% theta
  } else {
    eta_init <- X %*% beta
  }
  
  if (!requireNamespace("BayesLogit", quietly = TRUE)) {
    stop("Package 'BayesLogit' required for Polya-Gamma sampling")
  }
  omega_pg <- BayesLogit::rpg(n, 1, as.vector(eta_init))
  kappa_pg <- y - 0.5
  
  # MOVE COUNTER
  resMove <- move_tree(A_curr, M_curr, M_fix, N, move = "delete")
  numAdd_curr <- resMove[[5]]
  numDelete_curr <- resMove[[6]]
  
  # 2. SETUP STORAGE
  num_samples <- as.integer(floor((num_iter - num_burn) / num_thin))
  model_out <- model_out_leaf <- beta_out <- matrix(0, nrow = num_samples, ncol = d)
  if (p_theta > 0) {
    theta_out <- matrix(0, nrow = num_samples, ncol = p_theta)  # num_samples x p_theta
  } else {
    theta_out <- NULL
  }
  pAccept <- cMove <- meanD <- rep(0, num_samples)
  pi_out <- matrix(0, nrow = num_samples, ncol = d)
  ll_out <- rep(0, num_iter)
  Asamples <- vector("list", num_samples)
  Msamples <- vector("list", num_samples)
  
  if ((!is.null(Xnew) && p_theta == 0) || (!is.null(Xnew) && !is.null(Znew) && p_theta > 0)) {
    ypred_out <- matrix(0, nrow = ntest, ncol = num_samples)
  } else {
    ypred_out <- NULL
  }
  
  it <- 1
  betmodelac <- 0
  start_time <- proc.time()
  
  # LOADING BAR
  pb <- progress::progress_bar$new(
    format = "  MCMC Progress [:bar] :percent | ETA: :eta | Elapsed: :elapsed | Acceptance: :acc_rate",
    total = num_iter, 
    clear = FALSE,
    width = 80,
    show_after = 0
  )
  
  # 3. MCMC LOOP ---------------------------------------------------------------
  for (t in 1:num_iter) {
    
    pb$tick(tokens = list(acc_rate = sprintf("%.3f", betmodelac/max(1, t-1))))
    
    # a. Update Polya-Gamma variables
    if (p_theta > 0) {
      eta <- as.vector(X %*% beta + Z %*% theta)
    } else {
      eta <- as.vector(X %*% beta)
    }
    omega_pg <- BayesLogit::rpg(n, 1, eta)
    kappa_pg <- y - 0.5
    
    # b. New move proposal for compositional covariates
    if (adaptive_moves) {
      pEnt <- adaptiveMoves(A_curr, M_curr, lambda)
      pAdd <- pEnt$add
      metric <- pEnt$metric
    } else {
      pAdd <- 0.5
      metric <- 0
    }
    
    if (numDelete_curr == 0) {
      move <- "add"
    } else if (numAdd_curr == 0) {
      move <- "delete"
    } else {
      move <- sample(c("add", "delete"), 1, prob = c(pAdd, 1 - pAdd))
    }
    
    resMove <- move_tree(A_curr, M_curr, M_fix, N, move = move)
    A_prop <- resMove[[1]]
    M_prop <- resMove[[2]]
    numAdd_prop <- resMove[[3]]
    numDelete_prop <- resMove[[4]]
    numAdd_curr <- resMove[[5]]
    numDelete_curr <- resMove[[6]]
    
    proposemodel <- current_model(A_prop)
    
    proposal_ratio <- 1
    if (adaptive_moves) {
      pEnt_prop <- adaptiveMoves(A_prop, M_prop, lambda)
      pAdd_prop <- pEnt_prop$add
      metric_prop <- pEnt_prop$metric
      if (pAdd > 0 && pAdd < 1 && pAdd_prop > 0 && pAdd_prop < 1) {
        if (move == "add") {
          proposal_ratio <- (1 - pAdd_prop) / pAdd
        } else {
          proposal_ratio <- pAdd_prop / (1 - pAdd)
        }
      }
    }
    
    # c. ratio proposal
    r12 <- 1 / ifelse(move == "add", numAdd_curr, numDelete_curr)
    r21 <- 1 / ifelse(move == "add", numDelete_prop, numAdd_prop)
    
    # d. Evaluation of the proposed model
    newl <- which(proposemodel == 1)
    newk <- length(newl)
    Xselnew <- X[, newl, drop = FALSE] 
    
    # e. ratio prior
    log_prior_old <- log_prior_tree(M_fix, model, pi)
    log_prior_new <- log_prior_tree(M_fix, proposemodel, pi)
    log_prior_ratio <- log_prior_new - log_prior_old
    
    # f. Acceptance log-probability
    acceptMH <- accept_binary(Xsel, Xselnew, c_beta, kappa_pg, 
                              r21, r12, log_prior_ratio, omega_pg)
    
    if (proposal_ratio != 1) {
      acceptMH <- acceptMH + log(proposal_ratio)
    }
    
    acceptMH <- exp(min(acceptMH, 700))
    u <- runif(1)
    
    # g. Accept or reject the new model
    if (u <= acceptMH) {
      model <- proposemodel 
      Xsel <- Xselnew
      k <- newk
      l <- newl
      betmodelac <- betmodelac + 1
      A_curr <- A_prop
      M_curr <- M_prop
      numDelete_curr <- numDelete_prop
      numAdd_curr <- numAdd_prop
      if (adaptive_moves) {
        metric <- metric_prop
      }
    }
    
    # h. Update pi
    for (i in 1:d) {
      parent <- get_parent(M_fix, i)
      if (length(parent) == 0) {
        alpha_post <- alpha0_vec[i] + model[i]
        beta_post <- beta0_vec[i] + 1 - model[i]
        pi[i] <- min(max(rbeta(1, alpha_post, beta_post), 0.001), 0.999)
      } else {
        if (model[parent] == 0) {
          alpha_post <- alpha0_vec[i] + 0
          beta_post <- beta0_vec[i] + 1
          pi[i] <- min(max(rbeta(1, alpha_post, beta_post), 0.001), 0.999)
        } else {
          alpha_post <- alpha0_vec[i] + model[i]
          beta_post <- beta0_vec[i] + 1 - model[i]
          pi[i] <- min(max(rbeta(1, alpha_post, beta_post), 0.001), 0.999)
        }
      }
    }
    
    # i. Update regression parameters (compositional and non-compositional)
    leaf_model <- current_model_leaf(A_curr)
    l_leaf <- which(leaf_model == 1)
    k_leaf <- length(l_leaf)
    
    if (k_leaf > 0) {
      Xsel_leaf <- X[, l_leaf, drop = FALSE]
      
      # Prepara residui
      if (p_theta > 0) {
        kappa_res_pg <- kappa_pg - as.vector(Z %*% theta)
      } else {
        kappa_res_pg <- kappa_pg
      }
      
      # Aggiornamento beta con gestione robusta
      tryCatch({
        X_scaled <- Xsel_leaf * sqrt(omega_pg)
        XtWX <- crossprod(X_scaled)
        Sigma_inv <- XtWX + diag(1 / c_beta, k_leaf)
        
        # Stabilizzazione
        diag(Sigma_inv) <- diag(Sigma_inv) + 1e-8
        
        # Decomposizione spettrale
        eigen_decomp <- eigen(Sigma_inv, symmetric = TRUE)
        eigen_values <- pmax(eigen_decomp$values, 1e-10)  # Forza positività
        eigen_vectors <- eigen_decomp$vectors
        
        # Calcola Sigma in modo stabile
        Sigma <- eigen_vectors %*% diag(1/eigen_values) %*% t(eigen_vectors)
        
        # Media posteriore
        XtK <- crossprod(Xsel_leaf, kappa_res_pg)
        mu <- Sigma %*% XtK
        
        # Campionamento: usa L = V * D^(-1/2) dove Sigma = V * D^(-1) * V'
        z <- rnorm(k_leaf)
        L <- eigen_vectors %*% diag(1/sqrt(eigen_values))
        beta_temp <- as.vector(mu) + as.vector(L %*% z) 
        
        # Assegna a beta completo
        beta <- matrix(0, nrow = d, ncol = 1) 
        beta[l_leaf, 1] <- beta_temp
        
      }, error = function(e) {
        # Fallback: usa il modello precedente
        warning(sprintf("Errore aggiornamento beta (iter %d): %s", t, e$message))
      })
      
      # Aggiornamento theta
      if (p_theta > 0) {
        kappa_theta_pg <- kappa_pg - as.vector(X %*% beta)
        
        tryCatch({
          ZtW <- t(Z) * omega_pg
          Sigma_inv_theta <- ZtW %*% Z + diag(1/c_theta, p_theta)
          Sigma_inv_theta <- Sigma_inv_theta + diag(1e-8, p_theta)
          Sigma_theta_chol <- chol(Sigma_inv_theta)
          Sigma_theta <- chol2inv(Sigma_theta_chol)
          mu_theta <- Sigma_theta %*% (t(Z) %*% kappa_theta_pg)
          theta <- matrix(mvtnorm::rmvnorm(1, drop(mu_theta), Sigma_theta), ncol = 1)
        }, error = function(e) {
          warning(sprintf("Errore aggiornamento theta (iter %d): %s", t, e$message))
        })
      }
      
      # Log-likelihood con clipping
      if (p_theta > 0) {
        eta <- as.vector(Xsel_leaf %*% beta_temp + Z %*% theta)
      } else {
        eta <- as.vector(Xsel_leaf %*% beta_temp)
      }
      eta <- pmin(pmax(eta, -20), 20)
      ll_out[t] <- sum(y * eta - log1p(exp(eta)))
      
    } else {
      # Modello vuoto per covariate composizionali
      beta <- matrix(0, nrow = d, ncol = 1)
      
      if (p_theta > 0) {
        kappa_theta_pg <- kappa_pg
        tryCatch({
          ZtW <- t(Z) * omega_pg
          Sigma_inv_theta <- ZtW %*% Z + diag(1/c_theta, p_theta)
          Sigma_inv_theta <- Sigma_inv_theta + diag(1e-8, p_theta)
          Sigma_theta_chol <- chol(Sigma_inv_theta)
          Sigma_theta <- chol2inv(Sigma_theta_chol)
          mu_theta <- Sigma_theta %*% (t(Z) %*% kappa_theta_pg)
          theta <- matrix(mvtnorm::rmvnorm(1, drop(mu_theta), Sigma_theta), ncol = 1)
          
          eta_theta <- as.vector(Z %*% theta)
          eta_theta <- pmin(pmax(eta_theta, -20), 20)
          ll_out[t] <- sum(y * eta_theta - log1p(exp(eta_theta)))
        }, error = function(e) {
          ll_out[t] <- sum(y * 0 - log1p(exp(0)))
        })
      } else {
        ll_out[t] <- sum(y * 0 - log1p(exp(0)))  # modello nullo
      }
    }
    
    # j. SAVE OUTPUT
    if ((t > num_burn) & (t %% num_thin) == 0) {
      model_out[it, ] <- model
      model_out_leaf[it, ] <- leaf_model
      beta_out[it, ] <- as.vector(beta)
      if (p_theta > 0) {
        theta_out[it, ] <- as.vector(theta)
      }
      pi_out[it, ] <- pi
      cMove[it] <- move
      pAccept[it] <- acceptMH
      Asamples[[it]] <- A_curr
      Msamples[[it]] <- M_curr
      if (adaptive_moves) {
        meanD[it] <- metric
      }
      
      # PREDICTION
      if ((!is.null(Xnew) && p_theta == 0) || (!is.null(Xnew) && !is.null(Znew) && p_theta > 0)) {
        if (k_leaf > 0) {
          if (p_theta > 0) {
            eta_pred <- as.vector(Xnew[, l_leaf, drop = FALSE] %*% beta_temp + Znew %*% theta)
          } else {
            eta_pred <- as.vector(Xnew[, l_leaf, drop = FALSE] %*% beta_temp)
          }
          eta_pred <- pmin(pmax(eta_pred, -10), 10)
          prob_pred <- 1 / (1 + exp(-eta_pred))
          ypred_out[, it] <- rbinom(ntest, 1, prob_pred)
        } else {
          if (p_theta > 0) {
            eta_pred <- as.vector(Znew %*% theta)
            prob_pred <- 1 / (1 + exp(-eta_pred))
            ypred_out[, it] <- rbinom(ntest, 1, prob_pred)
          } else {
            ypred_out[, it] <- rbinom(ntest, 1, 0.5)
          }
        }
      }
      it <- it + 1
    }
  }
  
  pb$terminate()
  
  # 4. FINAL RESULTS --------------------------------------------------------
  total_time <- proc.time() - start_time
  
  return(list(
    model = model_out, 
    model_leaf = model_out_leaf, 
    beta = beta_out, 
    theta = theta_out,
    pi = pi_out,
    ypred = ypred_out, 
    move = cMove, 
    accept = pAccept, 
    ll = ll_out,
    M = Msamples, 
    A = Asamples, 
    nacc = betmodelac, 
    time = total_time,
    acceptance_rate = betmodelac / num_iter,
    metric = meanD,
    tau = tau,
    family = "binomial"
  ))
}


# R-hat (Gelman-Rubin diagnostic)
compute_rhat <- function(samples) {
  if (!require(posterior)) install.packages("posterior")
  library(posterior)
  
  # Verifica che l'input sia una lista di catene
  if (!is.list(samples)) {
    stop("L'input deve essere una lista di catene MCMC!")
  }
  
  n_chains <- length(samples)
  n_iter <- length(samples[[1]])
  
  all_zero <- all(sapply(samples, function(chain) all(chain == 0)))
  if (all_zero) {
    return(1.0) 
  }
  
  samples_array <- array(NA, dim = c(n_iter, n_chains, 1))
  for (chain in 1:n_chains) {
    samples_array[, chain, 1] <- samples[[chain]]
  }
  
  draws <- as_draws_array(samples_array)
  rhat_values <- rhat(draws)

  if (is.nan(rhat_values[1]) || is.infinite(rhat_values[1])) {
    return(1.0)
  }
  
  return(rhat_values[1])
}

check_rhat_convergence <- function(rhat_vector, threshold = 1.1, target_proportion = 0.95, plot_hist = TRUE) {
  converged_proportion <- mean(rhat_vector < threshold, na.rm = TRUE)
  
  has_converged <- converged_proportion >= target_proportion
  
  cat("Proportion of R-hat values below threshold:", round(converged_proportion, 4), "\n")
  if (has_converged) {
    cat("✅ CONVERGENCE ACHIEVED: At least", target_proportion * 100, "% of R-hat values are below", threshold, "\n")
  } else {
    cat("❌ CONVERGENCE NOT ACHIEVED: Less than", target_proportion * 100, "% of R-hat values are below", threshold, "\n")
  }
  
  if (plot_hist) {
    hist(rhat_vector, 
         main = "Distribution of R-hat Values", 
         xlab = "R-hat", 
         col = "lightblue",
         breaks = 20)
    abline(v = threshold, col = "red", lwd = 2)
    legend("topright", 
           legend = paste("Threshold =", threshold), 
           col = "red", 
           lwd = 2, 
           bty = "n")
  }
  
  invisible(list(
    converged_proportion = converged_proportion,
    has_converged = has_converged,
    threshold = threshold,
    target_proportion = target_proportion
  ))
}


compute_ess <- function(samples, method = "bulk") {
  if (!require(posterior)) install.packages("posterior")
  library(posterior)
  
  if (!is.list(samples)) {
    stop("Input must be a list of MCMC chains!")
  }
  
  if (!all(sapply(samples, is.vector))) {
    stop("All chains must be vectors (scalar parameters)!")
  }
  
  all_zero <- all(sapply(samples, function(chain) all(chain == 0)))
  if (all_zero) {
    n_total <- sum(sapply(samples, length))
    warning("Parameter is constantly zero. Returning total sample size as ESS.")
    return(n_total)
  }
  
  n_chains <- length(samples)
  n_iter <- length(samples[[1]])
  
  samples_array <- array(NA, dim = c(n_iter, n_chains, 1))
  for (chain in 1:n_chains) {
    samples_array[, chain, 1] <- samples[[chain]]
  }
  
  draws <- as_draws_array(samples_array)
  
  tryCatch({
    if (method == "bulk") {
      ess_value <- ess_bulk(draws)
    } else if (method == "tail") {
      ess_value <- ess_tail(draws)
    } else if (method == "mean") {
      ess_value <- ess_mean(draws)
    } else if (method == "sd") {
      ess_value <- ess_sd(draws)
    } else if (method == "median") {
      ess_value <- ess_median(draws)
    } else {
      stop("Method must be one of: 'bulk', 'tail', 'mean', 'sd', 'median'")
    }
    
    if (is.nan(ess_value[1]) || is.infinite(ess_value[1])) {
      warning("ESS calculation resulted in NaN or Inf. Returning total sample size.")
      return(sum(sapply(samples, length)))
    }
    
    return(ess_value[1])
  }, error = function(e) {
    warning("Error in ESS calculation: ", e$message, ". Returning total sample size.")
    return(sum(sapply(samples, length)))
  })
}


run_multiple_chains_binary <- function(y, X, Z = NULL, A, M, priors = c(5.0, 10.0, 1.0, 1.0), 
                                       N = 1, num_iter, num_burn, num_thin, 
                                       adaptive_levels = TRUE, tau = 0.5, kappa = 1,
                                       adaptive_moves = TRUE, lambda = 1,
                                       n_chains, parallel = TRUE, train_prop = 0.8, seed = NULL) {
  
  # Load required libraries
  if (!require(posterior)) install.packages("posterior")
  if (!require(doParallel)) install.packages("doParallel")
  library(posterior)
  library(doParallel)
  
  # Set up parallel backend
  if (parallel) {
    n_cores <- min(n_chains, detectCores() - 1)
    cl <- makeCluster(n_cores)
    registerDoParallel(cl)
    on.exit(stopCluster(cl))
  }
  
  # Run chains in parallel or sequentially
  chains <- foreach(i = 1:n_chains, .combine = 'c', .packages = c("mvtnorm", "MASS")) %dopar% {
    
    if(file.exists("R/utils.R")) {
      source("R/utils.R")
    } else {
      warning("File 'R/utils.R' don't found in the expected path.")
    }
    
    # Set unique seed for each chain
    chain_seed <- if (!is.null(seed)) seed + i else NULL
    
    # Split data into training and test sets for this chain
    n <- nrow(X)
    n_train <- floor(train_prop * n)
    
    # Random split
    train_indices <- sample(1:n, n_train)
    
    X_train <- X[train_indices, , drop = FALSE]
    y_train <- y[train_indices]
    X_test <- X[-train_indices, , drop = FALSE]
    y_test <- y[-train_indices]
    
    if (!is.null(Z)) {
      Z_train <- Z[train_indices, , drop = FALSE]
      Z_test <- Z[-train_indices, , drop = FALSE]
    } else {
      Z_train <- NULL
      Z_test <- NULL
    }
    
    # Run MCMC on training data
    result <- relativeShiftMCMC_binary(
      y = y_train, X = X_train, Z = Z_train, A = A, M = M, priors = priors,
      N = N, num_iter = num_iter, num_burn = num_burn, num_thin = num_thin,
      Xnew = X_test, Znew = Z_test,
      seed = chain_seed,
      adaptive_levels = adaptive_levels, tau = tau, kappa = kappa,
      adaptive_moves = adaptive_moves, lambda = lambda
    )
    
    result$X_test <- X_test
    result$y_test <- y_test
    result$Z_test <- Z_test
    result$train_indices <- train_indices
    
    list(result)
  }
  
  extract_samples <- function(chain, param) {
    switch(param,
           "beta" = chain$beta,
           "theta" = chain$theta)
  }
  
  params <- c("beta", "theta")
  
  samples_list <- lapply(params, function(p) {
    lapply(chains, function(chain) extract_samples(chain, p))
  })
  names(samples_list) <- params
  
  compute_diagnostics <- function(samples) {
    if (is.null(samples) || length(samples) == 0 || is.null(samples[[1]])) return(NULL)
    
    chain1 <- as.matrix(samples[[1]]) 
    n_iter <- nrow(chain1)
    n_params <- ncol(chain1)
    
    samples_array <- array(NA, dim = c(n_iter, n_chains, n_params))
    
    for (chain_idx in 1:n_chains) {
      mat_data <- as.matrix(samples[[chain_idx]])
      
      if(nrow(mat_data) == n_iter && ncol(mat_data) == n_params) {
        samples_array[, chain_idx, ] <- mat_data
      }
    }
    
    param_names <- colnames(samples[[1]])
    if (is.null(param_names)) {
      param_names <- paste0("p", 1:n_params)
    }
    dimnames(samples_array)[[3]] <- param_names
    
    draws <- posterior::as_draws_array(samples_array)
    
    tryCatch({
      res <- posterior::summarise_draws(draws, 
                                        posterior::rhat, 
                                        posterior::ess_bulk,
                                        posterior::ess_tail)
      return(res)
    }, error = function(e) {
      return(list(error = paste("Diagnostics failed:", e$message)))
    })
  }
  
  diagnostics <- lapply(samples_list, compute_diagnostics)
  
  return(list(chains = chains, diagnostics = diagnostics))
}
