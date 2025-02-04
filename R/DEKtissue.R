
#' @title Cross-cell type differential analysis
#'
#' @param K the number of cell types
#' @param Y a raw data matrix of complex samples
#' @param Prop cell-type proportions for all subjects
#' @param design a design matrix representing the status of samples (or continuous features: need to be developed)
#              example 1 of design: (0,...,0,1,...,1), 0 for control and 1 for cases
#              example 2 of design: (0,...,0) for all subjects belong to one group
#' @param WhichPar a number chosen from 1 to 2K (which parameter will be tested?)
#' @param contrast.vec contrast.vec should be a 2K*1 vector to specify a contrast vector
#' @param sort whether the p value is sorted
#' @param var.threshold var.threshold controls how much you want to bound your variation estimation. Default value is 0.1
#' @param data.threshold whether to bound varBeta
#'
#' @return A dataframe containing p-value, q-value and statistics
#' @import matrixStats
#' @export
#'

DEKTissue <- function(K, Y, Prop, design, WhichPar=NULL, contrast.vec=NULL, sort=F, var.threshold=0.1, data.threshold=TRUE){

  N = dim(Y)[2]
  if(dim(Prop)[1]!=N | dim(Prop)[2]!=K){
    stop("Dimension of proportion input is not correct!")
  }

  Y = t(na.omit(Y))
  G = dim(Y)[2]

  if(!all(design==0)){
    W = cbind(Prop,Prop*design)
  }else{
    W = Prop
  }
  H = solve(t(W)%*%W)%*%t(W)
  coefs = H%*%Y
  Ypred = W%*%coefs
  resi = Y-Ypred

  s2.case = colSums(resi^2) / (N - ncol(W))
  varBeta = matrix(diag(solve(t(W)%*%W)),ncol(W),1)%*%s2.case

  ## bound varBeta a bit
  if(data.threshold) {
    var.threshold = quantile(c(varBeta),0.1)
  }else {
    varBeta[varBeta<var.threshold] = var.threshold
  }

  res.table = data.frame(t.stat=rep(0,G), t.pval=rep(0,G), t.fdr=rep(0,G))
  rownames(res.table) = colnames(Y)

  if(!all(design==0)){
    if(is.null(contrast.vec)){
      res.table = data.frame(beta=rep(0,G), mu=rep(0,G), effect_size=rep(0,G), t.stat=rep(0,G), t.pval=rep(0,G), t.fdr=rep(0,G))
      rownames(res.table) = colnames(Y)
      res.table$beta=coefs[WhichPar,]
      res.table$mu=coefs[WhichPar-K,]
      res.table$effect_size = res.table$beta/(res.table$mu + res.table$beta/2)
      res.table$t.stat=coefs[WhichPar,]/sqrt(varBeta[WhichPar,])
      res.table$t.pval=2*pt(-abs(res.table$t.stat),df=N-ncol(W))
      res.table$t.fdr=p.adjust(res.table$t.pval, method = 'fdr')
      res.table = data.frame(res.table,t(coefs))
    }else{
      res.table = data.frame(t.stat=rep(0,G), t.pval=rep(0,G), t.fdr=rep(0,G))
      rownames(res.table) = colnames(Y)
      res.table$t.stat=as.numeric(contrast.vec%*%coefs)/as.numeric(sqrt(abs(contrast.vec)%*%varBeta))
      res.table$t.pval=2*pt(-abs(res.table$t.stat),df=N-ncol(W))
      res.table$t.fdr=p.adjust(res.table$t.pval, method = 'fdr')
      res.table = data.frame(res.table,t(coefs))
    }
  }else{

    if(is.null(contrast.vec)){
      res.table = data.frame(t.stat=rep(0,G), t.pval=rep(0,G), t.fdr=rep(0,G))
      rownames(res.table) = colnames(Y)
      res.table$t.stat=coefs[WhichPar,]/sqrt(varBeta[WhichPar,])
      res.table$t.pval=2*pt(-abs(res.table$t.stat),df=N-ncol(W))
      res.table$t.fdr=p.adjust(res.table$t.pval, method = 'fdr')
      res.table = data.frame(res.table,t(coefs))
    }else{
      res.table = data.frame(t.stat=rep(0,G), t.pval=rep(0,G), t.fdr=rep(0,G))
      rownames(res.table) = colnames(Y)
      i=which(contrast.vec!=0)[1]
      j=which(contrast.vec!=0)[2]
      coefs[coefs<0] = 0
      # res.table$muA = coefs[i,]
      # res.table$muB = coefs[j,]
      res.table$t.stat=as.numeric(contrast.vec%*%coefs)/as.numeric(sqrt(abs(contrast.vec)%*%varBeta))
      res.table$t.pval=2*pt(-abs(res.table$t.stat),df=N-ncol(W))
      res.table$t.fdr=p.adjust(res.table$t.pval, method = 'fdr')
      res.table = data.frame(res.table,t(coefs))
    }
  }

  if(sort){
    return(res.table[order(res.table$t.pval),])
  }else{
    return(res.table)
  }
}
