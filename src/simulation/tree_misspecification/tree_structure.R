#!/usr/bin/env Rscript

library(data.tree)
leaf2node_special = function(leafval, tree_structure){
  # pre-order
  K = ncol(leafval) # taxa are cols
  if (tree_structure == 'left'){
    aggval = t(apply(leafval, 1, function(x){rev(cumsum(x))}))
    y = aggval[,1:(K-1)]
    yl = aggval[,2:K]
    return(list(nodeval = y, leftval = yl))
  }
  if (tree_structure == 'right'){
    y = t(apply(leafval, 1, function(x){rev(cumsum(rev(x)))}))[,1:(K-1)]
    yl = leafval[,1:(K-1)]
    return(list(nodeval = y, leftval = yl))
  }
  if (tree_structure == 'balanced'){
    height = log2(K) # add asserting height is integer here
    df = data.frame(1:K)
    curr_label = K
    curr_level_num_nodes = K/2
    for (i in 1:height){
      df = cbind(df,rep((curr_label+1):(curr_label+curr_level_num_nodes),each=2^i))
      curr_label = curr_label+curr_level_num_nodes
      curr_level_num_nodes = curr_level_num_nodes/2
    }
    colnames(df) = paste0('level',height:0)
    df$pathString = apply(df,1,function(x){paste(rev(x),collapse = '/')})
    tree = as.Node(df)
    f = function(leafval_i){
      tree$Set(val=leafval_i,filterFun = isLeaf, traversal = 'pre-order')
      # y = as.numeric(tree$Get(Aggregate, "val", sum, traversal = 'pre-order', filterFun = isNotLeaf))
      tree$Do(function(x) {
        leftchild=names(x$children[1])
        rightchild=names(x$children[2])
        xl=x[[leftchild]]
        xr=x[[rightchild]]
        x$val <- xl$val+xr$val
        x$yl <- xl$val
      }, traversal = "post-order",filterFun = data.tree::isNotLeaf)
      y = as.numeric(tree$Get('val',traversal = 'pre-order',filterFun = isNotLeaf))
      yl = as.numeric(tree$Get('yl',traversal = 'pre-order',filterFun = isNotLeaf))
      return(list(nodeval=y,leftval=yl))
    }
    vallist=apply(leafval, 1, f)
    nodeval = do.call(rbind,lapply(vallist, function(x)x$nodeval))
    leftval = do.call(rbind,lapply(vallist, function(x)x$leftval))
    return(list(nodeval = nodeval, leftval = leftval))
    }
}

node2leaf_special = function(nodevals, tree_structure){
  # pre-order
  nodeval = nodevals$nodeval
  leftval = nodevals$leftval
  K = ncol(nodeval) + 1 # taxa are cols
  if (tree_structure == 'left'){
    return(cbind(leftval[,K-1],t(apply(nodeval-leftval,1,rev))))
  }
  if (tree_structure == 'right'){
    return(cbind(leftval, nodeval[,K-1]-leftval[,K-1]))
  }
  if (tree_structure == 'balanced'){
    height = log2(K) # add asserting height is integer here
    df = data.frame(1:K)
    curr_label = K
    curr_level_num_nodes = K/2
    for (i in 1:height){
      df = cbind(df,rep((curr_label+1):(curr_label+curr_level_num_nodes),each=2^i))
      curr_label = curr_label+curr_level_num_nodes
      curr_level_num_nodes = curr_level_num_nodes/2
    }
    colnames(df) = paste0('level',height:0)
    df$pathString = apply(df,1,function(x){paste(rev(x),collapse = '/')})
    tree = as.Node(df)
    f = function(nodeval_i,leftval_i){
      tree$Set(val=nodeval_i,filterFun = isNotLeaf, traversal = 'pre-order')
      tree$Set(yl=leftval_i,filterFun = isNotLeaf, traversal = 'pre-order')
      tree$Do(function(x) {
        leftchild=names(x$children[1])
        rightchild=names(x$children[2])
        xl=x[[leftchild]]
        xr=x[[rightchild]]
        xl$val <- x$yl 
        x$yr <- x$val - xl$val
        xr$val <- x$yr
      }, traversal = "pre-order",filterFun = data.tree::isNotLeaf)
      leafval = as.numeric(tree$Get('val',traversal = 'pre-order',filterFun = isLeaf))
      return(leafval)
    }
    leafval = do.call(rbind,lapply(1:nrow(nodeval), function(i)f(nodeval[i,],leftval[i,])))
    return(leafval)
  }
}
  
# test
# leafval = matrix(1:10,2,5)
# identical(node2leaf_special(leaf2node_special(leafval,'left'),'left'),leafval)
# identical(node2leaf_special(leaf2node_special(leafval,'right'),'right'),leafval)

theta2nodevals_special = function(theta, tree_structure){
  # pre-order
  K = ncol(theta) + 1
  if (tree_structure == 'left'){
    prods = t(apply(cbind(1,theta), 1, cumprod))
    nodevals = list(nodeval = prods[,1:(K-1)], leftval = prods[,2:K])
  }
  if (tree_structure == 'right'){
    prods = t(apply(cbind(1,1-theta), 1, cumprod))
    nodevals = list(nodeval = prods[,1:(K-1)], leftval = prods[,1:(K-1)] - prods[,2:K])
  }
  if (tree_structure == 'balanced'){
    height = log2(K) # add asserting height is integer here
    df = data.frame(1:K)
    curr_label = K
    curr_level_num_nodes = K/2
    for (i in 1:height){
      df = cbind(df,rep((curr_label+1):(curr_label+curr_level_num_nodes),each=2^i))
      curr_label = curr_label+curr_level_num_nodes
      curr_level_num_nodes = curr_level_num_nodes/2
    }
    colnames(df) = paste0('level',height:0)
    df$pathString = apply(df,1,function(x){paste(rev(x),collapse = '/')})
    tree = as.Node(df)
    f = function(theta_i){
      tree$Set(theta=theta_i,filterFun = isNotLeaf, traversal = 'pre-order')
      tree$prob = 1.0
      tree$Do(function(x) {
        leftchild=names(x$children[1])
        rightchild=names(x$children[2])
        xl=x[[leftchild]]
        xr=x[[rightchild]]
        xl$prob <- x$theta * x$prob
        xr$prob <- (1 - x$theta) * x$prob
        x$lprob <- xl$prob
      }, traversal = "pre-order",filterFun = data.tree::isNotLeaf)
      nodeval = as.numeric(tree$Get('prob',traversal = 'pre-order',filterFun = isNotLeaf))
      leftval = as.numeric(tree$Get('lprob',traversal = 'pre-order',filterFun = isNotLeaf))
      return(list(nodeval = nodeval,leftval=leftval))
    }
    vallist = apply(theta, 1, f)
    nodeval = do.call(rbind,lapply(vallist,function(x){x$nodeval}))
    leftval = do.call(rbind,lapply(vallist,function(x){x$leftval}))
    nodevals = list(nodeval = nodeval, leftval = leftval)
  }
  return(nodevals)
}

# test
# theta = matrix(rep(0.1,6),2,3)
# theta2nodevals_special(theta,'left')
# theta2nodevals_special(theta,'right')

right2left = function(mu,sigma,nmc){
  psi = rmvnorm(n = nmc,mean = mu,sigma = sigma)
  # tlr_inv
  nodevals = theta2nodevals_special(plogis(psi),'right')
  prob = node2leaf_special(nodevals, 'right')
  # tlr
  nodeprob = leaf2node_special(prob, 'left')
  psi = qlogis(nodeprob$leftval/nodeprob$nodeval)
  mu1 = apply(psi,2,mean)
  sigma1 = cov(psi)
  return(list(mu1,sigma1))
}

balanced2left = function(mu,sigma,nmc){
  psi = rmvnorm(n = nmc,mean = mu,sigma = sigma)
  # tlr_inv
  nodevals = theta2nodevals_special(plogis(psi),'balanced')
  prob = node2leaf_special(nodevals, 'balanced')
  # tlr
  nodeprob = leaf2node_special(prob, 'left')
  psi = qlogis(nodeprob$leftval/nodeprob$nodeval)
  mu1 = apply(psi,2,mean)
  sigma1 = cov(psi)
  return(list(mu1,sigma1))
}

left2right = function(mu,sigma,nmc){
  psi = rmvnorm(n = nmc,mean = mu,sigma = sigma)
  # tlr_inv
  nodevals = theta2nodevals_special(plogis(psi),'left')
  prob = node2leaf_special(nodevals, 'left')
  # tlr
  nodeprob = leaf2node_special(prob, 'right')
  psi = qlogis(nodeprob$leftval/nodeprob$nodeval)
  mu1 = apply(psi,2,mean)
  sigma1 = cov(psi)
  return(list(mu1,sigma1))
}

left2balanced = function(mu,sigma,nmc){
  psi = rmvnorm(n = nmc,mean = mu,sigma = sigma)
  # tlr_inv
  nodevals = theta2nodevals_special(plogis(psi),'left')
  prob = node2leaf_special(nodevals, 'left')
  # tlr
  nodeprob = leaf2node_special(prob, 'balanced')
  psi = qlogis(nodeprob$leftval/nodeprob$nodeval)
  mu1 = apply(psi,2,mean)
  sigma1 = cov(psi)
  return(list(mu1,sigma1))
}

