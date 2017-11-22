
neglogLikelihood <-
  function(theta) {
    c1 <- theta[1]
    a1 <- theta[2]
    sigma1 <- theta[3]
    sigma2 <- theta[4]
    R<- c1*exp(-D/a1)
    diag(R) <- 1
    diag(S) <- (dat$x2<25)*sigma1+(dat$x2>25)*sigma2
    V <- t(S) %*% R %*% S
    V_inv <- chol2inv(chol(V))
    XV <- crossprod(X, V_inv)
    XVX <- XV %*% X
    XVX_inv <- chol2inv(chol(XVX))
    I <- diag(nrow(D))
    logDetV <- determinant(x = V, logarithm = TRUE)$modulus
    logDetXVX <- determinant(x = XVX, logarithm = TRUE)$modulus
    P= I - X %*% chol2inv(chol(crossprod(X,X)))%*%t(X)
    y.t <-t(P)%*%y
    Q <- X %*% XVX_inv %*% XV
    logLikelihood <- -0.5 * (logDetV + logDetXVX + crossprod(y.t, V_inv) %*%  (I - Q) %*% y)
    neglogLikelihood <- -1 * logLikelihood
    return(neglogLikelihood)
  }