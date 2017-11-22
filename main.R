# Simple example showing REML estimation of non-stationary variance parameters,
#and kriging with an external drift with non-stationary variance paarmeters

library(DEoptim)
library(variography)
library(gstat)
library(sp)

#simulate field

#define discretisation grid
x1<-seq(1:50)-0.5
x2<-x1
grid<-expand.grid(x1,x2)
names(grid)<-c("x1","x2")

#compute spatial trend; x1 is used as covariate z
grid$z <- grid$x1
b1 <- 2
grid$mu<-b1*grid$z

#define variogram for simulation of residuals
sill<-1
range<-10
vgm<-sill*Exp(range)

#compute matrix with distances between simulation nodes
distx<-outer(grid$x1,grid$x1,FUN="-")
disty<-outer(grid$x2,grid$x2,FUN="-")
dist<-sqrt(distx^2+disty^2)

#compute matrix with mean covariances
cvm<-as(vgm,"CovariogramStructure") #coerce variogram to covariance function
f<-as(cvm, "function")
C<-f(h=dist)

#simulate values for residuals by Cholesky decomposition
set.seed(31415)
Upper<-chol(C)
G<-rnorm(n=nrow(grid),0,1) #simulate random numbers from standard normal distribution
grid$residuals <- crossprod(Upper,G)

#multiply residuals in upper half by a constant so that residual variance becomes non-stationary
ids <- which(grid$x2>25)
grid$residuals[ids] <- grid$residuals[ids]*3

#add residuals to trend
grid$y <- grid$mu+grid$residuals

#select simple random sample from grid
idssam <- sample.int(nrow(grid),size=100,replace=F)
dat <- grid[idssam,]

#compute matrix with distances between sampling points
D <- as.matrix(dist(cbind(dat$x1,dat$x2)))

olsmodel <- lm(y~z,data=dat)
X <- model.matrix(olsmodel)
y <- dat$y
coordinates(dat) <- ~x1+x2
vg <- variogram(y~z,data=dat)
plot(vg)

dat <- as.data.frame(dat)

S <- diag(nrow=nrow(dat))

# load the negative log-likelihood
source('likelihood.R')

#NB c1 is partial sill parameter of correlogram, so must be < 1
lbound <- c(c1 = 0, a1 = 0.1, sigma1=0.1, sigma2=0.1)
ubound <- c(c1 = 1, a1 = 20, sigma1=10, sigma2=10)
optPars <- DEoptim(
  fn = neglogLikelihood,
  lower = lbound,
  upper = ubound,
  control = DEoptim.control(strategy =2, bs=F, NP=40, itermax=200, CR=0.5, F=0.8, trace=T)
)

result<-optPars$optim$bestmem
c1 <- result[1]
a1 <- result[2]
sigma1 <- result[3]
sigma2 <- result[4]

#Now estimate regression coefficients by GLS using REML estimates of variogram parameters
R<- c1*exp(-D/a1)
diag(R) <- 1
diag(S) <- (dat$x2<25)*sigma1+(dat$x2>25)*sigma2
V <- t(S) %*% R %*% S
V_inv <- chol2inv(chol(V))
XV <- crossprod(X, V_inv)
XVX <- XV %*% X
XVX_inv <- chol2inv(chol(XVX))
Vy <- crossprod(y, V_inv)
XVy<-crossprod(X,t(Vy))
(beta_GLS<-XVX_inv%*%XVy)

#compare with OLS estimates of regression coefficients
#(beta_OLS <- coef(olsmodel))

#Predict at nodes of grid by kriging with an external drift with non-stationary variances

#Compute matrix with correlation between sampling points and prediction node
library(fields)
D0 <- rdist(cbind(grid$x1,grid$x2),cbind(dat$x1,dat$x2))
R0<- c1*exp(-D0/a1)

#compute residuals at sampling points
resid <- dat$y-X%*%beta_GLS

#compute sigma at the sampling points
sigmadata <- (dat$x2<25)*sigma1+(dat$x2>25)*sigma2

#compute sigma at the prediction nodes
sigma0 <- (grid$x2<25)*sigma1+(grid$x2>25)*sigma2

pred <- varpred <- numeric(length=nrow(grid))

for (i in 1:nrow(grid)) {
  x0 <-c(1,grid$z[i])
  mu0 <- x0%*%beta_GLS
  r0 <- R0[i,]
  c0 <- r0*sigmadata*sigma0[i]
  
  Cc0 <- solve(V,c0)
  pred[i]<- mu0 + crossprod(Cc0, as.vector(resid))
  
  x_a <- x0 - crossprod(X,Cc0)
  varpred[i] <- c1*sigma0[i]^2 - crossprod(c0, Cc0) + crossprod(x_a,solve(XVX, x_a))
}

grid$pred<-pred
grid$krigvar <- varpred

hist(varpred)

library(ggplot2)
ggplot(data = grid) +
  geom_tile(mapping = aes(x = x1, y = x2, fill = pred)) +
  scale_x_continuous(name = "Easting (km)") +
  scale_y_continuous(ame = "Northing (km) \n") +
  scale_fill_continuous(name = "pred")+
  coord_equal(ratio = 1)

ggplot(data = grid) +
  geom_tile(mapping = aes(x = x1,y = x2,fill = krigvar)) +
  scale_x_continuous(name = "Easting (km)") +
  scale_y_continuous(name = "Northing (km) \n") +
  scale_fill_continuous(name = "var")+
  coord_equal(ratio = 1)

