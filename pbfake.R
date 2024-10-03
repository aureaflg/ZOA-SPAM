
pbfake <- function(x, df = NULL, max.df = NULL, 
               inter = 20, degree= 3, order = 2,quantiles=F) 
{

bbase <- function(x, xl, xr, ndx, deg, quantiles=FALSE)
{
  tpower <- function(x, t, p)
    # Truncated p-th power function
    (x - t) ^ p * (x > t)
  # DS xl= min, xr=max, ndx= number of points within 
  # Construct B-spline basis
  # if quantiles=TRUE use different bases
  dx <- (xr - xl) / ndx # DS increment 
  if (quantiles) # if true use splineDesign
  {  #  this is not working and should be taken out
    knots <-  sort(c(seq(xl-deg*dx, xl, dx),quantile(x, prob=seq(0, 1, length=ndx)), seq(xr, xr+deg*dx, dx))) 
    B <- splineDesign(knots, x = x, outer.ok = TRUE, ord=deg+1)
    n <- dim(B)[2]  
    attr(B, "knots") <- knots[-c(1:(deg-1), (n-(deg-2)):n)]
    return(B)    
  }
  else # if false use Paul's
  { 
    knots <- seq(xl - deg * dx, xr + deg * dx, by = dx)
    P <- outer(x, knots, tpower, deg)# calculate the power in the knots
    n <- dim(P)[2]
    D <- diff(diag(n), diff = deg + 1) / (gamma(deg + 1) * dx ^ deg) # 
    B <- (-1) ^ (deg + 1) * P %*% t(D)
    attr(B, "knots") <- knots[-c(1:(deg-1), (n-(deg-2)):n)]
    B 
  }
}

no.dist.val <-  length(table(x))
#if (is.matrix(x)) stop("x is a matric declare it as a vector")
#lx <- length(x)
#control$inter <- if (lx<99) 10 else control$inter # this is to prevent singularities when length(x) is small:change to 99 30-11-11 MS
inter <- if (no.dist.val<=inter)  no.dist.val else inter 
xl <- min(x)
xr <- max(x)
xmax <- xr + 0.01 * (xr - xl)
xmin <- xl - 0.01 * (xr - xl)  
##                create the basis
X <- bbase(x, xmin, xmax, inter, degree, quantiles) # 
#        getting the base out    assign("X", X, envir = .GlobalEnv) 
r <- ncol(X)
##                the penalty matrix
D <- if(order==0) diag(r) else diff(diag(r), diff=order)
## ------      if df are set                
if(!is.null(df)) # degrees of freedom
{
  if (df>(dim(X)[2]-2)) 
  {df <- 3}  
  #warning("The df's exceed the number of columns of the design matrix", "\n",  "   they are set to 3") }
  if (df < 0)  #warning("the extra df's are set to 0")   
  df <- if (df < 0)  2  else  df+2
}
## -------- check max.df   (new 7-2018 MS)  
if (is.null(max.df)) max.df <- dim(X)[2]
if (max.df>(dim(X)[2])) 
{
  max.df <- dim(X)[2]
  warning("The max.df's are set to",  dim(X)[2],  "\n")
}   
xvar <- x
attr(xvar, "D") <- D
attr(xvar, "X") <- X
attr(xvar, "df")<- df
xvar
}
