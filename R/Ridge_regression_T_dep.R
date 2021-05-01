
#' Ridge_T_Dep_Solution
#'
#' @param a_input parameter
#'
#' @importFrom data.table data.table
#' @importFrom stats lm
#' @importFrom stats coef
#'
#' @return Ridge_T_Dep_Solution
#'
Ridge_T_Dep_Solution <- function(a_input){
  method <- "N" # optimization method ("Nelder-Mead"),

  # temperature ranges, or "Sections"
  x <- 0:Calphad.globals$RT # Temp. range
  x1 <- 0:Calphad.globals$bp # Section 1 (S1)
  x2 <- Calphad.globals$bp:Calphad.globals$RT # Section 2 (S2)

  #-------------------------------------------------------
  # Calculating Cp298 and S298 using Debye model

  cp298 <- Calphad.globals$NOA*DebyeVarLin_1(Calphad.globals$RT,a_input)

  int_value_1 <- function(a){integrate(DebyeVarLinOT_1, 0,Calphad.globals$RT,a)$val}
  s298 <- int_value_1(a_input) * Calphad.globals$NOA
  #-------------------------------------------------------
  # Creat a table (dt) for T and Cp values along the Temp. range
  dt <- data.table(x = x, y = Calphad.globals$NOA*DebyeVarLin_1(x,a_input))
  # Calculate S values and add them to the table
  for (i in seq(along=x))
    #dt$z[i] <- (integrate(Debye_constant_DOT, 0.0001,x[i])$value)
    dt$z[i] <- integrate(DebyeVarLinOT_1, 0.0001,x[i],a_input)$val * Calphad.globals$NOA

  # Creat a table (dt_cpot) for T and Cp/T values
  dt_cpot <- data.table(x = x, y = Calphad.globals$NOA*DebyeVarLin_1(x,a_input)/x)

  # Divide tables using the break point

  # First table (0-40 K)
  dat1 <- dt_cpot[x <= Calphad.globals$bp, ]
  dat1 <- dat1[-c(1), ]
  dat1_DF <- data.frame("x"=dat1$x,"y"=dat1$y)

  # Add origin point
  dat1 <- rbind(dat1_DF, c("x"=c(0), "y"=c(0)))

  # Arrange the table in ascending Order
  dat1 <- dat1[order(dat1$x),]

  # Second table (40-298.15 K)
  dat2 <- dt[x >= Calphad.globals$bp, ]
  #------------------------------------------------------
  # Fit the first section data (T vs Cp/T) using a fifth order polynomial to
  # provide initial values for the Ridge regression

  lm_poly1 <- lm(I(y) ~ x + I(x^2) +I(x^3)+I(x^4)+I(x^5)  , data = dat1)

  # S1 Initial values
  a1 <- 0
  b1 <- coef(lm_poly1)[1]
  c1 <- coef(lm_poly1)[2]
  d1 <- coef(lm_poly1)[3]
  e1 <- coef(lm_poly1)[4]
  f1 <- coef(lm_poly1)[5]
  g1 <- coef(lm_poly1)[6]

  # Ridge regression -  S1
  #-----------------------------------------------------
  # Obtain an approximation for lambda, while insuring that the (b1) value is above zero to
  # avoid having negative Cp values in the final result

  if (b1 < 0) {
    # Creat empty dataframe
    lambda_DF1 <- data.frame("Lambda"=c(0),"R2"=c(0))
    # potential Lambda values to loop over
    lambda1 <- c(1e-15, 1e-10,1e-9, 1e-8,1e-7, 1e-4, 1e-3,1e-2, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10,11,12,13,14,15,16,17,18,19, 20)
    # Find Lambda value
    for (i in lambda1){
      min.OLS <- function(data, par) {
        with(data, sum((   par[1]     +
                             par[2] *  x +
                             par[3] * (x^2) +
                             par[4] * (x^3) +
                             par[5] * (x^4) +
                             par[6] * (x^5) +
                             - y )^2+ i* abs(sum(par[1]+par[2]+par[3]+par[4]+par[5]+par[6])))
        )
      }

      result.opt <- optim(par = c(b1,c1,d1,e1,f1,g1),
                          min.OLS,
                          data = dat1,
                          lower=c(0,0,-5,-5,-5,-5),
                          upper=c(0.0000000001,1,5,5,5,5),
                          method=method
      )
      lambda_DF1 <- rbind(lambda_DF1, c("Lambda"=i, "R2"=result.opt$value))


    }
    lambda_DF1 <- lambda_DF1[-c(1), ]

    # select the lambda with minimum R2 corresponding value
    R2_min <- match(min(lambda_DF1$R2), lambda_DF1$R2)

    # Final Lambda value
    lambda_final1 <- lambda_DF1$Lambda[R2_min]
    lambda_final1 <- as.double(lambda_final1)

    # Apply Ridge regression using the obtained Lambda value
    min.OLS <- function(data, par) {
      with(data, sum((   par[1]     +
                           par[2] *  x +
                           par[3] * (x^2) +
                           par[4] * (x^3) +
                           par[5] * (x^4) +
                           par[6] * (x^5) +
                           - y )^2+ lambda_final1* abs(sum(par[1]+par[2]+par[3]+par[4]+par[5]+par[6])))
      )
    }

    result.opt <- optim(par = c(b1,c1,d1,e1,f1,g1),
                        min.OLS,
                        data = dat1,
                        lower=c(0,0,-5,-5,-5,-5),
                        upper=c(0.0000000001,1,5,5,5,5),
                        method=method
    )

    # Final S1 Coifficients
    Calphad.globals$T1S1 <- 0
    Calphad.globals$T2S1 <- result.opt$par[1]
    Calphad.globals$T3S1 <- result.opt$par[2]
    Calphad.globals$T4S1 <- result.opt$par[3]
    Calphad.globals$T5S1 <- result.opt$par[4]
    Calphad.globals$T6S1 <- result.opt$par[5]
    Calphad.globals$T7S1 <- result.opt$par[6]

  }
  else{
    Calphad.globals$T1S1 <- a1
    Calphad.globals$T2S1 <- b1
    Calphad.globals$T3S1 <- c1
    Calphad.globals$T4S1 <- d1
    Calphad.globals$T5S1 <- e1
    Calphad.globals$T6S1 <- f1
    Calphad.globals$T7S1 <- g1

  }
  #------------------------------------------------------
  # Fit the Second section data (T vs Cp/T) using a forth order polynomial + x^-1 + x^-2 to
  # provide initial values for the Ridge regression

  lm_poly2 <- lm(I(y) ~  x + I(x^2)+ I(x^3)+ I(x^4) +I(x^-1) +I(x^-2) , data = dat2)

  # S2 Initial values
  b2 <- coef(lm_poly2)[2]
  c2 <- coef(lm_poly2)[3]
  d2 <- coef(lm_poly2)[4]
  e2 <- coef(lm_poly2)[5]
  f2 <- coef(lm_poly2)[6]
  g2 <- coef(lm_poly2)[7]

  # a2 coefficient function, used to insure a continous curve at the break point
  a2 <-function(x,b,c,d,e,f,g) {(x*(b1-b)+x^2*(c1-c)+x^3*(d1-d)+x^4*(e1-e)+x^5*(f1)+x^6*(g1)- (f/x) - (g/(x^2)))}

  # a2 initial value
  a2p <- a2(Calphad.globals$bp,b2,c2,d2,e2,f2,g2)


  #------------------------------------------------------
  # Error function for Cp

  errorcp <- function(b,c,d,e,f,g){
    a2p1 <- a2(Calphad.globals$bp,b,c,d,e,f,g)
    error <- polyFunc_512pm(x2,a2p1,b,c,d,e,f,g)
    return(error)
  }

  # Entropy functions for S1 and S2
  S1 <- function(x,b,c,d,e,f,g){(b*x+1/2*c*x^2+1/3*d*x^3+1/4*e*x^4+1/5*f*x^5+1/6*g*x^6)}
  Sbp <- S1(Calphad.globals$bp,b1,c1,d1,e1,f1,g1)
  S2 <- function(x,a,b,c,d,e,f,g){(a*log(x/Calphad.globals$bp)+b*(x-Calphad.globals$bp)+1/2*c*(x^2-Calphad.globals$bp^2)+
                                     1/3*d*(x^3-Calphad.globals$bp^3)+1/4*e*(x^4-Calphad.globals$bp^4)-f*(((1/1)*x^(-1))-((1/1)*Calphad.globals$bp^(-1)))
                                   -g*(((1/2)*x^(-2))-((1/2)*Calphad.globals$bp^(-2)))+Sbp)}

  # Error function for S
  errorS <- function(b,c,d,e,f,g){
    a2p1 <- a2(Calphad.globals$bp,b,c,d,e,f,g)
    error <- S2(x2,a2p1,b,c,d,e,f,g)
    return(error)
  }

  #-------------------------------------------------------
  # Error function for Cp298

  errorcp298 <- function(b,c,d,e,f,g) {
    #f2p2 <- f2(bp,b,c,d,e)
    a2p2 <- a2(Calphad.globals$bp,b,c,d,e,f,g)
    error <- polyFunc_512pm(Calphad.globals$RT,a2p2,b,c,d,e,f,g)
    return(error)
  }

  # Error function for S298
  errorS298 <- function(b,c,d,e,f,g) {
    #f2p2 <- f2(bp,b,c,d,e)
    a2p2 <- a2(Calphad.globals$bp,b,c,d,e,f,g)
    error <- S2(Calphad.globals$RT,a2p2,b,c,d,e,f,g)
    return(error)
  }

  #-----------------------------------------------
  # Ridge regression -  S2
  #-----------------------------------------------
  # Obtain an approximation for lambda

  lambda_DF <- data.frame("Lambda"=c(0),"R2"=c(0))
  lambda <- c(1e-45,1e-35,1e-25,1e-15, 1e-10, 1e-8, 1e-4, 1e-3,1e-2, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10,11,12,13,14,15,16,17,18,19, 20,30,40,50,60,70)

  # optimization of parameters
  for (i in lambda){
    minimize <- function(data, par) {
      with(data, sum((errorcp(par[1],par[2],par[3],par[4],par[5],par[6])-dat2$y)^2) +
             sum((errorS(par[1],par[2],par[3],par[4],par[5],par[6])-dat2$z)^2) +
             10*sqrt((errorcp298(par[1],par[2],par[3],par[4],par[5],par[6])-cp298)^2 +
                       (i* sum(abs(par[1])+abs(par[2])+abs(par[3])))) +
             30*sqrt((errorS298(par[1],par[2],par[3],par[4],par[5],par[6])-s298)^2)
           + (i* sum(abs(par[1])+abs(par[2])+abs(par[3])+abs(par[4])+abs(par[5])+abs(par[6])))
      )
    }

    result.opt1 <- optim(par = c(b2,c2,d2,e2,f2,g2),
                         minimize,
                         data = dat2
                         ,method=method
                         #,lower=c(-50,-50,0,0),
                         #upper=c(50,50,0.00000000000000001,0.00000000000000001)
                         ,control=list(trace=1, maxit=1000000000)
    )

    lambda_DF <- rbind(lambda_DF, c("Lambda"=i, "R2"=result.opt1$value))


  }
  lambda_DF <- lambda_DF[-c(1), ]

  # select the lambda with minimum R2 corresponding value
  R2_min <- match(min(lambda_DF$R2), lambda_DF$R2)

  # Final Lambda value
  lambda_final <- lambda_DF$Lambda[R2_min]
  lambda_final <- as.double(lambda_final)

  # Apply Ridge regression using the obtained Lambda value
  minimize <- function(data, par) {
    with(data, sum((errorcp(par[1],par[2],par[3],par[4],par[5],par[6])-dat2$y)^2) +
           sum((errorS(par[1],par[2],par[3],par[4],par[5],par[6])-dat2$z)^2) +
           10*sqrt((errorcp298(par[1],par[2],par[3],par[4],par[5],par[6])-cp298)^2) +
           30*sqrt((errorS298(par[1],par[2],par[3],par[4],par[5],par[6])-s298)^2)
         + (lambda_final* sum(abs(par[1])+abs(par[2])+abs(par[3])+abs(par[4])+abs(par[5])+abs(par[6])))
    )
  }

  result.opt1 <- optim(par = c(b2,c2,d2,e2,f2,g2),
                       minimize,
                       data = dat2
                       ,method=method
                       #,lower=c(-50,-50,0,0),
                       #upper=c(50,50,0.00000000000000001,0.00000000000000001)
                       ,control=list(trace=1, maxit=1000000000)
  )

  # Final S2 Coefficients
  M22op <- result.opt1$par[1]
  M23op <- result.opt1$par[2]
  M24op <- result.opt1$par[3]
  M25op <- result.opt1$par[4]
  M26op <- result.opt1$par[5]
  M27op <- result.opt1$par[6]
  M21op <- a2(Calphad.globals$bp,M22op,M23op,M24op,M25op,M26op,M27op)

  Calphad.globals$T1S2 <- result.opt1$par[1]
  Calphad.globals$T2S2 <- result.opt1$par[2]
  Calphad.globals$T3S2 <- result.opt1$par[3]
  Calphad.globals$T4S2 <- result.opt1$par[4]
  Calphad.globals$T5S2 <- result.opt1$par[5]
  Calphad.globals$T6S2 <- result.opt1$par[6]
  Calphad.globals$T0S2 <- a2(Calphad.globals$bp,M22op,M23op,M24op,M25op,M26op,M27op)

  # Calculate the MTB Error

  MTB_Error <- sum((errorcp(Calphad.globals$T0S2,Calphad.globals$T1S2,Calphad.globals$T2S2,
                            Calphad.globals$T3S2,Calphad.globals$T4S2,Calphad.globals$T5S2)-dat2$y)^2) +
    sum((errorS(Calphad.globals$T0S2,Calphad.globals$T1S2,Calphad.globals$T2S2,
                Calphad.globals$T3S2,Calphad.globals$T4S2,Calphad.globals$T5S2)-dat2$z)^2) +
    10*sqrt((errorcp298(Calphad.globals$T0S2,Calphad.globals$T1S2,Calphad.globals$T2S2,
                        Calphad.globals$T3S2,Calphad.globals$T4S2,Calphad.globals$T5S2)-cp298)^2) +
    30*sqrt((errorS298(Calphad.globals$T0S2,Calphad.globals$T1S2,Calphad.globals$T2S2,
                       Calphad.globals$T3S2,Calphad.globals$T4S2,Calphad.globals$T5S2)-s298)^2)

  Calphad.globals$MTB_Error <- MTB_Error

}
