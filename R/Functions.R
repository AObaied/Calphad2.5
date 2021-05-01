
# This file contains all the functions used in this model

#=====================================================
# polynomial functions
#=====================================================
#' 6th order polynomial
#'
#' @param x Temperature
#' @param a fitting parameter
#' @param b fitting parameter
#' @param c fitting parameter
#' @param d fitting parameter
#' @param e fitting parameter
#' @param f fitting parameter
#' @param g fitting parameter
#'
#' @return 6th order polynomial
#'
polyFunc_6 <- function(x,a,b,c,d,e,f,g){(a + b*x + c*x^2 + d*x^3+ e*x^4+ f*x^5+ g*x^6 )}
#' 5th order polynomial with x^-1 and x^-2
#'
#' @param x Temperature
#' @param a fitting parameter
#' @param b fitting parameter
#' @param c fitting parameter
#' @param d fitting parameter
#' @param e fitting parameter
#' @param f fitting parameter
#' @param g fitting parameter
#'
#' @return 5th order polynomial with x^-1 and x^-2
#'
polyFunc_512pm <- function(x,a,b,c,d,e,f,g){(a + b*x + c*x^2 + d*x^3+ e*x^4 +f*x^(-1) +g*x^(-2) )}

#=====================================================
# Debye model
#=====================================================

#' Integral in Debye model
#'
#' @param x temperature
#'
#' @return Integral in Debye model
#'
intDebye <- function(x) { x^4*exp(-x)/(1-exp(-x))^2 }

#' Debye model
#'
#' @param x temperature
#' @param Td Debye temperature
#'
#' @importFrom stats integrate
#'
#' @return Debye model
#' @export
#'
#' @examples
#' Debye_constant(500,250)
Debye_constant <- function(x, Td) {
  #
  Xd    <-  Td/x
  DInt  <-  vector(length=length(Xd), mode="numeric")
  koef  <-  9*8.314  # 9R - physical constant
  #
  for (i in seq(along=Xd))
    DInt[i] <- as.numeric(integrate(intDebye, 0, Xd[i])$value)
  return((DInt*koef/(Xd)^3))}

#' Debye model - with linear term
#'
#' @param x temperature
#' @param Td Debye temperature
#' @param b linear term
#'
#' @importFrom stats integrate
#'
#' @return Debye model - with linear term
Debye_constant_b <- function(x, Td,b) {
  #
  Xd    <-  Td/x
  DInt  <-  vector(length=length(Xd), mode="numeric")
  koef  <-  9*8.314  # 9R - physical constant
  #
  for (i in seq(along=Xd))
    DInt[i] <- as.numeric(integrate(intDebye, 0, Xd[i])$value)
  return((DInt*koef/(Xd)^3)+b*x)}



# Debye model divided by T. This is integrated to find the Entropy of the Debye model

#' Debye model divided by T
#'
#' @param x Temperature
#'
#' @return Debye model divided by T
#'
Debye_constant_DOT <- function(x){
  CP <- vector(length = length(x), mode = "numeric")
  for(i in seq(along = x))
    CP[i] <- Calphad.globals$NOA*Debye_constant(x[i],Calphad.globals$Td)/x[i]
  return(CP)
}

#' Debye model divided by T - with linear term
#'
#' @param x Temperature
#'
#' @return Debye model divided by T - with linear term
#'
Debye_constant_DOT_b <- function(x){
  CP <- vector(length = length(x), mode = "numeric")
  for(i in seq(along = x))
    CP[i] <- Calphad.globals$NOA*Debye_constant_b(x[i],Calphad.globals$Td,Calphad.globals$b_sol_coef)/x[i]
  return(CP)
}

#=====================================================
# T_Dep_solution functions
#=====================================================

#' Temp dependent debye model
#'
#' @param x Temperature
#' @param a constant parameter
#'
#' @importFrom stats integrate
#'
#' @return Temp dependent debye model
#'
DebyeVarLin_1 <- function(x, a) {
  #
  b <- (Calphad.globals$Td - a)/Calphad.globals$RT
  Td_temp    <-  a+b*x
  Xd    <-  Td_temp/x
  DInt  <-  vector(length=length(Xd), mode="numeric")
  koef  <-  9*8.314  # 9R - physical const
  #
  for (i in seq(along=Xd))
    DInt[i] <- as.numeric(integrate(intDebye, 0, Xd[i])$value)
  return(DInt*koef/(Td_temp/x)^3)
}

#' Temp dependent debye model divided by T. used to find the entropy
#'
#' @param x Temperature
#' @param a constant parameter
#'
#' @return Temp dependent debye model divided by T. used to find the entropy
#'
DebyeVarLinOT_1 <- function(x, a) {
  #
  b <- (Calphad.globals$Td - a)/Calphad.globals$RT
  Td_temp    <-  a+b*x
  Xd    <-  Td_temp/x
  DInt  <-  vector(length=length(Xd), mode="numeric")
  koef  <-  9*8.314  # 9R - physical const
  #
  for (i in seq(along=Xd))
    DInt[i] <- as.numeric(integrate(intDebye, 0, Xd[i])$value)
  return( (DInt*koef/(Td_temp/x)^3)/x )
}

#' Debye_model Temp dependent solution
#'
#' @param x Temperature
#'
#' @return Debye_model Temp dependent solution
#'
Debye_model_above <- function(x){
  f <- vector(length = length(x), mode = "numeric")
  for(i in seq(along = x))
    if(x[i] >= 0 && x[i] < Calphad.globals$bp)
      f[i] <- S1_above(x[i])
  else if(x[i] >= Calphad.globals$bp && x[i] <= Calphad.globals$RT)
    f[i] <- S2_above(x[i])
  return(f)
}

#' polynomial function for the first section of the fitted debye model (0-40)
#'
#' @param x #' @param x Temperature
#'
#' @return polynomial function for the first section of the fitted debye model (0-40)
#'
#' @export
#'
S1_above <- function(x){
  f <- vector(length = length(x), mode = "numeric")
  for(i in seq(along = x))
    f[i] <- polyFunc_6(x[i],Calphad.globals$T1S1,Calphad.globals$T2S1,Calphad.globals$T3S1,
                       Calphad.globals$T4S1,Calphad.globals$T5S1,Calphad.globals$T6S1,
                       Calphad.globals$T7S1)
  return(f)
}

#' polynomial function for the second section of the fitted debye model (40-298.15)
#'
#' @param x Temperature
#'
#' @return polynomial function for the second section of the fitted debye model (40-298.15)
#' @export
#'
S2_above <- function(x){
  f <- vector(length = length(x), mode = "numeric")
  for(i in seq(along = x))
    f[i] <- polyFunc_512pm(x[i],Calphad.globals$T0S2,Calphad.globals$T1S2,
                            Calphad.globals$T2S2,Calphad.globals$T3S2,
                            Calphad.globals$T4S2,Calphad.globals$T5S2,
                            Calphad.globals$T6S2)
  return(f)
}

#' polynomial function for the first section of the fitted debye model (0-40) divided by T - Entropy
#'
#' @param x Temperature
#'
#' @return polynomial function for the first section of the fitted debye model (0-40) divided by T - Entropy
#' @export
#'
entr_1_above <- function(x){
  CP <- vector(length = length(x), mode = "numeric")
  for(i in seq(along = x))
    CP[i] <- S1_above(x[i])/x[i]
  return(CP)
}

#' polynomial function for the second section of the fitted debye model (40-298.15) divided by T - Entropy
#'
#' @param x Temperature
#'
#' @return polynomial function for the second section of the fitted debye model (40-298.15) divided by T - Entropy
#'
entr_2_above <- function(x){
  CP <- vector(length = length(x), mode = "numeric")
  for(i in seq(along = x))
    CP[i] <- S2_above(x[i])/x[i]
  return(CP)
}

# #====================================================================
# # b_Solution
# #=====================================================
# # polynomial function for the fitted debye model plus the b term
#
#' Debye_model_b_solution
#'
#' @param x Temp
#'
#' @return Debye_model_b_solution
#'
Debye_model_below <- function(x){
  f <- vector(length = length(x), mode = "numeric")
  for(i in seq(along = x))
    if(x[i] >= 0 && x[i] < Calphad.globals$bp)
      f[i] <- S1_below(x[i])
  else if(x[i] >= Calphad.globals$bp && x[i] <= Calphad.globals$RT)
    f[i] <- S2_below(x[i])
  return(f)
}

#' polynomial function for the first section of the fitted debye model (0-40)
#'
#' @param x Temp
#'
#' @return polynomial function for the first section of the fitted debye model (0-40)
#'
S1_below <- function(x){
  f <- vector(length = length(x), mode = "numeric")
  for(i in seq(along = x))
    f[i] <- polyFunc_6(x[i],Calphad.globals$T1S1,Calphad.globals$T2S1,Calphad.globals$T3S1,
                       Calphad.globals$T4S1,Calphad.globals$T5S1,Calphad.globals$T6S1,
                       Calphad.globals$T7S1)
  return(f+ Calphad.globals$b_sol_coef*x)
}

#' polynomial function for the second section of the fitted debye model (40-298.15)
#'
#' @param x Temp
#'
#' @return polynomial function for the second section of the fitted debye model (40-298.15)
#'
S2_below <- function(x){
  f <- vector(length = length(x), mode = "numeric")
  for(i in seq(along = x))
    f[i] <- polyFunc_512pm(x[i],Calphad.globals$T0S2,Calphad.globals$T1S2,
                            Calphad.globals$T2S2,Calphad.globals$T3S2,
                            Calphad.globals$T4S2,Calphad.globals$T5S2,
                            Calphad.globals$T6S2)
  return(f+ Calphad.globals$b_sol_coef*x)
}

#' polynomial function for the first section of the fitted debye model (0-40) divided by T
#'
#' @param x Temp
#'
#' @return polynomial function for the first section of the fitted debye model (0-40) divided by T
#'
entr_1_below <- function(x){
  CP <- vector(length = length(x), mode = "numeric")
  for(i in seq(along = x))
    CP[i] <- S1_below(x[i])/x[i]
  return(CP)
}


#' polynomial function for the second section of the fitted debye model (40-298.15) divided by T
#'
#' @param x Temp
#'
#' @return polynomial function for the second section of the fitted debye model (40-298.15) divided by T
#'
entr_2_below <- function(x){
  CP <- vector(length = length(x), mode = "numeric")
  for(i in seq(along = x))
    CP[i] <- S2_below(x[i])/x[i]
  return(CP)
}
