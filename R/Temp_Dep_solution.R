
# T_Dep_Solution
#=================================================
# Find Td that corresponds to the input Cp298 value

#' Find root for Td - Temp. Dep. Solution
#'
#' @param Td Debye temp
#'
#' @return Find root for Td - Temp. Dep. Solution
#'
g <- function(Td){Debye_constant(Calphad.globals$RT,Td) - Calphad.globals$CP298/Calphad.globals$NOA}

#' calculate Td
#'
#' @param x limit
#' @param y limit
#'
#' @importFrom stats uniroot
#'
#' @return calculate Td
#'
calc_TD <- function(x,y){
  Calphad.globals$Td <- uniroot(g,c(x,y))$root
}
#=================================================
# optimize the fitting parameter (a) in (Td = a + b*T) based on the input S298 value

#' sum_square_error
#'
#' @param a parameter
#'
#' @return returns the sum square error
#'
#'
sum_square_error <- function(a) {
  sum(round(integrate(DebyeVarLinOT_1, 0,Calphad.globals$RT,a)$val - Calphad.globals$S298/Calphad.globals$NOA,2)^2)
}

#' optim_test
#'
#' @param opt optimized parameters
#'
#' @return optim test results
#'
#'
#'
optim_test <- function(opt){
  sum_square_error(opt[1])
}

# Optimize the function (int_value_opt) and if no solution could be found then save the
# element in a file called (elements_with_problems.Csv)

#' Optimize_TDep_soution
#'
#' @param init initial value to calculate parameter a1
#' @importFrom stats optim
#'
#' @return calculate parameters a1 and b1 for Temp. Dep. solution
#'
Optimize_TDep_soution <- function(init){
skip_to_next <- FALSE
tryCatch({

    result.opt_1 <-  optim(init,optim_test,method="N")

  }, error=function(e){skip_to_next <<- TRUE})
  if(skip_to_next) {
  elements_with_problems <- rbind(elements_with_problems, c("Symbol"=Calphad.globals$ele))
  print("Unique solution couldn't be evaluated using input parameters")
  }

  Calphad.globals$a1 <- result.opt_1$par[1]

  #=================================================
  # calculate the coefficient (b) in (Td = a + b*T)
  Calphad.globals$b1 <- (Calphad.globals$Td-Calphad.globals$a1)/Calphad.globals$RT

  # # Values to add to the final results dataframe
  # a1_coeff <- a1
  # b1_coeff <- b1
  # b_sol_coef <- NA

}

# #=======================================================================================================
# calculate results

# Cp function
#Debye_model <- Debye_model_above(x,Calphad.globals$Td)

#' Temp_Dep_Results
#'
#' @param T1 Temp
#' @param T2 Temp
#' @param T3 Temp
#'
#' @return Temp_Dep_Results
#'
Temp_Dep_Results <- function(T1,T2,T3){
  # Cp298
  Calphad.globals$cp_d <- Debye_model_above(T1)
  # Cp200
  Calphad.globals$cp_d200 <- Debye_model_above(T2)
  # Cp100
  Calphad.globals$cp_d100 <- Debye_model_above(T3)

  # S function for S298
  S_1_above <- integrate(entr_1_above, 0,Calphad.globals$bp)$val
  S_2_above <- integrate(entr_2_above, Calphad.globals$bp,T1)$val
  S_all_above <-  S_1_above + S_2_above
  # S298
  Calphad.globals$s_d = round(S_all_above,5)

  #------------------
  # S function for S200
  S_1_above200 <- integrate(entr_1_above, 0,Calphad.globals$bp)$val
  S_2_above200 <- integrate(entr_2_above, Calphad.globals$bp,T2)$val
  S_all_above200 <-  S_1_above200 + S_2_above200
  # S200
  Calphad.globals$s_d200 = round(S_all_above200,5)

  # S function for S100
  S_1_above100 <- integrate(entr_1_above, 0,Calphad.globals$bp)$val
  S_2_above100 <- integrate(entr_2_above, Calphad.globals$bp,T3)$val
  S_all_above100 <-  S_1_above100 + S_2_above100
  # S100
  Calphad.globals$s_d100 = round(S_all_above100,5)
  #=============================================
  # S up to RT
  # xrt <- seq(1,Calphad.globals$RT,1)
  # S_DF <- data.table(x = xrt)
  #
  # for (i in seq(along=xrt)){
  #   if(x[i]<Calphad.globals$bp){
  #     S <- integrate(entr_1_above, 0,xrt[i])$val
  #     S_DF$y[i] <-  S
  #   }
  #   if(x[i]>=40){
  #     S1 <- integrate(entr_1_above, 0,Calphad.globals$bp)$val
  #     S2 <- integrate(entr_2_above, Calphad.globals$bp,xrt[i])$val
  #     S_DF$y[i] <-  S1 + S2
  #   }
  # }
}

#------------------


