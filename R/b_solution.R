
# b_Solution

# read debye model fitting coefficients from csv file (Cp_S_Coefficients.csv)
# This Csv file contains the coefficient obtained from fitting the debye model results, the (Cp298 vs. Td) and
# (S298 vs. Td) curves, for a Td range of (40 - 1800).
# Solving these two systems of equations will allow us to obtain the Td and b values



#' Obtain Td and b for b_solution
#'
#' @param L1 parameter
#' @param L2 parameter
#'
#' @importFrom stats uniroot
#' @importFrom utils read.csv
#'
#' @return  Solve the 2 systems of equations to obtain Td and b if there are no solution obtained,
#' then the element with the problem will be saved into a Csv file named (elements_with_problems.Csv)
#'
calc_TD_b <- function(L1,L2){
  D_Coefficients <- system.file("Subsidiary_Csv_files","Cp_S_Coefficients.csv", package="Calphad2.5")
  D_Coefficients <- read.csv(D_Coefficients, header=F,skip=1)$V1

  eq_1 <- function(x){(Calphad.globals$CP298/Calphad.globals$NOA - D_Coefficients[1] - D_Coefficients[2]*x - D_Coefficients[3]*x^2 - D_Coefficients[4]*x^3-D_Coefficients[5]*x^4-D_Coefficients[6]*x^5-D_Coefficients[7]*x^6-D_Coefficients[8]*x^7-D_Coefficients[9]*x^8-D_Coefficients[10]*x^9-D_Coefficients[11]*x^10)/Calphad.globals$RT}
  eq_2 <- function(x){(Calphad.globals$S298/Calphad.globals$NOA - D_Coefficients[12] - D_Coefficients[13]*x - D_Coefficients[14]*x^2 - D_Coefficients[15]*x^3 - D_Coefficients[16]*x^4- D_Coefficients[17]*x^5- D_Coefficients[18]*x^6 - D_Coefficients[19]*x^7- D_Coefficients[20]*x^8- D_Coefficients[21]*x^9- D_Coefficients[22]*x^10- D_Coefficients[23]*x^11- D_Coefficients[24]*x^12)/Calphad.globals$RT}
  #--------------------------------------------------
  skip_to_next <- FALSE
  tryCatch({
    Calphad.globals$Td <- round(uniroot(function(x) eq_2(x)-eq_1(x),c(L1,L2),
                                        extendInt = "yes")$root,1)
  }, error=function(e){skip_to_next <<- TRUE})
  if(skip_to_next) {
    elements_with_problems <- rbind(elements_with_problems, c("Symbol"=Calphad.globals$ele))
    print("Unique solution couldn't be evaluated using input parameters")
  }
  Calphad.globals$b_sol_coef <- eq_1(Calphad.globals$Td)*Calphad.globals$NOA
}

#---------------------------------------------------
# calculate results



#' b_solution_Results
#'
#' @param T1 Temp
#' @param T2 Temp
#' @param T3 Temp
#'
#' @return b_solution_Results
#'
b_solution_Results <- function(T1,T2,T3){
# Cp298
Calphad.globals$cp_d = Debye_model_below(T1)

# Cp200
Calphad.globals$cp_d200 = Debye_model_below(T2)
# Cp100
Calphad.globals$cp_d100 = Debye_model_below(T3)
#------------------
# S function for S298
S_1_below <- integrate(entr_1_below, 0,Calphad.globals$bp)$val
S_2_below <- integrate(entr_2_below, Calphad.globals$bp,T1)$val
S_all_below <-  S_1_below + S_2_below
# S298
Calphad.globals$s_d = round(S_all_below,5)
#------------------
# S function for S200
S_1_below200 <- integrate(entr_1_below, 0,Calphad.globals$bp)$val
S_2_below200 <- integrate(entr_2_below, Calphad.globals$bp,T2)$val
S_all_below200 <-  S_1_below200 + S_2_below200
# S200
Calphad.globals$s_d200 = round(S_all_below200,5)

# S function for S100
S_1_below100 <- integrate(entr_1_below, 0,Calphad.globals$bp)$val
S_2_below100 <- integrate(entr_2_below, Calphad.globals$bp,T3)$val
S_all_below100 <-  S_1_below100 + S_2_below100
# S100
Calphad.globals$s_d100 = round(S_all_below100,5)
#=============================================
# S upto RT
# xrt <- seq(1,298.15,1)
# S_DF <- data.table(x = xrt)
#
# for (i in seq(along=xrt)){
#   if(x[i]<40){
#     S <- integrate(entr_1_below, 0,xrt[i])$val
#     S_DF$y[i] <-  S
#   }
#   if(x[i]>=40){
#     S1 <- integrate(entr_1_below, 0,40)$val
#     S2 <- integrate(entr_2_below, 40,xrt[i])$val
#     S_DF$y[i] <-  S1 + S2
#   }
# }
}

