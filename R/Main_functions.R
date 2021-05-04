Calphad.globals <- new.env()

Calphad.globals$RT <- 298.15
Calphad.globals$bp <- 40
Calphad.globals$ele <- NULL
Calphad.globals$S_diff <- NULL
Calphad.globals$CP100 <- NULL
Calphad.globals$CP200 <- NULL
Calphad.globals$CP298 <- NULL
Calphad.globals$S100 <- NULL
Calphad.globals$S200 <- NULL
Calphad.globals$S298 <- NULL
Calphad.globals$NOA <- NULL
Calphad.globals$b_sol_coef <- NULL
Calphad.globals$a1 <- NULL
Calphad.globals$b1 <- NULL
Calphad.globals$MTB_Error <- NULL


Calphad.globals$T1S1 <- NULL
Calphad.globals$T2S1 <- NULL
Calphad.globals$T3S1 <- NULL
Calphad.globals$T4S1 <- NULL
Calphad.globals$T5S1 <- NULL
Calphad.globals$T6S1 <- NULL
Calphad.globals$T7S1 <- NULL
Calphad.globals$T1S1 <- NULL
Calphad.globals$T0S2 <- NULL
Calphad.globals$T1S2 <- NULL
Calphad.globals$T2S2 <- NULL
Calphad.globals$T3S2 <- NULL
Calphad.globals$T4S2 <- NULL
Calphad.globals$T5S2 <- NULL
Calphad.globals$T6S2 <- NULL

Calphad.globals$cp_d <- NULL
Calphad.globals$cp_d200 <- NULL
Calphad.globals$cp_d100 <- NULL
Calphad.globals$s_d <- NULL
Calphad.globals$s_d200 <- NULL
Calphad.globals$s_d100 <- NULL


#' Calphad_2.5
#'
#' @param ele input parameter
#' @param CP100 input parameter
#' @param CP200 input parameter
#' @param CP298 input parameter
#' @param S100 input parameter
#' @param S200 input parameter
#' @param S298 input parameter
#'
#' @importFrom utils read.table
#'
#' @return Result
#'
#' @export
#'
Calphad_2.5 <- function(ele,CP100,CP200,CP298,S100,S200,S298, env = Calphad.globals){

  # This dataframe will contain the final results
  Results_DB <- NA

  # This dataframe will contain a list with the elements/compounds that
  # had problems and the model wasn't able to assess correctly
  elements_with_problems <- NA

  # Starting the for loop over the list of elements/compounds listed in the input file
  #for(i in Data_Input$Element){
  Calphad.globals$ele <- ele
    #=======================================
    # Read the data entries for each element/compound from the input file

    Calphad.globals$CP298 <- CP298
    Calphad.globals$S298 <- S298

    if(missing(CP100)) {
      Calphad.globals$CP100 <- "-"
    } else {
      Calphad.globals$CP100 <- CP100
    }

    if(missing(CP200)) {
      Calphad.globals$CP200 <- "-"
    } else {
      Calphad.globals$CP200 <- CP200
    }

    if(missing(S100)) {
      Calphad.globals$S100 <- "-"
    } else {
      Calphad.globals$S100 <- S100
    }

    if(missing(S200)) {
      Calphad.globals$S200 <- "-"
    } else {
      Calphad.globals$S200 <- S200
    }


    #======================================
    # Count the number of atoms

    ele_number <- strsplit(x=ele, split=" ")
    ele_number <- ele_number[[1]][1]
    Calphad.globals$NOA <- count(ele_number) # NOA = Number Of Atoms

    #================================
    # Choose a solutions
    #
    # (24.942) is the Cp value of the Debye model if the Debye temperature (Td) is equal to (0.01)
    # at the room temperature. If the Cp value is higher than this value (24.942),
    # then the element's (Cp,S) position is below the cp298 Vs. S298 line created by plotting
    # the results (Cp and S) of the debye model for the Td range (0.1 - 1800)
    # ==> if the Cp value is above (24.942), then the element/compound should be assessed
    # using the (Debye model + bx) solution (will be referred to as [b_solution] from now on).
    #------------------------------
    # If the Cp value is higher than (24.942), we check the S298 value of the Debye model
    # corresponding to the input Cp value from the file named (Tabulated_Debye_model.Csv),
    # which contains tabulated results (Cp and S) for a Td range (40-1800).
    # if the Difference between the S298_Debye and the S298_input is negative then we keep using
    # the b_Solution, but if it is positive then we use (temperature dependent_debye,
    # denoted as the T_Dep_solution).
    # This solution is used if the element's (Cp,S) position is above the cp298 Vs. S298 line.

    if (Calphad.globals$CP298/Calphad.globals$NOA > 24.942){ # b solution ==> entropy difference is set to (-1) to select the b_solution
      Calphad.globals$S_diff <- -1}
    if (Calphad.globals$CP298/Calphad.globals$NOA < 24.942){

      # Read (Tabulated_Debye_model.csv) file
      File <- system.file("Subsidiary_Csv_files","Tabulated_Debye_model.csv", package="Calphad2.5")
      Data   <-  read.table(File, header = F, dec = ".", col.names=c("Td", "CP298","S298"), stringsAsFactors = F, skip=1, sep=",")
      CPDF <- data.frame("Td"=Data$Td,"CP298_D"=Data$CP298,"S298_D"=Data$S298)

      # Find corresponding S298 for the Debye model
      S298_D <- mean(CPDF[round(CPDF$CP298_D,1) == round(Calphad.globals$CP298/Calphad.globals$NOA,1), "S298_D"]) # used mean because the D_table in not accurate
      Calphad.globals$S_diff <- Calphad.globals$S298/Calphad.globals$NOA - S298_D

    }

    #------------------------------
    # Choose a solution based on the difference in entropy (S_diff)

    if (Calphad.globals$S_diff > 0){ # Choose T_Dep_solution
      calc_TD(-10,2500)
      hush(Optimize_TDep_soution(100))
      hush(Ridge_T_Dep_Solution(Calphad.globals$a1))
      parametersS2 <- c(Calphad.globals$T0S2,Calphad.globals$T1S2,Calphad.globals$T2S2,
                        Calphad.globals$T3S2,Calphad.globals$T4S2,Calphad.globals$T5S2,
                        Calphad.globals$T6S2)
      parametersS1 <- c(Calphad.globals$T1S1,Calphad.globals$T2S1,Calphad.globals$T3S1,
                      Calphad.globals$T4S1,Calphad.globals$T5S1,Calphad.globals$T6S1,
                      Calphad.globals$T7S1)
      hush(Temp_Dep_Results(Calphad.globals$RT,200,100))
      cps <- c(Calphad.globals$cp_d,Calphad.globals$cp_d100,Calphad.globals$cp_d200)
      #print(cps)
      ss <- c(Calphad.globals$s_d,Calphad.globals$s_d100,Calphad.globals$s_d200)
      #print(ss)


    }

    if (Calphad.globals$S_diff < 0){ # Choose b_Solution
      calc_TD_b(0,1000)
      hush(Ridge_b_Solution(Calphad.globals$Td))
      hush(b_solution_Results(Calphad.globals$RT,200,100))
    }
    #=====================================================================
    # plotting the results and saving them as a png figure in the (Results) folder

    # Uncomment to save as png
    #png(paste("Results/",ele,".png", sep = ""), width = 5, height = 5, units = 'in', res = 700)

    # Uncomment to view result
    #source(paste(getwd(),"/Subsidiary_R_files/Plotting.R",sep = ""),local = F)
    #source(paste(getwd(),"/Subsidiary_R_files/Plotting_paper_S_CP.R",sep = ""),local = F)
    #source(paste(getwd(),"/Subsidiary_R_files/Plotting_paper_S_CP_Debye_only.R",sep = ""),local = F)
    #source(paste(getwd(),"/Subsidiary_R_files/Plotting_paper_S_CP_Debye_100_250_500.R",sep = ""),local = F)

    # Uncomment to save as png
    #dev.off()

    #=====================================================================
    # Saving the results in a dataframe
    # Results_DB <-rbind(Results_DB, c("Element" = Calphad.globals$ele,
    #                                  "CP100_model"=Calphad.globals$cp_d100,"CP100_Input"=Calphad.globals$CP100,
    #                                  "CP200_model"=Calphad.globals$cp_d200,"CP200_Input"=Calphad.globals$CP200,
    #                                  "CP298_model"=Calphad.globals$cp_d,"CP298_Input"=Calphad.globals$CP298,
    #                                  "S100_model"=Calphad.globals$s_d100,"S100_Input"=Calphad.globals$S100,
    #                                  "S200_model"=Calphad.globals$s_d200,"S200_Input"=Calphad.globals$S200,
    #                                  "S298_model"=Calphad.globals$s_d,"S298_Input"=Calphad.globals$S298,
    #                                  "Td"=Calphad.globals$Td, "a1_T_dep_sol"=Calphad.globals$a1,"b1_T_dep_sol"=Calphad.globals$b1,
    #                                  "b_sol_coef" = Calphad.globals$b_sol_coef,'NOA'=Calphad.globals$NOA,'S-diff'=Calphad.globals$S_diff,
    #                                  "a_S1" = as.numeric(Calphad.globals$T1S1), "b_S1" = as.numeric(Calphad.globals$T2S1), "c_S1" = as.numeric(Calphad.globals$T3S1),
    #                                  "d_S1" = as.numeric(Calphad.globals$T4S1), "e_S1" = as.numeric(Calphad.globals$T5S1), "f_S1" = as.numeric(Calphad.globals$T6S1), "g_S1" = as.numeric(Calphad.globals$T7S1),
    #                                  "a_S2" = as.numeric(Calphad.globals$T0S2), "b_S2" = as.numeric(Calphad.globals$T1S2), "c_S2" = as.numeric(Calphad.globals$T2S2), "d_S2" = as.numeric(Calphad.globals$T3S2),
    #                                  "e_S2" = as.numeric(Calphad.globals$T4S2), "f_S2" = as.numeric(Calphad.globals$T5S2), "g_S2" = as.numeric(Calphad.globals$T6S2)))
    # Results_DB <- Results_DB[-1,]
    #write.csv(Results_DB, file = "Results",Calphad.globals$ele,".csv", row.names=FALSE)
    # Results_input <- data.frame("Element" = Calphad.globals$ele, "CP100_Input"=Calphad.globals$CP100, "CP200_Input"=Calphad.globals$CP200, "CP298_Input"=Calphad.globals$CP298,
    #                             "S100_Input"=Calphad.globals$S100, "S200_Input"=Calphad.globals$S200, "S298_Input"=Calphad.globals$S298)
    # Results_output_parameters <- NA
    # Results_output_parameters <- data.frame("Td"=Calphad.globals$Td, "a1_T_dep_sol"=Calphad.globals$a1,"b1_T_dep_sol"=Calphad.globals$b1,
    #                                         "b_sol_coef" = Calphad.globals$b_sol_coef,'Number of atoms'=Calphad.globals$NOA)
    # Results_output_poly <- data.frame( "a_S1" = as.numeric(Calphad.globals$T1S1), "b_S1" = as.numeric(Calphad.globals$T2S1), "c_S1" = as.numeric(Calphad.globals$T3S1),
    #                                    "d_S1" = as.numeric(Calphad.globals$T4S1), "e_S1" = as.numeric(Calphad.globals$T5S1), "f_S1" = as.numeric(Calphad.globals$T6S1), "g_S1" = as.numeric(Calphad.globals$T7S1),
    #                                    "a_S2" = as.numeric(Calphad.globals$T0S2), "b_S2" = as.numeric(Calphad.globals$T1S2), "c_S2" = as.numeric(Calphad.globals$T2S2), "d_S2" = as.numeric(Calphad.globals$T3S2),
    #                                    "e_S2" = as.numeric(Calphad.globals$T4S2), "f_S2" = as.numeric(Calphad.globals$T5S2), "g_S2" = as.numeric(Calphad.globals$T6S2))
    # Results <- c(Results_input, Results_output_parameters, Results_output_poly)
    # results_table <- knitr::kable(Results_DB, digits = 2, caption = "Calphad 2.5 model results",
    #                               col.names =c("parameter","value"),"simple")

    DF_Id <- c("CP_100","CP_200","CP_298.15","S_100","S_200","S_298.15")
    DF_input <- c(Calphad.globals$CP100,Calphad.globals$CP200,Calphad.globals$CP298,
                  Calphad.globals$S100,Calphad.globals$S200,Calphad.globals$S298)
    DF_model <- c(Calphad.globals$cp_d100,Calphad.globals$cp_d200,Calphad.globals$cp_d,
                  Calphad.globals$s_d100,Calphad.globals$s_d200,Calphad.globals$s_d)
    results_Cp_S <- data.frame(DF_Id,DF_input,DF_model)

    DF_para_names <-c("Td","a1_T_dep_sol","b1_T_dep_sol", "b_sol_coef", 'Number of atoms')
    DF_para_values <- c(Calphad.globals$Td, Calphad.globals$a1,Calphad.globals$b1,
                        Calphad.globals$b_sol_coef,Calphad.globals$NOA)
    Results_output_parameters <- data.frame(DF_para_names, DF_para_values)

    knitr::kables(
      list(
        knitr::kable(results_Cp_S, col.names = c('Quantity', 'Input value', 'Model Result'), valign = 't',format = '"html"'),
        knitr::kable(Results_output_parameters, col.names = c('Parameter', 'Value'),digits = 5, valign = 't',format = '"html"')
        #, knitr::kable(Results_output_parameters, col.names = c('Parameter', 'Value'),digits = 0, valign = 't')
      ),
      caption = 'Calphad 2.5 model results'
    )

  #}
    #return(Results_DB)
}
