# function for concentration_dependence, 7 parameters. 
library(rootSolve)

pr_prm_Intramol_protein_withDNA <- function(total_protein, deltadeltaF1, deltadeltadimer, deltadeltacoop ,deltadeltaB1, deltadeltaB2, whichOperator){
  ####  for singles, deltadeltaB1=CI binding to DNA, deltadeltaB2= on DNA (operator), 
  ### whichOperator options for singles, 0= not affecting operator sequences, 1= operator 1, 2= operator 2, 3= operator 3. 
  ### for doubles: whichOperator options: 0= All togehter from the protein,1 = pr+OR1, 2=pr+OR2, 3=pr+OR3, 4= OR1+OR2, 5= OR2+OR3, 6= OR1+OR3
  
  #1. define parameter values 
  
  R= 1.98*10^(-3) # kcal/mol
  Temp= 310.15
  Ka= 5 *10^7 # ka for dimmer formation 
  DNA= 10^(-9)
  
  deltadimer= -R*Temp*log(Ka) 
  Ka = exp((-(deltadimer+ deltadeltadimer))/(R*Temp))  
  
  wt_deltaG=  -1 
  deltaG_folding1= wt_deltaG+ deltadeltaF1

  wt_BdeltaGs<- c() 
  if (whichOperator==1) { # or1
    wt_BdeltaGs<- c(-11.7 + deltadeltaB2 + deltadeltaB1, -10.1+ deltadeltaB1, -10.1+ deltadeltaB1)  
  } else if (whichOperator==2) { # or2
    wt_BdeltaGs<- c(-11.7 + deltadeltaB1 , -10.1+ deltadeltaB1+deltadeltaB2, -10.1+ deltadeltaB1)  
  } else if (whichOperator==3) { # or3
    wt_BdeltaGs<- c(-11.7  + deltadeltaB1, -10.1+ deltadeltaB1, -10.1+ deltadeltaB1 + deltadeltaB2)  
  } else if (whichOperator==4){ # or1+ or2
    wt_BdeltaGs<- c(-11.7  + deltadeltaB1 , -10.1+ deltadeltaB2, -10.1)  
  } else if (whichOperator==5) { # or2+ or3
    wt_BdeltaGs<- c(-11.7  , -10.1+ deltadeltaB1, -10.1+ deltadeltaB2)  
  } else if (whichOperator==6) {# or1+ or3
    wt_BdeltaGs<- c(-11.7  + deltadeltaB1 , -10.1, -10.1+ deltadeltaB2)  
  }  else { # protein+ protein
    wt_BdeltaGs<- c(-11.7 + deltadeltaB1+deltadeltaB2 , -10.1+ deltadeltaB1+deltadeltaB2, -10.1+ deltadeltaB1+ deltadeltaB2) 
  } 
  
  
  coop= -2 + deltadeltacoop  # OR1 and OR2 ; or OR2 and OR3
  
  
  #####
  #2. CI-OR thermodynamics model with the updated parameters, empirically search free-dimer vs total concentraiton relationship.  
  
  relations=list(deltaGs=c(0, wt_BdeltaGs[1], wt_BdeltaGs[2], wt_BdeltaGs[3],  wt_BdeltaGs[1]+ wt_BdeltaGs[2]+ coop,  wt_BdeltaGs[1]+ wt_BdeltaGs[3],  wt_BdeltaGs[3]+ wt_BdeltaGs[2]+ coop, sum(wt_BdeltaGs)+ coop),
                 dimer_number=c(0,1,1,1,2,2,2,3)) # thermo_parameters
  
  
  configaration<-c() 
  free_dimer<- c() # free dimer ones 
  for (i in 2:50) { 
    free_dimer[1]<- 10^(-40) 
    free_dimer[i]<- free_dimer[i-1]*5.5
  }  
  Total_protein <- c()
  Propor_sup<- c()
  sum_config<- c()
  r_formd<- c()
  
  library(rootSolve)
  
  for (j in 1:50) { 
    for (i in 1:8) {
      configaration[i]<- exp(-relations[[1]][i]/(R*Temp))*free_dimer[j]^relations[[2]][i]
      sum_config<- sum(configaration)
      config_prob<- configaration/sum_config 
      r_formd[j]<- sum(config_prob[i]*relations[[1]][i])*2*DNA
    }
    Total_protein[j]= (free_dimer[j]/Ka)^0.5+ 2*free_dimer[j]+ r_formd[j]
  }
  wt_trial<- data.frame(Total_protein, free_dimer)
  loess_linear<- loess(wt_trial$free_dimer~ wt_trial$Total_protein, span=0.3) 
  
  fraction_folded1<- exp(-deltaG_folding1/(R*Temp))/(1+exp(-deltaG_folding1/(R*Temp)))
  
  function_protein<- total_protein*fraction_folded1

  free_dimer= predict(loess_linear, function_protein) 
  
  ##
  ### 3.calculate the PR and PRM activity and return 
  configaration<-c() 
  sum_config<- c()
  r_formd<- c() 
  for (i in 1:8) {
    configaration[i]<- exp(-relations[[1]][i]/(R*Temp))*free_dimer^relations[[2]][i]
    sum_config<- sum(configaration)
    config_prob<- configaration/sum_config 
    r_formd<- sum(config_prob[i]*relations[[2]][i])*2*DNA
  }
  
  Propor_sup= config_prob[2]+ config_prob[3]+ config_prob[5]+ config_prob[6]+ config_prob[7]+ config_prob[8]
  Propor_activation=  config_prob[3]+ config_prob[5] # OR2 bound or OR1, OR2 bound
  return(
    list( pr = log2(2^11.761-(2^11.761-23.24125)*Propor_sup), prm= log2(32+ (1024-32)*Propor_activation ))) 
} 
