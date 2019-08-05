#### making data files 
# toy 
source('~/directory/Function_thermodynamics_model.R')
toy<- seq(-0.5, 5, by= 0.1)
toys<- data.frame(expand.grid(toy, toy))
colnames(toys)<- c('delta_deltaA', 'delta_deltaB')
mywts<-c(8.35e-07, 5.52e-08) #  high, low as in the earlier. 
mywts_out<- list(pr= c(), prm=c())

for(j in 1:2){
  for( k in 1:2){
    x<- pr_prm_Intramol_protein_withDNA(mywts[k],0,0,0,0,0,0)
    mywts_out[[j]][k]= x[[j]]
  }
  print(mywts_out[[j]])
}

mysings<- list(
  pr= list(
    folding = list(h= data.frame(), l= data.frame()), 
    dimerization= list(h= data.frame(), l= data.frame()), 
    tetramerization= list(h= data.frame(), l= data.frame()), 
    binding to DNA = list(h= data.frame(), l= data.frame()), 
    OR1= list(h= data.frame(), l= data.frame()), 
    OR2= list(h= data.frame(), l= data.frame()), 
    OR3= list(h= data.frame(), l= data.frame())
  ), 
  prm= list(
    folding = list(h= data.frame(), l= data.frame()), 
    dimerization= list(h= data.frame(), l= data.frame()), 
    tetramerization= list(h= data.frame(), l= data.frame()), 
    binding to DNA = list(h= data.frame(), l= data.frame()), 
    OR1= list(h= data.frame(), l= data.frame()), 
    OR2= list(h= data.frame(), l= data.frame()), 
    OR3= list(h= data.frame(), l= data.frame())
  )
)

for(i in 1:2){
  for(j in 1:7){
    for (k in 1:2) {
      mysings[[i]][[j]][[k]] = data.frame(toy)
      for (l in 1:length(toy)){ 
        tmp<- list(
          x_f= pr_prm_Intramol_protein_withDNA(mywts[k], toy[l],0,0,0,0,0) ,
          x_dim= pr_prm_Intramol_protein_withDNA(mywts[k], 0, toy[l],0,0,0,0) ,
          x_tetra= pr_prm_Intramol_protein_withDNA(mywts[k], 0, 0, toy[l],0,0,0) ,
          x_b= pr_prm_Intramol_protein_withDNA(mywts[k], 0, 0,0,toy[l],0,0),
          x_OR1 = pr_prm_Intramol_protein_withDNA(mywts[k], 0, 0,0,0, toy[l],1),
          x_OR2 = pr_prm_Intramol_protein_withDNA(mywts[k], 0, 0,0,0, toy[l],2),  
          x_OR3 = pr_prm_Intramol_protein_withDNA(mywts[k], 0, 0,0,0, toy[l],3)
        )
        mysings[[i]][[j]][[k]][l, 'Output'] = tmp[[j]][[i]]
      }
      colnames(mysings[[i]][[j]][[k]])[1] = 'delta_deltaG'
    }
  }
}



for (i in 1:2){
  for(j in 1:7){
    for(k in 1:2) {
      mysings[[i]][[j]][[k]]$amount= names(mysings[[i]][[j]])[k]
      mysings[[i]][[j]][[k]]$mutation= names(mysings[[i]])[j]
    }
  }
}
sing_in_rows<- list(
  pr= rbind(mysings[[1]][[1]][[1]], mysings[[1]][[1]][[2]]), 
  prm= rbind(mysings[[2]][[1]][[1]], mysings[[2]][[1]][[2]])
)

for (i in 1:2) {
  for(j in 1:7){
    sing_in_rows[[i]] = rbind(sing_in_rows[[i]], mysings[[i]][[j]][[1]],  mysings[[i]][[j]][[2]])
  } 
  
  sing_in_rows[[i]]$mutation= factor(sing_in_rows[[i]]$mutation, 
                                     levels = c('folding', 'dimerization', 'binding to DNA','tetramerization',  'OR1', 'OR2', 'OR3') )
  
  sing_in_rows[[i]]$amount= factor(sing_in_rows[[i]]$amount, 
                                   levels = c('l', 'h') )
}

xx<- list(
  pr= list(), 
  prm=list()
)
for ( i in 1: 2) { 
  for (k in 1:7) {
    h= mysings_ci[[1]][[i]][[k]][[1]] 
    colnames(h)[2]= 'High'
    l= mysings_ci[[1]][[i]][[k]][[2]] 
    colnames(l)[2]= 'Low'
    xx[[i]][[k]]=  merge(l[, c(1,2,4)], h[,c(1,2)])  
  }
}
hl_comp<- list(
)

for (i in 1:2){ 
  hl_comp[[i]]= rbind(xx[[i]][[1]], xx[[i]][[2]],xx[[i]][[3]],xx[[i]][[4]],xx[[i]][[5]],xx[[i]][[6]], xx[[i]][[7]])
}
names(hl_comp)<- c('pr', 'prm')
mysings_ci[[3]]<- hl_comp
names(mysings_ci)[3] = 'hl_comp'
for (i in 1:2){
  mysings_ci[[3]][[i]]$mutation = factor(mysings_ci[[3]][[i]]$mutation, levels=levels(mysings_ci[[2]][[1]]$mutation))
}

save(mysings_ci, file='concentr_dep_mut_effect_epis/ci_sings_7_params_list. RData')
########
######## making doubles 

mydoubs<- list(
  pr= list(
    fold_fold = list(h= data.frame(), l= data.frame()), 
    dimer_dimer= list(h= data.frame(), l= data.frame()), 
    tetramer_tetramer= list(h= data.frame(), l= data.frame()), 
    bind_bind = list(h= data.frame(), l= data.frame()), 
    OR1_OR1= list(h= data.frame(), l= data.frame()), 
    OR2_OR2= list(h= data.frame(), l= data.frame()),
    OR3_OR3= list(
      h= data.frame(), l= data.frame()),
    fold_OR3 = list(
      h= data.frame(), l= data.frame()), 
    dimer_OR3 = list(
      h= data.frame(), l= data.frame()), 
    tetramer_OR3 = list(
      h= data.frame(), l= data.frame()),
    bind_OR3 = list(
      h= data.frame(), l= data.frame()), 
    
    fold_dimer = list(h= data.frame(), l= data.frame()), 
    fold_tetramer= list(h= data.frame(), l= data.frame()), 
    fold_bind = list(h= data.frame(), l= data.frame()), 
    fold_OR1 = list(h= data.frame(), l= data.frame()), 
    fold_OR2 = list(h= data.frame(), l= data.frame()), 
    
    dimer_tetramer= list(h= data.frame(), l= data.frame()), 
    dimer_bind = list(h= data.frame(), l= data.frame()), 
    dimer_OR1 = list(h= data.frame(), l= data.frame()), 
    dimer_OR2 = list(h= data.frame(), l= data.frame()), 
    
    tetramer_bind = list(h= data.frame(), l= data.frame()), 
    tetramer_OR1 = list(h= data.frame(), l= data.frame()), 
    tetramer_OR2 = list(h= data.frame(), l= data.frame()),
    
    bind_OR1 = list(h= data.frame(), l= data.frame()), 
    bind_OR2 = list(h= data.frame(), l= data.frame()) ,
    OR1_OR2= list(h= data.frame(), l= data.frame()), 
    OR2_OR3= list(h= data.frame(), l= data.frame()),
    OR1_OR3= list(h= data.frame(), l= data.frame())
    
    # don't do OR1 and OR2, since still we are interested involving protein. 
    
  ), 
  
  prm= list(
    fold_fold = list(h= data.frame(), l= data.frame()),  
    dimer_dimer= list(h= data.frame(), l= data.frame()), 
    tetramer_tetramer= list(h= data.frame(), l= data.frame()), 
    bind_bind = list(h= data.frame(), l= data.frame()), 
    OR1_OR1= list(h= data.frame(), l= data.frame()), 
    OR2_OR2= list(h= data.frame(), l= data.frame()),
    OR3_OR3= list(
      h= data.frame(), l= data.frame()),
    fold_OR3 = list(
      h= data.frame(), l= data.frame()), 
    dimer_OR3 = list(
      h= data.frame(), l= data.frame()), 
    tetramer_OR3 = list(
      h= data.frame(), l= data.frame()),
    bind_OR3 = list(
      h= data.frame(), l= data.frame()), 
    
    fold_dimer = list(h= data.frame(), l= data.frame()), 
    fold_tetramer= list(h= data.frame(), l= data.frame()), 
    fold_bind = list(h= data.frame(), l= data.frame()), 
    fold_OR1 = list(h= data.frame(), l= data.frame()), 
    fold_OR2 = list(h= data.frame(), l= data.frame()), 
    
    dimer_tetramer= list(h= data.frame(), l= data.frame()), 
    dimer_bind = list(h= data.frame(), l= data.frame()), 
    dimer_OR1 = list(h= data.frame(), l= data.frame()), 
    dimer_OR2 = list(h= data.frame(), l= data.frame()), 
    
    tetramer_bind = list(h= data.frame(), l= data.frame()), 
    tetramer_OR1 = list(h= data.frame(), l= data.frame()), 
    tetramer_OR2 = list(h= data.frame(), l= data.frame()),
    
    bind_OR1 = list(h= data.frame(), l= data.frame()), 
    bind_OR2 = list(h= data.frame(), l= data.frame()),
    OR1_OR2= list(h= data.frame(), l= data.frame()), 
    OR2_OR3= list(h= data.frame(), l= data.frame()),
    OR1_OR3= list(h= data.frame(), l= data.frame())
    
  )
)
for(i in 1:2){
  for(j in 1:28){
    for (k in 1:2) {
      mydoubs[[i]][[j]][[k]] = data.frame(toys)
      for (l in 1:length(rownames(toys))){ 
        tmp<- list(
          f_f= pr_prm_Intramol_protein_withDNA(mywts[k], toys[l,1]+toys[l,2] ,0,0,0,0,0) ,
          dim_dim= pr_prm_Intramol_protein_withDNA(mywts[k], 0, toys[l,1]+toys[l,2],0,0,0,0) ,
          tetra_tetra= pr_prm_Intramol_protein_withDNA(mywts[k], 0, 0, toys[l,1]+toys[l,2],0,0,0) ,
          b_b= pr_prm_Intramol_protein_withDNA(mywts[k], 0, 0,0,toys[l,1]+toys[l,2],0,0),
          OR1_OR1 = pr_prm_Intramol_protein_withDNA(mywts[k], 0, 0,0,0, toys[l,1]+toys[l,2],1),
          OR2_OR2 = pr_prm_Intramol_protein_withDNA(mywts[k], 0, 0,0,0, toys[l,1]+toys[l,2],2),  
          
          OR3_OR3 = pr_prm_Intramol_protein_withDNA(mywts[k], 0, 0,0,0, toys[l,1]+toys[l,2],3), 
          f_OR3 = pr_prm_Intramol_protein_withDNA(mywts[k], toys[l,1], 0,0,0, toys[l,2],3), 
          dim_OR3 = pr_prm_Intramol_protein_withDNA(mywts[k], 0, toys[l,1], 0,0, toys[l,2],3),  
          tetra_OR3 = pr_prm_Intramol_protein_withDNA(mywts[k], 0, 0, toys[l,1], 0, toys[l,2],3), 
          b_OR3 = pr_prm_Intramol_protein_withDNA(mywts[k], 0, 0, 0, toys[l,1], toys[l,2],3), 
          
          f_dim= pr_prm_Intramol_protein_withDNA(mywts[k], toys[l,1] ,toys[l,2] ,0,0,0,0) ,
          f_tetra= pr_prm_Intramol_protein_withDNA(mywts[k], toys[l,1], 0, toys[l,2],0,0,0) ,
          f_b= pr_prm_Intramol_protein_withDNA(mywts[k], toys[l,1], 0,0,toys[l,2],0,0),
          f_OR1 = pr_prm_Intramol_protein_withDNA(mywts[k], toys[l,1] , 0,0,0, toys[l,2],1),
          f_OR2 = pr_prm_Intramol_protein_withDNA(mywts[k], toys[l,1], 0,0,0, toys[l,2],2), 
          
          dim_tetra= pr_prm_Intramol_protein_withDNA(mywts[k],  0,toys[l,1], toys[l,2],0,0,0) ,
          dim_b= pr_prm_Intramol_protein_withDNA(mywts[k], 0,toys[l,1], 0,toys[l,2],0,0),
          dim_OR1 = pr_prm_Intramol_protein_withDNA(mywts[k], 0, toys[l,1] , 0,0, toys[l,2],1),
          dim_OR2 = pr_prm_Intramol_protein_withDNA(mywts[k], 0, toys[l,1], 0,0, toys[l,2],2),  
          
          tetra_b= pr_prm_Intramol_protein_withDNA(mywts[k], 0,0,toys[l,1], toys[l,2],0,0),
          tetra_OR1 = pr_prm_Intramol_protein_withDNA(mywts[k], 0, 0, toys[l,1] ,0, toys[l,2],1),
          tetra_OR2 = pr_prm_Intramol_protein_withDNA(mywts[k], 0, 0, toys[l,1], 0, toys[l,2],2), 
          
          b_OR1 = pr_prm_Intramol_protein_withDNA(mywts[k], 0, 0, 0,toys[l,1] , toys[l,2],1),
          b_OR2 = pr_prm_Intramol_protein_withDNA(mywts[k], 0, 0, 0, toys[l,1], toys[l,2],2),
          
          OR1_OR2 = pr_prm_Intramol_protein_withDNA(mywts[k], 0, 0, 0, tmp[[i]][[j]][[k]][l,1], tmp[[i]][[j]][[k]][l,2],4),
          OR2_OR3 = pr_prm_Intramol_protein_withDNA(mywts[k], 0, 0, 0, tmp[[i]][[j]][[k]][l,1], tmp[[i]][[j]][[k]][l,2],5),
          OR1_OR3 = pr_prm_Intramol_protein_withDNA(mywts[k], 0, 0, 0, tmp[[i]][[j]][[k]][l,1], tmp[[i]][[j]][[k]][l,2],6)
        )
        mydoubs[[i]][[j]][[k]][l, 'thermodynamic'] = tmp[[j]][[i]]
      }
    }
  }
}
for(i in 1:2){
  for(j in 1:25){
    for (k in 1:2) {
      mydoubs[[i]][[j]][[k]]$amount = names(mydoubs[[i]][[j]])[k] 
      mydoubs[[i]][[j]][[k]]$mutation = names(mydoubs[[i]])[j]
      mydoubs[[i]][[j]][[k]]$promoter = names(mydoubs)[i]
    }
  }
}

library(stringr)
whichwhich<- str_split(names(mydoubs[[1]]), pattern = "_")

for(i in 1:2){
  for(j in 1:25){
    for (k in 1:2) {
      x1<- mysings_ci[[1]][[i]][[whichwhich[[j]][1]]][[k]]$Output
      x2<- mysings_ci[[1]][[i]][[whichwhich[[j]][2]]][[k]]$Output
      xx<- expand.grid(x1, x2)
      colnames(xx) = c('path1', 'path2')
      mydoubs[[i]][[j]][[k]] = cbind(mydoubs[[i]][[j]][[k]], xx)
    }
  }
}

for (i in 1:2){
  for(j in 1:28){
    for (k in 1:2) {
      
      mydoubs[[i]][[j]][[k]]$additive =  mydoubs[[2]][[j]][[k]]$path1 +   mydoubs[[2]][[j]][[k]]$path2 - mywts_out[[2]][k]
      mydoubs[[i]][[j]][[k]]$epis =  mydoubs[[2]][[j]][[k]]$obs -  mydoubs[[2]][[j]][[k]]$exp
    }
  }
}

specify_decimal <- function(x, k) trimws(format(round(x, k), nsmall=k))

for (i in 1:2){
  for(j in 1:28){
    for(k in 1:2){
      mydoubs[[i]][[j]][[k]]$epis_deci= as.numeric(specify_decimal(mydoubs[[i]][[j]][[k]]$epis, 2))
      mydoubs[[i]][[j]][[k]]$thermo_deci= as.numeric(specify_decimal( mydoubs[[i]][[j]][[k]]$thermodynamic, 2))
      mydoubs[[i]][[j]][[k]]$path1_deci= as.numeric(specify_decimal( mydoubs[[i]][[j]][[k]]$path1, 2))
      mydoubs[[i]][[j]][[k]]$path2_deci= as.numeric(specify_decimal( mydoubs[[i]][[j]][[k]]$path2, 2))
      
      mydoubs[[i]][[j]][[k]][((mydoubs[[i]][[j]][[k]]$thermo_deci== mydoubs[[i]][[j]][[k]]$path1_deci)| 
                                       (doub_list[[1]][[i]][[j]][[k]]$thermo_deci== mydoubs[[i]][[j]][[k]]$path2_deci)) & 
                                      doub_list[[1]][[i]][[j]][[k]]$epis_deci!=0, 
                                    "epis_cla" ] = 'masking'
      
      mydoubs[[i]][[j]][[k]][doub_list[[1]][[i]][[j]][[k]]$epis_deci==0,'epis_cla']= 'none'  
      
      
    }
  }
}

for (i in 1:2){
  for(j in 1:28){
    for(k in 1:2){
      for (l in 1:length(rownames(mydoubs[[i]][[j]][[k]]))){
        
        if ((mydoubs[[i]][[j]][[k]][l, 'path1']> mywts_out[[i]][k] & mydoubs[[i]][[j]][[k]][l, 'path2']> mywts_out[[i]][k] & 
             mydoubs[[i]][[j]][[k]][l, 'thermo_deci'] > min(mydoubs[[i]][[j]][[k]][l, 'path1_deci'], mydoubs[[i]][[j]][[k]][l, 'path2_deci']) & 
             mydoubs[[i]][[j]][[k]][l, 'thermo_deci'] < max(mydoubs[[i]][[j]][[k]][l, 'path1_deci'], mydoubs[[i]][[j]][[k]][l, 'path2_deci']))|
            
            (mydoubs[[i]][[j]][[k]][l, 'path1']< mywts_out[[i]][k] & mydoubs[[i]][[j]][[k]][l, 'path2']< mywts_out[[i]][k] & 
             mydoubs[[i]][[j]][[k]][l, 'thermo_deci'] > min(mydoubs[[i]][[j]][[k]][l, 'path1_deci'], mydoubs[[i]][[j]][[k]][l, 'path2_deci']) & 
             mydoubs[[i]][[j]][[k]][l, 'thermo_deci'] < max(mydoubs[[i]][[j]][[k]][l, 'path1_deci'], mydoubs[[i]][[j]][[k]][l, 'path2_deci']))) {
          mydoubs[[i]][[j]][[k]][l, 'epis_cla'] = 'sign' 
        } else if (
          (mydoubs[[i]][[j]][[k]][l, 'path1']> mywts_out[[i]][k] & mydoubs[[i]][[j]][[k]][l, 'path2']> mywts_out[[i]][k] & 
           mydoubs[[i]][[j]][[k]][l, 'thermo_deci'] < min(mydoubs[[i]][[j]][[k]][l, 'path1_deci'], mydoubs[[i]][[j]][[k]][l, 'path2_deci']) )|
          
          (mydoubs[[i]][[j]][[k]][l, 'path1']< mywts_out[[i]][k] & mydoubs[[i]][[j]][[k]][l, 'path2']< mywts_out[[i]][k] & 
           mydoubs[[i]][[j]][[k]][l, 'thermo_deci'] > max(mydoubs[[i]][[j]][[k]][l, 'path1_deci'], mydoubs[[i]][[j]][[k]][l, 'path2_deci'])) | 
          (min(mydoubs[[i]][[j]][[k]][l, 'path1_deci'], mydoubs[[i]][[j]][[k]][l, 'path2_deci'])< mywts_out[[i]][k] & 
           max(mydoubs[[i]][[j]][[k]][l, 'path1_deci'], mydoubs[[i]][[j]][[k]][l, 'path2_deci'])> mywts_out[[i]][k] & 
           (mydoubs[[i]][[j]][[k]][l, 'thermo_deci'] > max(mydoubs[[i]][[j]][[k]][l, 'path1_deci'], mydoubs[[i]][[j]][[k]][l, 'path2_deci'])| 
            mydoubs[[i]][[j]][[k]][l, 'thermo_deci'] < min (mydoubs[[i]][[j]][[k]][l, 'path1_deci'], mydoubs[[i]][[j]][[k]][l, 'path2_deci'])) )) {
          mydoubs[[i]][[j]][[k]][l, 'epis_cla'] = 'reciprocal sign' 
        }
        else if (mydoubs[[i]][[j]][[k]][l, 'epis_deci']!= 0 & mydoubs[[i]][[j]][[k]][l, 'epis_cla']!='reciprocal sign' & 
                 mydoubs[[i]][[j]][[k]][l, 'epis_cla']!='sign' & mydoubs[[i]][[j]][[k]][l, 'epis_cla']!='masking' ) { 
          mydoubs[[i]][[j]][[k]][l, 'epis_cla'] = 'magnitude' 
        }
      }
    }
  }
}

mydoub_hl<- list( 
  pr= list(
    thermodynamic= list(), 
    additive= list(), 
    epis= list() 
  ), 
  prm=list(
    thermodynamic= list(), 
    additive= list(), 
    epis= list()

  )
)

for(i in 1:2){
  for(j in 1:28){
    mydoub_hl[[i]][[1]][[j]] = data.frame(High=  mydoubs[[i]][[j]][[1]]$thermodynamic, Low= mydoubs[[i]][[j]][[2]]$thermodynamic )
    mydoub_hl[[i]][[2]][[j]] = data.frame(High=  mydoubs[[i]][[j]][[1]]$additive, Low= mydoubs[[i]][[j]][[2]]$additive )
    mydoub_hl[[i]][[3]][[j]] = data.frame(High=  mydoubs[[i]][[j]][[1]]$epis, Low= mydoubs[[i]][[j]][[2]]$epis )

  }
}
for(i in 1:2) {
  for(k in 1:3) {
    names(mydoub_hl[[i]][[k]]) = names(mydoubs[[1]])
  }
}

for(i in 1:2) {
  for(k in 1:3) {
    for (l in 1:28) {
      mydoub_hl[[i]][[k]][[l]]$combination= names(mydoub_hl[[i]][[k]])[l]
      mydoub_hl[[i]][[k]][[l]]$whichmeasure = names(mydoub_hl[[i]])[k]
      print(head(mydoub_hl[[i]][[k]][[l]]))
    }
  }
}
mydoub_inrows <- list(
  pr= rbind(mydoub_hl[[1]][[1]][[1]], mydoub_hl[[1]][[2]][[1]],mydoub_hl[[1]][[3]][[1]] ), 
  prm=rbind(mydoub_hl[[2]][[1]][[1]], mydoub_hl[[2]][[2]][[1]],mydoub_hl[[2]][[3]][[1]])
)
for(i in 1:2) {
  for (j in 2:28) {
    mydoub_inrows[[i]]= rbind(mydoub_inrows[[i]],mydoub_hl[[i]][[1]][[j]], mydoub_hl[[i]][[2]][[j]],mydoub_hl[[i]][[3]][[j]])
  }
}
for(i in 1:2) {
  
  mydoub_inrows[[i]]$combination= factor(
    mydoub_inrows[[i]]$combination, levels= 
      c("fold_fold" ,   "dimer_dimer" ,   "bind_bind" ,   "tetramer_tetramer" ,        "OR1_OR1"     ,      "OR2_OR2"      ,'OR3_OR3', 
        "fold_dimer"  ,      "fold_tetramer"   , "fold_bind" ,"fold_OR1"   ,       "fold_OR2"    ,   'fold_OR3', 
        "dimer_tetramer" ,  "dimer_bind",  "tetramer_bind"   , 'tetramer_OR3', 
        "dimer_OR1"   ,      "dimer_OR2"    ,  'dimer_OR3',    
        "tetramer_OR1"  ,    "tetramer_OR2"  ,    "bind_OR1"    ,      "bind_OR2"    , 'bind_OR3' , 'OR1_OR2','OR2_OR3','OR1_OR3'
      )
  )
  
}
doub_list<- list(indepdently=mydoubs, high_low= mydoub_hl, inrows= mydoub_inrows)

doub_list[[4]]<- list(
  pr=data.frame(), 
  prm= data.frame())

for(i in 1:2){
  doub_list[[4]][[i]]= doub_list[[3]][[i]][, c(1,3:8)]
  colnames(doub_list[[4]][[i]])[1] = 'val' 
  doub_list[[4]][[i]]$amount = 'High'
  
  x= doub_list[[3]][[i]][, c(2:8)]
  colnames(x)[1] = 'val' 
  x$amount = 'Low'
  
  doub_list[[4]][[i]]= rbind(doub_list[[4]][[i]], x)
}

names(doub_list)[4] = 'all_rows_for_tile_plot'
for (i in 1:2){
  doub_list[[4]][[i]]$whichmeasure = factor(doub_list[[4]][[i]]$whichmeasure, levels= c('additive', 'thermodynamic', 'epis'))
  doub_list[[4]][[i]]$amount = factor(doub_list[[4]][[i]]$amount, levels= c('Low', 'High'))
}
# epis only, classes. 
epis<- list(
  pr= doub_list[[3]][[1]][doub_list[[3]][[1]]$whichmeasure=='epis',],
  prm= doub_list[[3]][[2]][doub_list[[3]][[2]]$whichmeasure=='epis',]
)
specify_decimal <- function(x, k) trimws(format(round(x, k), nsmall=k))

for (i in 1:2){
  epis[[i]]$H_min_L = as.numeric(epis[[i]]$High) - as.numeric(epis[[i]]$Low)
  epis[[i]]$H_deci =as.numeric(specify_decimal(as.numeric(epis[[i]]$High), 1))
  epis[[i]]$L_deci =as.numeric(specify_decimal(as.numeric(epis[[i]]$Low), 1))
  
  epis[[i]]$LH_class = 'same' 
  epis[[i]][(epis[[i]]$H_deci==0 & epis[[i]]$L_deci!=0)| (epis[[i]]$H_deci!=0 & epis[[i]]$L_deci==0) ,'LH_class'] = 'one_0' 
  
  epis[[i]][(epis[[i]]$H_deci>0 & epis[[i]]$L_deci<0) ,'LH_class'] = 'L<0,H>0' 
  epis[[i]][(epis[[i]]$H_deci<0 & epis[[i]]$L_deci>0) ,'LH_class'] = 'L>0,H<0' 
  
}
for (i in 1:2){
  epis[[i]]$LH_class = factor(epis[[i]]$LH_class, levels= c('same' , 'one_0' , 'L<0,H>0', 'L>0,H<0'))
}

doub_list[[5]]<- epis
names(doub_list)[5] = 'epis_only_for_h_l'
tmp<- doub_list[[3]]
for(i in 1:2){
  names(tmp[[i]])= c("High" ,  "Low", "combination",  "whichmeasure" ,"path2", "path1" , "delta_deltaB", "delta_deltaA", 'p2','p1','do' 
  )
}

duplicated_rows<- list(
  pr= rbind(doub_list[[3]][[1]], tmp[[1]]), 
  prm= rbind(doub_list[[3]][[2]], tmp[[2]])
)
doub_list[[6]]<- duplicated_rows
names(doub_list)[6] = 'duplicated_rows_for_inrow'
doub_list[[7]]<-  list()

names(doub_list)[7]= 'low_high_exp_obs_epis_togehter'

for (i in 1:2) {
  
  doub_list[[7]][[i]] = rbind(doub_list[[1]][[i]][[1]][[1]], doub_list[[1]][[i]][[1]][[2]])
  
  for (l in 2:28){
    doub_list[[7]][[i]] = rbind( doub_list[[7]][[i]],doub_list[[1]][[i]][[l]][[1]], doub_list[[1]][[i]][[l]][[2]])
  }
}
names(doub_list[[7]])= c('pr', 'prm')
for (i in 1:2){
  doub_list[[7]][[i]][doub_list[[7]][[i]]$amount=='l','amount'] = 'Low'
  doub_list[[7]][[i]][doub_list[[7]][[i]]$amount=='h','amount'] = 'High' 
  
  doub_list[[7]][[i]]$amount = factor(doub_list[[7]][[i]]$amount, levels= c('Low', 'High')) 
  
  doub_list[[7]][[i]][doub_list[[7]][[i]]$mutation=='fold_fold'  , 'mutation'] = 'folding' 
  doub_list[[7]][[i]][doub_list[[7]][[i]]$mutation=='dimer_dimer'  , 'mutation'] = 'dimerization' 
  doub_list[[7]][[i]][doub_list[[7]][[i]]$mutation=='bind_bind'  , 'mutation'] = 'binding to DNA' 
  doub_list[[7]][[i]][doub_list[[7]][[i]]$mutation=='tetramer_tetramer'  , 'mutation'] = 'tetramerization' 
  doub_list[[7]][[i]][doub_list[[7]][[i]]$mutation=='OR1_OR1'  , 'mutation'] = 'OR1' 
  doub_list[[7]][[i]][doub_list[[7]][[i]]$mutation=='OR2_OR2'  , 'mutation'] = 'OR2' 
  doub_list[[7]][[i]][doub_list[[7]][[i]]$mutation=='OR3_OR3'  , 'mutation'] = 'OR3' 
  doub_list[[7]][[i]][doub_list[[7]][[i]]$mutation=='fold_dimer'  , 'mutation'] = 'F+D' 
  doub_list[[7]][[i]][doub_list[[7]][[i]]$mutation=='fold_tetramer'  , 'mutation'] = 'F+T' 
  doub_list[[7]][[i]][doub_list[[7]][[i]]$mutation=='fold_bind'  , 'mutation'] = 'F+B' 
  doub_list[[7]][[i]][doub_list[[7]][[i]]$mutation=='fold_OR1'  , 'mutation'] = 'F+OR1' 
  doub_list[[7]][[i]][doub_list[[7]][[i]]$mutation=='fold_OR2'  , 'mutation'] = 'F+OR2' 
  doub_list[[7]][[i]][doub_list[[7]][[i]]$mutation=='fold_OR3'  , 'mutation'] = 'F+OR3' 
  doub_list[[7]][[i]][doub_list[[7]][[i]]$mutation=='dimer_tetramer'  , 'mutation'] = 'D+T' 
  doub_list[[7]][[i]][doub_list[[7]][[i]]$mutation=='dimer_bind'  , 'mutation'] = 'D+B' 
  doub_list[[7]][[i]][doub_list[[7]][[i]]$mutation=='tetramer_bind'  , 'mutation'] = 'T+B' 
  doub_list[[7]][[i]][doub_list[[7]][[i]]$mutation=='tetramer_OR3'  , 'mutation'] = 'T+OR3' 
  doub_list[[7]][[i]][doub_list[[7]][[i]]$mutation=='dimer_OR1'  , 'mutation'] = 'D+OR1' 
  doub_list[[7]][[i]][doub_list[[7]][[i]]$mutation=='dimer_OR2'  , 'mutation'] = 'D+OR2' 
  doub_list[[7]][[i]][doub_list[[7]][[i]]$mutation=='dimer_OR3'  , 'mutation'] = 'D+OR3' 
  doub_list[[7]][[i]][doub_list[[7]][[i]]$mutation=='tetramer_OR1'  , 'mutation'] = 'T+OR1' 
  doub_list[[7]][[i]][doub_list[[7]][[i]]$mutation=='tetramer_OR2'  , 'mutation'] = 'T+OR2' 
  doub_list[[7]][[i]][doub_list[[7]][[i]]$mutation=='bind_OR1'  , 'mutation'] = 'B+OR1' 
  doub_list[[7]][[i]][doub_list[[7]][[i]]$mutation=='bind_OR2'  , 'mutation'] = 'B+OR2' 
  doub_list[[7]][[i]][doub_list[[7]][[i]]$mutation=='bind_OR3'  , 'mutation'] = 'B+OR3' 
  doub_list[[7]][[i]][doub_list[[7]][[i]]$mutation=='OR1_OR2'  , 'mutation'] = 'OR1+OR2' 
  doub_list[[7]][[i]][doub_list[[7]][[i]]$mutation=='OR2_OR3'  , 'mutation'] = 'OR2+OR3' 
  doub_list[[7]][[i]][doub_list[[7]][[i]]$mutation=='OR1_OR3'  , 'mutation'] = 'OR1+OR3' 
  
  doub_list[[7]][[i]]$mutation = factor(doub_list[[7]][[i]]$mutation , 
                                        levels= c('folding', 'dimerization', 'binding to DNA','tetramerization','OR1','OR2','OR3', 
                                                  'F+D', 'F+B', 'F+OR1', 'D+B', 'D+OR1', 'B+OR1', 'OR1+OR2',
                                                  'F+OR2', 'D+OR2', 'B+OR2', 'F+OR3', 'D+OR3', 'B+OR3', 'OR2+OR3',
                                                  'OR1+OR3','F+T', 'D+T','T+B', 'T+OR1','T+OR2', 'T+OR3')) 
  
}
for_distribution<- list()
mypath<- list()
for(i in 1:2) {
  pathH= doub_list[[1]][[i]][[1]][[1]][,c("path1" , "path2" )]
  pathL= doub_list[[1]][[i]][[1]][[2]][,c("path1" , "path2" )]
  
  for (j in 2:28) {
    
    pathH= rbind(pathH, doub_list[[1]][[i]][[j]][[1]][,c("path1" , "path2" )])
    pathL= rbind(pathL, doub_list[[1]][[i]][[j]][[2]][,c("path1" , "path2" )])
  }
  
  names(pathH) = c('H:path1', 'H:paht2') 
  names(pathL) = c('L:path1', 'L:paht2') 
  
  mypath[[i]]= cbind(pathH, pathL)
  
}

for(i in 1:2){
  for_distribution[[i]]= doub_list[[3]][[i]][doub_list[[3]][[i]]$whichmeasure=='thermodynamic', ]
   for_distribution[[i]]= cbind(for_distribution[[i]], mypath[[i]])
  
  names(for_distribution[[i]])[12:15] = c("H_path1" ,     "H_path_2",      "L_path1",     
                                          "L_path2"     )
}

dens<- list()
for (i in 1:2){
  h12<- for_distribution[[i]][,c('combination', 'path1',   'path2', 'delta_deltaA', 'delta_deltaB', 
                                 'High')]
  colnames(h12)[6] = 'val'
  h12$whichmeasure = 'do'
  h1<- for_distribution[[i]][,c('combination', 'path1',   'path2', 'delta_deltaA', 'delta_deltaB', 
                                'H_path1')]
  colnames(h1)[6] = 'val'
  h1$whichmeasure = 'p1'
  h2<- for_distribution[[i]][,c('combination', 'path1',   'path2', 'delta_deltaA', 'delta_deltaB', 
                                'H_path_2')]
  h2$whichmeasure = 'p2'
  
  colnames(h2)[6] = 'val'
  
  h<- rbind(h12, h1, h2) 
  h$amount = 'High'
  
  l12<- for_distribution[[i]][,c('combination', 'path1',   'path2', 'delta_deltaA', 'delta_deltaB', 
                                 'Low')]
  colnames(l12)[6] = 'val'
  l12$whichmeasure = 'do'
  
  l1<- for_distribution[[i]][,c('combination', 'path1',   'path2', 'delta_deltaA', 'delta_deltaB', 
                                'L_path1')]
  colnames(l1)[6] = 'val'
  l1$whichmeasure = 'p1'
  l2<- for_distribution[[i]][,c('combination', 'path1',   'path2', 'delta_deltaA', 'delta_deltaB', 
                                'L_path2')]
  colnames(l2)[6] = 'val'
  l2$whichmeasure = 'p2'
  l<- rbind(l12, l1, l2) 
  l$amount = 'Low'
  
  dens[[i]]= rbind(h,l)
  
}
for (i in 1:2){
  dens[[i]]$whichmeasure== factor(dens[[i]]$whichmeasure, levels= c('do', 'p1','p2'))
}
doub_list[[8]]=dens
names(doub_list)[8] = 'for_density_of_dme'

save(doub_list, file='~/path/ci_doubs_7_params_list28_comb.RData')



