## plotting 

library(ggplot2)
library(ggpubr)

# The palette with grey:
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") # first gray, second near yellow, third near blue
# The palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
## mypalettes
mygreenPalette<- c('#C3F780', '#8BC34A', '#76A63F', '#648C35', '#486627', '#364D1D', '#243313')
mygreyPalette<- c('#EEEEEE', '#DDDDDD', '#CCCCCC', '#BBBBBB', '#AAAAAA', '#999999', '#777777', '#555555', '#333333','#111111')
myredPalette<- c('#FFEBEE', '#FFCDD2', '#EF9A9A', '#E57373', '#EF5350', '#F44336', '#E53935', '#D32F2F', '#C62828', '#B71C1C')

grey_to_green<- c( '#243313','#364D1D','#486627','#648C35','#76A63F','#DDDDDD')

#####
##### 
load('~/directory/mysings_ci.RData')
load('~/directory/doub_list.RData')

########
## 1. energy vs output 
ggplot(data=mysings_ci[[2]][[1]]) + geom_line(aes(x= delta_deltaG, y= Output, col= amount)) + 
  facet_grid(.~ mutation) + scale_color_manual(values= c('blue', 'red')) + 
  theme_classic() + geom_vline(xintercept = 0, linetype= 2) + labs(y= 'log2(PR.AU)')

ggplot(data=mysings_ci[[2]][[2]]) + geom_line(aes(x= delta_deltaG, y= Output, col= amount)) + 
  facet_grid(.~ mutation) + scale_color_manual(values= c('blue', 'red')) + 
  theme_classic() + geom_vline(xintercept = 0, linetype= 2) + labs(y= 'log2(PRM.AU)')
#1.1, tetramerization & or2 looks very similar. Plotting two 
ggplot() + geom_line(data=sing_in_rows[[1]][sing_in_rows[[1]]$mutation=='tetramer', ], aes(x= delta_deltaG, y= Output, col= amount)) + 
  geom_line(data=sing_in_rows[[1]][sing_in_rows[[1]]$mutation=='OR2', ], aes(x= delta_deltaG, y= Output, col= amount), lty=2) +   scale_color_manual(values= c('blue', 'red')) + 
  theme_classic() + geom_vline(xintercept = 0, linetype= 2) + labs(y= 'log2(PR.AU)')

#####
## 2. High vs low 


#2.1, comparing phenotypes 
ggplot(data=mysings_ci[[3]][[1]]) + 
  geom_point(aes(x=as.numeric(Low), y=as.numeric(High), col= delta_deltaG), size=2, shape=15)+ 
  theme_classic()+ 
  geom_abline(intercept = 0, slope=1, lty=2, col='gray') + 
  geom_point(aes(x= 7.069273, y=4.564443), col='red', shape=4) + 
  scale_color_gradientn(colors = mygreyPalette) + facet_grid(.~mutation)+ 
  labs(y= 'High: PR (log2 AU)', x= 'Low: PR (log2 AU)') +
  scale_y_continuous(breaks=seq(0,18, 3)) + scale_x_continuous(breaks=seq(0,18,3))  

ggplot(data=mysings_ci[[3]][[2]]) + 
  geom_point(aes(x=as.numeric(Low), y=as.numeric(High), col= delta_deltaG), size=2, shape=15)+ 
  theme_classic()+ 
  geom_abline(intercept = 0, slope=1, lty=2, col='gray') + 
  geom_point(aes(x= 9.38, y=7.64), col='red', shape=4) + 
  scale_color_gradientn(colors = mygreyPalette) + 
  labs(y= 'High: PRM (log2 AU)', x= 'Low: PRM (log2 AU)') +facet_grid(.~mutation)+ 
  scale_y_continuous(breaks=seq(0,18, 3)) + scale_x_continuous(breaks=seq(0,18,3))  
#2.2, differences 
ggplot(data= mysings_ci[[3]][[2]]) + 
  geom_point(aes(x=Low, y=High-Low), alpha= 0.6, col='gray', shape= 1) + 
  geom_point(aes(x= mywts_out[[2]][2], y= -1.731983), col='red', shape=4) + theme_classic()+ 
  labs(x= "Low: PRM (log2 AU)", y= 'High-Low : PRM(log2 AU)') + geom_hline(yintercept = 0, lty=2, col='gray') + 
  facet_grid(.~mutation)
ggplot(data= mysings_ci[[3]][[2]]) + 
  geom_point(aes(x=High, y=High-Low), alpha= 0.6, col='gray', shape= 1) + 
  geom_point(aes(x= mywts_out[[2]][1], y= -1.731983), col='red', shape=4) + theme_classic()+ 
  labs(x= "High: PRM (log2 AU)", y= 'High-Low : PRM(log2 AU)') + geom_hline(yintercept = 0, lty=2, col='gray') + 
  facet_grid(.~mutation)
ggplot(data= mysings_ci[[3]][[2]]) + 
  geom_line(aes(x=delta_deltaG, y=High-Low)) + 
  geom_point(aes(x= 0, y= -1.731983), col='red', shape=4) + theme_classic()+ 
  labs(x= "", y= 'High-Low : PRM(log2 AU)') + geom_hline(yintercept = 0, lty=2, col='gray') + 
  facet_grid(.~mutation)+ geom_vline(xintercept = 0, lty=2, col='gray')
ggplot(data= mysings_ci[[3]][[1]]) + 
  geom_line(aes(x=delta_deltaG, y=High-Low)) + 
  geom_point(aes(x= 0, y= -2.50483), col='red', shape=4) + theme_classic()+ 
  labs(x= "", y= 'High-Low : PR(log2 AU)') + geom_hline(yintercept = 0, lty=2, col='gray') + 
  facet_grid(.~mutation) + geom_vline(xintercept = 0, lty=2, col='gray')

###### how doubles combine 
# low vs high 
######  25 combinations (5*5, excluding tetramerization and OR3), low vs high 

ggplot(data=doub_list[[6]][[1]][doub_list[[6]][[1]]$whichmeasure=='thermodynamic' & doub_list[[6]][[1]]$path1%in% c('folding','dimerization', 'binding to DNA', 'OR1', 'OR2') & 
                                  doub_list[[6]][[1]]$path2%in% c('folding','dimerization', 'binding to DNA', 'OR1', 'OR2') , ]) + 
  geom_point(aes(x=as.numeric(Low), y=as.numeric(High), col= delta_deltaA), size=1.6, shape=15)+ 
  theme(panel.background = element_rect(fill = "white", colour = "grey50"), legend.position = 'none',strip.text.x = element_text(size = 8),strip.text.y = element_text(size = 8))+ 
  geom_abline(intercept = 0, slope=1, lty=2, col='gray') + 
  geom_point(aes(y=4.564443, x=7.069273), col='red', shape=4) + 
  scale_color_gradientn(colors = mygreyPalette, name='') + facet_grid(path1~path2) + 
  labs(y= 'High: PR (log2 AU)', x= 'Low: PR (log2 AU)')  +
  scale_y_continuous(breaks=seq(0,18, 3)) + scale_x_continuous(breaks=seq(0,18,3))  

ggplot(data=doub_list[[6]][[2]][doub_list[[6]][[2]]$whichmeasure=='thermodynamic'& doub_list[[6]][[2]]$path1%in% c('folding','dimerization', 'binding to DNA', 'OR1', 'OR2') & 
                                  doub_list[[6]][[2]]$path2%in% c('folding','dimerization', 'binding to DNA', 'OR1', 'OR2'), ]) + 
  geom_point(aes(x=as.numeric(Low), y=as.numeric(High), col= delta_deltaA), size=1.6, shape=15)+ 
  theme(panel.background = element_rect(fill = "white", colour = "grey50"), legend.position = 'none',strip.text.x = element_text(size = 8),strip.text.y = element_text(size = 8))+ 
  geom_abline(intercept = 0, slope=1, lty=2, col='gray') + 
  geom_point(aes(x= 9.38, y=7.64), col='red', shape=4) + 
  scale_color_gradientn(colors = mygreyPalette) + facet_grid(path1~path2)+ 
  labs(y= 'High: PRM (log2 AU)', x= 'Low: PRM (log2 AU)') +
  scale_y_continuous(breaks=seq(0,18, 3)) + scale_x_continuous(breaks=seq(0,18,3))  

x1<- ggplot(data=doub_list[[6]][[1]][doub_list[[6]][[1]]$whichmeasure=='thermodynamic' & doub_list[[6]][[1]]$path1%in% c('tetramerization', 'OR3'), ]) + 
  geom_point(aes(x=as.numeric(Low), y=as.numeric(High), col= delta_deltaA), size=1.6, shape=15)+ 
  theme(panel.background = element_rect(fill = "white", colour = "grey50"), legend.position = 'none',strip.text.x = element_text(size = 8),strip.text.y = element_text(size = 8))+ 
  geom_abline(intercept = 0, slope=1, lty=2, col='gray') + 
  geom_point(aes(y=4.564443, x=7.069273), col='red', shape=4) + 
  scale_color_gradientn(colors = mygreyPalette, name='') + facet_grid(path1~path2) + 
  labs(y= 'High: PR (log2 AU)', x= 'Low: PR (log2 AU)') +
  scale_y_continuous(breaks=seq(0,18, 3)) + scale_x_continuous(breaks=seq(0,18,3))  

x2<- ggplot(data=doub_list[[6]][[1]][doub_list[[6]][[1]]$whichmeasure=='thermodynamic' & doub_list[[6]][[1]]$path1%in% c('tetramerization', 'OR3'), ]) + 
  geom_point(aes(x=as.numeric(Low), y=as.numeric(High), col= delta_deltaB), size=1.6, shape=15)+ 
  theme(panel.background = element_rect(fill = "white", colour = "grey50"), legend.position = 'none',strip.text.x = element_text(size = 8),strip.text.y = element_text(size = 8))+ 
  geom_abline(intercept = 0, slope=1, lty=2, col='gray') + 
  geom_point(aes(y=4.564443, x=7.069273), col='red', shape=4) + 
  scale_color_gradientn(colors = mygreyPalette, name='') + facet_grid(path1~path2) + 
  labs(y= 'High: PR (log2 AU)', x= 'Low: PR (log2 AU)') +
  scale_y_continuous(breaks=seq(0,18, 3)) + scale_x_continuous(breaks=seq(0,18,3))  
ggarrange(x1,x2, nrow=2, ncol=1)
# 8*5

x1<- ggplot(data=doub_list[[6]][[2]][doub_list[[6]][[2]]$whichmeasure=='thermodynamic'& doub_list[[6]][[1]]$path1%in% c('tetramerization', 'OR3'), ]) + 
  geom_point(aes(x=as.numeric(Low), y=as.numeric(High), col= delta_deltaA), size=1.6, shape=15)+ 
  theme(panel.background = element_rect(fill = "white", colour = "grey50"), legend.position = 'none',strip.text.x = element_text(size = 8),strip.text.y = element_text(size = 8))+ 
  geom_abline(intercept = 0, slope=1, lty=2, col='gray') + 
  geom_point(aes(x= 9.38, y=7.64), col='red', shape=4) + 
  scale_color_gradientn(colors = mygreyPalette) +  facet_grid(path1~path2)+ 
  labs(y= 'High: PRM (log2 AU)', x= 'Low: PRM (log2 AU)') +
  scale_y_continuous(breaks=seq(0,18, 3)) + scale_x_continuous(breaks=seq(0,18,3))  
x2<- ggplot(data=doub_list[[6]][[2]][doub_list[[6]][[2]]$whichmeasure=='thermodynamic'& doub_list[[6]][[1]]$path1%in% c('tetramerization', 'OR3'), ]) + 
  geom_point(aes(x=as.numeric(Low), y=as.numeric(High), col= delta_deltaB), size=1.6, shape=15)+ 
  theme(panel.background = element_rect(fill = "white", colour = "grey50"), legend.position = 'none',strip.text.x = element_text(size = 8),strip.text.y = element_text(size = 8))+ 
  geom_abline(intercept = 0, slope=1, lty=2, col='gray') + 
  geom_point(aes(x= 9.38, y=7.64), col='red', shape=4) + 
  scale_color_gradientn(colors = mygreyPalette) +  facet_grid(path1~path2)+ 
  labs(y= 'High: PRM (log2 AU)', x= 'Low: PRM (log2 AU)') +
  scale_y_continuous(breaks=seq(0,18, 3)) + scale_x_continuous(breaks=seq(0,18,3))  
ggarrange(x1,x2, nrow=2, ncol=1)
# 8*5

#####
##### 3. aditive vs thermodynamic 

#pr. 1
ggplot() + 
  geom_point(data=doub_list[[7]][[1]][(doub_list[[7]][[1]]$delta_deltaA<0 | doub_list[[7]][[1]]$delta_deltaB<0) &
                                        doub_list[[7]][[1]]$mutation%in% levels(doub_list[[7]][[1]]$mutation)[c(1:7)], ], 
             aes(x=additive, y=thermodynamic), col='cyan3', alpha=0.1, shape=19) + 
  geom_point(data=doub_list[[7]][[1]][doub_list[[7]][[1]]$delta_deltaA>=0 & doub_list[[7]][[1]]$delta_deltaB>=0 &
                                        doub_list[[7]][[1]]$mutation%in% levels(doub_list[[7]][[1]]$mutation)[c(1:7)], ], 
             aes(x=additive, y=thermodynamic), col='gray', alpha=0.2, shape=19)+
  stat_density2d(data=doub_list[[7]][[1]][doub_list[[7]][[1]]$delta_deltaA>=0 & doub_list[[7]][[1]]$delta_deltaB>=0 &
                                            doub_list[[7]][[1]]$mutation%in% levels(doub_list[[7]][[1]]$mutation)[c(1:7)], ], 
                 aes(x=additive, y=thermodynamic), col= 'red',alpha=0.5) + scale_fill_gradientn(colours = c(mygreyPalette[3],myredPalette[3:10])) + 
  theme(panel.background = element_rect(fill = "white", colour = "grey50"), strip.text.x = element_text(size = 10),strip.text.y = element_text(size = 10))+ 
  geom_abline(intercept = 0, slope=1, lty=2) +
  facet_grid(amount~mutation)  + xlim(2,20) + ylim(2,14) + labs(x= 'additive: PR (log2 AU)', y= 'thermodynamic: PR (log2 AU)')
# 8.5*2.28    

# prm.1
ggplot() + 
  geom_point(data=doub_list[[7]][[2]][(doub_list[[7]][[2]]$delta_deltaA<0 | doub_list[[7]][[2]]$delta_deltaB<0) &
                                        doub_list[[7]][[2]]$mutation%in% levels(doub_list[[7]][[2]]$mutation)[c(1:7)], ], 
             aes(x=additive, y=thermodynamic), col='cyan3', alpha=0.1, shape=19) + 
  geom_point(data=doub_list[[7]][[2]][doub_list[[7]][[2]]$delta_deltaA>=0 & doub_list[[7]][[2]]$delta_deltaB>=0 &
                                        doub_list[[7]][[2]]$mutation%in% levels(doub_list[[7]][[2]]$mutation)[c(1:7)], ], 
             aes(x=additive, y=thermodynamic), col='gray', alpha=0.2, shape=19)+
  stat_density2d(data=doub_list[[7]][[2]][doub_list[[7]][[2]]$delta_deltaA>=0 & doub_list[[7]][[2]]$delta_deltaB>=0 &
                                            doub_list[[7]][[2]]$mutation%in% levels(doub_list[[7]][[2]]$mutation)[c(1:7)], ], 
                 aes(x=additive, y=thermodynamic), col= 'red',alpha=0.5) + scale_fill_gradientn(colours = c(mygreyPalette[3],myredPalette[3:10])) + 
  theme(panel.background = element_rect(fill = "white", colour = "grey50"), strip.text.x = element_text(size = 9),strip.text.y = element_text(size = 9))+ 
  geom_abline(intercept = 0, slope=1, lty=2) +
  facet_grid(amount~mutation)  + xlim(0,14) + ylim(3,11) + labs(x= 'additive: PRM (log2 AU)', y= 'thermodynamic: PRM (log2 AU)')
#  

#pr. 2
ggplot() + 
  geom_point(data=doub_list[[7]][[1]][(doub_list[[7]][[1]]$delta_deltaA<0 | doub_list[[7]][[1]]$delta_deltaB<0) &
                                        doub_list[[7]][[1]]$mutation%in% levels(doub_list[[7]][[1]]$mutation)[c(8:14)], ], 
             aes(x=additive, y=thermodynamic), col='cyan3', alpha=0.1, shape=19) + 
  geom_point(data=doub_list[[7]][[1]][doub_list[[7]][[1]]$delta_deltaA>=0 & doub_list[[7]][[1]]$delta_deltaB>=0 &
                                        doub_list[[7]][[1]]$mutation%in% levels(doub_list[[7]][[1]]$mutation)[c(8:14)], ], 
             aes(x=additive, y=thermodynamic), col='gray', alpha=0.2, shape=19)+
  stat_density2d(data=doub_list[[7]][[1]][doub_list[[7]][[1]]$delta_deltaA>=0 & doub_list[[7]][[1]]$delta_deltaB>=0 &
                                            doub_list[[7]][[1]]$mutation%in% levels(doub_list[[7]][[1]]$mutation)[c(8:14)], ], 
                 aes(x=additive, y=thermodynamic), col= 'red',alpha=0.5) + scale_fill_gradientn(colours = c(mygreyPalette[3],myredPalette[3:10])) + 
  theme(panel.background = element_rect(fill = "white", colour = "grey50"), strip.text.x = element_text(size = 10),strip.text.y = element_text(size = 10))+ 
  geom_abline(intercept = 0, slope=1, lty=2) +
  facet_grid(amount~mutation)  + xlim(2,20) + ylim(2,14) + labs(x= 'additive: PR (log2 AU)', y= 'thermodynamic: PR (log2 AU)')
# 8.5*2.28    
#prm.2
ggplot() + 
  geom_point(data=doub_list[[7]][[2]][(doub_list[[7]][[2]]$delta_deltaA<0 | doub_list[[7]][[2]]$delta_deltaB<0) &
                                        doub_list[[7]][[2]]$mutation%in% levels(doub_list[[7]][[2]]$mutation)[c(8:14)], ], 
             aes(x=additive, y=thermodynamic), col='cyan3', alpha=0.1, shape=19) + 
  geom_point(data=doub_list[[7]][[2]][doub_list[[7]][[2]]$delta_deltaA>=0 & doub_list[[7]][[2]]$delta_deltaB>=0 &
                                        doub_list[[7]][[2]]$mutation%in% levels(doub_list[[7]][[2]]$mutation)[c(8:14)], ], 
             aes(x=additive, y=thermodynamic), col='gray', alpha=0.2, shape=19)+
  stat_density2d(data=doub_list[[7]][[2]][doub_list[[7]][[2]]$delta_deltaA>=0 & doub_list[[7]][[2]]$delta_deltaB>=0 &
                                            doub_list[[7]][[2]]$mutation%in% levels(doub_list[[7]][[2]]$mutation)[c(8:14)], ], 
                 aes(x=additive, y=thermodynamic), col= 'red',alpha=0.5) + scale_fill_gradientn(colours = c(mygreyPalette[3],myredPalette[3:10])) + 
  theme(panel.background = element_rect(fill = "white", colour = "grey50"), strip.text.x = element_text(size = 9),strip.text.y = element_text(size = 9))+ 
  geom_abline(intercept = 0, slope=1, lty=2) +
  facet_grid(amount~mutation)  + xlim(0,14) + ylim(3,11) + labs(x= 'additive: PRM (log2 AU)', y= 'thermodynamic: PRM (log2 AU)')
#  

#pr. 3+4
x1<-ggplot() + 
  geom_point(data=doub_list[[7]][[1]][(doub_list[[7]][[1]]$delta_deltaA<0 | doub_list[[7]][[1]]$delta_deltaB<0) &
                                        doub_list[[7]][[1]]$mutation%in% levels(doub_list[[7]][[1]]$mutation)[c(15:21)], ], 
             aes(x=additive, y=thermodynamic), col='cyan3', alpha=0.1, shape=19) + 
  geom_point(data=doub_list[[7]][[1]][doub_list[[7]][[1]]$delta_deltaA>=0 & doub_list[[7]][[1]]$delta_deltaB>=0 &
                                        doub_list[[7]][[1]]$mutation%in% levels(doub_list[[7]][[1]]$mutation)[c(15:21)], ], 
             aes(x=additive, y=thermodynamic), col='gray', alpha=0.2, shape=19)+
  stat_density2d(data=doub_list[[7]][[1]][doub_list[[7]][[1]]$delta_deltaA>=0 & doub_list[[7]][[1]]$delta_deltaB>=0 &
                                            doub_list[[7]][[1]]$mutation%in% levels(doub_list[[7]][[1]]$mutation)[c(15:21)], ], 
                 aes(x=additive, y=thermodynamic), col= 'red',alpha=0.5) + scale_fill_gradientn(colours = c(mygreyPalette[3],myredPalette[3:10])) + 
  theme(panel.background = element_rect(fill = "white", colour = "grey50"), strip.text.x = element_text(size = 10),strip.text.y = element_text(size = 10))+ 
  geom_abline(intercept = 0, slope=1, lty=2) +
  facet_grid(amount~mutation)  + xlim(2,20) + ylim(2,14) + labs(x= 'additive: PR (log2 AU)', y= 'thermodynamic: PR (log2 AU)')

x2<- ggplot() + 
  geom_point(data=doub_list[[7]][[1]][(doub_list[[7]][[1]]$delta_deltaA<0 | doub_list[[7]][[1]]$delta_deltaB<0) &
                                        doub_list[[7]][[1]]$mutation%in% levels(doub_list[[7]][[1]]$mutation)[c(22:28)], ], 
             aes(x=additive, y=thermodynamic), col='cyan3', alpha=0.1, shape=19) + 
  geom_point(data=doub_list[[7]][[1]][doub_list[[7]][[1]]$delta_deltaA>=0 & doub_list[[7]][[1]]$delta_deltaB>=0 &
                                        doub_list[[7]][[1]]$mutation%in% levels(doub_list[[7]][[1]]$mutation)[c(22:28)], ], 
             aes(x=additive, y=thermodynamic), col='gray', alpha=0.2, shape=19)+
  stat_density2d(data=doub_list[[7]][[1]][doub_list[[7]][[1]]$delta_deltaA>=0 & doub_list[[7]][[1]]$delta_deltaB>=0 &
                                            doub_list[[7]][[1]]$mutation%in% levels(doub_list[[7]][[1]]$mutation)[c(22:28)], ], 
                 aes(x=additive, y=thermodynamic), col= 'red',alpha=0.5) + scale_fill_gradientn(colours = c(mygreyPalette[3],myredPalette[3:10])) + 
  theme(panel.background = element_rect(fill = "white", colour = "grey50"), strip.text.x = element_text(size = 10),strip.text.y = element_text(size = 10))+ 
  geom_abline(intercept = 0, slope=1, lty=2) +
  facet_grid(amount~mutation)  + xlim(2,20) + ylim(2,14) + labs(x= 'additive: PR (log2 AU)', y= 'thermodynamic: PR (log2 AU)')
ggarrange(x1, x2, nrow=2, ncol=1, common.legend = F)
# 8.5*5 

x1<- ggplot() + 
  geom_point(data=doub_list[[7]][[2]][(doub_list[[7]][[2]]$delta_deltaA<0 | doub_list[[7]][[2]]$delta_deltaB<0) &
                                        doub_list[[7]][[2]]$mutation%in% levels(doub_list[[7]][[2]]$mutation)[c(15:21)], ], 
             aes(x=additive, y=thermodynamic), col='cyan3', alpha=0.1, shape=19) + 
  geom_point(data=doub_list[[7]][[2]][doub_list[[7]][[2]]$delta_deltaA>=0 & doub_list[[7]][[2]]$delta_deltaB>=0 &
                                        doub_list[[7]][[2]]$mutation%in% levels(doub_list[[7]][[2]]$mutation)[c(15:21)], ], 
             aes(x=additive, y=thermodynamic), col='gray', alpha=0.2, shape=19)+
  stat_density2d(data=doub_list[[7]][[2]][doub_list[[7]][[2]]$delta_deltaA>=0 & doub_list[[7]][[2]]$delta_deltaB>=0 &
                                            doub_list[[7]][[2]]$mutation%in% levels(doub_list[[7]][[2]]$mutation)[c(15:21)], ], 
                 aes(x=additive, y=thermodynamic), col= 'red',alpha=0.5) + scale_fill_gradientn(colours = c(mygreyPalette[3],myredPalette[3:10])) + 
  theme(panel.background = element_rect(fill = "white", colour = "grey50"), strip.text.x = element_text(size = 9),strip.text.y = element_text(size = 9))+ 
  geom_abline(intercept = 0, slope=1, lty=2) +
  facet_grid(amount~mutation)  + xlim(0,14) + ylim(3,11) + labs(x= 'additive: PRM (log2 AU)', y= 'thermodynamic: PRM (log2 AU)')
x2<- ggplot() + 
  geom_point(data=doub_list[[7]][[2]][(doub_list[[7]][[2]]$delta_deltaA<0 | doub_list[[7]][[2]]$delta_deltaB<0) &
                                        doub_list[[7]][[2]]$mutation%in% levels(doub_list[[7]][[2]]$mutation)[c(22:28)], ], 
             aes(x=additive, y=thermodynamic), col='cyan3', alpha=0.1, shape=19) + 
  geom_point(data=doub_list[[7]][[2]][doub_list[[7]][[2]]$delta_deltaA>=0 & doub_list[[7]][[2]]$delta_deltaB>=0 &
                                        doub_list[[7]][[2]]$mutation%in% levels(doub_list[[7]][[2]]$mutation)[c(22:28)], ], 
             aes(x=additive, y=thermodynamic), col='gray', alpha=0.2, shape=19)+
  stat_density2d(data=doub_list[[7]][[2]][doub_list[[7]][[2]]$delta_deltaA>=0 & doub_list[[7]][[2]]$delta_deltaB>=0 &
                                            doub_list[[7]][[2]]$mutation%in% levels(doub_list[[7]][[2]]$mutation)[c(22:28)], ], 
                 aes(x=additive, y=thermodynamic), col= 'red',alpha=0.5) + scale_fill_gradientn(colours = c(mygreyPalette[3],myredPalette[3:10])) + 
  theme(panel.background = element_rect(fill = "white", colour = "grey50"), strip.text.x = element_text(size = 9),strip.text.y = element_text(size = 9))+ 
  geom_abline(intercept = 0, slope=1, lty=2) +
  facet_grid(amount~mutation)  + xlim(0,14) + ylim(3,11) + labs(x= 'additive: PRM (log2 AU)', y= 'thermodynamic: PRM (log2 AU)')
#  
###### 4. 
#####  epistasis low vs high
mygs<- list()
for (i in 1:4){
  mygs[[i]] = 
    ggplot(data=doub_list[[3]][[1]][
      doub_list[[3]][[1]]$whichmeasure=='epis' & doub_list[[3]][[1]]$combination%in% 
        levels( doub_list[[3]][[1]]$combination)[c(((i-1)*7+1) : (i*7))], ],
      aes(x=as.numeric(Low), y=as.numeric(High))) + 
    geom_vline(xintercept = 0, lty=2, col='gray') + 
    geom_hline(yintercept = 0, lty=2, col='gray') + 
    geom_point(shape=19, alpha =0.1, col='gray', size=2) +  
    stat_density2d(col= 'red',alpha=0.5, n=50) + 
    theme(panel.background = element_rect(fill = "white", colour = "grey50"), strip.text.x = element_text(size = 10),strip.text.y = element_text(size = 10))+ 
    geom_abline(intercept = 0, slope=1, lty=2, col='gray') + 
    facet_grid(.~combination) + xlim(-8,8) + ylim(-8,8)+ labs(y= "High: Epistasis", x= "Low: Epistasis")
}
mygs[[1]] # 8.5*1.67
mygs[[2]] # 8.5*1.67
myg2<- list()
for (i in 1:4){
  myg2[[i]] = 
    ggplot(data=doub_list[[3]][[2]][
      doub_list[[3]][[2]]$whichmeasure=='epis' & doub_list[[3]][[2]]$combination%in% 
        levels( doub_list[[3]][[2]]$combination)[c(((i-1)*7+1) : (i*7))], ],
      aes(x=as.numeric(Low), y=as.numeric(High))) + 
    geom_vline(xintercept = 0, lty=2, col='gray') + 
    geom_hline(yintercept = 0, lty=2, col='gray') + 
    geom_point(shape=19, alpha =0.1, col='gray', size=2) +  
    stat_density2d(col= 'red',alpha=0.5, n=50) + 
    theme(panel.background = element_rect(fill = "white", colour = "grey50"), strip.text.x = element_text(size = 10),strip.text.y = element_text(size = 10))+ 
    geom_abline(intercept = 0, slope=1, lty=2, col='gray') + 
    facet_grid(.~combination) + xlim(-7,7) + ylim(-7,7)+ labs(x= "Low: Epistasis", y= "High: Epistasis")
}
myg2[[2]]

ggarrange(mygs[[3]], mygs[[4]],nrow=2,ncol=1) # 8.5*3.5
ggarrange(myg2[[3]], myg2[[4]],nrow=2,ncol=1)
#####
#### 5. tile plots 
myggs_epis1<- list(pr=list(), prm=list()) # 

for(j in 1:2){
  for (i in 1:4){
    myggs_epis1[[j]][[i]] = 
      ggplot(data=doub_list[[4]][[j]][
        doub_list[[4]][[j]]$whichmeasure=='epis' & doub_list[[4]][[j]]$combination%in% 
          levels( doub_list[[4]][[j]]$combination)[c(((i-1)*7+1) : (i*7))], ]) + 
      geom_tile(aes(x= delta_deltaA, y=delta_deltaB, fill= as.numeric(val))) + 
      theme_classic() + 
      scale_fill_gradient2('Epistasis',low= 'dark green' , high='magenta', mid='grey88', midpoint = 0) + 
      facet_grid(amount~combination) + 
      geom_contour(aes(x= delta_deltaA, y= delta_deltaB, z =as.numeric(val) ,colour = ..level..),colour = "black", binwidth=1 )  + 
      geom_vline(xintercept = 0, lty=2,col='grey') + 
      geom_hline(yintercept = 0, lty=2,col='grey') + 
      geom_abline(intercept = 0, slope =1, lty=2,col='grey') + labs(x='', y='')
  }
  
}

#####
### 6. sign epis tile plot
for(i in 1:2){
  doub_list[[7]][[i]]$epis_cla = factor(doub_list[[7]][[i]]$epis_cla, levels= c('none', 'magnitude','masking','sign','reciprocal sign'))
}

myggs_epis2<- list(pr=list(), prm=list()) # 

for(j in 1:2){
  for (i in 1:4){
    myggs_epis2[[j]][[i]] = 
      ggplot(data=doub_list[[7]][[j]][
        doub_list[[7]][[j]]$mutation%in% 
          levels( doub_list[[7]][[j]]$mutation)[c(((i-1)*7+1) : (i*7))], ]) + 
      geom_tile(aes(x= delta_deltaA, y=delta_deltaB, fill= epis_cla)) + 
      theme_classic() + 
      scale_fill_manual('Epistasis',values= c(mygreyPalette[c(1,3,5)], '#fff0db','#d9b99b')) + 
      facet_grid(amount~mutation) + 
      geom_vline(xintercept = 0, lty=2,col='grey') + 
      geom_hline(yintercept = 0, lty=2,col='grey') + 
      geom_abline(intercept = 0, slope =1, lty=2,col='grey') + labs(x='', y='')
  }
  
}
###
### delta epis & epis classes

delta<- list(pr=list(), prm=list())
for(i in 1:2) {
  for (j in 1:4){
    delta[[i]][[j]] = 
      ggplot(data=doub_list[[5]][[i]][
        doub_list[[5]][[i]]$combination%in% 
          levels( doub_list[[5]][[i]]$combination)[c(((j-1)*7+1) : (j*7))], ]) + 
      geom_tile(aes(x=delta_deltaA, y=delta_deltaB, fill=  as.numeric(H_min_L))) + 
      scale_fill_gradient2('Epistasis: High- Low',low= 'blue' , high='red', mid='grey88', midpoint = 0) + facet_grid(.~combination) + 
      theme_classic() + 
      geom_contour(aes(x= delta_deltaA, y= delta_deltaB, z =as.numeric(H_min_L) ,colour = ..level..),colour = "black", binwidth=1 )  + 
      geom_vline(xintercept = 0, lty=2,col='grey') + 
      geom_hline(yintercept = 0, lty=2,col='grey') + 
      geom_abline(intercept = 0, slope =1, lty=2,col='grey') + 
      labs(x='', y='') 
  }
  
}

myclass<- list(pr=list(), prm=list())
for(i in 1:2) {
  for (j in 1:4){
    myclass[[i]][[j]] = 
      ggplot(data=doub_list[[5]][[i]][
        doub_list[[5]][[i]]$combination%in% 
          levels( doub_list[[5]][[i]]$combination)[c(((j-1)*7+1) : (j*7))], ]) + 
      geom_tile(aes(x=delta_deltaA, y=delta_deltaB, fill=  LH_class)) + 
      scale_fill_manual('Epistasis: High-Low',values= c('grey', 'white', 'yellow','tan1'))+ facet_grid(.~combination) + 
      theme_classic() + 
      geom_contour(aes(x= delta_deltaA, y= delta_deltaB, z =as.numeric(H_min_L) ,colour = ..level..),colour = "black", binwidth=1 )  + 
      geom_vline(xintercept = 0, lty=2,col='grey') + 
      geom_hline(yintercept = 0, lty=2,col='grey') + 
      geom_abline(intercept = 0, slope =1, lty=2,col='grey') + 
      labs(x='', y='')
  }
  
}

ggarrange(delta[[1]][[1]], myclass[[1]][[1]], nrow=2, ncol=1)

######
#### 6. additive vs observed 
add_therm<- list(PR=list(), PRM=list())
mywts<- c(7, 7.6)

for (i in 1:2) {
  for(j in 1:28) {
    add_therm[[i]][[j]] = 
      ggplot(data= doub_list[[4]][[i]][
        doub_list[[4]][[i]]$whichmeasure%in% c('additive', 'thermodynamic') &
          doub_list[[4]][[i]]$combination==levels(doub_list[[4]][[i]]$combination)[j],]) + 
      geom_tile(aes(x= delta_deltaA, y=delta_deltaB, fill= as.numeric(val))) + 
      theme_classic() + labs(x= '', y='', title= levels(doub_list[[4]][[i]]$combination)[j], fill=  paste(names(add_therm)[i], 'log2 (AU)', sep=' '))+
      scale_fill_gradient(low='gray', high='dark green',
                          breaks= c(min(as.numeric(doub_list[[4]][[i]][doub_list[[4]][[i]]$whichmeasure=='thermodynamic', 'val'])), 
                                    mywts[i],
                                    max(as.numeric(doub_list[[4]][[i]][doub_list[[4]][[i]]$whichmeasure=='additive', 'val']))), 
                          labels=c("","","")) + facet_grid(amount~whichmeasure) + 
      geom_contour(aes(x= delta_deltaA, y= delta_deltaB, z = as.numeric(val) ,colour = ..level..),colour = "black", binwidth=1 )  + 
      geom_vline(xintercept = 0, lty=2,col='grey') + 
      geom_hline(yintercept = 0, lty=2,col='grey') + 
      geom_abline(intercept = 0, slope =1, lty=2,col='grey') 
  }
}


ggarrange(add_therm[[1]][[1]], add_therm[[1]][[2]], add_therm[[1]][[3]], add_therm[[1]][[4]],
          add_therm[[1]][[5]], add_therm[[1]][[6]], add_therm[[1]][[7]], add_therm[[1]][[8]],
          add_therm[[1]][[9]], add_therm[[1]][[10]], add_therm[[1]][[11]], add_therm[[1]][[12]],
          add_therm[[1]][[13]], add_therm[[1]][[14]], add_therm[[1]][[15]], add_therm[[1]][[16]],
          add_therm[[1]][[17]], add_therm[[1]][[18]], add_therm[[1]][[19]], add_therm[[1]][[20]],
          add_therm[[1]][[21]], add_therm[[1]][[22]], add_therm[[1]][[23]], add_therm[[1]][[24]],
          add_therm[[1]][[25]], add_therm[[1]][[26]], add_therm[[1]][[27]], add_therm[[1]][[28]],
          nrow=7,ncol=4, common.legend = T) # 12*24

ggarrange(add_therm[[2]][[1]], add_therm[[2]][[2]], add_therm[[2]][[3]], add_therm[[2]][[4]],
          add_therm[[2]][[5]], add_therm[[2]][[6]], add_therm[[2]][[7]], add_therm[[2]][[8]],
          add_therm[[2]][[9]], add_therm[[2]][[10]], add_therm[[2]][[11]], add_therm[[2]][[12]],
          add_therm[[2]][[13]], add_therm[[2]][[14]], add_therm[[2]][[15]], add_therm[[2]][[16]],
          add_therm[[2]][[17]], add_therm[[2]][[18]], add_therm[[2]][[19]], add_therm[[2]][[20]],
          add_therm[[2]][[21]], add_therm[[2]][[22]], add_therm[[2]][[23]], add_therm[[2]][[24]],
          add_therm[[2]][[25]], add_therm[[2]][[26]], add_therm[[2]][[27]], add_therm[[2]][[28]],
          nrow=7,ncol=4, common.legend = T) # 12*24
##### 7. 
#### distribution of mutational effects 

mygs_high<- list()
for(i in 1:28){
  mygs_high[[i]] = ggplot(data= doub_list[[8]][[1]][ doub_list[[8]][[1]]$amount=='High' & doub_list[[8]][[1]]$combination== levels(doub_list[[8]][[1]]$combination)[i], ]) + 
    geom_line(aes(x= as.numeric(val), linetype= whichmeasure, col=whichmeasure),stat='density') +
    theme_classic() + labs(title= levels(doub_list[[8]][[1]]$combination)[i], x= 'High: PR (log2 AU)')  + 
    scale_colour_manual(values= c('black', 'dark green',  cbbPalette[2]))+ xlim(4, 12)
}

ggarrange(mygs_high[[1]],mygs_high[[2]],mygs_high[[3]],mygs_high[[4]],mygs_high[[5]],mygs_high[[6]],mygs_high[[7]],
          mygs_high[[8]],mygs_high[[9]],mygs_high[[10]],mygs_high[[11]],mygs_high[[12]],mygs_high[[13]],mygs_high[[14]], 
          mygs_high[[15]],mygs_high[[16]],mygs_high[[17]],mygs_high[[18]],mygs_high[[19]],mygs_high[[20]],mygs_high[[21]],
          mygs_high[[22]],mygs_high[[23]],mygs_high[[24]],mygs_high[[25]],mygs_high[[26]],mygs_high[[27]],mygs_high[[28]], 
          nrow=4, ncol=7, common.legend = T)
# 15 & 8 

mygs_low<- list()
for(i in 1:28){
  mygs_low[[i]] = ggplot(data= doub_list[[8]][[1]][ doub_list[[8]][[1]]$amount=='Low' & doub_list[[8]][[1]]$combination== levels(doub_list[[8]][[1]]$combination)[i], ]) + 
    geom_line(aes(x= as.numeric(val), linetype= whichmeasure, col=whichmeasure),stat='density') +
    theme_classic() + labs(title= levels(doub_list[[8]][[1]]$combination)[i], x= 'Low: PR (log2 AU)')  + 
    scale_colour_manual(values= c('black', 'dark green',  cbbPalette[2]))+ xlim(4, 12)
}


ggarrange(mygs_low[[1]],mygs_low[[2]],mygs_low[[3]],mygs_low[[4]],mygs_low[[5]],mygs_low[[6]],mygs_low[[7]],
          mygs_low[[8]],mygs_low[[9]],mygs_low[[10]],mygs_low[[11]],mygs_low[[12]],mygs_low[[13]],mygs_low[[14]], 
          mygs_low[[15]],mygs_low[[16]],mygs_low[[17]],mygs_low[[18]],mygs_low[[19]],mygs_low[[20]],mygs_low[[21]],
          mygs_low[[22]],mygs_low[[23]],mygs_low[[24]],mygs_low[[25]],mygs_low[[26]],mygs_low[[27]],mygs_low[[28]], 
          nrow=4, ncol=7, common.legend = T)
####### with prm 

mygs_high<- list()
for(i in 1:28){
  mygs_high[[i]] = ggplot(data= doub_list[[8]][[2]][ doub_list[[8]][[2]]$amount=='High' & doub_list[[8]][[2]]$combination== levels(doub_list[[8]][[2]]$combination)[i], ]) + 
    geom_line(aes(x= as.numeric(val), linetype= whichmeasure, col=whichmeasure),stat='density') +
    theme_classic() + labs(title= levels(doub_list[[8]][[2]]$combination)[i], x= 'High: PRM (log2 AU)')  + 
    scale_colour_manual(values= c('black', 'dark green',  cbbPalette[2]))+ xlim(4, 10.5)
}

mygs_high[[10]]
ggarrange(mygs_high[[1]],mygs_high[[2]],mygs_high[[3]],mygs_high[[4]],mygs_high[[5]],mygs_high[[6]],mygs_high[[7]],
          mygs_high[[8]],mygs_high[[9]],mygs_high[[10]],mygs_high[[11]],mygs_high[[12]],mygs_high[[13]],mygs_high[[14]], 
          mygs_high[[15]],mygs_high[[16]],mygs_high[[17]],mygs_high[[18]],mygs_high[[19]],mygs_high[[20]],mygs_high[[21]],
          mygs_high[[22]],mygs_high[[23]],mygs_high[[24]],mygs_high[[25]],mygs_high[[26]],mygs_high[[27]],mygs_high[[28]], 
          nrow=4, ncol=7, common.legend = T)
# 15 & 8 

mygs_low<- list()
for(i in 1:28){
  mygs_low[[i]] = ggplot(data= doub_list[[8]][[2]][ doub_list[[8]][[2]]$amount=='Low' & doub_list[[8]][[2]]$combination== levels(doub_list[[8]][[2]]$combination)[i], ]) + 
    geom_line(aes(x= as.numeric(val), linetype= whichmeasure, col=whichmeasure),stat='density') +
    theme_classic() + labs(title= levels(doub_list[[8]][[1]]$combination)[i], x= 'Low: PRM (log2 AU)')  + 
    scale_colour_manual(values= c('black', 'dark green',  cbbPalette[2]))+ xlim(4, 10.5)
}

mygs_low[[10]]
ggarrange(mygs_low[[1]],mygs_low[[2]],mygs_low[[3]],mygs_low[[4]],mygs_low[[5]],mygs_low[[6]],mygs_low[[7]],
          mygs_low[[8]],mygs_low[[9]],mygs_low[[10]],mygs_low[[11]],mygs_low[[12]],mygs_low[[13]],mygs_low[[14]], 
          mygs_low[[15]],mygs_low[[16]],mygs_low[[17]],mygs_low[[18]],mygs_low[[19]],mygs_low[[20]],mygs_low[[21]],
          mygs_low[[22]],mygs_low[[23]],mygs_low[[24]],mygs_low[[25]],mygs_low[[26]],mygs_low[[27]],mygs_low[[28]], 
          nrow=4, ncol=7, common.legend = T)

