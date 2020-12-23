
library(data.table)
library(ggplot2)
library(dplyr)
library(cowplot)

crand1 = fread("gwas/complex_dem/test/prs/grandom/complexdem_prs_grandom.1.c.sscore")
crand2 = fread("gwas/complex_dem/test/prs/grandom/complexdem_prs_grandom.1.c.p.sscore")
crand3 = fread("gwas/complex_dem/test/prs/grandom/complexdem_prs_grandom.1.nc.sscore")

colnames(crand1)= colnames(crand2) = colnames(crand3) = c("IID","dosage","pcs0","cm","re","cmre")

pop.test = fread("gwas/complex_dem/genos_complex_l1e7_ss500_m0.07.test.pop")
pop.test$IID = paste("tsk_",seq(1,17999,2),sep="")

plt_all = function(crand){
crand$rep = rep(1:20,each=9000)
crand = merge(crand,
               pop.test[,c("IID","longitude","latitude")],
               by="IID")

mcrand = melt(crand,
               id.vars = c("IID","dosage","rep","longitude","latitude"))

mcrand.r = mcrand%>%
  group_by(variable,rep)%>%
  summarize(rlat = cor(value,latitude),
            rlong = cor(value,longitude))

plts.rand = lapply(seq(1,20),function(i){
  df = mcrand%>%
    filter(rep == i & variable =="pcs0")%>%
    group_by(variable,longitude,latitude)%>%
    summarize(gvalue = mean(value))
  scale_mean = mean(df$gvalue)
  return(ggplot(df)+
           geom_tile(aes(longitude,latitude,fill=gvalue))+
           scale_fill_gradient2(high = "#fc8d59", 
                                mid = "#ffffbf", 
                                low = "#91bfdb",
                                midpoint = scale_mean)+
           theme_void()+
           theme(legend.position="none"))
  
})

return(list(plot_grid(plotlist = plts.rand), mcrand.r))
}

cmplx1 = plt_all(crand1)
cmplx2 = plt_all(crand2)
cmplx3 = plt_all(crand3)
