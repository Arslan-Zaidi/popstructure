
library(readr)
library(here)
library(ggplot2)
library(patchwork)
library(data.table)

rds=list.files(here("gwas/grid/analysis/ldsc/"),pattern="*.rds")

rnd_ge=c("ldsc_feffects_ge_pcs0_grandom.rds",
                   "ldsc_feffects_ge_cm_grandom.rds",
                   "ldsc_feffects_ge_re_grandom.rds")

rnd_noge=c("ldsc_feffects_noge_pcs0_random.rds",
                   "ldsc_feffects_noge_cm_random.rds",
                   "ldsc_feffects_noge_re_random.rds")

sm_ge=c("ldsc_feffects_ge_pcs0_smooth.rds",
                   "ldsc_feffects_ge_cm_smooth.rds",
                   "ldsc_feffects_ge_re_smooth.rds")

sm_noge=c("ldsc_feffects_noge_pcs0_smooth.rds",
                   "ldsc_feffects_noge_cm_smooth.rds",
                   "ldsc_feffects_noge_re_smooth.rds")


shp_ge=c("ldsc_feffects_ge_pcs0_sharp.rds",
                   "ldsc_feffects_ge_cm_sharp.rds",
                   "ldsc_feffects_ge_re_sharp.rds")

shp_noge=c("ldsc_feffects_noge_pcs0_sharp.rds",
                    "ldsc_feffects_noge_cm_sharp.rds",
                    "ldsc_feffects_noge_re_sharp.rds")

rndge_plts=lapply(rnd_ge,function(x){
  plt=read_rds(here(
    paste("gwas/grid/analysis/ldsc/",x,sep="")
  ))
  return(plt)
})


smge_plts=lapply(sm_ge,function(x){
  plt=read_rds(here(
    paste("gwas/grid/analysis/ldsc/",x,sep="")
  ))
  return(plt)
})

shpge_plts=lapply(shp_ge,function(x){
  plt=read_rds(here(
    paste("gwas/grid/analysis/ldsc/",x,sep="")
    ))
  return(plt)
})

#function to get x and y limits from plot
get_lims=function(plt){
  #xlims=layer_scales(plt)$x$range$range
  ylims=layer_scales(plt)$y$range$range
  return(ylims)
}

#get x and y lims from all plots and choose the max


fplt=function(lplt){
  lims=sapply(lplt,get_lims)
  max.y=max(lims[,2])
  patches=( lplt[[1]] +
              theme(axis.title.x = element_blank(),
                    plot.title = element_text(hjust = 0.5,
                                              size=11),
                    axis.title.y = element_text(size=11),
                    axis.text = element_text(size=10)) +
              labs(title = "No correction") +
              ylim( c(1,max.y) ) ) +

    ( lplt[[2]] +
        theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x=element_text(size=10),
        axis.title.x = element_text(size=11),
        plot.title = element_text(hjust = 0.5,
                                  size=11)) +
        labs(title = "Common-PCA") +
        ylim(c(1,max.y)) ) +

    (lplt[[3]] +
        theme(axis.title = element_blank(),
              axis.text.y = element_blank(),
              axis.text.x = element_text(size=11),
              plot.title = element_text(hjust = 0.5,
                                        size=11)) +
        labs(title="Rare-PCA")+
       ylim(c(1,max.y)))

return(patches)
}

rndge_plt_comb=fplt(rndge_plts)
rndge_plt_comb[[2]] = rndge_plt_comb[[2]] + theme(axis.title.x=element_blank())

smge_plt_comb=fplt(smge_plts)
smge_plt_comb[[2]] = smge_plt_comb[[2]] + theme(axis.title.x=element_blank())
smge_plt_comb[[1]] = smge_plt_comb[[1]] + theme(plot.title=element_blank())
smge_plt_comb[[2]] = smge_plt_comb[[2]] + theme(plot.title=element_blank())
smge_plt_comb[[3]] = smge_plt_comb[[3]] + theme(plot.title=element_blank())

shpge_plt_comb=fplt(shpge_plts)
shpge_plt_comb[[1]] = shpge_plt_comb[[1]] + theme(plot.title=element_blank())
shpge_plt_comb[[2]] = shpge_plt_comb[[2]] + theme(plot.title=element_blank())
shpge_plt_comb[[3]] = shpge_plt_comb[[3]] + theme(plot.title=element_blank())


ge_plt=rndge_plt_comb / smge_plt_comb / shpge_plt_comb


##### no genetic variance

rnd_noge_plts=lapply(rnd_noge,function(x){
  plt=read_rds(here(
    paste("gwas/grid/analysis/ldsc/",x,sep="")
  ))
  return(plt)
})

sm_noge_plts=lapply(sm_noge,function(x){
  plt=read_rds(here(
    paste("gwas/grid/analysis/ldsc/",x,sep="")
  ))
  return(plt)
})

shp_noge_plts=lapply(shp_noge,function(x){
  plt=read_rds(here(
    paste("gwas/grid/analysis/ldsc/",x,sep="")
  ))
  return(plt)
})

rnd_noge_plt_comb=fplt(rnd_noge_plts)
rnd_noge_plt_comb[[2]] = rnd_noge_plt_comb[[2]] + theme(axis.title.x=element_blank())

sm_noge_plt_comb=fplt(sm_noge_plts)
sm_noge_plt_comb[[2]] = sm_noge_plt_comb[[2]] + theme(axis.title.x=element_blank())
sm_noge_plt_comb[[1]] = sm_noge_plt_comb[[1]] + theme(plot.title=element_blank())
sm_noge_plt_comb[[2]] = sm_noge_plt_comb[[2]] + theme(plot.title=element_blank())
sm_noge_plt_comb[[3]] = sm_noge_plt_comb[[3]] + theme(plot.title=element_blank())

shp_noge_plt_comb=fplt(shp_noge_plts)
shp_noge_plt_comb[[1]] = shp_noge_plt_comb[[1]] + theme(plot.title=element_blank())
shp_noge_plt_comb[[2]] = shp_noge_plt_comb[[2]] + theme(plot.title=element_blank())
shp_noge_plt_comb[[3]] = shp_noge_plt_comb[[3]] + theme(plot.title=element_blank())


noge_plt=rnd_noge_plt_comb / sm_noge_plt_comb / shp_noge_plt_comb
noge_plt


ggsave(here("analyses/ldsc/ldsc_grid_tau100_train_feffects_noge.png"),
       noge_plt,
       height=10,
       width=14,
       units="cm")

ggsave(here("analyses/ldsc/ldsc_grid_tau100_train_feffects_ge.png"),
       ge_plt,
       height=10,
       width=14,
       units="cm")
