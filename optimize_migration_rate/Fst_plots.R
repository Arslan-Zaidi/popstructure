
library(data.table)
library(ggplot2)
library(dplyr)
library(here)
library(ggrepel)

fst100 = fread(here('optimize_migration_rate/grid/tau100/fst/genos_grid_t100_ss500.mall.fst'))

colnames(fst100) = c("m","chrom","rsid","position","count","fst")

fst100.mean= fst100%>%
  group_by(m)%>%
  summarize(fst = mean(fst))

fst9 = fread(here('optimize_migration_rate/grid/tau-9/fst/genos_grid_t9_ss250.mall.fst'))

colnames(fst9) = c("m","chrom","rsid","position","count","fst")

fst9.mean= fst9%>%
  group_by(m)%>%
  summarize(fst = mean(fst))

ggplot(fst100)+
  geom_density(aes(fst,fill=m))+
  facet_grid(m~.)+
  geom_vline(data=fst100.mean,
             aes(xintercept=fst),
             color="black")+
  geom_text(data=fst100.mean,
             aes(x=0.007,
                 y=1000,
                 label=paste("Fst :",round(fst,4),sep="")),
                 hjust=0,
             color="black")+
  theme_minimal()+
  theme(axis.text.y=element_blank(),
        panel.grid = element_blank(),
        strip.text = element_blank(),
        strip.background = element_blank(),
        )+
  labs(x="Fst", y ="Density",fill='Migration\n rate')+
  scale_fill_gradient(low="orange",high="white")

ggplot(fst9)+
  geom_density(aes(fst,fill=m))+
  facet_grid(m~.)+
  geom_vline(data=fst9.mean,
             aes(xintercept=fst),
             color="black")+
  geom_text(data=fst9.mean,
            aes(x=0.0025,
                y=1000,
                label=paste("Fst :",round(fst,4),sep="")),
            hjust=0,
            color="black")+
  theme_minimal()+
  theme(axis.text.y=element_blank(),
        panel.grid = element_blank(),
        strip.text = element_blank(),
        strip.background = element_blank(),
  )+
  labs(x="Fst", y ="Density",fill='Migration\n rate')+
  scale_fill_gradient(low="orange",high="white")

