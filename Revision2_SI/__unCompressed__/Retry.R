library(ggforce)

DataAS %>% 
  arrange(target,ligand,QAS) %>% 
  head()

DataAS %>% 
  filter(target=="csf1r") %>% 
  ggplot(aes(x=factor(QAS),y=ASVal,
             color=factor(type))) +
  #geom_sina() +
  #geom_violin() +
  geom_boxplot() +
  theme_bw() +
  coord_cartesian(ylim=c(-20,0))

DataAS %>% 
  group_by(target,QAS) %>% 
  distinct(type) %>% 
  summarise(count=n()) %>% 
  ungroup() %>% 
  distinct(count)

DataAS_doROC <- DataAS %>% 
  group_by(target,QAS) %>% 
  do(doROC=roc(.$type,.$ASVal,direction=">"))

DataAS_AUC <- DataAS_doROC %>% 
  mutate(ROC=unlist(doROC,recursive=F)$auc) %>% 
  select(-doROC) %>% 
  ungroup()

DataAS_AUC_mean <- DataAS_AUC %>%
  group_by(QAS) %>% 
  summarise(ROC=mean(ROC)) %>% 
  arrange(desc(ROC))
DataAS_AUC_mean


#####
Plot_Self_Means_sepRuns_Full <- DataOTD_Self_sepRuns %>% 
  filter(QRunSurf%in%c(0,0.3,0.5,0.7,1)) %>% 
  ggplot(aes(x=QSurf,y=ROC,color=factor(QRunSurf))) +
  geom_hline(data=DataAS_AUC_mean %>% mutate(SF="dG"),
             aes(yintercept=ROC),color="grey10",size=0.5) +
  #geom_smooth(aes(group=paste0(QAS,QRunSurf)),size=0.5,level=0.75,se=F) +
  geom_line(aes(group=paste0(QAS,QRunSurf)),size=0.5) +
  scale_color_manual(values=c("#b15928","#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", 
                              "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a"),
                     name="on-top runs\nscores quantile") +
  theme_bw() +
  facet_grid(target~QAS) +
  theme(legend.position="bottom",        
        legend.box.spacing = unit(0,"mm"),
        plot.margin = margin(1,5,0,1,unit="mm"),
        legend.box.margin =  margin(0,0,0,0,unit="mm"),
        legend.key.height = unit(1,"mm"),
        legend.key.width = unit(10,"mm")) +
  scale_x_continuous(expand=c(0,0)) +
  ylab("AUROC") + xlab("On-top sites scores quantile") 
Plot_Self_Means_sepRuns_Full
ggsave(filename="Plot_Self_Quantiles_SI_Prod_35prot-2.png",Plot_Self_Means_sepRuns_Full, 
       width = 16, height = 55, units = "cm", dpi=1000)
