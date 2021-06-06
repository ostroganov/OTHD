#Other metrics

QNames2 <- c(
  `auroc` = "AUROC",
  `bedroc_5` = "BEDROC, alpha=5",
  `ef_5` = "Enrichment factor at 5%"
)

DataAS_othMetr <- DataAS %>% 
  group_by(target,QAS) %>% 
  summarise(auroc=enrichvs::auc(ASVal,type,decreasing=F),
            bedroc_5=bedroc(ASVal,type,decreasing=F,alpha=5),
            ef_5=enrichment_factor(ASVal,type,decreasing=F)) %>% 
  gather("Metrics","Accuracy",bedroc_5,ef_5,auroc)

DataOTD_Self_sepRuns_othMetr <- 
  left_join(DataAS,DataSurf_sepRuns) %>% 
  filter(QAS==0) %>% 
  group_by(target,QAS,QSurf,QRunSurf) %>% 
  summarise(auroc=enrichvs::auc(ASVal-SurfVal,type,decreasing=F),
            bedroc_5=bedroc(ASVal-SurfVal,type,decreasing=F,alpha=5),
            ef_5=enrichment_factor(ASVal-SurfVal,type,decreasing=F)) %>% 
  gather("Metrics","Accuracy",bedroc_5,ef_5,auroc)

Plot_Self_Means_sepRuns_othMetr <- DataOTD_Self_sepRuns_othMetr %>% 
  filter(!target%in%c(SkipPept),
         QRunSurf!="Arb") %>% 
  ggplot(aes(x=QSurf,y=Accuracy,color=factor(QRunSurf))) +
  geom_hline(data=DataAS_othMetr %>% mutate(SF="dG") %>% 
               filter(QAS==0) %>% 
               mutate(QAS=case_when(QAS==0~"Min",
                                    QAS==0.5~"Med",
                                    QAS==1~"Max")) %>% group_by(Metrics) %>% 
               filter(QAS=="Min") %>% summarise(Accuracy=mean(Accuracy)),
             aes(yintercept=Accuracy),color="grey10",size=0.5) +
  geom_smooth(aes(group=paste0(QAS,QRunSurf)),size=0.5,level=0.75,se=F) +
  #geom_line(aes(group=paste0(QAS,QRunSurf))) +
  #geom_point() +
  scale_color_manual(values=c("#b15928","#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", 
                              "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a"),
                     name="On-top runs\nscores quantile") +
  theme_bw() +
  facet_grid(Metrics~.,scales="free_y",labeller = as_labeller(QNames2)) +
  theme(legend.position="bottom",        
        legend.box.spacing = unit(0,"mm"),
        plot.margin = margin(1,5,0,1,unit="mm"),
        legend.box.margin =  margin(0,0,0,0,unit="mm"),
        legend.key.height = unit(1,"mm"),
        legend.key.width = unit(10,"mm")) +
  scale_x_continuous(expand=c(0,0)) +
  ylab("Average accuracy (see metrics on the right)") + xlab("On-top sites scores quantile") 
Plot_Self_Means_sepRuns_othMetr
ggsave(filename="Plot_Self_Quantiles_SI_allMetrics.png",Plot_Self_Means_sepRuns_othMetr, width = 16, height = 19, units = "cm", dpi=2000)


DataOTD_grp1_sepRuns_othMetr <- bind_rows(DataOTD_Self_sepRuns_othMetr %>% 
                                            mutate(QRunSurf=case_when(QRunSurf==1~"Max",
                                                                      QRunSurf==0.5~"Med",
                                                                      QRunSurf==0~"Min",
                                                                      TRUE~QRunSurf),
                                                   Baseline="Self") %>% 
                                            filter(QRunSurf%in%c("Max","Med","Min","Arb")),
                                          DataAS_othMetr %>% 
                                            crossing(QRunSurf=c("Max","Med","Min","Arb")) %>% 
                                            mutate(Baseline="None"))

#GainTab_allMetr <- 
DataOTD_grp1_sepRuns_othMetr %>% 
  filter(Baseline%in%c("Self","None"),
         QSurf%in%c(NA,0.5),
         QAS%in%c(NA,0),
         QRunSurf%in%c(NA,"Max")) %>% 
  mutate(Baseline=factor(Baseline,levels=c("None","Self"))) %>%
  left_join(TestSet %>% select(ID,Type,`Full name`,OpenAS) %>% 
              rename(target=ID,Target=`Full name`) %>% 
              mutate(OpenAS=ifelse(OpenAS==1,"yes","no"))) %>% 
  left_join(TestSet2 %>% select(ID,TotCh) %>% rename(target=ID)) %>% 
  ungroup() %>% 
  select(target,Target,Type,OpenAS,TotCh,Metrics,Baseline,Accuracy) %>% 
  spread(Baseline,Accuracy) %>% 
  mutate(Gain=Self-None) %>% 
  gather("key","val",None,Self,Gain) %>% 
  mutate(key=factor(key,levels=c("None","Self","Gain"))) %>% 
  dcast(Target+target+Type+OpenAS+TotCh~Metrics+key,value.var = "val") %>% 
  arrange(auroc_None) %>% 
  write_delim("OTD_GainTable_allMetr_Prod.csv", ";")

