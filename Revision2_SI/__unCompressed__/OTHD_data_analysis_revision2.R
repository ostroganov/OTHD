library(tidyverse)
library(reshape2)
library(pROC)
library(readxl)
library(multidplyr) # available at https://github.com/tidyverse/multidplyr
library(enrichvs) # available at https://github.com/cran/enrichvs
#library(doSNOW) # for parallelization on Windows
library(doParallel) # for parallelization on Linux
#library(plotly) #required for ggplotly interactive plot

cl <- new_cluster(10) # Creates a 10-threads cluster. Reduce if you do not have 10 threads.
cluster_library(cl, c("tidyverse","pROC"))

#####
#Load data
#####

#Raw data:
#
Data10sampAS <- read_delim(file="DataAS_10runs.csv",";")
Data10sampOTS <- read_delim(file="DataOTS_10runs.csv",";")
Data3sampBasOTS <- read_delim(file="DataBasOTS_3runs.csv",";")

#Semiproduct data:
#
DataAS <- read_delim("DataAS_Prod.csv",",")
DataSurf_sepRuns <- read_delim("DataSurf_sepRuns_Prod.csv",",",
                               col_types = cols(
                                 target = col_character(),
                                 ligand = col_character(),
                                 QSurf = col_double(),
                                 QRunSurf = col_character(),
                                 SurfVal = col_double()
                               ))
DataOTD_grp1_sepRuns <- read.csv("DataOTD_grp1_sepRuns_Prod.csv")
DataOTD_grp2_sepRuns <- read.csv("DataOTD_grp2_sepRuns_Prod.csv")

TestSet <- read_xlsx("TestSet.xlsx")

SkipPept <- c("")
QNames <- c(
  `1` = "Maximum",
  `0.5` = "Median",
  `0` = "Minimum",
  `Arb` = "Arbitrary"
)
scaleFUN <- function(x) sprintf("%.2f", x)

#Check
Data10sampAS %>% filter(type==0) %>% distinct(target,ligand) %>% group_by(target) %>% summarise(count=n()) %>% View()
Data10sampOTS %>% distinct(target,ligand) %>% group_by(target) %>% summarise(count=n()) %>% View()
Data3sampBasOTS %>% distinct(target,ligand) %>% group_by(target) %>% summarise(count=n()) %>% View()

Data10sampOTS %>% group_by(target,grid) %>% 
  distinct(ligand) %>% 
  summarise(countLig=n()) %>% 
  View()

#####
#Active site
#####

DataAS <- Data10sampAS %>% 
  group_by(target,ligand) %>% 
  crossing(.,data_frame(QAS=c(0,0.5,1))) %>%
  group_by(target,ligand,QAS,type) %>% 
  partition(cluster=cl) %>% 
  summarize(ASVal=quantile(top.dG,probs=unique(QAS))) %>% 
  collect()

colnames(DataAS)
DataAS %>% 
  write.csv("DataAS_Prod.csv", row.names = F, quote = F)

#Continue here if you have loaded data from DataAS_Prod.csv

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

DataAS_AUC %>% 
  dcast(target~QAS)

DataAS_AUC_mean <- DataAS_AUC %>%
  group_by(QAS) %>% 
  summarise(ROC=mean(ROC)) %>% 
  arrange(desc(ROC))
DataAS_AUC_mean

DataAS_ROC <- DataAS_doROC %>% 
  group_by(target,QAS) %>%
  mutate(sens=list(unlist(doROC,recursive=F)$sensitivities),
         spec=list(unlist(doROC,recursive=F)$specificities)) %>% 
  select(-doROC) %>% 
  unnest(cols = c(sens, spec))

DataAS %>%
  filter(QAS==0) %>% 
  group_by(target) %>% 
  summarise(NumLig=sum(type==1),
            NumDec=sum(type==0)) %>% 
  View()

#####
#On-top Docking
#####

GetSurfRun <- compiler::cmpfun(function(QRS,Energies){
  Quant=unique(QRS)
  if(Quant=="Arb"){
    sample(Energies,1) %>% return()
  }
  else{
    quantile(Energies,probs=as.numeric(Quant)) %>% return()
  }
})

GetSurfRunVal <- compiler::cmpfun(function(df){
  df %>%
    crossing(data_frame(QRunSurf=c(seq(0,1,0.1),"Arb"))) %>%
    group_by(QRunSurf,ligand,grid) %>% 
    partition(cluster=cl) %>% 
    summarize(SurfRunVal=GetSurfRun(QRunSurf,top.dG)) %>% 
    collect()
})

cluster_copy(cl,"GetSurfRun")
DataSurf_sepRuns_0 <- Data10sampOTS %>% 
  group_by(target) %>% 
  do(SurfRunVal0=GetSurfRunVal(.)) %>% 
  ungroup() %>% 
  unnest(cols = c(SurfRunVal0))

GetSurfVal <- compiler::cmpfun(function(df){
  df %>% 
    crossing(data_frame(QSurf=seq(0,1,0.05))) %>%
    group_by(ligand,QSurf,QRunSurf) %>% 
    partition(cluster=cl) %>%
    summarize(SurfVal=quantile(SurfRunVal,probs=unique(QSurf))) %>% 
    collect()
})

DataSurf_sepRuns <- DataSurf_sepRuns_0 %>% 
  group_by(target) %>%
  do(SurfVal0=GetSurfVal(.)) %>% 
  ungroup() %>% 
  unnest(cols = c(SurfVal0))

colnames(DataSurf_sepRuns)
DataSurf_sepRuns %>% 
  write.csv("DataSurf_sepRuns_Prod.csv", row.names = F, quote = F)

#Continue here if you have loaded data from DataSurf_sepRuns_Prod.csv and DataAS_Prod.csv

DataOTD_Self_sepRuns_doROC <- 
  left_join(DataAS,DataSurf_sepRuns) %>% 
  #filter(QAS==0) %>% 
  group_by(target,QAS,QSurf,QRunSurf) %>% 
  do(doROC=roc(.$type,.$ASVal-.$SurfVal,direction=">")) 

DataOTD_Self_sepRuns <- DataOTD_Self_sepRuns_doROC %>% 
  mutate(ROC=unlist(doROC,recursive=F)$auc) %>% 
  select(-doROC) %>% 
  ungroup() %>% 
  mutate(Baseline="Self") 

Plot_Self_Means_sepRuns_Full <- DataOTD_Self_sepRuns %>% 
  filter(QRunSurf!="Arb") %>% 
  ggplot(aes(x=QSurf,y=ROC,color=factor(QRunSurf))) +
  geom_hline(data=DataAS_AUC_mean %>% mutate(SF="dG") %>% 
               filter(QAS==0) %>% 
               mutate(QAS=case_when(QAS==0~"Min",
                                    QAS==0.5~"Med",
                                    QAS==1~"Max")) %>% filter(QAS=="Min"),
             aes(yintercept=ROC),color="grey10",size=0.5) +
  geom_smooth(aes(group=paste0(QAS,QRunSurf)),size=0.5,level=0.75,se=F) +
  scale_color_manual(values=c("#b15928","#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", 
                              "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a"),
                     name="on-top runs\nscores quantile") +
  theme_bw() +
  theme(legend.position="bottom",        
        legend.box.spacing = unit(0,"mm"),
        plot.margin = margin(1,5,0,1,unit="mm"),
        legend.box.margin =  margin(0,0,0,0,unit="mm"),
        legend.key.height = unit(1,"mm"),
        legend.key.width = unit(10,"mm")) +
  scale_x_continuous(expand=c(0,0)) +
  ylab("Average AUROC") + xlab("On-top sites scores quantile") 
Plot_Self_Means_sepRuns_Full
ggsave(filename="Plot_Self_Quantiles_SI_Prod.png",Plot_Self_Means_sepRuns_Full, width = 16, height = 10, units = "cm", dpi=2000)

Plot_Self_Means_sepRuns <- DataOTD_Self_sepRuns %>% 
  filter(QRunSurf%in%c(0,0.5,1,"Arb")) %>% 
  mutate(QRunSurf=factor(QRunSurf,levels=c("Arb",0,0.5,1))) %>% 
  ggplot(aes(x=QSurf,y=ROC,color=factor(QRunSurf))) +
  geom_hline(data=DataAS_AUC_mean %>% mutate(SF="dG") %>% 
               filter(QAS==0) %>% 
               mutate(QAS=case_when(QAS==0~"Min",
                                    QAS==0.5~"Med",
                                    QAS==1~"Max")) %>% filter(QAS=="Min"),
             aes(yintercept=ROC),color="grey10",size=0.5) +
  geom_smooth(aes(group=paste0(QAS,QRunSurf)),size=0.5,level=0.75) +
  scale_color_manual(values=c("#984ea3","#4daf4a","#377eb8","#e41a1c","#ff7f00"),
                     guide = guide_legend(reverse = TRUE), 
                     name="Quantile of the\non-top runs used") +
  theme_bw() +
  guides(color=F) +
  theme(legend.position="bottom",     
        axis.title.y = element_text(vjust = 0.25),
        text=element_text(size=7),
        legend.box.spacing = unit(0,"mm"),
        plot.margin = margin(1,2,0,1,unit="mm"),
        legend.box.margin =  margin(0,0,0,0,unit="mm"),
        legend.key.size = unit(1,"mm")) +
  scale_x_continuous(expand=c(0,0)) +
  ylab("Average AUROC") + xlab("On-top scores quantile") 
Plot_Self_Means_sepRuns
ggsave(filename="Plot_Self_Quantiles_Prod.png",Plot_Self_Means_sepRuns, width = 4.3, height = 3.0, units = "cm", dpi=2000)

#####
#ROC curves comparison
#####

AS_ord <- DataAS_AUC %>% 
  filter(QAS==0) %>%
  arrange(ROC) %>% 
  .$target

DataOTD_Self_sepRuns_ROC <- DataOTD_Self_sepRuns_doROC %>% 
  group_by(target,QAS,QSurf,QRunSurf) %>% 
  mutate(sens=list(unlist(doROC,recursive=F)$sensitivities),
         spec=list(unlist(doROC,recursive=F)$specificities)) %>% 
  select(-doROC) %>% 
  unnest(cols = c(sens, spec)) %>% 
  mutate(Baseline="Self")

ROC_curves_Plot <- bind_rows(DataAS_ROC %>% mutate(Baseline="None"),
          DataOTD_Self_sepRuns_ROC) %>% 
  ungroup() %>% 
  mutate(PDB=factor(target,levels=rev(AS_ord))) %>%  
  filter(QAS==0,QRunSurf%in%c(NA,1),QSurf%in%c(NA,0.5)) %>% 
  ggplot(aes(y=sens,x=1-spec,color=Baseline)) +
  scale_color_brewer(palette="Set1",
                     breaks=c("Self","None"),
                     labels=c("On-top docking","Conventional docking")) +
  scale_x_continuous(name="False positive rate",breaks=c(0,0.5,1)) +
  scale_y_continuous(name="True positive rate",breaks=c(0,0.5,1)) +
  geom_line() +
  facet_wrap(~PDB,ncol=5) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.direction = "vertical",
        legend.title = element_blank())
ROC_curves_Plot
ggsave(filename="ROC_curves_Plot.png",ROC_curves_Plot, width = 16.5, height = 20, units = "cm", dpi=2000)

#####
#Compare rankings by conventional and On-Top docking
#####

RankTab <- left_join(DataAS,DataSurf_sepRuns) %>% 
  filter(QAS==0,QSurf==0.5,QRunSurf==1) %>% 
  group_by(target,QAS,QSurf,QRunSurf) %>% 
  mutate(OTVal=ASVal-SurfVal,
         ASrank=rank(ASVal),
         OTrank=rank(OTVal)) %>% 
  filter(type==1)

RankPlot <- RankTab %>% 
  mutate(PDB=factor(target,levels=rev(AS_ord))) %>%
  ggplot(aes(x=ASrank,y=OTrank)) +
  geom_abline(slope=1,intercept=0,color="grey30") +
  geom_point(color="red") +
  facet_wrap(~PDB,ncol=5) +
  theme_bw() +
  scale_x_continuous(breaks=seq(0,10000,250),name="Ranking according to the conventional docking", limits=c(0,1100)) +
  scale_y_continuous(breaks=seq(0,10000,250),name="Ranking according to the On-Top docking", limits=c(0,1100)) +
  theme(legend.position = "bottom",
        legend.direction = "vertical",
        axis.text.x = element_text(angle = -45, vjust=1.0, hjust=0.0),
        legend.title = element_blank())
RankPlot
ggsave(filename="Ranking_Plot.png", RankPlot, width = 16.5, height = 20, units = "cm", dpi=2000)

#####
#Gain Plot
#####

DataOTD_grp1_sepRuns <- bind_rows(DataOTD_Self_sepRuns %>% 
                                     mutate(QRunSurf=case_when(QRunSurf==1~"Max",
                                                               QRunSurf==0.5~"Med",
                                                               QRunSurf==0~"Min",
                                                               TRUE~QRunSurf)) %>% 
                                     filter(QRunSurf%in%c("Max","Med","Min","Arb")),
                                   DataAS_AUC %>% 
                                     crossing(QRunSurf=c("Max","Med","Min","Arb")) %>% 
                                     mutate(Baseline="None"))

DataOTD_grp1_sepRuns %>% head()
DataOTD_grp1_sepRuns %>% 
  write.csv("DataOTD_grp1_sepRuns_Prod.csv", row.names = F, quote = F)

#Continue here if you have loaded data from DataOTD_grp1_sepRuns_Prod.csv

OTD_GainPlot_Prod <- DataOTD_grp1_sepRuns %>% filter(Baseline%in%c("Self","None"),
                                                       QSurf%in%c(NA,0.5),
                                                       QAS%in%c(NA,0),
                                                       QRunSurf%in%c(NA,"Max")) %>% 
  mutate(Baseline=factor(Baseline,levels=c("None","Self"))) %>%
  left_join(TestSet %>% select(ID,SType) %>% rename(target=ID)) %>% 
  ggplot(aes(y=ROC,x=Baseline,color=SType)) +
  geom_line(aes(group=target),arrow=arrow(angle=15,type="open",length=unit(1,"mm")),alpha=0.8) +
  guides(color=guide_legend(ncol=2)) +
  scale_color_manual(values=c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", 
                              "#c994c7", "#a65628", "#f781bf", "#999999")) +
  scale_y_continuous(breaks=seq(0,1,0.05), name="AUROC", labels=scaleFUN) + 
  scale_x_discrete(limits=c("None","Self"),
                   labels=c("ASD","OTD"),
                   expand = c(0.1,0.1),
                   name="Docking") +
  theme_bw() +
  theme(legend.position="bottom",
        text=element_text(size=7),
        axis.title.y = element_text(vjust = 0.25),
        legend.title = element_blank(),
        plot.margin = margin(1,1,0,1,unit="mm"),
        axis.title.x = element_blank(),
        panel.spacing = unit(0,"mm"),
        legend.text = element_text(size=5),
        legend.margin=margin(0,0,0,0),
        legend.box.spacing = unit(0.0,"mm"),
        legend.box.margin =  margin(-0.5,5,0,0,unit="mm"),
        legend.key.size = unit(1,"mm"))
OTD_GainPlot_Prod
ggsave(filename="OTD_GainPlot_Prod.png",OTD_GainPlot_Prod, width = 2.9, height = 6, units = "cm", dpi=2000)

GainTab <- DataOTD_grp1_sepRuns %>% filter(Baseline%in%c("Self","None"),
                                      QSurf%in%c(NA,0.5),
                                      QAS%in%c(NA,0),
                                      QRunSurf%in%c(NA,"Max")) %>% 
  mutate(Baseline=factor(Baseline,levels=c("None","Self"))) %>%
  left_join(TestSet %>% select(ID,Type,`Full name`) %>% 
              rename(target=ID,Target=`Full name`)) %>% 
  dcast(Target+target+Type~Baseline,value.var = "ROC") %>% 
  mutate(Gain=Self-None)

GainTab %>% 
  mutate(Gain=sprintf("%.1f", 100*Gain),
         Self=sprintf("%.1f", 100*Self),
         None=sprintf("%.1f", 100*None)) %>% 
  arrange(as.numeric(None)) %>% 
  format(digits=5) %>% 
  write.csv2("OTD_GainTable_Prod.csv",row.names = F, quote = F)

GainTab %>% 
  mutate(Gain=sprintf("%.1f", 100*Gain),
         Self=sprintf("%.1f", 100*Self),
         None=sprintf("%.1f", 100*None)) %>% 
  arrange(as.numeric(None)) %>% 
  format(digits=5) %>% 
  left_join(TestSet %>% rename(target=ID) %>% select(target,TotCharge)) %>% 
  arrange(None) %>% 
  write_delim("ProtList_wCharges.csv",";")

DataOTD_grp1_sepRuns %>% filter(Baseline%in%c("Self","None"),
                                 QSurf%in%c(NA,0.5),
                                 QAS%in%c(NA,0),
                                 QRunSurf%in%c(NA,"Max")) %>% 
  mutate(Baseline=factor(Baseline,levels=c("None","Self"))) %>%
  left_join(TestSet %>% select(ID,Type,`Full name`) %>% 
              rename(target=ID,Target=`Full name`)) %>% 
  dcast(Target+target+Type~Baseline,value.var = "ROC") %>% 
  summarise(Gain=sprintf("%.1f", 100*mean(Self-None)),
         Self=sprintf("%.1f", 100*mean(Self)),
         None=sprintf("%.1f", 100*mean(None))) %>% 
  select(None,Self,Gain)

#####
#Analyse vPAINs (frequent hitters)
#requires DataAS and DataSurf_sepRuns
#####

NumDecoys <- DataAS %>%
  filter(QAS==0) %>% 
  group_by(target) %>% 
  summarise(NumDec=sum(type==0)) %>% 
  .$NumDec %>% unique() %>% unlist()

GetTop <- function(ntop,Val,ligand,type) {
  NTOP <- unique(ntop)
  data_frame(ligand=ligand,Val=Val,type=type) %>% 
    top_n(NTOP,-Val) %>% 
    filter(type==0) %>% 
    .$ligand
}

FreqHitAS <- DataAS %>% 
  ungroup() %>% 
  filter(QAS==0) %>%   
  mutate(ligand=as.character(ligand)) %>% 
  crossing(ntop=seq(10,NumDecoys,10)) %>% 
  group_by(ntop,target) %>% 
  do(ligand = GetTop(.$ntop,.$ASVal,.$ligand,.$type)) %>% 
  unnest(cols = c(ligand)) %>% 
  group_by(ntop,ligand) %>% 
  summarise(count=n()) %>% 
  arrange(desc(count))

FreqHitAS %>% filter(ntop==20) %>%
  ungroup() %>% 
  select(ligand,count) %>% 
  write.csv2("vPAINs_AS.csv", row.names = F, quote = F)

FreqHitSurf <- left_join(DataAS,DataSurf_sepRuns) %>% 
  filter(QAS==0,QSurf==0.5,QRunSurf==1) %>% 
  mutate(Val=ASVal-SurfVal) %>% 
  crossing(ntop=seq(10,NumDecoys,10)) %>% 
  group_by(ntop,target) %>% 
  do(ligand = GetTop(.$ntop,.$Val,.$ligand,.$type)) %>% 
  unnest(cols = c(ligand)) %>% 
  group_by(ntop,ligand) %>% 
  summarise(count=n()) %>% 
  arrange(desc(count))

FreqHitSurf %>% filter(ntop==20) %>%
  ungroup() %>% 
  select(ligand,count) %>% 
  write.csv2("vPAINs_OTD.csv", row.names = F, quote = F)

FreqHitDec1p <- FreqHitSurf %>% filter(ntop==20, count>=10) %>% .$ligand
FreqHitAS %>% filter(ntop==20,ligand%in%c(FreqHitDec1p))
FreqHitSurf %>% filter(ntop==20,ligand%in%c(FreqHitDec1p)) %>% 
  mutate(ligand=factor(ligand,levels=c(FreqHitDec1p))) %>% 
  arrange(ligand)

#####
#Analyse improvement sources
#####

ImprovementPlot <- GainTab %>% 
  rename(ID=target) %>% 
  left_join(TestSet %>% 
              filter(!ID%in%c("None","Self"))) %>%
  mutate(OPEN=ifelse(OpenAS==1,"Open AS","Closed AS"),
         UnitCharge=as.numeric(TotCharge)/as.numeric(NumStAl)) %>% 
  ggplot(aes(x=as.numeric(TotCharge),y=100*Gain,color=OPEN,group=OPEN)) +
  scale_color_brewer(palette="Set1", name="") +
  scale_x_continuous(name="Protein total charge at pH 7") +
  scale_y_continuous(name="Improvement (%)") +
  geom_smooth(span=3,level=0.75,aes(group=NA),color="grey25") +
  geom_point(size=3,aes(group=ID),stroke=1.5,shape=1) +
  theme_bw()
ImprovementPlot
ggsave("ImprovementPlot_TotCharge.png",ImprovementPlot,height=9,width=18,dpi=1000,units="cm")
#ggplotly(ImprovementPlot)

ImprovementPlot2 <- GainTab %>% 
  rename(ID=target) %>% 
  left_join(TestSet %>% 
              filter(!ID%in%c("None","Self"))) %>%
  mutate(OPEN=ifelse(OpenAS==1,"Open AS","Closed AS"),
         UnitCh=as.numeric(TotCharge)/as.numeric(NumStAl)) %>%
  ggplot(aes(x=as.numeric(TotCharge),y=100*Gain,shape=OPEN,color=None)) +
  scale_color_gradientn(colors=c("dark red","red","orange","yellow","dark green"), 
                        name="Conventional docking accuracy",
                        labels = scales::percent_format(accuracy = 1)) +
  scale_shape_manual(values=c(16,1),name="") +
  scale_x_continuous(name="Protein total charge at pH 7") +
  scale_y_continuous(name="Improvement (%)",
                     position = "left") +
  geom_smooth(span=3,level=0.75,aes(group=NA),color="grey25") +
  geom_point(size=3,aes(group=ID),stroke=1.2,alpha=0.8) +
  theme_bw() +
  geom_text(label="Ribonuclease A",x=11,y=46,size=3,hjust=1,color="#a23333") +
  geom_text(label="Thymidylate synthase",x=-17,y=38,size=3,hjust=0,color="#fc3333") +
  geom_text(label="Ribonuclease T1",x=-27,y=32,size=3,hjust=0,color="#ff7733") +
  geom_text(label="Poly(ADP-ribose) polymerase",x=15,y=24.58,size=3,hjust=1,color="#b03335") +
  theme(legend.position = "bottom",
        legend.box.margin = margin(-3.5,0,-2,-10,"mm"),
        legend.direction = "horizontal",
        legend.box = "horizontal") +
  guides(color=guide_colorbar(title.position = "top",
                              title.theme = element_text(size = 10,color="grey20"),
                              barwidth = unit(50,"mm"),order=2,
                              barheight = unit(2,"mm"),nbin=100,ticks.linewidth=1.2,
                              label.theme = element_text(size=8,color="grey30")),
         shape=guide_legend(nrow=2,order=1))
ImprovementPlot2
ggsave("ImprovementPlot_TotCharge_v2.png",ImprovementPlot2,height=8,width=9,dpi=1000,units="cm")

ROCsTab <- GainTab %>% 
  rename(ID=target) %>% 
  left_join(TestSet %>% 
              left_join(TestSet %>% 
                          filter(!ID%in%c("None","Self")))) %>%
  mutate(OPEN=ifelse(OpenAS==1,"Open","Closed"),
         #TotCh=TotCh/as.numeric(NumStAl),
         None_plot=ifelse(abs(Gain)>0.031,None+0.015,None),
         Self_plot=ifelse(abs(Gain)>0.031,Self-0.015,Self))
ROCsPlot <- ROCsTab %>% 
  ggplot(aes(group=OPEN)) +
  scale_color_manual(values=c("#377eb8","#e41a1c"), name="Docking type") +
  scale_x_continuous(name="Protein total charge at pH 7") +
  scale_y_continuous(name="VS accuracy (%)") +
  #scale_linetype(values=c("solid","dashed"),name="Active site type") +
  scale_alpha_manual(values=c(0.8,0.4),name="Active site type") +
  geom_point(data=. %>% rename(ASD=None,OTD=Self) %>% 
               gather("key","value",ASD,OTD),
               aes(x=as.numeric(TotCharge),y=100*value,color=key),size=2) +
  geom_segment(aes(x=as.numeric(TotCharge),xend=as.numeric(TotCharge),y=100*None,yend=100*Self,alpha=factor(OPEN)),
               #arrow=arrow(length = unit(0.14, "cm"),angle=15),
               lineend = "round",
               color="grey10",size=0.4) +
  theme_bw() +
  #geom_text(label="Ribonuclease A",x=0,y=47,size=1.7,hjust=1) +
  #geom_text(label="Thymidylate synthase",x=-34.5,y=60.5,size=1.7,hjust=0) +
  #geom_text(label="Ribonuclease T1",x=-35.5,y=65.2,size=1.7,hjust=0) +
  #geom_text(label="Poly(ADP-ribose) polymerase",x=-8,y=43.7,size=1.7,hjust=1) +
  theme(legend.position = "bottom",
        legend.box = "vertical",
        legend.key.width = unit(4,"mm"),
        legend.text = element_text(margin = margin(0,0,0,-2,"mm")),
        legend.box.margin = margin(2,0,0,0,"mm"),
        legend.margin = margin(-5,0,0,0,"mm"))
ROCsPlot
ggsave("ROCchangePlot_TotCharge.png",ROCsPlot,height=7,width=9,dpi=1000,units="cm")
#ggplotly(ImprovementPlot)

#####
#Study properties influence on improvement
#####

GainTab %>% 
  rename(ID=target) %>% 
  left_join(TestSet %>% 
              left_join(TestSet %>% 
                          filter(!ID%in%c("None","Self")))) %>%
  mutate(OPEN=ifelse(OpenAS==1,"Open AS","Closed AS")) %>% 
  filter(as.numeric(TotCharge)<=-15|as.numeric(TotCharge)>=0) %>% 
  summarise(MeanNone=mean(None),
            MeanSelf=mean(Self),
            MeanGain=mean(Gain),
            n=n())

ROCsTab %>% 
  filter(as.numeric(TotCharge)<=-14|as.numeric(TotCharge)>=0) %>% 
  filter(OPEN=="Open") %>%
  #filter(None<0.71) %>% 
  group_by(OPEN) %>% 
  summarise(None=mean(None)*100,
            Self=mean(Self)*100,
            Gain=mean(Gain)*100,
            count=n())

#####
#BEDROC and Enrichment factor analyses
#####

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

#####
#Each-on-each with separate treating of runs
#####

GetSurfRunVal_p2p <- compiler::cmpfun(function(df){
  df %>%
    crossing(QRunSurf=c(1)) %>%
    group_by(QRunSurf,ligand,grid) %>% 
    #partition(cluster=cl) %>% 
    summarize(SurfRunVal=GetSurfRun(QRunSurf,top.dG)) #%>% 
    #collect()
})

cluster_copy(c("GetSurfRunVal_p2p","GetSurfRun","Data3sampBasOTS"),cluster=cl)
DataSurf_sepRuns_p2p_0 <- Data3sampBasOTS %>% 
  group_by(target) %>% 
  do(SurfRunVal0=GetSurfRunVal_p2p(.)) %>% 
  ungroup() %>% 
  unnest(cols = c(SurfRunVal0))

GetSurfVal_p2p <- compiler::cmpfun(function(df){
  df %>% 
    crossing(QSurf=c(0.5)) %>%
    group_by(ligand,QSurf,QRunSurf) %>% 
    #partition(cluster=cl) %>%
    summarize(SurfVal=quantile(SurfRunVal,probs=unique(QSurf))) #%>% 
    #collect()
})

DataSurf_sepRuns_p2p_1 <- DataSurf_sepRuns_p2p_0 %>% 
  group_by(target) %>%
  do(SurfVal0=GetSurfVal_p2p(.)) %>% 
  ungroup() %>% 
  unnest(cols = c(SurfVal0))

left_join(DataAS,DataSurf_sepRuns_p2p_1 %>% rename(Baseline=target)) %>% 
  filter(QAS==0,QSurf==0.5) %>% 
  group_by(target,Baseline,QSurf,QAS,QRunSurf) %>% 
  summarise(rat=sum(str_detect(ligand,"decoy"))/n()) %>% 
  View()

DataOTD_p2p_sepRuns_doROC <- left_join(DataAS,DataSurf_sepRuns_p2p_1 %>% rename(Baseline=target)) %>% 
  filter(QAS==0,QSurf==0.5) %>% 
  group_by(target,Baseline,QSurf,QAS,QRunSurf) %>% 
  do(doROC=roc(.$type,.$ASVal-.$SurfVal,direction=">")) 

DataOTD_p2p_sepRuns <- DataOTD_p2p_sepRuns_doROC %>% 
  mutate(ROC=unlist(doROC,recursive=F)$auc) %>% 
  select(-doROC)
  
DataOTD_p2p_sepRuns %>% head(3)
DataOTD_Self_sepRuns %>% head(3)
DataAS_AUC %>% head(3)

DataOTD_grp2_sepRuns <- bind_rows(DataOTD_p2p_sepRuns %>% 
                                    mutate(QRunSurf=as.character(QRunSurf)),
                                   DataOTD_Self_sepRuns %>% 
                                    filter(QRunSurf==1),
                                   DataAS_AUC %>% 
                                    crossing(QRunSurf=c(1)) %>% 
                                    mutate(QRunSurf=as.character(QRunSurf),
                                           Baseline="None"))
DataOTD_grp2_sepRuns %>% 
  write.csv("DataOTD_grp2_sepRuns_Prod.csv", row.names = F, quote = F)

DataOTD_grp2_sepRuns_mean <- DataOTD_grp2_sepRuns %>% 
  group_by(QSurf,QAS,QRunSurf,Baseline) %>% 
  summarise(MeanROC=mean(ROC),
            MinROC=min(ROC),
            MedROC=median(ROC)) %>% 
  arrange(desc(MeanROC))

DataOTD_grp2_sepRuns_mean %>%
  filter(QSurf%in%c(0.5,NA),QAS==0,QRunSurf==1) %>% 
  arrange(desc(MeanROC)) %>% 
  View()

DOTD_grp2_sepRuns_Ord <- DataOTD_grp2_sepRuns_mean %>%
  filter(QSurf%in%c(0.5,NA),QAS==0,QRunSurf%in%c(1,NA)) %>% 
  mutate(Baseline=as.character(Baseline)) %>% 
  arrange(desc(MeanROC)) %>% 
  .$Baseline

DataOTD_plot_sepRuns_SI <- DataOTD_grp2_sepRuns %>% 
  filter(QSurf%in%c(0.5,NA),QAS==0,QRunSurf%in%c(1,0.5,0,"Arb")) %>% 
  ungroup() %>% 
  mutate(Baseline=factor(Baseline,levels=c(DOTD_grp2_sepRuns_Ord)),
         QRunSurf=factor(QRunSurf,levels=c(1,0.5,0,"Arb"))) %>%
  ggplot(aes(y=ROC,x=Baseline,color=target)) +
  geom_point(alpha=0.7) +
  geom_boxplot(notch = T,aes(group=Baseline),fill=NA,outlier.alpha=0) +
  geom_point(data=DataOTD_grp2_sepRuns_mean %>% 
               filter(QSurf%in%c(0.5,NA),QAS==0,QRunSurf%in%c(1,0.5,0,"Arb")) %>%
               ungroup() %>% 
               mutate(Baseline=factor(Baseline,levels=c(DOTD_grp2_sepRuns_Ord)),
                      QRunSurf=factor(QRunSurf,levels=c("Arb",0,0.5,1))) %>%
               rename(ROC=MeanROC),shape="—",size=5,color="red") +
  #facet_grid(QRunSurf~.,labeller = as_labeller(QNames)) +
  theme_bw() +
  scale_y_continuous(breaks=seq(0,1,0.1),name="Virtual screening accuracy (AUROC)") + 
  xlab("PDB ID of the protein, which surface was used in on-top docking") +
  theme(axis.text.x = element_text(angle = -45, vjust=1.0, hjust=0.0,
                                   colour=ifelse(unique(DataOTD_grp2_sepRuns_mean %>% 
                                                          filter(QSurf%in%c(0.5,NA),QAS==0,QRunSurf%in%c(1,NA)) %>% 
                                                          arrange(desc(MeanROC)) %>% 
                                                          .$Baseline)%in%c("None","Self","United"),"red","black")),
        legend.position="none",
        plot.margin = margin(5,10,5,5))
DataOTD_plot_sepRuns_SI
ggsave(filename="Data_p2p_sepRuns_Prod_SI.png",DataOTD_plot_sepRuns_SI, width = 18, height = 8.5, units = "cm", dpi=2000)
#ggplotly(DataOTD_plot_sepRuns_SI)

#####
# VS-score in Active Site
#####

DataVS <- read_delim("VSscore_data.csv",";")

DataVS_AS <- DataVS %>% 
  filter(grid==0) 

DataVS_Surf <- DataVS %>% 
  filter(grid!=0) %>% 
  select(-type)

DataVS_AS %>% 
  group_by(target,ligand,grid) %>% 
  summarise(c=n()) %>% 
  ungroup() %>% 
  distinct(c)

DataVS_Surf %>% 
  group_by(target,ligand,grid) %>% 
  summarise(c=n()) %>% 
  ungroup() %>% 
  distinct(c)

DataVS_AS2 <- DataVS_AS %>%   
  crossing(QAS=seq(0,1,0.5)) %>% 
  group_by(target,ligand,QAS,type) %>% 
  summarize(ASVal=quantile(top.VS,probs=unique(QAS)))

colnames(DataVS_AS2)
DataVS_AS2 %>% 
  write.csv("DataAS_VS_Prod.csv", row.names = F, quote = F)

DataVS_AS2 %>% 
  group_by(target,QAS) %>% 
  distinct(type) %>% 
  summarise(count=n()) %>% 
  arrange(desc(count))

DataVS_AS_AUC <- DataVS_AS2 %>%
  ungroup() %>% 
  #filter(target=="1efy",QAS=="0") %>% 
  group_by(target,QAS) %>% 
  summarise(ROC=as.numeric(roc(type,ASVal,direction=">")$auc)) 

DataVS_AS_AUC_mean <- DataVS_AS_AUC %>%
  group_by(QAS) %>% 
  summarise(ROC=mean(ROC)) %>% 
  arrange(desc(ROC))
DataVS_AS_AUC_mean

DataVS_AS_AUC_SEPmean <- DataVS_AS_AUC %>%
  left_join(TestSet %>% select(ID,VStrain) %>% 
              rename(target=ID)) %>% 
  group_by(QAS,VStrain) %>% 
  summarise(ROC=mean(ROC)) %>% 
  arrange(desc(ROC))
DataVS_AS_AUC_SEPmean


#####
# On-Top Docking with separate treating of runs using VS-score
#####

GetSurfRun <- function(QRS,Energies){
  Quant=unique(QRS)
  if(Quant=="Arb"){
    sample(Energies,1) %>% return()
  }
  else{
    quantile(Energies,probs=as.numeric(Quant)) %>% return()
  }
}

DataVS_Surf_sepRuns <- DataVS_Surf %>%
  crossing(QRunSurf=c(seq(0,1,0.5),"Arb")) %>%
  group_by(target,ligand,QRunSurf,grid) %>% 
  summarize(SurfRunVal=GetSurfRun(QRunSurf,top.VS)) %>% 
  ungroup() %>% 
  mutate(QRunSurf=case_when(QRunSurf==0~"Min",
                            QRunSurf==0.5~"Med",
                            QRunSurf==1~"Max",
                            QRunSurf=="Arb"~"Arb")) %>% 
  crossing(QSurf=seq(0,1,0.05)) %>% 
  group_by(target,ligand,QSurf,QRunSurf) %>%
  summarize(SurfVal=quantile(SurfRunVal,probs=unique(QSurf)))

colnames(DataVS_Surf_sepRuns)
DataVS_Surf_sepRuns %>% 
  write.csv("DataSurf_VS_sepRuns_Prod.csv", row.names = F, quote = F)

DataVS_Surf_sepRuns <- read.csv("DataSurf_VS_sepRuns_Prod.csv")

DataOTD_VS_Self_sepRuns <- left_join(DataVS_AS2,DataVS_Surf_sepRuns) %>% 
  filter(QAS==0) %>% 
  group_by(target,QAS,QSurf,QRunSurf) %>% 
  summarise(ROC=as.numeric(roc(type,ASVal-SurfVal,direction=">")$auc)) %>% 
  mutate(Baseline="Self", SF="VS")

DataOTD_VS_Self_sepRuns %>% 
  #filter(QAS==0,QSurf==0.5,QRunSurf=="Max") %>% 
  left_join(TestSet %>% select(ID,VStrain) %>% 
              rename(target=ID)) %>% 
  group_by(VStrain,QAS,QSurf,QRunSurf) %>% 
  summarise(ROC=mean(ROC)) %>% 
  arrange(desc(ROC)) %>% View()

PlotVS_Self_Means_sepRuns <- DataOTD_VS_Self_sepRuns %>% 
  mutate(QRunSurf=factor(QRunSurf,levels=c("Arb","Min","Med","Max"))) %>% 
  ggplot(aes(x=QSurf,y=ROC,color=factor(QRunSurf))) +
  geom_smooth(aes(group=paste0(QAS,QRunSurf)),size=0.5,level=0.75) +
  #scale_color_brewer(palette = "Set1", name="Quantile of the\nover-the-hood\nruns used") +
  scale_color_manual(values=c("#984ea3","#4daf4a","#377eb8","#e41a1c","#ff7f00"),
                     guide = guide_legend(reverse = TRUE), 
                     name="Quantile of the\nover-the-hood\nruns used") +
  #scale_linetype_discrete(name="Quantile of the\nunder-the-hood\nruns used") +
  #labs(title="Dependense of the areas Under the ROC Curves (AUC)",
  #     subtitle="on the statistical procedure used for the Over-The-Hood-corrected docking",
  #     caption="Gray areas denote 95% confidence intervals of mean curves") +
  theme_bw() +
  guides(color=F) +
  theme(legend.position="bottom",        
        text=element_text(size=7),
        #axis.title.y = element_text(hjust=1,size=7.2),
        axis.title.x = element_text(hjust=0.9),
        legend.box.spacing = unit(0,"mm"),
        legend.box.margin =  margin(0,0,0,0,unit="mm"),
        legend.key.size = unit(1,"mm")) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(breaks=seq(0,1,0.05)) +
  ylab("AUC value") + xlab("Quantile of the over-the-hood grids energies") 
PlotVS_Self_Means_sepRuns
ggsave(filename="PlotVS_Self_Quantiles_Prod.png",PlotVS_Self_Means_sepRuns, width = 5.2, height = 3.7, units = "cm", dpi=2000)

#####
#Compare dG-score and VS-score
#####

PlotBoth_lab <- c(
  `1` = "Targets used for\nVS-score parameterization",
  `0` = "Targets not used for\nVS-score parameterization",
  `dG` = "dG-score",
  `VS` = "VS-score"
)

DataDG_AS_AUC_SEPmean <- DataAS_AUC %>%
  left_join(TestSet %>% select(ID,VStrain) %>% 
              rename(target=ID)) %>% 
  group_by(QAS,VStrain) %>% 
  summarise(ROC=mean(ROC)) %>% 
  arrange(desc(ROC))
DataDG_AS_AUC_SEPmean

PlotBoth_Self_Means_sepRuns <- bind_rows(DataOTD_VS_Self_sepRuns %>% 
                                           mutate(QRunSurf=case_when(QRunSurf=="Min"~"Minimum",
                                                                     QRunSurf=="Med"~"Median",
                                                                     QRunSurf=="Max"~"Maximum",
                                                                     QRunSurf=="Arb"~"Arbitrary")),
                                         DataOTD_Self_sepRuns %>% 
                                           mutate(SF="dG") %>% 
                                           filter(QRunSurf%in%c(0,0.5,1,"Abs")) %>% 
                                           mutate(QRunSurf=case_when(QRunSurf==0~"Minimum",
                                                                     QRunSurf==0.5~"Median",
                                                                     QRunSurf==1~"Maximum",
                                                                     QRunSurf=="Arb"~"Arbitrary"))) %>% 
  mutate(QRunSurf=factor(QRunSurf,levels=c("Arbitrary","Minimum","Median","Maximum"))) %>% 
  left_join(TestSet %>% select(ID,VStrain) %>% 
              rename(target=ID)) %>% 
  ggplot(aes(x=QSurf,y=ROC,color=factor(QRunSurf))) +
  geom_hline(data=bind_rows(DataVS_AS_AUC_SEPmean %>% mutate(SF="VS"),
                            DataDG_AS_AUC_SEPmean %>% mutate(SF="dG")) %>%
               ungroup() %>% 
               mutate(QAS=case_when(QAS==0~"Minimum",
                                    QAS==0.5~"Median",
                                    QAS==1~"Maximum",
                                    QAS=="Arb"~"Arbitrary")) %>% 
               filter(QAS=="Minimum"),
             aes(yintercept=ROC),color="black") +
  geom_smooth(aes(group=paste0(QAS,QRunSurf)),size=0.5,level=0.75) +
  scale_color_manual(values=c("#984ea3","#4daf4a","#377eb8","#e41a1c","#ff7f00"),
                     guide = guide_legend(reverse = TRUE), 
                     name="Quantile of the on-top runs used") +
  theme_bw() +
  facet_grid(SF~VStrain,labeller = as_labeller(PlotBoth_lab)) +
  guides(color=F) +
  theme(legend.position="bottom",        
        legend.direction = "vertical",
        text=element_text(size=7),
        #axis.title.y = element_text(hjust=1,size=7.2),
        #axis.title.x = element_text(hjust=0.9),
        #panel.spacing.x = unit(6,"mm"),
        legend.box.spacing = unit(0,"mm"),
        legend.box.margin =  margin(0,0,0,0,unit="mm"),
        legend.key.size = unit(3,"mm")) +
  scale_x_continuous(expand=c(0,0),name="On-top scores quantile",breaks=c(0,0.25,0.5,0.75,1),labels=c(0,0.25,0.5,0.75,1)) +
  scale_y_continuous(breaks=seq(0,1,0.05),name="Average AUROC") +
  coord_cartesian(ylim=c(0.8,NA)) +
  ylab("Average AUC") + xlab("Over-the-hood scores quantile") 
PlotBoth_Self_Means_sepRuns
ggsave(filename="VS-vs-dG_Prod.png",PlotBoth_Self_Means_sepRuns, width = 8.5, height = 8, units = "cm", dpi=2000)

#####
# Separate proteins

PlotBoth_Self_Means_sepRuns_eachProt <- bind_rows(DataOTD_VS_Self_sepRuns %>% 
                                           mutate(QRunSurf=case_when(QRunSurf=="Min"~"Minimum",
                                                                     QRunSurf=="Med"~"Median",
                                                                     QRunSurf=="Max"~"Maximum",
                                                                     QRunSurf=="Arb"~"Arbitrary")),
                                         DataOTD_Self_sepRuns %>% 
                                           mutate(SF="dG") %>% 
                                           filter(QRunSurf%in%c(0,0.5,1,"Abs")) %>% 
                                           mutate(QRunSurf=case_when(QRunSurf==0~"Minimum",
                                                                     QRunSurf==0.5~"Median",
                                                                     QRunSurf==1~"Maximum",
                                                                     QRunSurf=="Arb"~"Arbitrary"))) %>% 
  mutate(QRunSurf=factor(QRunSurf,levels=c("Arbitrary","Minimum","Median","Maximum"))) %>% 
  left_join(TestSet %>% select(ID,VStrain) %>% 
              rename(target=ID)) %>% 
  filter(QRunSurf=="Maximum") %>% 
  ggplot(aes(x=QSurf,y=ROC,color=target)) +
  geom_hline(data=bind_rows(DataVS_AS_AUC_SEPmean %>% mutate(SF="VS"),
                            DataDG_AS_AUC_SEPmean %>% mutate(SF="dG")) %>%
               ungroup() %>% 
               mutate(QAS=case_when(QAS==0~"Minimum",
                                    QAS==0.5~"Median",
                                    QAS==1~"Maximum",
                                    QAS=="Arb"~"Arbitrary")) %>% 
               filter(QAS=="Minimum"),
             aes(yintercept=ROC),color="black") +
  geom_line(size=0.5,alpha=0.8) +
  scale_color_hue(name="Target") +
  theme_bw() +
  facet_grid(SF~VStrain,labeller = as_labeller(PlotBoth_lab)) +
  theme(legend.position="right",        
        legend.direction = "vertical",
        text=element_text(size=7),
        #axis.title.y = element_text(hjust=1,size=7.2),
        #axis.title.x = element_text(hjust=0.9),
        #panel.spacing.x = unit(6,"mm"),
        legend.box.spacing = unit(0,"mm"),
        legend.box.margin =  margin(0,0,0,0,unit="mm"),
        legend.key.size = unit(3,"mm")) +
  scale_x_continuous(expand=c(0,0),name="On-top scores quantile",breaks=c(0,0.25,0.5,0.75,1),labels=c(0,0.25,0.5,0.75,1)) +
  scale_y_continuous(breaks=seq(0,1,0.05),name="Average AUROC") +
  coord_cartesian(ylim=c(0.6,NA)) +
  ylab("Average AUC") + xlab("Over-the-hood scores quantile") 
PlotBoth_Self_Means_sepRuns_eachProt
ggsave(filename="VS-vs-dG_eachProt.png",PlotBoth_Self_Means_sepRuns_eachProt, width = 12, height = 8, units = "cm", dpi=2000)


#####
#VS-score vPAINS
#####

FreqHitAS_VS <- DataVS_AS2 %>% 
  ungroup() %>% 
  filter(QAS==0) %>%   
  mutate(ligand=as.character(ligand)) %>% 
  crossing(ntop=seq(10,NumDecoys,10)) %>% 
  group_by(ntop,target) %>% 
  do(ligand = GetTop(.$ntop,.$ASVal,.$ligand,.$type)) %>% 
  unnest(cols = c(ligand)) %>% 
  group_by(ntop,ligand) %>% 
  summarise(count=n()) %>% 
  arrange(desc(count))

FreqHitSurf_VS <- left_join(DataVS_AS2,DataVS_Surf_sepRuns) %>% 
  filter(QAS==0,QSurf==0.5,QRunSurf=="Max") %>% 
  mutate(Val=ASVal-SurfVal) %>% 
  crossing(ntop=seq(10,NumDecoys,10)) %>% 
  group_by(ntop,target) %>% 
  do(ligand = GetTop(.$ntop,.$Val,.$ligand,.$type)) %>% 
  unnest(cols = c(ligand)) %>% 
  group_by(ntop,ligand) %>% 
  summarise(count=n()) %>% 
  arrange(desc(count))

FreqHitDec1p_VS <- FreqHitAS_VS %>% filter(ntop==20, count>=10) %>% .$ligand
FreqHitAS_VS %>% filter(ntop==20,ligand%in%c(FreqHitDec1p_VS))
FreqHitSurf_VS %>% filter(ntop==20,ligand%in%c(FreqHitDec1p_VS)) %>% 
  mutate(ligand=factor(ligand,levels=c(FreqHitDec1p_VS))) %>% 
  arrange(ligand)

