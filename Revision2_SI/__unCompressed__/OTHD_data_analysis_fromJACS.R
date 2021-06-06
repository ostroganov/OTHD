library(tidyverse)
library(reshape2)
library(pROC)
library(readxl)
library(multidplyr)
library(doSNOW)
library(rPref)
#library(plotly) #required for ggplotly interactive plot

cl <- makeCluster(10)

#####
#Load data
#####

Data10sampAS <- read_delim(file="DataAS_10runs.csv",";")
Data10sampOTH <- read_delim(file="DataOTH_10runs.csv",";")
Data3sampBasOTH <- read_delim(file="DataBasOTH_3runs.csv",";")
TestSet <- read_xlsx("TestSet.xlsx")

SkipPept <- c("")
QNames <- c(
  `1` = "Maximum",
  `0.5` = "Median",
  `0` = "Minimum",
  `Arb` = "Arbitrary"
)

#Optionally loaded data:
#
#DataAS <- read_delim("DataAS_Prod.csv",",")
#DataSurf_sepRuns <- read_delim("DataSurf_sepRuns_Prod.csv",",")
#DataOTHD_grp1_sepRuns <- read.csv("DataOTHD_grp1_sepRuns_Prod.csv")
#DataOTHD_grp2_sepRuns <- read.csv("DataOTHD_grp2_sepRuns_Prod.csv")

scaleFUN <- function(x) sprintf("%.2f", x)

#####
#Active site
#####

DataAS <- Data10sampAS %>% 
  group_by(target,ligand) %>% 
  crossing(.,data_frame(QAS=c(0))) %>% 
  partition(target,ligand,QAS,type,cluster=cl) %>% 
  summarize(ASVal=quantile(top.dG,probs=unique(QAS))) %>% 
  collect()

colnames(DataAS)
DataAS %>% 
  write.csv("DataAS_Prod.csv", row.names = F, quote = F)

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

DataAS_ROC <- DataAS_doROC %>% 
  group_by(target,QAS) %>%
  mutate(sens=list(unlist(doROC,recursive=F)$sensitivities),
         spec=list(unlist(doROC,recursive=F)$specificities)) %>% 
  select(-doROC) %>% 
  unnest()

DataAS %>%
  filter(QAS==0) %>% 
  group_by(target) %>% 
  summarise(NumLig=sum(type==1),
            NumDec=sum(type==0)) %>% 
  View()

#####
#Over-the-Hood Docking
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
    partition(QRunSurf,ligand,grid, cluster=cl) %>% 
    summarize(SurfRunVal=GetSurfRun(QRunSurf,top.dG)) %>% 
    collect()
})

clusterExport(cl=cl, list("group_by","GetSurfRun","str_detect","paste0","ifelse","unique","sample","quantile","%>%"), envir=environment())
DataSurf_sepRuns_0 <- Data10sampOTH %>% 
  group_by(target) %>% 
  do(SurfRunVal0=GetSurfRunVal(.)) %>% 
  ungroup() %>% 
  unnest()

GetSurfVal <- compiler::cmpfun(function(df){
  df %>% 
    crossing(data_frame(QSurf=seq(0,1,0.05))) %>%
    partition(ligand,QSurf,QRunSurf, cluster=cl) %>%
    summarize(SurfVal=quantile(SurfRunVal,probs=unique(QSurf))) %>% 
    collect()
})

DataSurf_sepRuns <- DataSurf_sepRuns_0 %>% 
  group_by(target) %>%
  do(SurfVal0=GetSurfVal(.)) %>% 
  ungroup() %>% 
  unnest()

colnames(DataSurf_sepRuns)
DataSurf_sepRuns %>% 
  write.csv("DataSurf_sepRuns_Prod.csv", row.names = F, quote = F)

DataOTHD_Self_sepRuns_doROC <- left_join(DataAS,DataSurf_sepRuns) %>% 
  filter(QAS==0) %>% 
  group_by(target,QAS,QSurf,QRunSurf) %>% 
  do(doROC=roc(.$type,.$ASVal-.$SurfVal,direction=">")) 

DataOTHD_Self_sepRuns <- DataOTHD_Self_sepRuns_doROC %>% 
  mutate(ROC=unlist(doROC,recursive=F)$auc) %>% 
  select(-doROC) %>% 
  ungroup() %>% 
  mutate(Baseline="Self")

Plot_Self_Means_sepRuns_Full <- DataOTHD_Self_sepRuns %>% 
  filter(!target%in%c(SkipPept),
         QRunSurf!="Arb") %>% 
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
                     name="OTH-runs\nscores quantile") +
  theme_bw() +
  theme(legend.position="bottom",        
        legend.box.spacing = unit(0,"mm"),
        plot.margin = margin(1,5,0,1,unit="mm"),
        legend.box.margin =  margin(0,0,0,0,unit="mm"),
        legend.key.height = unit(1,"mm"),
        legend.key.width = unit(10,"mm")) +
  scale_x_continuous(expand=c(0,0)) +
  ylab("Average AUC") + xlab("OTH-sites scores quantile") 
Plot_Self_Means_sepRuns_Full
ggsave(filename="Plot_Self_Quantiles_SI_Prod.png",Plot_Self_Means_sepRuns_Full, width = 16, height = 10, units = "cm", dpi=2000)

Plot_Self_Means_sepRuns <- DataOTHD_Self_sepRuns %>% 
  filter(QRunSurf%in%c(0,0.5,1,"Arb")) %>% 
  mutate(QRunSurf=factor(QRunSurf,levels=c("Arb",0,0.5,1))) %>% 
  filter(!target%in%c(SkipPept)) %>% 
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
                     name="Quantile of the\nover-the-hood\nruns used") +
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
  ylab("Average AUC") + xlab("OTH scores quantile") 
Plot_Self_Means_sepRuns
ggsave(filename="Plot_Self_Quantiles_Prod.png",Plot_Self_Means_sepRuns, width = 4.3, height = 3.0, units = "cm", dpi=2000)

#####
#ROC curves comparison
#####

AS_ord <- DataAS_AUC %>% 
  filter(QAS==0) %>%
  arrange(ROC) %>% 
  .$target

DataOTHD_Self_sepRuns_ROC <- DataOTHD_Self_sepRuns_doROC %>% 
  group_by(target,QAS,QSurf,QRunSurf) %>% 
  mutate(sens=list(unlist(doROC,recursive=F)$sensitivities),
         spec=list(unlist(doROC,recursive=F)$specificities)) %>% 
  select(-doROC) %>% 
  unnest() %>% 
  mutate(Baseline="Self")

ROC_curves_Plot <- bind_rows(DataAS_ROC %>% mutate(Baseline="None"),
          DataOTHD_Self_sepRuns_ROC) %>% 
  ungroup() %>% 
  mutate(PDB=factor(target,levels=rev(AS_ord))) %>%  
  filter(QAS==0,QRunSurf%in%c(NA,1),QSurf%in%c(NA,0.5)) %>% 
  ggplot(aes(y=sens,x=1-spec,color=Baseline)) +
  scale_color_brewer(palette="Set1",
                     breaks=c("Self","None"),
                     labels=c("OTH-docking","AS-docking")) +
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
#Gain Plot
#####

DataOTHD_grp1_sepRuns <- bind_rows(DataOTHD_Self_sepRuns %>% 
                                     mutate(QRunSurf=case_when(QRunSurf==1~"Max",
                                                               QRunSurf==0.5~"Med",
                                                               QRunSurf==0~"Min",
                                                               TRUE~QRunSurf)) %>% 
                                     filter(QRunSurf%in%c("Max","Med","Min","Arb")),
                                   DataAS_AUC %>% 
                                     crossing(QRunSurf=c("Max","Med","Min","Arb")) %>% 
                                     mutate(Baseline="None"))

DataOTHD_grp1_sepRuns %>% head()
DataOTHD_grp1_sepRuns %>% 
  write.csv("DataOTHD_grp1_sepRuns_Prod.csv", row.names = F, quote = F)

OTHD_GainPlot_Prod <- DataOTHD_grp1_sepRuns %>% filter(Baseline%in%c("Self","None"),
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
  scale_y_continuous(breaks=seq(0,1,0.05), name="AUC", labels=scaleFUN) + 
  scale_x_discrete(limits=c("None","Self"),
                   labels=c("ASD","OTHD"),
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
OTHD_GainPlot_Prod
ggsave(filename="OTHD_GainPlot_Prod.png",OTHD_GainPlot_Prod, width = 2.9, height = 6, units = "cm", dpi=2000)

GainTab <- DataOTHD_grp1_sepRuns %>% filter(Baseline%in%c("Self","None"),
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
  arrange(desc(as.numeric(None))) %>% 
  format(digits=5) %>% 
  write.csv2("OTHD_GainTable_Prod.csv",row.names = F, quote = F)

DataOTHD_grp1_sepRuns %>% filter(Baseline%in%c("Self","None"),
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
#Analyse frequent hitters
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
  unnest() %>% 
  group_by(ntop,ligand) %>% 
  summarise(count=n()) %>% 
  arrange(desc(count))

FreqHitAS %>% filter(ntop==20) %>%
  ungroup() %>% 
  select(ligand,count) %>% 
  write.csv2("FrequentHitters_AS.csv", row.names = F, quote = F)

FreqHitSurf <- left_join(DataAS,DataSurf_sepRuns) %>% 
  filter(QAS==0,QSurf==0.5,QRunSurf==1) %>% 
  mutate(Val=ASVal-SurfVal) %>% 
  crossing(ntop=seq(10,NumDecoys,10)) %>% 
  group_by(ntop,target) %>% 
  do(ligand = GetTop(.$ntop,.$Val,.$ligand,.$type)) %>% 
  unnest() %>% 
  group_by(ntop,ligand) %>% 
  summarise(count=n()) %>% 
  arrange(desc(count))

FreqHitSurf %>% filter(ntop==20) %>%
  ungroup() %>% 
  select(ligand,count) %>% 
  write.csv2("FrequentHitters_OTHD.csv", row.names = F, quote = F)

FreqHitDec1p <- FreqHitAS %>% filter(ntop==10, count>=10) %>% .$ligand
FreqHitAS %>% filter(ntop==10,ligand%in%c(FreqHitDec1p))
FreqHitSurf %>% filter(ntop==10,ligand%in%c(FreqHitDec1p)) %>% 
  mutate(ligand=factor(ligand,levels=c(FreqHitDec1p))) %>% 
  arrange(ligand)

FreqHitIdeal <- DataAS %>%   
  ungroup() %>% 
  mutate(ligand=as.character(ligand)) %>% 
  filter(QAS==0) %>% 
  select(target,ligand,type) %>% 
  sample_n(nrow(.)) %>% 
  arrange(desc(type)) %>% 
  mutate(nr=row_number()) %>% 
  crossing(ntop=seq(10,NumDecoys,10)) %>% 
  group_by(ntop,target) %>% 
  do(ligand = GetTop(.$ntop,.$nr,.$ligand,.$type)) %>% 
  unnest() %>% 
  group_by(ntop,ligand) %>% 
  summarise(count=n()) %>% 
  arrange(desc(count))

FreqHitPlot <- bind_rows(FreqHitAS %>% mutate(Proc="ASD"),
                         FreqHitSurf %>% mutate(Proc="OTHD"),
                         FreqHitIdeal %>% mutate(Proc="Ideal")) %>% 
  mutate(Proc=factor(Proc,levels=c("ASD","OTHD","Ideal"))) %>% 
  filter(count>30) %>% 
  group_by(ntop,Proc) %>% 
  summarise(c=n()) %>% 
  ungroup() %>% 
  ggplot(aes(x=100*ntop/NumDecoys,y=100*c/NumDecoys,color=Proc)) +
  geom_line() +
  scale_x_continuous(name="Top-list size (% of ligands)") +
  scale_y_continuous(name="Percent of decoys, residing in the top\nof more than 30 targets (%)",
                     limits=c(NA,NA)) +
  scale_color_brewer(palette = "Set1", name="Procedure") +
  theme_bw()
FreqHitPlot
ggsave(filename="FrequentHittersPlot.png",FreqHitPlot, width = 16, height = 10, units = "cm", dpi=2000)

#####
#Analyse improvement sources
#####

ResidTab <- data_frame(Type=c("Asp", "Glu", "Cys", "Tyr", "His", "Lys", "Arg", "TerN", "TerC"),
                       Prot=c("a", "a", "a", "a", "c", "c", "c", "c", "a"))

TargFils <- list.files("./targets")
system("mkdir CalcDir")
for(i in TargFils) {
  system(paste0("cp .\\targets\\",i,"\\protein.pdb .\\CalcDir\\",i,".pdb"),
         intern=T)
}

CalcTotCh <- function(pH,file,resid=seq(1,10000)) {
  pH <- unique(pH)
  file <- unique(file)
  if(!file.exists(paste0(file,"-",pH,".log"))){
    system(paste0("build_model.exe -f ",file,".pdb -omm ",file,"-",pH,"-prep.pdb -mode normal -noligand -olog ",file,"-",pH,".log -pH ",pH),
                           intern=T)
    }
  lines <- readLines(paste0(file,"-",pH,".log"))
  lines[which(str_detect(lines,"Group\tType")):
          (which(str_detect(lines,"Build model"))-2)] %>%
    paste0("\n", collapse = "") %>% 
    read_delim("\t",trim_ws=T) %>% 
    rowwise() %>% 
    mutate(Type=ifelse(Type=="Ligand",str_to_title(str_trunc(ResName,3,"right",ellipsis = "")),Type)) %>% 
    ungroup() %>% 
    filter(Type%in%c(ResidTab %>% filter(!Type%in%c("TerN","TerC")) %>% 
                       distinct(Type) %>% unlist()),
           ResNum%in%c(resid)) %>% 
    left_join(ResidTab,by = "Type") %>% 
    mutate(Q=ifelse(Prot=="a",
                    -1/(1+10**(`pKa(calc)`-pH)),
                    1/(1+10**(pH-`pKa(calc)`)))) %>% 
    summarise(NQ=sum(Q)) %>% 
    .$NQ
}

clusterExport(cl=cl, list("mutate","CalcTotCh","system",
                          "paste0","ifelse","readLines",
                          "str_detect","str_to_title","%>%",
                          "str_trunc","read_delim","rowwise",
                          "filter","distinct","unlist","left_join",
                          "summarise","%in%","ungroup", "ResidTab"),
              envir=environment())

TestSet2 <- TestSet %>% 
  filter(!ID%in%c("None","Self")) %>% 
  crossing(pH=c(7.0)) %>% 
  partition(ID,pH,cluster = cl) %>% 
  do(TotCh=tryCatch(CalcTotCh(.$pH,paste0("./CalcDir/",.$ID)),finally = NA)) %>% 
  collect() %>% 
  unnest()

ImprovementPlot <- GainTab %>% 
  rename(ID=target) %>% 
  left_join(TestSet2 %>% 
              left_join(TestSet %>% 
                          select(ID,OpenAS,NumStAl))) %>%
  mutate(OPEN=ifelse(OpenAS==1,"Open AS","Closed AS"),
         TotCh=TotCh/as.numeric(NumStAl)) %>% 
  filter(pH==7) %>% 
  ggplot(aes(x=TotCh,y=100*Gain,color=OPEN,group=OPEN)) +
  scale_color_brewer(palette="Set1", name="") +
  scale_x_continuous(name="Protein total charge at pH 7") +
  scale_y_continuous(name="Improvement (%)") +
  geom_smooth(span=3,level=0.75,aes(group=NA),color="grey25") +
  geom_point(size=3,aes(group=ID),stroke=1.5,shape=1) +
  theme_bw()
ImprovementPlot
ggsave("ImprovementPlot_TotCharge.png",ImprovementPlot,height=9,width=18,dpi=1000,units="cm")
#ggplotly(ImprovementPlot)

ROCsPlot <- GainTab %>% 
  rename(ID=target) %>% 
  left_join(TestSet2 %>% 
              left_join(TestSet %>% 
                          select(ID,OpenAS,NumStAl))) %>%
  mutate(OPEN=ifelse(OpenAS==1,"Open AS","Closed AS"),
         TotCh=TotCh/as.numeric(NumStAl)) %>% 
  filter(pH==7) %>% 
  ggplot(aes(x=TotCh,xend=TotCh,y=100*None,yend=100*Self,color=factor(OPEN),group=OPEN)) +
  scale_color_brewer(palette="Set1", name="") +
  scale_x_continuous(name="Protein total charge at pH 7") +
  scale_y_continuous(name="VS performance (ASD-to-OTHD) (%)") +
  geom_segment(arrow=arrow(length = unit(0.17, "cm")),alpha=0.8,size=0.8) +
  theme_bw()
ROCsPlot
ggsave("ROCchangePlot_TotCharge.png",ROCsPlot,height=10,width=18,dpi=1000,units="cm")
#ggplotly(ImprovementPlot)

GainTab %>% 
  rename(ID=target) %>% 
  left_join(TestSet2 %>% 
              left_join(TestSet %>% 
                          select(ID,OpenAS,NumStAl))) %>%
  mutate(OPEN=ifelse(OpenAS==1,"Open AS","Closed AS"),
         TotCh=TotCh/as.numeric(NumStAl)) %>% 
  filter(pH==7) %>% 
  filter(TotCh<=-25|TotCh>=0) %>% 
  summarise(MeanNone=mean(None),
            MeanSelf=mean(Self),
            MeanROC=mean(Gain),
            n=n())

#####
#Each-on-each with separate treating of runs
#####

GetSurfRunVal_p2p <- compiler::cmpfun(function(df){
  df %>%
    crossing(data_frame(QRunSurf=c(seq(0,1,0.5),"Arb"))) %>%
    partition(QRunSurf,ligand,grid, cluster=cl) %>% 
    summarize(SurfRunVal=GetSurfRun(QRunSurf,top.dG)) %>% 
    collect()
})

clusterExport(cl=cl, list("group_by","GetSurfRun","str_detect","paste0","ifelse","unique","sample","quantile","%>%"), envir=environment())
DataSurf_sepRuns_p2p_0 <- Data3sampBasOTH %>% 
  group_by(target) %>% 
  do(SurfRunVal0=GetSurfRunVal_p2p(.)) %>% 
  ungroup() %>% 
  unnest()

GetSurfVal_p2p <- compiler::cmpfun(function(df){
  df %>% 
    crossing(data_frame(QSurf=c(0.25,0.5,0.75))) %>%
    partition(ligand,QSurf,QRunSurf, cluster=cl) %>%
    summarize(SurfVal=quantile(SurfRunVal,probs=unique(QSurf))) %>% 
    collect()
})

DataSurf_sepRuns_p2p_1 <- DataSurf_sepRuns_p2p_0 %>% 
  group_by(target) %>%
  do(SurfVal0=GetSurfVal_p2p(.)) %>% 
  ungroup() %>% 
  unnest()

DataOTHD_p2p_sepRuns <- left_join(DataAS,DataSurf_sepRuns_p2p_1 %>% rename(Baseline=target)) %>% 
  filter(QAS==0,QSurf==0.5) %>% 
  group_by(target,Baseline,QSurf,QAS,QRunSurf) %>% 
  summarise(ROC=roc(type,ASVal-SurfVal,direction=">")$auc)

DataOTHD_p2p_sepRuns %>% head(3)
DataOTHD_Self_sepRuns %>% head(3)
DataAS_AUC %>% head(3)

DataOTHD_grp2_sepRuns <- bind_rows(DataOTHD_p2p_sepRuns,
                                   DataOTHD_Self_sepRuns,
                                   DataAS_AUC %>% 
                                     crossing(QRunSurf=c(1,0.5,0,"Arb")) %>% 
                                     mutate(Baseline="None"))
DataOTHD_grp2_sepRuns %>% 
  write.csv("DataOTHD_grp2_sepRuns_Prod.csv", row.names = F, quote = F)

DataOTHD_grp2_sepRuns_mean <- DataOTHD_grp2_sepRuns %>% 
  filter(!target%in%c(SkipPept)) %>%
  group_by(QSurf,QAS,QRunSurf,Baseline) %>% 
  summarise(MeanROC=mean(ROC),
            MinROC=min(ROC),
            MedROC=median(ROC)) %>% 
  arrange(desc(MeanROC))

DataOTHD_grp2_sepRuns_mean %>%
  filter(QSurf%in%c(0.5,NA),QAS==0,QRunSurf==1) %>% 
  arrange(desc(MeanROC)) %>% 
  View()

DOTHD_grp2_sepRuns_Ord <- DataOTHD_grp2_sepRuns_mean %>%
  filter(QSurf%in%c(0.5,NA),QAS==0,QRunSurf%in%c(1,NA)) %>% 
  mutate(Baseline=as.character(Baseline)) %>% 
  arrange(desc(MeanROC)) %>% 
  .$Baseline

DataOTHD_plot_sepRuns_SI <- DataOTHD_grp2_sepRuns %>% 
  filter(!target%in%c(SkipPept)) %>%
  filter(QSurf%in%c(0.5,NA),QAS==0,QRunSurf%in%c(1,0.5,0,"Arb")) %>% 
  ungroup() %>% 
  mutate(Baseline=factor(Baseline,levels=c(DOTHD_grp2_sepRuns_Ord)),
         QRunSurf=factor(QRunSurf,levels=c(1,0.5,0,"Arb"))) %>%
  ggplot(aes(y=ROC,x=Baseline,color=target)) +
  geom_point(alpha=0.7) +
  geom_boxplot(notch = T,aes(group=Baseline),fill=NA,outlier.alpha=0) +
  geom_point(data=DataOTHD_grp2_sepRuns_mean %>% 
               filter(QSurf%in%c(0.5,NA),QAS==0,QRunSurf%in%c(1,0.5,0,"Arb")) %>%
               ungroup() %>% 
               mutate(Baseline=factor(Baseline,levels=c(DOTHD_grp2_sepRuns_Ord)),
                      QRunSurf=factor(QRunSurf,levels=c("Arb",0,0.5,1))) %>%
               rename(ROC=MeanROC),shape="—",size=5,color="red") +
  facet_grid(QRunSurf~.,labeller = as_labeller(QNames)) +
  theme_bw() +
  scale_y_continuous(breaks=seq(0,1,0.1),name="AUC") + 
  xlab("Protein, which exteriority was used in OTH-docking") +
  theme(axis.text.x = element_text(angle = -45, vjust=1.0, hjust=0.0,
                                   colour=ifelse(unique(DataOTHD_grp2_sepRuns_mean %>% 
                                                          filter(QSurf%in%c(0.5,NA),QAS==0,QRunSurf%in%c(1,NA)) %>% 
                                                          arrange(desc(MeanROC)) %>% 
                                                          .$Baseline)%in%c("None","Self","United"),"red","black")),
        legend.position="none")
DataOTHD_plot_sepRuns_SI
ggsave(filename="Data_p2p_sepRuns_Prod_SI.png",DataOTHD_plot_sepRuns_SI, width = 16.5, height = 20, units = "cm", dpi=2000)
#ggplotly(DataOTHD_plot_sepRuns_SI)

DrawTarg <- c("Self",
              DOTHD_grp2_sepRuns_Ord[2],
              DOTHD_grp2_sepRuns_Ord[length(DOTHD_grp2_sepRuns_Ord)-1],
              "None")

DataOTHD_grp2_sepRuns_Plot <- DataOTHD_grp2_sepRuns %>% 
  filter(QSurf%in%c(0.5,NA),QAS==0) %>% 
  filter(Baseline%in%c(DrawTarg)) %>%
  ungroup() %>%
  mutate(Baseline=case_when(Baseline=="Self"~"Self",
                            Baseline==DrawTarg[2]~"Best",
                            Baseline==DrawTarg[3]~"Worst",
                            Baseline=="None"~"None")) %>% 
  mutate(Baseline=factor(Baseline,levels=c("Self","Best","Worst","None")),
         QRunSurf=factor(QRunSurf,levels=c(1,0.5,0,"Arb")))

DataOTHD_grp2_sepRuns_mean_Plot <- DataOTHD_grp2_sepRuns_mean %>% 
  filter(QSurf%in%c(0.5,NA),QAS==0) %>%
  filter(Baseline%in%c(DrawTarg)) %>%
  ungroup() %>%
  mutate(Baseline=case_when(Baseline=="Self"~"Self",
                            Baseline==DrawTarg[2]~"Best",
                            Baseline==DrawTarg[3]~"Worst",
                            Baseline=="None"~"None")) %>% 
  mutate(Baseline=factor(Baseline,levels=c("Self","Best","Worst","None")),
         QRunSurf=factor(QRunSurf,levels=c("Arb",0,0.5,1))) %>%
  rename(ROC=MeanROC)

DataOTHD_plot_sepRuns <- ggplot(data=arrange(distinct(DataOTHD_grp2_sepRuns_Plot,Baseline),Baseline),
                                mapping=aes(y=ROC,x=Baseline)) +
  geom_jitter(data=filter(DataOTHD_grp2_sepRuns_Plot,Baseline!="None",QRunSurf==1),
              color="#e41a1c",
              width=0.3, alpha=0.3, size=0.5) +
  geom_jitter(data=filter(DataOTHD_grp2_sepRuns_Plot,Baseline=="None",QRunSurf=="Arb"),
              color="grey10",
              width=0.3, alpha=0.3, size=0.5) +
  geom_point(data=filter(DataOTHD_grp2_sepRuns_mean_Plot,Baseline!="None",QRunSurf==1),
             color="#e41a1c",
             shape="—",size=7) +
  geom_point(data=filter(DataOTHD_grp2_sepRuns_mean_Plot,Baseline=="None",QRunSurf=="Arb"),
             color="grey10",
             shape="—",size=7) +
  theme_bw() +
  #scale_color_manual(values=c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00"),
  #                   name="Quantile of the\nover-the-hood\nruns used") +
  scale_y_continuous(breaks=seq(0,1,0.2), name="AUC", labels=scaleFUN) + 
  scale_x_discrete(limits=c( "Self","Best","Worst","None"), 
                   name="Protein, which exteriority was used in OTH-docking") +
  theme(legend.position="none",
        text=element_text(size=7),
        plot.margin = margin(1,1,0,1,unit="mm"),
        panel.spacing = unit(0.5,"mm"),
        legend.box.spacing = unit(0,"mm"),
        legend.box.margin =  margin(0,0,0,0,unit="mm"),
        legend.key.size = unit(1,"mm"))
DataOTHD_plot_sepRuns
ggsave(filename="Data_p2p_sepRuns_Prod.png",DataOTHD_plot_sepRuns, width = 6.7, height = 5.1, units = "cm", dpi=2000)


#####
#Analysing AS procedure-optimization bootstrap results
#####

POB_AS <- read_delim("AS_ProcOpt_Bootstrap.csv",",") %>% 
  melt(id.vars=c("ASRuns","Try")) %>% 
  rename(MeanROC=value,QAS=variable) %>% 
  mutate(QAS=case_when(QAS=="Mean_ROC_min" ~ "Minimum",
                       QAS=="Mean_ROC_median" ~ "Median")) %>% 
  mutate(Baseline="None")

POB_AS %>% head()

POB_AS %>% distinct(Try)

POB_AS_meanSD <- POB_AS %>% 
  group_by(ASRuns,QAS,Baseline) %>% 
  summarise(ROCmean=mean(MeanROC),
            ROCsd=sd(MeanROC),
            ROCse=ROCsd/sqrt(n()),
            lQ=quantile(MeanROC,probs=0.25),
            hQ=quantile(MeanROC,probs=0.75)) %>% 
  ungroup()
POB_AS_meanSD %>% 
  filter(Baseline=="None") %>% 
  head()

POB_AS_PlotSD <- POB_AS_meanSD %>%
  filter(!is.na(QAS)) %>% 
  bind_rows(POB_AS_meanSD %>% 
              filter(is.na(QAS)) %>% 
              select(-QAS) %>% 
              crossing(distinct(POB_AS_meanSD,QAS))) %>%
  ggplot(aes(y=ROCmean,x=ASRuns,color=factor(QAS))) +
  scale_color_brewer(palette="Accent",name="AS score",
                     na.value="black",limits=c("Minimum","Median")) +
  scale_linetype_discrete(name="AS score",limits=c("Minimum","Median")) +
  geom_line(aes(group=QAS,
                linetype=factor(QAS)), color="grey40",size=0.2) +
  scale_x_continuous(name="Number of docking runs in AS",breaks=seq(1,101,2)) +
  #scale_x_log10(name="Number of docking runs in AS",breaks=2**seq(1,101,2)) +
  scale_y_continuous(breaks=seq(0,1,0.05),name="Average AUC", labels=scaleFUN) +
  coord_cartesian(ylim = c(0.77,0.92)) +
  #geom_errorbar(aes(ymin=ROCmean-ROCse, ymax=ROCmean+ROCse), width=.55, size=0.2, show.legend = F) +
  geom_errorbar(aes(ymin=ROCmean-ROCsd, ymax=ROCmean+ROCsd), width=.9, size=0.3, show.legend = F) +
  #geom_errorbar(aes(ymin=lQ, ymax=hQ), width=.55, size=0.2, show.legend = F) +
  geom_point(size=0.8) +
  theme_bw() +
  theme(text=element_text(size=7),
        plot.margin = margin(1,1,0,1,unit="mm"),
        axis.title.x = element_text(hjust = 0.75),
        axis.title.y = element_text(vjust = 0.25),
        legend.justification = c(0.1,0.65), 
        legend.position = c(0,1),
        legend.spacing = unit(2,"mm"),
        legend.box.just = "bottom",
        legend.title.align = 0.5,
        legend.margin = margin(0,0,0,0,unit="mm"),
        legend.box.spacing = unit(0,"mm"),
        legend.box.margin =  margin(5,0,0,2,unit="mm"),
        legend.key.width = unit(5,"mm"),
        legend.key.height = unit(0.05,"mm"),
        legend.text = element_text(size=5),
        legend.title = element_text(size=5),
        legend.background = element_rect(fill="white",
                                         size=0.5),
        legend.box="horizontal")
POB_AS_PlotSD
ggsave(filename="ProcOptAS_PlotSD_Prod.png",POB_AS_PlotSD, width = 4, height = 4, units = "cm", dpi=2000)

#####
#Analysing OTH ProcOpt bootstrap results
#####

POB_OTH <- read_delim("OTH_ProcOpt_Bootstrap.csv",",") %>% 
  mutate(Baseline="Self")

POB_OTH %>% colnames()

POB_OTH %>% distinct(Try)

POB_OTH_meanSD <- POB_OTH %>% 
  group_by(ASRuns,OTHRuns,OTHGrids) %>% 
  summarise(ROCmean=mean(MeanROC),
            ROCsd=sd(MeanROC),
            ROCse=ROCsd/sqrt(n()),
            lQ=quantile(MeanROC,probs=0.25),
            hQ=quantile(MeanROC,probs=0.75)) %>% 
  ungroup() %>% 
  bind_rows(POB_AS_meanSD %>% 
              filter(QAS%in%c(NA,"Minimum")) %>% 
              select(-QAS)) %>% 
  mutate(DockRuns=ASRuns+ifelse(is.na(OTHRuns),0,(OTHRuns*OTHGrids)))
POB_OTH_meanSD %>% 
  arrange(desc(ROCmean))

POB_OTH_meanSD %>% 
  arrange(desc(ROCmean)) %>%
  rename(MeanROC=ROCmean) %>% 
  write.table("OptProc_fullTable.csv",quote=F,sep=";",row.names=F)

POB_OTH_meanSD %>% group_by(Baseline) %>% 
  summarise(MinMROC=min(ROCmean),
            MaxMROC=max(ROCmean))

POB_pareto <- POB_OTH_meanSD %>% 
  psel(high(ROCmean)*low(DockRuns))

POB_pareto %>% 
  arrange(desc(ROCmean)) %>%
  rename(MeanROC=ROCmean) %>% 
  write.table("OptProc_ParetoFront.csv",quote=F,sep=";",row.names=F)

POB_PlotSD <- 
  bind_rows(filter(POB_pareto,
                   ASRuns%in%c(1,2,3,4,5,7,9,10),
                   OTHRuns%in%c(NA,1,3,5,9),
                   OTHGrids%in%c(NA,3,5,7,9,11,15)),
            filter(POB_pareto,
                   ASRuns%in%c(10),
                   OTHRuns%in%c(10),
                   OTHGrids%in%c(25))) %>% 
  distinct() %>% 
  ggplot(aes(y=ROCmean,x=DockRuns,color=factor(OTHGrids),shape=factor(OTHRuns),fill=factor(ASRuns))) +
  scale_color_brewer(palette="Dark2",name="# of OTH sites",na.value="black") +
  scale_shape_manual(name="# of runs in\neach OTH site",values=c(21,24,23,25,22),na.value=10) +
  scale_fill_brewer(palette="Dark2",name="# of runs in\nactive site",na.value=NA) +
  #scale_x_continuous(name="Total number of docking runs",breaks=seq(0,1000,10)) +
  scale_x_log10(name="Total number of docking runs",breaks=2**seq(0,1000,1)) +
  scale_y_continuous(breaks=seq(0,1,0.05),name="Average AUC", labels=scaleFUN) +
  coord_cartesian(ylim = c(0.77,0.92)) +
  #geom_errorbar(aes(ymin=ROCmean-ROCse, ymax=ROCmean+ROCse), width=.5, size=0.2, show.legend = F) +
  geom_errorbar(aes(ymin=ROCmean-ROCsd, ymax=ROCmean+ROCsd), width=.05, size=0.3, show.legend = F) +
  #geom_errorbar(aes(ymin=lQ, ymax=hQ), width=.05, size=0.2, show.legend = F) +
  geom_point(size=1.5) +
  theme_bw() +
  guides(color = guide_legend(order=2),
         fill = guide_legend(order=3,override.aes=list(shape=21)),
         shape = guide_legend(order=1),
         errorbar=F) +
  theme(text=element_text(size=7),
        plot.margin = margin(1,1,0,1,unit="mm"),
        axis.title.y = element_text(vjust = 0.25),
        legend.justification = c(0.99, 0.00), 
        legend.position = c(0.99, 0.01),
        panel.spacing = unit(0.5,"mm"),
        legend.spacing = unit(2,"mm"),
        legend.box.just = "bottom",
        legend.title.align = 0.5,
        legend.margin = margin(0,0,0,0,unit="mm"),
        legend.box.spacing = unit(0,"mm"),
        legend.box.margin =  margin(5,0,0,2,unit="mm"),
        legend.key.width = unit(5,"mm"),
        legend.key.height = unit(0.05,"mm"),
        legend.text = element_text(size=5),
        legend.title = element_text(size=5),
        legend.background = element_rect(fill="white",
                                         size=0.5),
        legend.box="horizontal")
POB_PlotSD
ggsave(filename="POB_PlotSD_Prod.png",POB_PlotSD, width = 13.81, height = 4, units = "cm", dpi=2000)
