library(tidyverse)
library(reshape2)
library(pROC)
library(readxl)
library(multidplyr)
library(doSNOW)
#library(plotly) #required for ggplotly interactive plot

cl <- makeCluster(10)

#####
#Load data
#####

Data10sampAS <- read_delim(file="DataAS_10runs.csv",";")
Data10sampOTS <- read_delim(file="DataOTS_10runs.csv",";")
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
#DataOTD_grp1_sepRuns <- read.csv("DataOTD_grp1_sepRuns_Prod.csv")

scaleFUN <- function(x) sprintf("%.2f", x)

#####
#Active site
#####

DataAS <- Data10sampAS %>%
  group_by(target,ligand,type) %>% 
  summarise(top.dG=min(top.dG)) %>% 
  group_by(target) %>% 
  mutate(C0=exp(-1000*top.dG/(300*1.987)),
         C=C0/(sum(C0)+1)) %>%
  group_by(target,ligand,type) %>% 
  summarise(C=sum(C)) %>% 
  arrange(target,desc(C))

DataAS_doROC <- DataAS %>% 
  group_by(target) %>% 
  do(doROC=roc(.$type,.$C,direction="<"))

DataAS_AUC <- DataAS_doROC %>% 
  mutate(ROC=unlist(doROC,recursive=F)$auc) %>% 
  select(-doROC) %>% 
  ungroup()

DataAS_AUC_mean <- DataAS_AUC %>%
  summarise(ROC=mean(ROC)) %>% 
  arrange(desc(ROC))
DataAS_AUC_mean

DataAS_ROC <- DataAS_doROC %>% 
  group_by(target) %>%
  mutate(sens=list(unlist(doROC,recursive=F)$sensitivities),
         spec=list(unlist(doROC,recursive=F)$specificities)) %>% 
  select(-doROC) %>% 
  unnest()

DataAS %>%
  group_by(target) %>% 
  summarise(NumLig=sum(type==1),
            NumDec=sum(type==0)) %>% 
  View()

#####
#On-top Docking
#####

GetSurfRun <- compiler::cmpfun(function(QRS,Energies){
  Quant=unique(QRS)
  if(Quant=="All"){
    return(Energies)
  }
  else{
    return(quantile(Energies,probs=as.numeric(Quant)))
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
DataSurf_sepRuns_0 <- Data10sampOTS %>% View()
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

DataSurf_sepRuns <- read.csv("DataSurf_sepRuns_Prod.csv")

DataSurf_sepRuns %>% head()
DataSurf_sepRuns %>% distinct(QSurf)
DataSurf_sepRuns %>% distinct(QRunSurf)

#####
#LOOK
#####

DataOTD_Self_sepRuns_BT <- left_join(DataAS,DataSurf_sepRuns) %>% 
  filter(QAS==0)

Plot1 <- DataOTD_Self_sepRuns_BT %>% 
  ungroup() %>% 
  filter(QRunSurf!="Arb",
         ASVal < 5,
         SurfVal < 5) %>% 
  ggplot(aes(y=ASVal,x=SurfVal,color=factor(type))) +
  facet_grid(factor(QSurf)~factor(QRunSurf)) +
  scale_color_brewer(palette="Set1") +
  #coord_cartesian(xlim=c(NA,5),ylim=c(NA,5)) +
  #geom_point(size=0.8) +
  geom_density2d() +
  theme_bw()
Plot1
ggsave(filename="Plot1.png",Plot1, width = 30, height = 30, units = "cm", dpi=1000)

DataOTD_Self_sepRuns_BT2 <- DataOTD_Self_sepRuns_BT %>% 
  ungroup() %>% 
  filter(QSurf==0.5) %>% 
  crossing(Int=seq(-2,2,0.1)) %>% 
  #filter(target=="1b5j",ligand=="1b5j/active-10",QAS==0,QSurf==0,QRunSurf==0.1) %>% 
  group_by(target,QAS,QSurf,QRunSurf,Int) %>% 
  mutate(C0_AS=exp(-1000*ASVal/(300*1.987)),
         C0_Surf=exp(-1000*(SurfVal-Int)/(300*1.987)),
         C_AS=C0_AS/(C0_AS+C0_Surf+1),
         Contrib=C_AS/sum(C_AS))

####

DataOTD_Self_sepRuns_doROC <- DataOTD_Self_sepRuns_BT2 %>% 
  group_by(target,QAS,QSurf,QRunSurf,Int) %>% 
  do(doROC=roc(.$type,.$Contrib,direction="<")) 

DataOTD_Self_sepRuns <- DataOTD_Self_sepRuns_doROC %>% 
  mutate(ROC=unlist(doROC,recursive=F)$auc) %>% 
  select(-doROC) %>% 
  ungroup() %>% 
  mutate(Baseline="Self") 

Plot_Self_Means_sepRuns_Full <- DataOTD_Self_sepRuns %>% 
  filter(!target%in%c(SkipPept),
         QRunSurf!="Arb") %>% 
  ggplot(aes(x=Int,y=ROC,color=factor(QRunSurf))) +
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
  scale_y_continuous(breaks=c(0.70,0.80,0.90,0.92,0.93,0.94,0.95,0.96)) +
  ylab("Average AUC") + xlab("Value added to the on-top binding energy") 
Plot_Self_Means_sepRuns_Full
ggsave(filename="Plot_Self_Quantiles-Int_SI_Prod.png",Plot_Self_Means_sepRuns_Full, width = 16, height = 10, units = "cm", dpi=2000)

Plot_Self_sepRuns_Full <- DataOTD_Self_sepRuns %>% 
  filter(!target%in%c(SkipPept),
         QRunSurf!="Arb",
         QRunSurf==1) %>% 
  ggplot(aes(x=Int,y=ROC,color=factor(target))) +
  geom_hline(data=DataAS_AUC_mean %>% mutate(SF="dG") %>% 
               filter(QAS==0) %>% 
               mutate(QAS=case_when(QAS==0~"Min",
                                    QAS==0.5~"Med",
                                    QAS==1~"Max")) %>% filter(QAS=="Min"),
             aes(yintercept=ROC),color="grey10",size=0.5) +
  #geom_smooth(aes(group=paste0(QAS,QRunSurf)),size=0.5,level=0.75,se=F) +
  geom_line(aes(group=paste0(QAS,QRunSurf,target)),size=0.5) +
  #scale_color_manual(values=c("#b15928","#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", 
  #                            "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a"),
  #                   name="on-top runs\nscores quantile") +
  theme_bw() +
  theme(legend.position="bottom",        
        legend.box.spacing = unit(0,"mm"),
        plot.margin = margin(1,5,0,1,unit="mm"),
        legend.box.margin =  margin(0,0,0,0,unit="mm"),
        legend.key.height = unit(1,"mm"),
        legend.key.width = unit(10,"mm")) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(breaks=c(0.70,0.80,0.90,0.92,0.93,0.94,0.95,0.96)) +
  ylab("Average AUC") + xlab("Value added to the on-top binding energy") 
Plot_Self_sepRuns_Full
ggsave(filename="Plot_Self_Quantiles-Int2_SI_Prod.png",Plot_Self_sepRuns_Full, width = 16, height = 10, units = "cm", dpi=2000)


Plot_Self_Means_sepRuns <- DataOTD_Self_sepRuns %>% 
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
  ylab("Average AUC") + xlab("On-top scores quantile") 
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
  group_by(target,QAS,QSurf,QRunSurf,Int) %>% 
  mutate(sens=list(unlist(doROC,recursive=F)$sensitivities),
         spec=list(unlist(doROC,recursive=F)$specificities)) %>% 
  select(-doROC) %>% 
  unnest() %>% 
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
#Gain Plot
#####

DataOTD_grp1_sepRuns <- bind_rows(DataOTD_Self_sepRuns %>% 
                                     mutate(QRunSurf=case_when(QRunSurf==1~"Max",
                                                               QRunSurf==0.5~"Med",
                                                               QRunSurf==0~"Min")) %>% 
                                     filter(QRunSurf%in%c("Max","Med","Min","Arb"),
                                            Int==0),
                                   DataAS_AUC %>% 
                                     crossing(QRunSurf=c("Max","Med","Min","Arb")) %>% 
                                     mutate(Baseline="None"))

DataOTD_grp1_sepRuns %>% head()
DataOTD_grp1_sepRuns %>% 
  write.csv("DataOTD_grp1_sepRuns_Prod.csv", row.names = F, quote = F)

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
  scale_y_continuous(breaks=seq(0,1,0.05), name="AUC", labels=scaleFUN) + 
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
  write.csv2("FrequentHitters_OTD.csv", row.names = F, quote = F)

FreqHitDec1p <- FreqHitAS %>% filter(ntop==20, count>=10) %>% .$ligand
FreqHitAS %>% filter(ntop==20,ligand%in%c(FreqHitDec1p))
FreqHitSurf %>% filter(ntop==20,ligand%in%c(FreqHitDec1p)) %>% 
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
                         FreqHitSurf %>% mutate(Proc="OTD"),
                         FreqHitIdeal %>% mutate(Proc="Ideal")) %>% 
  mutate(Proc=factor(Proc,levels=c("ASD","OTD","Ideal"))) %>% 
  filter(count>30) %>% 
  group_by(ntop,Proc) %>% 
  summarise(c=n()) %>% 
  ungroup() %>% 
  ggplot(aes(x=100*ntop/NumDecoys,y=100*c/NumDecoys,color=Proc)) +
  geom_line() +
  scale_x_continuous(name="Top-list size (% of ligands)") +
  scale_y_continuous(name="Percent of decoys, residing in the\ntop of more than 30 targets (%)",
                     limits=c(NA,NA)) +
  scale_color_brewer(palette = "Set1", name="Procedure",labels=c("Conventional docking","On-top docking","Ideal case")) +
  theme_bw() +
  theme(legend.position = c(0.3,0.75))
FreqHitPlot
ggsave(filename="FrequentHittersPlot.png",FreqHitPlot, width = 9, height = 7, units = "cm", dpi=2000)

FreqHitPlot_SI <- bind_rows(FreqHitAS %>% mutate(Proc="ASD"),
                         FreqHitSurf %>% mutate(Proc="OTD"),
                         FreqHitIdeal %>% mutate(Proc="Ideal")) %>% 
  mutate(Proc=factor(Proc,levels=c("ASD","OTD","Ideal"))) %>% 
  filter(count>30) %>% 
  group_by(ntop,Proc) %>% 
  summarise(c=n()) %>% 
  ungroup() %>% 
  ggplot(aes(x=100*ntop/NumDecoys,y=100*c/NumDecoys,color=Proc)) +
  geom_line(size=1.5,lineend="round") +
  scale_x_continuous(name="Top-list size (% of ligands)") +
  scale_y_continuous(name="Percent of decoys, residing in the\ntop of more than 30 targets (%)",
                     limits=c(NA,NA)) +
  scale_color_brewer(palette = "Set1", name="Procedure",labels=c("Conventional docking","On-top docking","Ideal case")) +
  theme_bw() +
  theme(legend.position = c(0.3,0.75))
FreqHitPlot_SI
ggsave(filename="FrequentHittersPlot_SI.png",FreqHitPlot_SI, width = 11, height = 12, units = "cm", dpi=2000)


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

ImprovementPlot2 <- GainTab %>% 
  rename(ID=target) %>% 
  left_join(TestSet2 %>% 
              left_join(TestSet %>% 
                          select(ID,OpenAS,NumStAl))) %>%
  mutate(OPEN=ifelse(OpenAS==1,"Open AS","Closed AS"),
         TotCh=TotCh/as.numeric(NumStAl)) %>% 
  filter(pH==7.0) %>% 
  ggplot(aes(x=TotCh,y=100*Gain,shape=OPEN,color=None)) +
  scale_color_gradientn(colors=c("dark red","red","orange","yellow","dark green"), 
                        name="Conventional docking accuracy",
                        labels = scales::percent_format(accuracy = 1)) +
  scale_shape_manual(values=c(16,1),name="") +
  scale_x_continuous(name="Protein total charge at pH 7") +
  scale_y_continuous(name="Improvement (%)",
                     #breaks=c(0,10,20,30,40,50),
                     #labels=c(0,10,"",30,"",50),
                     position = "left") +
  geom_smooth(span=3,level=0.75,aes(group=NA),color="grey25") +
  geom_point(size=3,aes(group=ID),stroke=1.2,alpha=0.8) +
  theme_bw() +
  geom_text(label="Ribonuclease A",x=-0.5,y=46,size=3,hjust=1,color="#a23333") +
  geom_text(label="Thymidylate synthase",x=-34,y=38,size=3,hjust=0,color="#fc3333") +
  geom_text(label="Ribonuclease T1",x=-35.5,y=32,size=3,hjust=0,color="#ff7733") +
  geom_text(label="Poly(ADP-ribose) polymerase",x=-0.2,y=24.58,size=3,hjust=1,color="#b03335") +
  theme(legend.position = "bottom",
        #legend.margin = margin(0,0,-5,0,"mm"),
        #axis.text.y = element_text(color=c("grey30","grey50","grey30","grey30","grey30","grey30")),
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
  left_join(TestSet2 %>% 
              left_join(TestSet %>% 
                          select(ID,OpenAS,NumStAl))) %>%
  mutate(OPEN=ifelse(OpenAS==1,"Open","Closed"),
         TotCh=TotCh/as.numeric(NumStAl),
         None_plot=ifelse(abs(Gain)>0.031,None+0.015,None),
         Self_plot=ifelse(abs(Gain)>0.031,Self-0.015,Self)) %>% 
  filter(pH==7)
ROCsPlot <- ROCsTab %>% 
  ggplot(aes(group=OPEN)) +
  scale_color_manual(values=c("#377eb8","#e41a1c"), name="Docking type") +
  scale_x_continuous(name="Protein total charge at pH 7") +
  scale_y_continuous(name="VS accuracy (%)") +
  #scale_linetype(values=c("solid","dashed"),name="Active site type") +
  scale_alpha_manual(values=c(0.8,0.4),name="Active site type") +
  geom_point(data=. %>% rename(ASD=None,OTD=Self) %>% 
               gather("key","value",ASD,OTD),
               aes(x=TotCh,y=100*value,color=key),size=2) +
  geom_segment(aes(x=TotCh,xend=TotCh,y=100*None,yend=100*Self,alpha=factor(OPEN)),
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
  left_join(TestSet2 %>% 
              left_join(TestSet %>% 
                          select(ID,OpenAS,NumStAl))) %>%
  mutate(OPEN=ifelse(OpenAS==1,"Open AS","Closed AS"),
         TotCh=TotCh/as.numeric(NumStAl)) %>% 
  filter(pH==7) %>% 
  filter(TotCh<=-25|TotCh>=0) %>% 
  summarise(MeanNone=mean(None),
            MeanSelf=mean(Self),
            MeanGain=mean(Gain),
            n=n())

ROCsTab %>% 
  filter(TotCh > -15|TotCh < -25) %>% 
  filter(OPEN=="Open",None<0.71) %>% 
  group_by(OPEN) %>% 
  summarise(None=mean(None)*100,
            Self=mean(Self)*100,
            Gain=mean(Gain)*100)
