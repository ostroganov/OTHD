library(tidyverse)
library(reshape2)
library(pROC)
library(multidplyr)
library(hablar)

Program <- "GOLD"
TSet <- "DUDEp_OrigLigs"
mainSF="plp"

CurrentDir <- "GOLD_DUDEp_OrigLigs_runs123"
PreparedDataFile <- "GOLD_DUDEp_OrigLigs_runs123.csv"

##>>Data preparation for GOLD----

read_gold <- function(filename, dir ){
  metadata <- filename %>% str_sub(1,-5) %>% str_split("_") %>% unlist()
  path <- paste0(dir, "/", filename)
  read_csv(path, col_types="cdddddddddd") %>%
    mutate(target = metadata[1],
           SF = metadata[2],
           Grid = ifelse(metadata[4]=="ref",0,as.numeric(metadata[4])+1),
           ligand_name = paste(target, ligand_name, sep = "_")) %>%
    rename(Ligand= ligand_name) %>%
    gather(Run, Score, starts_with("dock"))
}
read_gold("1b5j_plp_ac_ref.csv",paste0(CurrentDir,"/actives"))

read_gold_dir <- function(dir){
  dir.act <- paste0(dir, "/actives")
  dir.dec <- paste0(dir, "/decoys")
  files.act <- dir(dir.act, "[.]csv$")
  files.dec <- dir(dir.dec, "[.]csv$")
  actives <- map(files.act, read_gold, dir = dir.act) %>%
    bind_rows() %>%
    mutate(Type = 1,
           Ligand = paste0(Ligand, "_ac"))
  decoys <- map(files.dec, read_gold, dir = dir.dec) %>%
    bind_rows() %>%
    mutate(Type = 0,
           Ligand = paste0(Ligand, "_dec"))
  bind_rows(actives, decoys)
}

full_data <- read_gold_dir(CurrentDir)

full_data_prep <- full_data %>% 
  filter(!str_detect(Ligand,"ref")) %>% 
  mutate(Run=as.numeric(str_extract(Run,"[0-9]+")),
         Score=-Score) %>% 
  spread(SF,Score)

full_data_prep %>% View()

full_data_prep %>% write_delim(PreparedDataFile,";")

full_data_prep %>% colnames()

##>>Data preparation for LeadFinder----

ResDir <- paste0(CurrentDir,"/results")
targets <- dir(ResDir)
compl_targ_l <- map(paste0(ResDir,"/",targets),dir) %>% map_int(length)
compl_targ_l
compl_targ_l[compl_targ_l>500] %>% unique()
compl_targ <- targets[compl_targ_l>500]

#targets <- c("ace","ada")

read_tsf <- function(fname,targ) {
  print(paste(fname,targ,sep = " "))
  res <- read_delim(paste0(ResDir,"/",targ,"/",fname,".tsv"),"\t",trim_ws = T,skip = 4,
                    col_types = c(
                      Ligand = col_character(),
                      `SDF molecule id` = col_integer(),
                      Rank = col_character(),
                      `dG, kcal/mol` = col_double(),
                      `VS score` = col_double(),
                      `Rank score` = col_double(),
                      `RMS, A` = col_double(),
                      `Docking time` = col_double(),
                      `Ligand Efficiency (dG/Nheavy)` = col_double(),
                      `Molar weight` = col_double(),
                      `N(heavy atoms)` = col_integer(),
                      `N(FRBs)` = col_integer(),
                      `E(LJ)` = col_double(),
                      `E(metal)` = col_double(),
                      `E(me-penalty)` = col_double(),
                      `E(sol)` = col_double(),
                      `E(surf,polar,prot)` = col_double(),
                      `E(surf,polar,solv)` = col_double(),
                      `E(surf,nonpolar,prot)` = col_double(),
                      `E(surf,nonpolar,solv)` = col_double(),
                      `E(HB-polar)` = col_double(),
                      `E(HB-nonpolar)` = col_double(),
                      `E(HB-penalty)` = col_double(),
                      `E(HB-const)` = col_double(),
                      `E(HB-extra)` = col_double(),
                      `E(coul-buried)` = col_double(),
                      `E(coul-intermediate)` = col_double(),
                      `E(coul-surface)` = col_double(),
                      `E(Born-ligand)` = col_double(),
                      `E(Born-ligand-groups)` = col_double(),
                      `E(Born-protein)` = col_double(),
                      `E(internal attraction)` = col_double(),
                      `E(internal repulsion)` = col_double(),
                      `E(internal water)` = col_double(),
                      `E(internal-14)` = col_double(),
                      `E(covalent-bond)` = col_double(),
                      `E(covalent-angle)` = col_double(),
                      `E(tors)` = col_double(),
                      `E(dihedral)` = col_double(),
                      `E(constaints)` = col_double(),
                      `E(restraints)` = col_double(),
                      `E(penalty)` = col_double(),
                      `mol header` = col_character(),
                      `mol comment` = col_character())) 
  dim(res)
  #print(res)
  res2 <- res %>% 
    select(-`SDF molecule id`) %>% 
    filter(Rank==0) %>%  
    mutate(Rank=as.numeric(Rank),
           `dG, kcal/mol`=as.numeric(`dG, kcal/mol`),
           `VS score`=as.numeric(`VS score`),
           tsf_name=fname)
  print(dim(res2)[1])
  #print(res2)
  res2 %>% 
    convert(num(`E(LJ)`,`E(metal)`,`E(me-penalty)`,`E(sol)`,`E(surf,polar,prot)`,`E(surf,polar,solv)`,`E(surf,nonpolar,prot)`,`E(surf,nonpolar,solv)`,`E(HB-polar)`,`E(HB-nonpolar)`,`E(HB-penalty)`,`E(HB-const)`,`E(HB-extra)`,`E(coul-buried)`,`E(coul-intermediate)`,`E(coul-surface)`,`E(Born-ligand)`,`E(Born-ligand-groups)`,`E(Born-protein)`,`E(internal attraction)`,`E(internal repulsion)`,`E(internal water)`,`E(internal-14)`,`E(covalent-bond)`,`E(covalent-angle)`,`E(tors)`,`E(dihedral)`,`E(constaints)`,`E(restraints)`,`E(penalty)`,`dG, kcal/mol`,`VS score`,`Rank score`,`RMS, A`,`Docking time`,`Ligand Efficiency (dG/Nheavy)`,`Molar weight`),
            chr(Ligand,Rank,`mol header`,`mol comment`),
            int(`N(heavy atoms)`,`N(FRBs)`))
  #return(res2)
}

read_tsf("asite-actives-docked-3","1kim")

read_dir_tsf <- function(targ) {
  files <- dir(paste0(ResDir,"/",targ), pattern="*.tsv")
  #print(files)
  bases <- str_split(files, pattern="\\.tsv", simplify=TRUE)[,1] %>% unique() %>% sort()
  print(bases)
  bases %>% lapply(read_tsf,targ) %>% bind_rows %>% 
    mutate(target=targ)
}

read_dir_tsf("1kim")

full_data <- compl_targ %>% lapply(read_dir_tsf) %>% bind_rows()

full_data %>% filter(target=="1kim",tsf_name=="asite-actives-docked-3") %>% View()

full_data %>% distinct(tsf_name,target) %>% 
  filter(str_detect(tsf_name,"decoys_10")) 
head()

full_data_prep_noSel <- full_data %>%
  #slice_sample(n=1000) %>% 
  #filter(target%in%compl_targ) %>%
  mutate(Type=ifelse(str_detect(tsf_name,"actives"),1,0),
         Grid=ifelse(str_starts(tsf_name,"surface"),str_extract_all(tsf_name,"[0-9]+",simplify=T)[,1],0),
         Run=ifelse(str_starts(tsf_name,"surface"),str_extract_all(tsf_name,"[0-9]+",simplify=T)[,2],
                    ifelse(str_starts(tsf_name,"asite"),str_extract_all(tsf_name,"[0-9]+",simplify=T)[,1],"11"))) %>%
  convert(num(Grid,Run)) %>% 
  rename(dG=`dG, kcal/mol`,
         VS=`VS score`,
         RS=`Rank score`)

full_data_prep <- full_data_prep_noSel %>% 
  select(Ligand,target,Grid,Run,Type,dG,VS,RS)

full_data_prep %>% View()

full_data_prep_noSel %>% write_rds(paste0(PreparedData,".rds"))
full_data_prep %>% write_delim(PreparedData,";")

full_data_prep %>% colnames()

##Fast analysis----
##>Conventional docking----

full_data_prep <- read_delim(PreparedDataFile,";")

DATA_AStab <- full_data_prep %>% filter(Grid==0) %>% 
  select(-Grid) %>% 
  melt(id.vars=c("Ligand","target","Run","Type")) %>% 
  rename(SF=variable,ScoreAS=value) 
  
#separate(Ligand,sep="_",into=c("a","b","c")) %>% 
#  mutate(b=ifelse(Type==0,sprintf("%04d",as.numeric(b)+1),as.numeric(b)+1)) %>% 
#  mutate(Ligand=paste(a,b,c,sep = "_")) %>% 
#  select(-a,-b,-c)

DATA_OTtab <- full_data_prep %>% filter(Grid!=0) %>% 
  melt(id.vars=c("Ligand","target","Grid","Run","Type")) %>% 
  rename(SF=variable,ScoreOT=value)

AbsentLigands <- bind_rows(anti_join(DATA_AStab %>% distinct(Ligand), DATA_OTtab %>% distinct(Ligand)), anti_join(DATA_OTtab %>% distinct(Ligand),DATA_AStab %>% distinct(Ligand)))$Ligand

AbsentLigands

DATA_AS <- DATA_AStab %>% filter(!Ligand%in%c(AbsentLigands)) %>% 
  group_by(target,Ligand,Type,SF) %>% 
  summarise(RunAS=Run[which(ScoreAS==min(ScoreAS))],
            ScoreAS=min(ScoreAS))

DATA_AS_doROC <- DATA_AS %>% 
  group_by(target,SF) %>% 
  do(doROC=roc(.$Type,.$ScoreAS,direction=">"))

DATA_AS_AUC <- DATA_AS_doROC %>% 
  mutate(ROC=unlist(doROC,recursive=F)$auc) %>% 
  select(-doROC) %>% 
  ungroup()
DATA_AS_AUC

DATA_AS_AUC %>% 
  mutate(ROC=100*ROC) %>% 
  spread(SF,ROC) %>% 
  write_delim(paste0(CurrentDir,"/",Program,"_",TSet,"_ASD.csv"),";")

AS_ord <- DATA_AS_AUC %>% 
  filter(SF==mainSF) %>% 
  arrange(ROC) %>% 
  .$target

DATA_AS_ROC <- DATA_AS_doROC %>% 
  group_by(target,SF) %>%
  mutate(sens=list(unlist(doROC,recursive=F)$sensitivities),
         spec=list(unlist(doROC,recursive=F)$specificities)) %>% 
  select(-doROC) %>% 
  unnest(cols = c(sens, spec))

##>On-Top docking----

DATA_OTD <- DATA_OTtab %>% filter(!Ligand%in%c(AbsentLigands)) %>% 
  group_by(target,Ligand,Grid,SF) %>% 
  summarise(RunOT=Run[which(ScoreOT==max(ScoreOT,na.rm=T))],
            ScoreOT0=max(ScoreOT,na.rm=T),
            ScoreOT0=ifelse(ScoreOT0==-Inf,NA,ScoreOT0)) %>%
  filter(!is.na(ScoreOT0)) %>% 
  group_by(target,Ligand,SF) %>% 
  summarise(RunOT=RunOT[which.min(abs(ScoreOT0-median(ScoreOT0,na.rm=T)))],
            GridOT=Grid[which.min(abs(ScoreOT0-median(ScoreOT0,na.rm=T)))],
            ScoreOT=median(ScoreOT0,na.rm=T),
            ScoreOT=ifelse(ScoreOT==-Inf,NA,ScoreOT)) %>%
  left_join(DATA_AS) %>% 
  mutate(ScoreOTD=ScoreAS-ScoreOT)

DATA_OTD %>% 
  ungroup() %>% 
  filter(is.na(ScoreOT))

DATA_OTD %>% 
  ungroup() %>% 
  slice_min(ScoreAS,n=10)

DATA_OTD %>% 
  ungroup() %>% 
  slice_min(ScoreOTD,n=10)

DATA_OTD_doROC <- DATA_OTD %>% 
  group_by(target,SF) %>% 
  do(doROC=roc(.$Type,.$ScoreOTD,direction=">"))

DATA_OTD_AUC <- DATA_OTD_doROC %>% 
  mutate(ROC=unlist(doROC,recursive=F)$auc) %>% 
  select(-doROC) %>% 
  ungroup()
DATA_OTD_AUC

DATA_OTD_AUC %>% 
  mutate(ROC=100*ROC) %>% 
  spread(SF,ROC) %>% 
  write_delim(paste0(CurrentDir,"/",Program,"_",TSet,"_OTD.csv"),";")

DATA_OTD_ROC <- DATA_OTD_doROC %>% 
  group_by(target,SF) %>%
  mutate(sens=list(unlist(doROC,recursive=F)$sensitivities),
         spec=list(unlist(doROC,recursive=F)$specificities)) %>% 
  select(-doROC) %>% 
  unnest(cols = c(sens, spec))

left_join(DATA_AS_AUC %>% rename(AS_ROC=ROC),
          DATA_OTD_AUC %>% rename(OTD_ROC=ROC)) %>% 
  mutate(AS_ROC=round(AS_ROC,2),
         OTD_ROC=round(OTD_ROC,2))

left_join(DATA_AS_AUC %>% rename(AS_ROC=ROC),
          DATA_OTD_AUC %>% rename(OTD_ROC=ROC)) %>% 
  mutate(Gain=100*(OTD_ROC-AS_ROC)) %>% 
  select(-AS_ROC,-OTD_ROC) %>%
  #filter(abs(Gain)>2) %>% 
  spread(SF,Gain) %>% 
  bind_rows(.,left_join(DATA_AS_AUC %>% rename(AS_ROC=ROC),
                        DATA_OTD_AUC %>% rename(OTD_ROC=ROC)) %>% 
              mutate(Gain=100*(OTD_ROC-AS_ROC)) %>% 
              select(-AS_ROC,-OTD_ROC) %>% 
              spread(SF,Gain) %>% 
              summarise(plp=mean(plp)) %>% 
              mutate(target="Average")) %>% View()

left_join(DATA_AS_AUC %>% rename(AS_ROC=ROC),
          DATA_OTD_AUC %>% rename(OTD_ROC=ROC)) %>% 
  mutate(Gain=100*(OTD_ROC-AS_ROC)) %>% 
  select(-AS_ROC,-OTD_ROC) %>%
  #filter(abs(Gain)>2) %>% 
  spread(SF,Gain) %>% 
  bind_rows(.,left_join(DATA_AS_AUC %>% rename(AS_ROC=ROC),
                        DATA_OTD_AUC %>% rename(OTD_ROC=ROC)) %>% 
              mutate(Gain=100*(OTD_ROC-AS_ROC)) %>% 
              select(-AS_ROC,-OTD_ROC) %>% 
              spread(SF,Gain) %>% 
              summarise(plp=mean(plp)) %>% 
              mutate(target="Average")) %>%
  write_delim(paste0(CurrentDir,"/",Program,"_",TSet,"_Gains.csv"),";")

##>Analyse Scores disctributions----

ASRunsActNRJplot <- DATA_AStab %>% filter(Type==1) %>% 
  filter(SF==mainSF) %>%
  arrange(ScoreAS) %>% 
  group_by(Ligand,target) %>% 
  mutate(Run=row_number()) %>% 
  ggplot(aes(x=Run,y=ScoreAS, color=Ligand)) +
  geom_line(show.legend = F) +
  geom_smooth(color="black",size=1.5) +
  ylim(NA,20) +
  theme_bw() +
  facet_wrap(~target,ncol=5)
ggsave(paste0(CurrentDir,"/",Program,"_",TSet,"_ASRuns-ActNRJ.pdf"),ASRunsActNRJplot,units="cm",
       width=26,height = 40,dpi=2000)

OTRunsNRJplot <- DATA_OTtab %>% #filter(Type==1) %>%
  filter(SF==mainSF) %>%
  group_by(Ligand,target,Grid,Type) %>% 
  summarise(ScoreOT=max(ScoreOT,na.rm=T),
            ScoreOT=ifelse(ScoreOT==-Inf,NA,ScoreOT)) %>%
  filter(!is.na(ScoreOT)) %>% 
  arrange(ScoreOT) %>% 
  group_by(Ligand,target) %>% 
  mutate(Grid=percent_rank(row_number())) %>% 
  ungroup() %>% 
  arrange(Type) %>% 
  ggplot(aes(x=Grid,y=ScoreOT,group=Ligand, 
             color=factor(Type,levels=c(1,0)), 
             alpha=factor(Type,levels=c(1,0)))) +
  geom_line(show.legend = F) +
  scale_alpha_manual(values=c(0.5,0.15),name="Active") +
  scale_color_brewer(palette="Set1",name="Active") +
  #geom_smooth(alpha=1,size=1.5) +
  geom_line(size=0.5) +
  ylim(NA,20) +
  theme_bw() +
  facet_wrap(~target,ncol=5)
ggsave(paste0(CurrentDir,"/",Program,"_",TSet,"_OTRunsNRJ.pdf"),OTRunsNRJplot,units="cm",
       width=26,height = 40,dpi=2000)

ASandOTH_ActNRJplot <- DATA_OTtab %>% filter(Type==1) %>% 
  filter(SF==mainSF) %>%
  group_by(Ligand,target,Grid) %>% 
  summarise(Score=max(ScoreOT,na.rm=T),
            Score=ifelse(Score==-Inf,NA,Score)) %>% 
  arrange(Score) %>% 
  filter(!is.na(Score)) %>% 
  group_by(Ligand,target) %>% 
  mutate(Grid=percent_rank(row_number()),
         Site="OT") %>% 
  ungroup() %>% 
  bind_rows(DATA_AStab %>% filter(Type==1) %>%
              filter(SF==mainSF) %>%
              rename(Score=ScoreAS) %>% 
              arrange(Score) %>% 
              filter(!is.na(Score)) %>% 
              group_by(Ligand,target) %>% 
              mutate(Grid=percent_rank(row_number()),
                     Site="AS")) %>% 
  arrange(Site) %>% 
  ggplot(aes(x=Grid,y=Score,group=paste0(Ligand,Site),color=factor(Site,levels=c("AS","OT")), 
             alpha=factor(Site,levels=c("AS","OT")))) +
  geom_line() +
  scale_alpha_manual(values=c(0.15,0.15),name="Site") +
  scale_color_brewer(palette="Set1",name="Site") +
  geom_smooth(alpha=1,size=1.5,aes(group=Site),show.legend=F,se=F) +
  ylim(NA,20) +
  theme_bw() +
  facet_wrap(~target,ncol=5)
ggsave(paste0(CurrentDir,"/",Program,"_",TSet,"_AS-and-OTH_ActNRJ.pdf"),ASandOTH_ActNRJplot,units="cm",
       width=26,height = 40,dpi=2000)

ASandOTH_NRJdata <- DATA_OTtab %>% #filter(Type==1) %>% 
  filter(SF==mainSF) %>% 
  group_by(Ligand,Type,target,Grid) %>% 
  summarise(Score=max(ScoreOT,na.rm=T),
            Score=ifelse(ScoreOT==-Inf,NA,Score)) %>% 
  arrange(Score) %>% 
  filter(!is.na(Score)) %>% 
  group_by(Ligand,Type,target) %>% 
  mutate(Grid=percent_rank(row_number()),
         Site="OT") %>% 
  ungroup() %>% 
  bind_rows(DATA_AStab %>% #filter(Type==1) %>%
              filter(SF==mainSF) %>%
              rename(Score=ScoreAS) %>% 
              arrange(Score) %>% 
              filter(!is.na(Score)) %>% 
              group_by(Ligand,Type,target) %>% 
              mutate(Grid=percent_rank(row_number()),
                     Site="AS")) %>% 
  arrange(Site)

ASandOTH_NRJdataP <- ASandOTH_NRJdata %>% filter(target=="1qbo",Ligand=="1qbo_0443_dec") %>% 
  filter(Ligand%in%c(DATA_AStab %>% filter(target=="1qbo") %>% distinct(Ligand) %>% slice_sample(n=50) %>% .$Ligand))
  
ASandOTH_NRJplot <- ASandOTH_NRJdata %>% 
  filter(Score<100) %>% 
  ggplot(aes(x=Grid,y=Score,group=paste0(Ligand,Site),color=factor(Site,levels=c("AS","OT")), 
             alpha=factor(Site,levels=c("AS","OT")))) +
  #geom_line(size=0.05) +
  scale_alpha_manual(values=c(0.15,0.15),name="Site") +
  scale_color_brewer(palette="Set1",name="Site") +
  geom_smooth(alpha=1,size=1.5,aes(group=Site),se=F,span=1.25) +
  #ylim(NA,20) +
  coord_cartesian(ylim=c(-125,20)) +
  theme_bw() +
  facet_wrap(~target,ncol=5)
ggsave(paste0(CurrentDir,"/",Program,"_",TSet,"_AS-and-OTH_allNRJ.pdf"),ASandOTH_NRJplot,units="cm",
       width=26,height = 40,dpi=2000)

ASandOTH_NRJdata %>% filter(is.na(Type))

ASandOTH_NRJplot2 <- ASandOTH_NRJdata %>% 
  filter(Score<100) %>% 
  ggplot(aes(x=Grid,y=Score,group=paste0(Ligand,Site),color=factor(Site,levels=c("AS","OT")), 
             alpha=factor(Site,levels=c("AS","OT")))) +
  #geom_line(size=0.05) +
  scale_alpha_manual(values=c(0.15,0.15),name="Site") +
  scale_color_brewer(palette="Set1",name="Site") +
  geom_smooth(alpha=1,size=1.5,aes(group=paste0(Site,Type),linetype=factor(Type,levels=c(1,0))),se=F,span=1.25) +
  scale_linetype(name="Is active?",labels=c("Yes","No")) +
  #ylim(NA,20) +
  coord_cartesian(ylim=c(-125,20)) +
  theme_bw() +
  facet_wrap(~target,ncol=5)
ggsave(paste0(CurrentDir,"/",Program,"_",TSet,"_AS-and-OTH_Act-and-Dec-NRJ.pdf"),ASandOTH_NRJplot2,units="cm",
       width=26,height = 40,dpi=2000)

ASandOTH_NRJplotQuant <- ASandOTH_NRJdata %>% 
  filter(Score<100) %>% 
  ggplot(aes(x=Grid,y=Score,group=paste0(Ligand,Site),color=factor(Site,levels=c("AS","OT")), 
             alpha=factor(Site,levels=c("AS","OT")))) +
  #geom_line(size=0.05) +
  scale_alpha_manual(values=c(0.15,0.15),name="Site") +
  scale_color_brewer(palette="Set1",name="Site") +
  geom_quantile(quantiles=c(0.1),alpha=0.5,size=0.3,
                aes(group=paste0(Site,Type),
                    linetype=factor(Type,levels=c(1,0)))) +
  geom_quantile(quantiles=c(0.5),alpha=1,size=0.5,
                aes(group=paste0(Site,Type),
                    linetype=factor(Type,levels=c(1,0)))) +
  scale_linetype(name="Is active?",labels=c("Yes","No")) +
  #ylim(NA,20) +
  coord_cartesian(ylim=c(-125,20)) +
  theme_bw() +
  facet_wrap(~target,ncol=5)
ggsave(paste0(CurrentDir,"/",Program,"_",TSet,"_AS-and-OTH_Act-and-Dec-NRJ_Median.pdf"),
       ASandOTH_NRJplotQuant,units="cm",
       width=26,height = 40,dpi=2000)

##>Analyse ligand movements----

EnPlot <- DATA_OTD %>%
  ungroup() %>% 
  arrange(Type) %>% 
  #filter(!target%in%c("fpps")) %>% 
  ggplot(aes(x=ScoreAS,y=ScoreOTD,color=factor(Type,levels=c(1,0)),alpha=factor(Type,levels=c(1,0)))) +
  geom_point(size=1) +
  geom_smooth(method="lm",aes(group=paste0(target,SF)),color="black",size=0.2,show.legend=F) +
  scale_alpha_manual(values=c(0.9,0.15),name="Active") +
  #facet_grid(target~SF) +
  facet_wrap(~target,ncol=5) +
  scale_color_brewer(palette="Set1",name="Active") +
  theme_bw()
ggsave(paste0(CurrentDir,"/",Program,"_",TSet,"_EnergyPlot.pdf"),EnPlot,units="cm",
       width=26,height = 40,dpi=2000)

EnDensPlot <- DATA_OTD %>%
  ungroup() %>% 
  arrange(Type) %>% 
  #filter(!target%in%c("fpps")) %>% 
  ggplot(aes(x=ScoreOTD,color=factor(Type,levels=c(1,0)))) +
  geom_density(alpha=0.7) +
  xlim(NA,20) +
  #facet_grid(target~SF) +
  facet_wrap(~target,ncol=5) +
  scale_color_brewer(palette="Set1",name="Active") +
  theme_bw() +
  theme(legend.position="bottom")
ggsave(paste0(CurrentDir,"/",Program,"_",TSet,"_EnDensPlot.pdf"),EnDensPlot,units="cm",
       width=26,height = 40,dpi=2000)

DATA_Ranks <- DATA_OTD %>%
  group_by(target,SF) %>% 
  mutate(RankAS=rank(ScoreAS),
         RankOTD=rank(ScoreOTD),
         RankDiff=RankOTD-RankAS) 

DATA_Ranks %>% 
  group_by(target,SF) %>% 
  filter(Type==1) %>% 
  slice_max(RankDiff,n=10) %>% 
  arrange(target,SF,desc(RankDiff)) %>% 
  write_delim(paste0(CurrentDir,"/",Program,"_",TSet,"_DegradingActives.csv"),";")

DATA_Ranks %>% 
  group_by(target,SF) %>% 
  filter(Type==0) %>% 
  slice_min(RankDiff,n=10) %>% 
  arrange(target,SF,RankDiff) %>% 
  write_delim(paste0(CurrentDir,"/",Program,"_",TSet,"_FlorishingDecoys.csv"),";")

RankPlot <- DATA_Ranks %>% 
  ungroup() %>% 
  filter(SF==mainSF) %>%
  arrange(Type) %>% 
  #filter(Type==1) %>% 
  mutate(PDB=factor(target,levels=rev(AS_ord))) %>%
  ggplot(aes(x=RankAS,y=RankOTD)) +
  geom_abline(slope=1,intercept=0,color="grey30") +
  geom_point(aes(color=factor(Type,levels=c(1,0)),alpha=factor(Type,levels=c(1,0)))) +
  scale_alpha_manual(values=c(0.9,0.15),name="Active") +
  scale_color_brewer(palette="Set1",name="Active") +
  facet_wrap(~PDB,ncol=5) +
  theme_bw() +
  scale_x_continuous(breaks=seq(0,10000,250),name="Ranking according to the conventional docking", limits=c(0,1100)) +
  scale_y_continuous(breaks=seq(0,10000,250),name="Ranking according to the On-Top docking", limits=c(0,1100)) +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        axis.text.x = element_text(angle = -45, vjust=1.0, hjust=0.0))
RankPlot
ggsave(filename=paste0(CurrentDir,"/",Program,"_",TSet,"_RankingPlot.pdf"), RankPlot, width = 16.5, height = 20, units = "cm", dpi=2000)


##>ROC curves comparison----

ROC_curves_Plot <- bind_rows(DATA_AS_ROC %>% mutate(Dock="ASD"),
                             DATA_OTD_ROC %>% mutate(Dock="OTD")) %>% 
  ungroup() %>% 
  filter(SF==mainSF) %>% 
  mutate(PDB=factor(target,levels=rev(AS_ord))) %>%  
  ggplot(aes(y=sens,x=1-spec,color=Dock)) +
  scale_color_brewer(palette="Set1",
                     breaks=c("OTD","ASD"),
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
ggsave(filename=paste0(CurrentDir,"/",Program,"_",TSet,"_ROCcurves.pdf"),ROC_curves_Plot, width = 16.5, height = 20, units = "cm", dpi=2000)


##Repeat full analysis----

cl <- new_cluster(10) # Creates a 10-threads cluster. Reduce if you do not have 10 threads.
cluster_library(cl, c("tidyverse","pROC"))

##>FullAn: Active site----

DataAS <- DATA_AStab %>% 
  filter(SF==mainSF,
         !is.na(ScoreAS)) %>% 
  select(-SF) %>% 
  rename(ligand=Ligand,
         type=Type) %>% 
  group_by(target,ligand) %>% 
  crossing(.,tibble(QAS=c(0,0.5,1.0))) %>%
  group_by(target,ligand,QAS,type) %>% 
  partition(cluster=cl) %>% 
  summarize(ASVal=quantile(ScoreAS,probs=unique(QAS)),
            count=n()) %>% 
  collect()

colnames(DataAS)
DataAS %>% group_by(count) %>% 
  summarise(n=n())

DataAS %>% 
  write.csv(paste0(CurrentDir,"/",Program,"_",TSet,"_DataAS_Prod.csv"), row.names = F, quote = F)

#Continue here if you have loaded data from DataAS_Prod.csv

DataAS %>% 
  #filter(count==10) %>%
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

##>FullAn: On-Top docking----

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
    crossing(tibble(QRunSurf=c(0,0.25,0.5,0.75,1.0))) %>%
    group_by(QRunSurf,ligand,grid) %>% 
    partition(cluster=cl) %>% 
    summarize(SurfRunVal=GetSurfRun(QRunSurf,ScoreOT)) %>% 
    collect()
})

cluster_copy(cl,"GetSurfRun")
DataSurf_sepRuns_0 <- DATA_OTtab %>% 
  filter(SF==mainSF,
         !is.na(ScoreOT)) %>% 
  select(-SF) %>% 
  rename(ligand=Ligand,
         type=Type,
         grid=Grid) %>% 
  group_by(target) %>% 
  do(SurfRunVal0=GetSurfRunVal(.)) %>% 
  ungroup() %>% 
  unnest(cols = c(SurfRunVal0))

GetSurfVal <- compiler::cmpfun(function(df){
  df %>% 
    crossing(tibble(QSurf=seq(0,1,0.1))) %>%
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
  write.csv(paste0(CurrentDir,"/",Program,"_",TSet,"_DataSurf_sepRuns_Prod.csv"), row.names = F, quote = F)

#Continue here if you have loaded data from DataSurf_sepRuns_Prod.csv and DataAS_Prod.csv

left_join(DataAS,DataSurf_sepRuns) %>% 
  group_by(target,QAS,QSurf,QRunSurf) %>% 
  mutate(count=length(unique(type))) %>% View()

DataOTD_Self_sepRuns_doROC <- 
  left_join(DataAS,DataSurf_sepRuns) %>% 
  filter(!ligand%in%c(AbsentLigands)) %>% 
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
  geom_hline(data=DataAS_AUC %>% mutate(SF=mainSF),
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
ggsave(filename=paste0(CurrentDir,"/",Program,"_",TSet,"_Self_Quantiles_SI_Prod.pdf"),Plot_Self_Means_sepRuns_Full, 
       width = 16, height = 130, units = "cm", dpi=1000, limitsize = FALSE)
ggsave(filename=paste0(CurrentDir,"/",Program,"_",TSet,"_Self_Quantiles_SI_Prod.png"),Plot_Self_Means_sepRuns_Full, 
       width = 16, height = 18, units = "cm", dpi=1000)

Plot_Self_Means_sepRuns <- DataOTD_Self_sepRuns %>% 
  filter(QRunSurf%in%c(0,0.5,1,"Arb")) %>% 
  mutate(QRunSurf=factor(QRunSurf,levels=c("Arb",0,0.5,1))) %>% 
  ggplot(aes(x=QSurf,y=ROC,color=factor(QRunSurf))) +
  geom_hline(data=DataAS_AUC_mean %>% mutate(SF=mainSF),
             aes(yintercept=ROC),color="grey10",size=0.5) +
  geom_smooth(aes(group=paste0(QAS,QRunSurf)),size=0.5,level=0.75,se=F) +
  scale_color_manual(values=c("#984ea3","#4daf4a","#377eb8","#e41a1c","#ff7f00"),
                     guide = guide_legend(reverse = TRUE), 
                     name="Quantile of the\non-top runs used") +
  theme_bw() +
  #guides(color=F) +
  facet_grid(~QAS) +
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
ggsave(filename=paste0(CurrentDir,"/",Program,"_",TSet,"_Self_Quantiles_Prod.png"),
       Plot_Self_Means_sepRuns, width = 12, height = 6.0, units = "cm", dpi=2000)
