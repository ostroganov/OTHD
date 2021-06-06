
Data10sampOTS %>% head()

DataSurf_sepRuns %>% head()

AllMaxOTS <- Data10sampOTS %>% 
  group_by(target,ligand,grid) %>% 
  summarise(max.dG=max(top.dG))

AllMaxOTS2 <- AllMaxOTS %>% 
  left_join(DataSurf_sepRuns %>% 
              filter(QSurf==0.5,QRunSurf==1) %>% 
              select(target,ligand,SurfVal)) %>% 
  mutate(diffSV=max.dG-SurfVal)

BestSurfSites <- AllMaxOTS2 %>% 
  #filter(abs(diffSV)<0.5) %>% 
  group_by(target,grid) %>% 
  summarise(ratio=sum(abs(diffSV)<0.5)/n()) %>%
  group_by(target) %>% 
  filter(ratio==max(ratio))
BestSurfSites %>% write_delim("BestSites.csv",";")

#Calc new OTD----

DataSurf_sepRuns %>% head()

DataOTD_BSS_sepRuns_doROC <- 
  left_join(DataAS %>% filter(QAS==0) %>% select(-QAS),
            AllMaxOTS %>% semi_join(BestSurfSites) %>% select(-grid)) %>% 
  group_by(target) %>% 
  do(doROC=roc(.$type,.$ASVal-.$max.dG,direction=">")) 

DataOTD_BSS_sepRuns <- DataOTD_BSS_sepRuns_doROC %>% 
  mutate(ROC=unlist(doROC,recursive=F)$auc) %>% 
  select(-doROC) %>% 
  ungroup() %>% 
  mutate(Baseline="BestSurfSite") 

DataOTD_BSS_sepRuns

DataOTD_Self_sepRuns %>% 
  filter(QAS==0,QSurf==0.5,QRunSurf==1) %>% 
  rename(ROC_OTD=ROC) %>% 
  select(target,ROC_OTD) %>% 
  left_join(DataOTD_BSS_sepRuns %>% 
              rename(ROC_BSS=ROC) %>% 
              select(-Baseline)) %>% 
  mutate(diff=100*(ROC_BSS-ROC_OTD)) %>% write_delim("BestSites_performance.csv",";")

