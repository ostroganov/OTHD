Data10sampAS <- read_delim(file="./old/DataAS_10runs.csv",";")
Data10sampOTS <- read_delim(file="./old/DataOTS_10runs.csv",";")

Data10sampAS %>% head()

Data10sampOTS %>% head()

Data10sampOTS %>% distinct(grid)

Data10sampOTS %>% 
  filter(target=="1efy",grid==9) %>% 
  mutate(type=ifelse(str_detect(ligand,"active"),1,0),
         grid=0) 

Data10sampAS_new <- Data10sampAS %>% 
  filter(target!="1efy") %>% 
  bind_rows(Data10sampOTS %>% 
              filter(target=="1efy",grid==9) %>% 
              mutate(type=ifelse(str_detect(ligand,"active"),1,0),
                     grid=0)) %>% 
  arrange(target,ligand,grid,run)

Data10sampAS %>% 
  filter(target=="1efy") %>% 
  select(-type) %>% 
  mutate(grid=9) 

Data10sampOTS_new <- Data10sampOTS %>% 
  filter(!(target=="1efy"&grid==9)) %>% 
  bind_rows(Data10sampAS %>% 
              filter(target=="1efy") %>% 
              select(-type) %>% 
              mutate(grid=9)) %>% 
    arrange(target,ligand,grid,run)

Data10sampAS_new %>% 
  write.table(file="DataAS_10runs.csv",quote=F,sep=";",dec=".",row.names=F)

Data10sampOTS_new %>% 
  write.table(file="DataOTS_10runs.csv",quote=F,sep=";",dec=".",row.names=F)
