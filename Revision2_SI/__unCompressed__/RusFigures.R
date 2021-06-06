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
  scale_y_continuous(breaks=seq(0,1,0.05), name="Площадь под кривой", labels=scaleFUN) + 
  scale_x_discrete(limits=c("None","Self"),
                   labels=c("Докинг","Докинг по\nповерхности"),
                   expand = c(0.1,0.1),
                   name="Docking") +
  theme_bw() +
  theme(legend.position="bottom",
        text=element_text(size=7),
        axis.title.y = element_text(vjust = 0.25),
        legend.title = element_blank(),
        plot.margin = margin(2,2,0,1,unit="mm"),
        axis.title.x = element_blank(),
        panel.spacing = unit(0,"mm"),
        legend.text = element_text(size=5),
        legend.margin=margin(0,0,0,0),
        legend.box.spacing = unit(0.0,"mm"),
        legend.box.margin =  margin(0.0,5,0.4,0,unit="mm"),
        legend.key.size = unit(1,"mm"))
OTD_GainPlot_Prod
ggsave(filename="OTD_GainPlot_Prod_Rus.png",OTD_GainPlot_Prod, width = 4.2, height = 6, units = "cm", dpi=2000)

ImprovementPlot2 <- GainTab %>% 
  rename(ID=target) %>% 
  left_join(TestSet2 %>% 
              left_join(TestSet %>% 
                          select(ID,OpenAS,NumStAl))) %>%
  mutate(OPEN=ifelse(OpenAS==1,"Открытый АЦ","Закрытый АЦ"),
         TotCh=TotCh/as.numeric(NumStAl)) %>% 
  filter(pH==7.0) %>% 
  ggplot(aes(x=TotCh,y=100*Gain,shape=OPEN,color=None)) +
  scale_color_gradientn(colors=c("dark red","red","orange","yellow","dark green"), 
                        name="Точность традиционного докинга",
                        labels = scales::percent_format(accuracy = 1)) +
  scale_shape_manual(values=c(16,1),name="") +
  scale_x_continuous(name="Заряд белка при pH 7") +
  scale_y_continuous(name="Улучшение (%)",
                     #breaks=c(0,10,20,30,40,50),
                     #labels=c(0,10,"",30,"",50),
                     position = "left") +
  geom_smooth(span=3,level=0.75,aes(group=NA),color="grey25") +
  geom_point(size=3,aes(group=ID),stroke=1.2,alpha=0.8) +
  theme_bw() +
  geom_text(label="Рибонуклеаза A",x=-0.5,y=46,size=3,hjust=1,color="#a23333") +
  geom_text(label="Тимидилатсинтаза",x=-34,y=38,size=3,hjust=0,color="#fc3333") +
  geom_text(label="Рибонуклеаза T1",x=-35.5,y=32,size=3,hjust=0,color="#ff7733") +
  geom_text(label="Поли(АДФ-рибоза)-полимераза",x=-0.2,y=24.58,size=3,hjust=1,color="#b03335") +
  theme(legend.position = "bottom",
        #legend.margin = margin(0,0,-5,0,"mm"),
        #axis.text.y = element_text(color=c("grey30","grey50","grey30","grey30","grey30","grey30")),
        legend.box.margin = margin(-3.5,0,-2,-10,"mm"),
        legend.direction = "horizontal",
        legend.box = "horizontal") +
  guides(color=guide_colorbar(title.position = "top",
                              title.theme = element_text(size = 10,color="grey20"),
                              barwidth = unit(55,"mm"),order=2,
                              barheight = unit(2,"mm"),nbin=100,ticks.linewidth=1.2,
                              label.theme = element_text(size=8,color="grey30")),
         shape=guide_legend(nrow=2,order=1))
ImprovementPlot2
ggsave("ImprovementPlot_TotCharge_v2_Rus.png",ImprovementPlot2,height=8.75,width=9.5,dpi=1000,units="cm")
