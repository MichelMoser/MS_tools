###############################################################################################
#
#   timsR    package to inspect raw timsTOF data 
#
###############################################################################################
#library(devtools)
# install bug fixed for R4.3.3
#install_github("sneumann/opentims", subdir="opentimsr", ref = "fix/issue25")
#install_github("MatteoLacki/timsr")
#devtools::install_github("daattali/ggExtra")

# setup 

path_to_bruker_dll <- "/home/momi/tools/timsdata/linux64/libtimsdata.so"
setup_bruker_so(path_to_bruker_dll) 
all_columns = c('frame','scan','tof','intensity','mz','inv_ion_mobility','retention_time')


library(timsr)
library(tidyverse)
library(gganimate)
library(viridis)
library(ggExtra)



#########################
#  choose run as input 
########################

#DIA medium 
path = '/home/momi/DATA/DIA_NN/3325_Abrin_pure_DIA_Slot2-10_1_10-3-2023_14-10-44.d'
#DIA short 
path = '/home/momi/DATA/DIA_NN/4464_0036-24-ST-0031_Slot1-22_1_2-13-2024_22-39-30.d/'
#DDA
path = '/home/momi/DATA/DIA_NN/4736_GT240313_c-7_Slot2-10_25min_AUR15_300nl_50deg_3-14-2024_17-46-54.d/'
path


D = TimsR(path) # get data handle

print(length(D)) # The number of peaks in the run
summary(D)

#########################################################################
#
# inspect tables of timsR object
#  
# query tdf sqlite database 
#########################################################################
 
# name of tables in the tdf database
tables_names(D)


# global data
table2dt(D, "GlobalMetadata")
# calibration used for the run
table2dt(D, "TimsCalibration")
table2dt(D, "MzCalibration")
table2dt(D, "CollisionEnergySweepingInfo")
table2dt(D, "Frames")                 # list of all the frames with MSMSType and overall intensity
table2dt(D, "FrameMsMsInfo")
table2dt(D, "DiaFrameMsMsInfo")
table2dt(D, "DiaFrameMsMsWindows")     # DIA window settings, start and end in mz and 
table2dt(D, "DiaFrameMsMsWindowGroups")
table2dt(D, "PropertyDefinitions") %>% names()
table2dt(D, "PropertyDefinitions") 


#check out Funnel 1 radio frequency
table2dt(D, "Properties") %>% filter(Properties.Property == 386) 

table2dt(D, "Properties") %>% names()
table2dt(D, "Properties") %>% filter(Properties.Property %in% c(3, 10, 12, 53, 54, 341,343, 362, 363, 441)) %>% view



# inspect all tables  
for (i in tables_names(D)){
  
  print(i)
  print(table2dt(D, i))
}





# MS1 frames versus MS2 frames

seq(1,17473,9)
table2dt(D, "Frames")$Frames.MsMsType[seq(1,17473,9)]  #pick only MS1

dim(table2dt(D, "Frames"))

D[1, all_columns] %>% view
table2dt(D, 1)

D

?plot_TIC
plot_TIC(D)
?rt_query

rt_query(D, 1001, 1011)

RTS = opentimsr::retention_times(D)
length(RTS)

# frame : a collection of scans taken over the DIA windows with identical retention time 

#####################################
# function to plot single frame 
#####################################

plot_frame <- function(N){
  frameN <- D[N, all_columns]
  
  frameN %>%   
    rename("mz_start" = "mz", "IM_start" = "inv_ion_mobility") %>% 
    
    ggplot(aes(mz_start, IM_start, colour = intensity)) + 
    ggtitle(paste("Frame", N)) +
    geom_point(size = 0.1) +
    coord_fixed(ratio = 500) +
    #geom_rect(data = dia_windows, mapping=aes(xmin=mz_start, xmax=mz_end, ymin=IM_start, ymax=IM_end, fill=factor(DiaFrameMsMsWindows.WindowGroup)), color="white", alpha=0.5) +
    #geom_rect(data = F1_slice1, mapping=aes(xmin=mz_start, xmax=mz_end, ymin=IM_start, ymax=IM_end, fill = factor(cycle_id)), color="white", alpha=0.3) +
    #geom_rect(data = F1_steep, mapping=aes(xmin=mz_start, xmax=mz_end, ymin=IM_start, ymax=IM_end, fill = factor(cycle_id)), color="white", alpha=0.3) +
    #geom_rect(data = F1_flat, mapping=aes(xmin=mz_start, xmax=mz_end, ymin=IM_start, ymax=IM_end, fill = factor(cycle_id)), color="white", alpha=0.3) +
    scale_fill_viridis(discrete = TRUE) +
    guides(fill=guide_legend(title="DIA windows"))
}



plot_frame(4645)
plot_frame(46)

# scan : single ion-mobility snapshot and all m/z thereof with identical ion-mobility and retention time 

#################################
# function to plot single scan 
################################

D[3, all_columns] %>% filter(scan == 925)

plot_scan <- function(F, S){
  frameS <- D[F, all_columns] %>% filter(scan == S)
  frameS
  dim(frameS)
  frameS %>%   
    rename("mz_start" = "mz", "IM_start" = "inv_ion_mobility") %>% 
    
    ggplot(aes(mz_start, IM_start, colour = intensity)) + 
    ggtitle(paste("Frame", F, " scan", S)) +
    geom_point(size = 0.1) +
    coord_fixed(ratio = 500) +
    #geom_rect(data = dia_windows, mapping=aes(xmin=mz_start, xmax=mz_end, ymin=IM_start, ymax=IM_end, fill=factor(DiaFrameMsMsWindows.WindowGroup)), color="white", alpha=0.5) +
    #geom_rect(data = F1_slice1, mapping=aes(xmin=mz_start, xmax=mz_end, ymin=IM_start, ymax=IM_end, fill = factor(cycle_id)), color="white", alpha=0.3) +
    #geom_rect(data = F1_steep, mapping=aes(xmin=mz_start, xmax=mz_end, ymin=IM_start, ymax=IM_end, fill = factor(cycle_id)), color="white", alpha=0.3) +
    #geom_rect(data = F1_flat, mapping=aes(xmin=mz_start, xmax=mz_end, ymin=IM_start, ymax=IM_end, fill = factor(cycle_id)), color="white", alpha=0.3) +
    scale_fill_viridis(discrete = TRUE) +
    guides(fill=guide_legend(title="DIA windows"))
}

plot_scan(5, 36)


# get raw data 
frameN <- D[4640, all_columns]

dim(frameN)

max(frameN$mz)
max(frameN$inv_ion_mobility)
min(frameN$inv_ion_mobility)

length(seq(0.6, 1.5, 0.005))
seq(0.6, 1.5, 0.02)

tm <-  frameN %>% mutate(new_mz = cut(mz, breaks=seq(100, 1500, 1)), 
                  new_mob = cut(inv_ion_mobility, breaks=seq(0.6, 1.5, 0.01))) %>% 
  group_by(new_mz, new_mob) %>% 
  summarise(int_bin = mean(intensity)) %>% 
pivot_wider(names_from = new_mob, values_from = int_bin) 

tm %>% mutate(new_mz)
tmdf <- as.data.frame(tm)
tmdf$new_mz <- NULL


as.matrix(tmdf)
tmdf[is.na(tmdf)] <- 0


gplots::heatmap.2((as.matrix(tmdf)),dendrogram='none', Rowv=FALSE, Colv=FALSE,trace='none', col = redgreen(250))

(as.matrix(tmdf))
?heatmap
row.names(tmdf)

####  3d plot
#plot_ly(z =~as.matrix(tmdf)) %>% add_surface() %>% 
plot_ly() %>%   layout(title = 'frame 4640, tims-TOF data 3D view', 
         scene = list(
           yaxis = list(title = 'm/z'),
           xaxis = list(title = 'inverse ion mobility'),
           zaxis = list(title = 'intensity'))) %>%  add_trace(z = ~as.matrix(tmdf),
            name = 'ZED',
            type = 'surface',
            colorbar=list(title='Intensity'))






plot_frame(4656)


frameN %>%   
  rename("mz_start" = "mz", "IM_start" = "inv_ion_mobility") %>% 
  
  ggplot(aes(mz_start, intensity, colour = intensity)) + 
  geom_point(size = 0.1) 
  coord_fixed(ratio = 500) 

volcano
library(plotly)



plot_ly(data = frameN, x=~mz, y = ~inv_ion_mobility, z =~intensity, type = "scatter3d", mode = "markers", color =~intensity) %>% add_surface()
plot_ly(z = ~volcano) %>% add_surface()



pf <- frameN %>%   
  rename("mz_start" = "mz", "IM_start" = "inv_ion_mobility") %>% 
  
  ggplot(aes(mz_start, IM_start, colour = intensity)) + 
  geom_point(size = 0.1) +
  coord_fixed(ratio = 500) 

ggMarginal(pf)  

length(intensity_per_frame(D))

tibble("int" = intensity_per_frame(D)[seq(1,length(intensity_per_frame(D)),9)], 
            "time" = retention_times(D)[seq(1,length(intensity_per_frame(D)),9)]) %>% 
  ggplot(aes(time, int)) + geom_line()

RTS


########################################
#   interactive plot of TICs along an MS run
#########################################
library(wesanderson)

# get data to plot TIC for complete MS run
TICplot <- tibble(TIC_intensity = intensity_per_frame(D),
frame = seq(1, length(RTS)), 
time = retention_times(D),
type = table2dt(D, "Frames")$Frames.MsMsType) %>% 
  mutate(type = ifelse(type == 0, "MS1", "MS2")) %>% 
  ggplot(aes(time, TIC_intensity, group = type, color = type, info = frame)) + geom_line() + 
    theme_bw(base_size = 18) + 
  ggtitle("TIC plot") + 
  ylab("TIC intensity") +
  xlab("RT (seconds)") + 
  
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_color_brewer(palette = "Set1")


plotly::ggplotly(TICplot, tooltip = c("info", "time", "TIC_intensity") )




plot_frame(4656)

%>% filter(DiaFrameMsMsWindows.WindowGroup == 1)
plot_frame(4640)

dia_windows %>% filter(DiaFrameMsMsWindows.WindowGroup == 1)

t <- ggMarginal(plot_frame(4645)+ ggtitle("MS1")+ylab("inverse ion mobility")+xlab("m/z")+
                  geom_rect(data = dia_windows , mapping=aes(xmin=mz_start, xmax=mz_end, ymin=IM_start, ymax=IM_end, fill=factor(DiaFrameMsMsWindows.WindowGroup)), color="white", alpha=0.5) 
                  , type = "histogram", xparams = list(binwidth = 1, fill = "orange"), yparams = list(binwidth = 0.01))
t1 <- ggMarginal(plot_frame(4646)+ 
                   ggtitle("MS2 DIA window 1")+
                   geom_rect(data = dia_windows  %>% filter(DiaFrameMsMsWindows.WindowGroup == 1), mapping=aes(xmin=mz_start, xmax=mz_end, ymin=IM_start, ymax=IM_end, fill=factor(DiaFrameMsMsWindows.WindowGroup)), color="white", alpha=0.5)+ 
                   ylab("inverse ion mobility")+xlab("m/z"), type = "histogram", xparams = list(binwidth = 1, fill = "orange"), yparams = list(binwidth = 0.01))
t2 <- ggMarginal(plot_frame(4647)+ 
                   ggtitle("MS2 DIAwindow 2")+
                   geom_rect(data = dia_windows  %>% filter(DiaFrameMsMsWindows.WindowGroup == 2), mapping=aes(xmin=mz_start, xmax=mz_end, ymin=IM_start, ymax=IM_end, fill=factor(DiaFrameMsMsWindows.WindowGroup)), color="white", alpha=0.5)+ 
                   ylab("inverse ion mobility")+
                   xlab("m/z"),
                 type = "histogram", xparams = list(binwidth = 1, fill = "orange"), yparams = list(binwidth = 0.01))
t3 <- ggMarginal(plot_frame(4648)+ ggtitle("MS2 DIAwindow 3")+
                   geom_rect(data = dia_windows  %>% filter(DiaFrameMsMsWindows.WindowGroup == 3), mapping=aes(xmin=mz_start, xmax=mz_end, ymin=IM_start, ymax=IM_end, fill=factor(DiaFrameMsMsWindows.WindowGroup)), color="white", alpha=0.5)+ 
                   ylab("inverse ion mobility")+
                   xlab("m/z"),
                 type = "histogram", xparams = list(binwidth = 1, fill = "orange"), yparams = list(binwidth = 0.01))
t4 <- ggMarginal(plot_frame(4649)+ ggtitle("MS2 DIAwindow 4")+
                   geom_rect(data = dia_windows  %>% filter(DiaFrameMsMsWindows.WindowGroup == 4), mapping=aes(xmin=mz_start, xmax=mz_end, ymin=IM_start, ymax=IM_end, fill=factor(DiaFrameMsMsWindows.WindowGroup)), color="white", alpha=0.5)+ 
                   ylab("inverse ion mobility")+
                   xlab("m/z"),
                 type = "histogram", xparams = list(binwidth = 1, fill = "orange"), yparams = list(binwidth = 0.01))
t5 <- ggMarginal(plot_frame(4650)+ ggtitle("MS2 DIAwindow 5")+
                   geom_rect(data = dia_windows  %>% filter(DiaFrameMsMsWindows.WindowGroup == 5), mapping=aes(xmin=mz_start, xmax=mz_end, ymin=IM_start, ymax=IM_end, fill=factor(DiaFrameMsMsWindows.WindowGroup)), color="white", alpha=0.5)+ 
                   ylab("inverse ion mobility")+
                   xlab("m/z"),
                 type = "histogram", xparams = list(binwidth = 1, fill = "orange"), yparams = list(binwidth = 0.01))
t6 <- ggMarginal(plot_frame(4651)+ ggtitle("MS2 DIAwindow 6")+
                   geom_rect(data = dia_windows  %>% filter(DiaFrameMsMsWindows.WindowGroup == 6), mapping=aes(xmin=mz_start, xmax=mz_end, ymin=IM_start, ymax=IM_end, fill=factor(DiaFrameMsMsWindows.WindowGroup)), color="white", alpha=0.5)+ 
                   ylab("inverse ion mobility")+
                   xlab("m/z"),
                 type = "histogram", xparams = list(binwidth = 1, fill = "orange"), yparams = list(binwidth = 0.01))
t7 <- ggMarginal(plot_frame(4652)+ ggtitle("MS2 DIAwindow 7")+
                   geom_rect(data = dia_windows  %>% filter(DiaFrameMsMsWindows.WindowGroup == 7), mapping=aes(xmin=mz_start, xmax=mz_end, ymin=IM_start, ymax=IM_end, fill=factor(DiaFrameMsMsWindows.WindowGroup)), color="white", alpha=0.5)+ 
                   ylab("inverse ion mobility")+
                   xlab("m/z"),
                 type = "histogram", xparams = list(binwidth = 1, fill = "orange"), yparams = list(binwidth = 0.01))
t8 <- ggMarginal(plot_frame(4653)+ ggtitle("MS2 DIAwindow 8")+
                   geom_rect(data = dia_windows  %>% filter(DiaFrameMsMsWindows.WindowGroup == 8), mapping=aes(xmin=mz_start, xmax=mz_end, ymin=IM_start, ymax=IM_end, fill=factor(DiaFrameMsMsWindows.WindowGroup)), color="white", alpha=0.5)+ 
                   ylab("inverse ion mobility")+
                   xlab("m/z"),
                 type = "histogram", xparams = list(binwidth = 1, fill = "orange"), yparams = list(binwidth = 0.01))

t
t1
t2
t3
t4
t5
t6
t7
t8

















table2dt(D, 'GlobalMetadata')
#table2dt(D, 'DiaFrameMsMsInfo') %>% view

table2dt(D, 'Properties')


# get ion mobility to  DIA windows information

#table2dt(D, 'DiaFrameMsMsWindows') %>% as_tibble() %>% view


#single frame 

frame90000 <- D[9*1000, all_columns]

frame1 <- D[1, all_columns]
frame1 <- as_tibble(frame1)

scan_to_im <- frame1 %>% select(scan, inv_ion_mobility) %>% unique()
scan_to_im


# MS 2 frame
frame2 <- D[2, all_columns]
frame2 <- as_tibble(frame2)
frame2

frame2 %>% select(scan, inv_ion_mobility) %>% unique()


frame90000





frameN <- D[9001, all_columns]

frameN %>% ggplot(aes(mz)) + geom_density()


plot_frame(9001)



60 * 18


### read timstof windows from timsControl file: 

timsControlWindows <- read_csv("~/Downloads/diaParameters_dia-PASEF-shortgradient_piAid3_BM.txt") %>%   filter(!row_number() %in% c(1)) %>% 
  mutate(across(c(3:6), as.numeric), 
         source = "timsControl") %>% 
  rename("IM_start" = "Start IM [1/K0]",  "IM_end" = "End IM [1/K0]", "mz_start" = "Start Mass [m/z]", "mz_end" = "End Mass [m/z]", "cycle_id" = "Cycle Id")  


timsControlWindows

timsControlWindows %>%  ggplot() +
  geom_rect(mapping=aes(xmin=mz_saatart, xmax=mz_end, ymin=IM_start, ymax=IM_end, fill=factor(cycle_id)), color="black", alpha=0.5) +
  xlab("m/z")+
  ylab("1/IM")


timsControlWindows <- timsControlWindows %>% select(-c(1,7))
timsControlWindows


dia_windows




# add ion mobility to  DIA windows information
dia_windows <- table2dt(D, 'DiaFrameMsMsWindows') %>% as_tibble() %>% 
  left_join(scan_to_im, by = c("DiaFrameMsMsWindows.ScanNumBegin" = "scan")) %>%  
  left_join(scan_to_im, by = c("DiaFrameMsMsWindows.ScanNumEnd" = "scan")) %>% 
  rename("IM_start" = "inv_ion_mobility.x" ,  "IM_end" = "inv_ion_mobility.y" ) %>% 
  mutate("mz_start" = DiaFrameMsMsWindows.IsolationMz - DiaFrameMsMsWindows.IsolationWidth/2, 
         "mz_end" = DiaFrameMsMsWindows.IsolationMz + DiaFrameMsMsWindows.IsolationWidth/2, 
         source = "tdf")  


dia_windows


dia_windows %>% ggplot() + 
  geom_rect(data = dia_windows, mapping=aes(xmin=mz_start, xmax=mz_end, ymin=IM_start, ymax=IM_end, fill=factor(DiaFrameMsMsWindows.WindowGroup)), color="black", alpha=0.5) 





tdfWindows<- dia_windows %>% select(DiaFrameMsMsWindows.WindowGroup, IM_start,  IM_end,  mz_start,  mz_end, source) %>% rename("cycle_id" = "DiaFrameMsMsWindows.WindowGroup")
tdfWindows


#rbind(timsControlWindows, tdfWindows) %>% 
tdfWindows %>% 
    #mutate(mz_end = ifelse(source == "tdf", mz_start, mz_end)) %>% 
  ggplot() +
  geom_rect( mapping=aes(xmin=mz_start, xmax=mz_end, ymin=IM_start, ymax=IM_end, fill=cycle_id), color="black", alpha=0.5) 




frame1 %>% select(scan, inv_ion_mobility) %>% unique() %>% mutate(Diff = lag(inv_ion_mobility) - inv_ion_mobility)
frame2 %>% select(scan, inv_ion_mobility) %>% unique() %>% mutate(Diff = lag(inv_ion_mobility) - inv_ion_mobility)


frame1 %>% select(scan, inv_ion_mobility) %>% unique() %>% tail()
frame2 %>% select(scan, inv_ion_mobility) %>% unique() %>% tail()
table2dt(D, 'Frames')





F1moser2 <- read_csv("/home/momi/DATA/SLICE_PASEF/1F_moser_method2.txt")
F1_slice1 <- read_csv("/home/momi/DATA/SLICE_PASEF/slice_test.txt")

names(F1_slice1) <- c("ms_type", "cycle_id", "IM_start", "IM_end", "mz_start", "mz_end", "CE")
names(F1moser2) <- c("ms_type", "cycle_id", "IM_start", "IM_end", "mz_start", "mz_end", "CE")

F1_slice1 <- F1_slice1 %>% mutate(IM_start = as.numeric(IM_start),
                                IM_end = as.numeric(IM_end),
                                mz_start = as.numeric(mz_start), 
                                mz_end = as.numeric(mz_end))  
F1_slice1 <- F1_slice1 %>% filter(cycle_id == 1)
F1_slice1 %>% view

# steep curve
F1_steep <- read_csv("/home/momi/DATA/SLICE_PASEF/slice_steep1.txt")
F1_steep <- read_csv("/media/momi/UBUNTU 20_0/TIMS_TOF/SLICE_PASEF/reverse_slice_steep.txt")
names(F1_steep) <- c("ms_type", "cycle_id", "IM_start", "IM_end", "mz_start", "mz_end", "CE")
F1_steep <- F1_steep %>% mutate(IM_start = as.numeric(IM_start),
                                  IM_end = as.numeric(IM_end),
                                  mz_start = as.numeric(mz_start), 
                                  mz_end = as.numeric(mz_end))  

F1_steep


F1_flat <- read_csv("/media/momi/UBUNTU 20_0/TIMS_TOF/SLICE_PASEF/reverse_slice_flat.txt")
names(F1_flat) <- c("ms_type", "cycle_id", "IM_start", "IM_end", "mz_start", "mz_end", "CE")
F1_flat <- F1_flat %>% mutate(IM_start = as.numeric(IM_start),
                                IM_end = as.numeric(IM_end),
                                mz_start = as.numeric(mz_start), 
                                mz_end = as.numeric(mz_end))  

F1_steep


F1moser2 <- F1moser2 %>% mutate(IM_start = as.numeric(IM_start),
       IM_end = as.numeric(IM_end),
       mz_start = as.numeric(mz_start), 
       mz_end = as.numeric(mz_end))  

F1moser2

F1_slice1
F1_steep <- F1_steep %>% mutate(cycle_id = 2)
F1_flat <- F1_flat %>% mutate(cycle_id = 3)

################

# PLOT single Frame


frameN %>%   
  rename("mz_start" = "mz", "IM_start" = "inv_ion_mobility") %>% 

ggplot(aes(mz_start, IM_start)) + 
  geom_point(size = 0.1) +
  coord_fixed(ratio = 500) +
  geom_rect(data = dia_windows, mapping=aes(xmin=mz_start, xmax=mz_end, ymin=IM_start, ymax=IM_end, fill=factor(DiaFrameMsMsWindows.WindowGroup)), color="white", alpha=0.5) +
  #geom_rect(data = F1_slice1, mapping=aes(xmin=mz_start, xmax=mz_end, ymin=IM_start, ymax=IM_end, fill = factor(cycle_id)), color="white", alpha=0.3) +
  #geom_rect(data = F1_steep, mapping=aes(xmin=mz_start, xmax=mz_end, ymin=IM_start, ymax=IM_end, fill = factor(cycle_id)), color="white", alpha=0.3) +
  #geom_rect(data = F1_flat, mapping=aes(xmin=mz_start, xmax=mz_end, ymin=IM_start, ymax=IM_end, fill = factor(cycle_id)), color="white", alpha=0.3) +
  scale_fill_viridis(discrete = TRUE) +
  guides(fill=guide_legend(title="cycles"))

  


  
# geom_density2d

  frame1 %>%   
    rename("mz_start" = "mz", "IM_start" = "inv_ion_mobility") %>% 
    
    ggplot(aes(mz_start, IM_start)) +
    geom_density_2d_filled(contour_var = "count")+
    #geom_rect(data = dia_windows, mapping=aes(xmin=mz_start, xmax=mz_end, ymin=IM_start, ymax=IM_end, fill=factor(DiaFrameMsMsWindows.WindowGroup)), color="white", alpha=0.5) +
    scale_fill_viridis(discrete = TRUE) 
  
  

  
  
  
      
    stat_density_2d(
      geom = "raster",
      aes(fill = after_stat(density)),
      contour = FALSE
    )+
    
    #geom_density_2d_filled() +
    #coord_fixed(ratio = 1000) +
    geom_rect(data = dia_windows, mapping=aes(xmin=mz_start, xmax=mz_end, ymin=IM_start, ymax=IM_end, fill=factor(DiaFrameMsMsWindows.WindowGroup)), color="white", alpha=0.5) 
    
  scale_fill_viridis(discrete = TRUE) 
  scale_fill_brewer(palette = "Spectral")
  







###################################################################
#   ANIMATION
# 12 frames for all MSMS windows

aw
  
out_t <- tibble()

frame_order <- seq(1,180,9)
frame_order <- c(1,2,1,3,1,4,5,1,6,1,7)

counter <- 0

for (fram in frame_order){
  counter <- counter + 1
#for (fram in seq(1,13)){
  frame <- D[fram, all_columns]
  frame <- frame %>% mutate("windowgroup" = fram) %>% 
    #filter(mz > 322.04 & mz < 322.05) %>% 
    mutate(frame_n = counter)
  print(fram)
  out_t <- rbind(out_t, frame)
}

dia_windows

out_t %>% left_join(dia_windows, by = c("windowgroup" = "DiaFrameMsMsWindows.WindowGroup"))


dia_windows

out_t %>% 
select(mz, frame_n, inv_ion_mobility, intensity, windowgroup) %>% 
  rename("mz_start" = "mz", "IM_start" = "inv_ion_mobility") %>% 

ggplot(aes(mz_start, IM_start, colour = intensity, group = interaction(mz_start, windowgroup))) + # interaction to avoid fading
  geom_point(size = 3) +
  #xlim(322.045, 322.05) +
  #coord_fixed(ratio = 1000) +
  #geom_rect(data = dia_windows, mapping=aes(xmin=mz_start, xmax=mz_end, ymin=IM_start, ymax=IM_end, fill=DiaFrameMsMsWindows.WindowGroup), color="black", alpha=0.5) +
  transition_states(factor(frame_n), 2,1)+
  enter_fade() + 
    exit_fade() +
    labs(title = 'MS1 frame {closest_state}')





p



anim <- ggplot(mtcars, aes(factor(gear), mpg)) +
  geom_boxplot() +
  transition_states(gear, 2, 1)

anim





frame2 %>% filter(mz > 322 & mz < 322.05) %>% view


names(frame2)

frame2 %>% ggplot(aes(mz, inv_ion_mobility, colour = intensity)) + 
  #geom_density_2d_filled(contour_var = "count")+
  geom_point() + 
  #geom_line()+ 
  xlim(322.045, 322.05)




#MS2

f2 <- frame2 %>% select(mz, inv_ion_mobility, intensity) %>% 
  rename("mz_start" = "mz", "inv_ion_mobility_start" = "inv_ion_mobility") %>% 
  ggplot(aes(mz_start, inv_ion_mobility_start, color = intensity)) + 
  geom_point(size = 0.1) +
   coord_fixed(ratio = 1000) +
  geom_rect(data = dia_windows, mapping=aes(xmin=mz_start, xmax=mz_end, ymin=inv_ion_mobility_start, ymax=inv_ion_mobility_end, fill=DiaFrameMsMsWindows.WindowGroup), color="black", alpha=0.5) 

  
f2

plotly::ggplotly(f2)  

#MS1  
frame2

f1 <- frame2 %>% 
    ggplot(aes(mz, inv_ion_mobility, color = intensity)) + 
    geom_point(size = 0.1) +
  ylim(c(1.42, 1.45))

f1
plotly::ggplotly(f1)


