###############################################################

## Toward evidence-based conservation of subterranean ecosystems
## Mammola S., Meierhofer M., ... Cardoso, P.

## ------------------------------------------------------------------------
# 'R script to reproduce the analyses'
## ------------------------------------------------------------------------

# Analysis performed with R (v. R 4.0.3) and R studio (v. 1.4.1103)
# Authors: Stefano Mammola & Melissa Meierhofer
# Location: Helsinki, May-Jul 2021

###############################################################

# clean the workspace -----------------------------------------------------

rm(list = ls())

# Working directory -------------------------------------------------------

setwd("/Users/stefanomammola/Desktop/Practical conservation review/analysis") #change with your working directory

# Loading R package -------------------------------------------------------

library("circlize")      
library("cowplot")
library("dplyr")         
library("ggplot2")       
library("ggpubr")
library("grid")
library("gridExtra")
library("maps")
library("parameters")
library("performance")
library("scatterpie")
library("tidyr")

# Loading useful functions ------------------------------------------------

# Custom function to get standard error
SE <- function(x) sd(x) / sqrt( length(x) )

# Custom function to predict logistic regressions
logisticline <- function(z,model) {
  eta <- model$coefficients[1] + model$coefficients[2]*z ;
  1 / (1 + exp(-eta))
}

logisticline_min <- function(z,model) {
  eta <- model$coefficients[1] + model$coefficients[2]*z - 1.96*summary(model)$coefficients[2] ;
  1 / (1 + exp(-eta))
}

logisticline_max <- function(z,model) {
  eta <- model$coefficients[1] + model$coefficients[2]*z + 1.96*summary(model)$coefficients[2] ;
  1 / (1 + exp(-eta))
}

# Custom function to split columns having semicolon as a separator
semi_colon_splitter <- function(input1, input2, names = c("input1","input2")){
  
  df        <- data.frame(input1,input2)  
  df$input1 <- as.factor(df$input1)
  df$input2 <- as.factor(df$input2)
  
  to_separate <- levels(df$input1)[grepl(";", levels(df$input1))]
  
  df_all <- df[df$input1 %in% to_separate ,]
  df     <- df[!df$input1 %in% to_separate,]
  df$input1 <- droplevels(df$input1)
  
  df_all$input1 <- as.character(df_all$input1)
  
  for(i in nrow(df_all)) {
    
    df_i <- df_all[i,]
    split   <- strsplit(df_all$input1, ";")[[1]]
    split   <- trimws(split, which = c("both"))
    
    df <- rbind(df,data.frame(input1  = split,
                              input2  = rep(df_i$input2, length(split))))
                              
  }
  
  colnames(df) <- names
  return(df)
}

# Parameters for plots ----------------------------------------------------

# Plot style. Modified from:
# https://ourcodingclub.github.io/tutorials/dataviz-beautification-synthesis/

theme_niwot <- function(){
  theme_bw() +
    theme(#text = element_text(family = "Arial"),
      axis.text = element_text(size = 10), 
      axis.title = element_text(size = 12),
      axis.line.x = element_line(color="black"), 
      axis.line.y = element_line(color="black"),
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),                                          
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.y = element_blank(),  
      plot.margin = unit(c(1, 1, 1, 1), units = , "cm"),
      plot.title = element_text(size = 18, vjust = 1, hjust = 0),
      legend.text = element_text(size = 12),          
      legend.title = element_blank(),                              
      legend.position = c(0.95, 0.15), 
      legend.key = element_blank(),
      legend.background = element_rect(color = "black", 
                                       fill = "transparent", 
                                       size = 2, linetype = "blank"))
}

# Colors for Figure 2
col_fig2 <- c("grey10","darkorchid3")

# Colors ad parameters for Figure 4
COL  <- c(rep("grey40",4),"darkgoldenrod2","brown4")

COL2 <- c("blue", "palevioletred4",
          rep("grey40",2), #n.s.
          "seagreen4", 
          rep("grey40",2), #n.s.
          "black")

alpha1 <- 0.2
alpha2 <- 0.1
alpha3 <- 0.6

###############################################################

## Data preparation:

###############################################################

# Loading the Database ----------------------------------------------------

db <- read.csv(file = "Data/Database_Practical_conservation.csv", sep = '\t', dec = '.', header = TRUE, as.is = FALSE)

# Information on the columns are in the file METADATA_Database_Practical_conservation.pdf  found in the same folder of the database): 

str(db)
summary(db)

#Database only with distinct paper
db_unique <- distinct(db, ID, .keep_all = TRUE) 

#Summary statistics
table(db_unique$Source) ; sum(table(db_unique$Source)) # N° of unique sources
table(db$Tested_statistically) ; table(db$Tested_statistically)[2] / sum(table(db$Tested_statistically)) #N° and % testing
mean(table(db$ID)) ; SE(table(db$ID)) # mean number of actions/paper

round(table(db$Tested_statistically == "yes",db$Conservation_Action)[2,]/table(db$Conservation_Action),2)#% testing by conservation action

#Redefining impact
db$Impact2 <- db$Impact

levels(db$Impact2)   <- c("Alien species\nPathogens",
                          "All",
                          "Climate\nchange",
                          "None",
                          "Alien species\nPathogens",
                          "Overexploitation\nPoaching",
                          "Pollution",
                          "Habitat change\n(subterranean)",
                          "Habitat change\n(surface)",
                          "Visitors")

levels(db$Taxon_Group)
levels(db$Tested_statistically)

#####################################
### FIGURE 2:::::::::::::::::::::::::
#####################################

# Map ---------------------------------------------------------------------

pie1 <- semi_colon_splitter(input1 = db$Higher_Geography,
                            input2 = db$Tested_statistically, 
                            names = c("Higher_Geography","Tested_statistically"))

#summary stats
table(pie1$Higher_Geography) # tot
table(pie1$Higher_Geography)/ sum(table(pie1$Higher_Geography)) # %

pie_1 <- data.frame(table(pie1$Higher_Geography,pie1$Tested_statistically))

radius <- data.frame(sqrt(table(pie1$Higher_Geography))) ; rownames(radius) <- NULL
n <- data.frame(table(pie1$Higher_Geography)) ; rownames(radius) <- NULL

pie <- data.frame(region = levels(pie_1$Var1),
                  n =  n$Freq,
                  long  = c(35, 140.0, -154.0, 82.0, -106.0, -57.0, 22.0),
                  lat  = c(2, -23.5, 21, 27.0, 49.0, -17.0, 49.0),
                  no = pie_1$Freq[1:nlevels(pie_1$Var1)],
                  yes = pie_1$Freq[c(nlevels(pie_1$Var1)+1):c(nrow(pie_1))],
                  radius = radius$Freq)

# Using map_data()

world <- map_data("world")

map1 <- ggplot() +
    geom_map(map = world, data = world,
             aes(map_id = region), 
             color = "gray65", fill = "gray65", size = 0.3) +
    labs(title = NULL) +
    
    #global
    annotate(geom="text", x=-154, y=35, label="Global",
             color="black")+
    
    #Palearctic
    annotate(geom="text", x=24, y=84, label="Palearctic",
             color="black")+
    
    #Nearctic
    annotate(geom="text", x=-108, y=84, label="Nearctic",
             color="black")+
    
    #Neotropical
    annotate(geom="text", x=-26.5, y=-28, label="Neotropical",
             color="black")+
    
    #Afrotropical
    annotate(geom="text", x=62, y=-3, label="Afrotropical",
             color="black")+
    
    #Afrotropical
    annotate(geom="text", x=86, y=39, label="Indomalaysian",
             color="black")+
    
    #Australasian
    annotate(geom="text", x=178, y=-29, label="Australasian",
             color="black")+
    
    theme_bw() +
    theme(
      axis.text.x  = element_blank(), 
      axis.text.y  = element_blank(),
      axis.title.y = element_blank(),
      axis.title.x = element_blank(), 
      axis.line.x = element_blank(), 
      axis.line.y = element_blank(),
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),                                          
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.y = element_blank(), 
      axis.ticks = element_blank(),
      plot.margin = unit(c(1,1,1,1), 'cm'),
      plot.title = element_text(size = 18, vjust = 1, hjust = 0),
      legend.text = element_text(size = 12),          
      legend.title = element_blank(),                              
      legend.position = c(0.1, 0.2), 
      legend.key = element_blank(),
      legend.background = element_rect(color = "black", 
                                       fill = "white", 
                                       size = 2, linetype = "blank"))

(map2 <- map1 + geom_scatterpie(aes(x=long, y=lat, group=region, r=radius),
                                data=pie, cols = c("no","yes"), color=NA, alpha=.9) +
    geom_scatterpie_legend(pie$radius, x=-150, y=-45, n = 2,
                           labeller = function (x) x=c(min(pie$n),max(pie$n)))+
    #scale_fill_manual(values=col_fig2)
    scale_fill_manual("",labels=c("Not tested", "Tested"),values=col_fig2)
  + theme(legend.position = "top"))

# Barplots ----------------------------------------------------------------

## Tested statistically by Impact
bar_1 <- data.frame(table(db$Tested_statistically,db$Impact2))

levels(bar_1$Var2)[2] <- "Multiple"

#summary stats
table(db$Impact2) # tot
table(db$Impact2)/ sum(table(db$Impact2)) # %

bar_1$Var2 <- factor(bar_1$Var2,levels = c("Alien species\nPathogens","Climate\nchange","Overexploitation\nPoaching",
                                           "Pollution","Habitat change\n(subterranean)","Habitat change\n(surface)",     
                                           "Visitors", "Multiple", "None"))  

(bar_p1 <-  ggplot(bar_1, aes(x=Var2,y=Freq, fill=Var1)) +
    
    geom_bar(stat="identity",position=position_dodge(), color = "grey20")+
    geom_text(aes(label=Freq), vjust=-1, color="black",
              position = position_dodge(0.9), size=2.5)+
    
    ylim(0,250)+
    scale_fill_manual("",labels=c("Not tested", "Tested"),values=col_fig2)+
    labs(title=NULL, subtitle = NULL,x=NULL, y = "Frequency")+
    theme_niwot()+
    theme(legend.position =  "none",
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          plot.margin = unit(c(0.2,0.2,0.2,0.2), 'cm')) 
)

## Tested statistically by System
bar_3 <- semi_colon_splitter(input1 = db$System,
                             input2 = db$Tested_statistically, 
                             names = c("System","Tested_statistically"))
#summary stats
table(bar_3$System) # tot
table(bar_3$System) / sum(table(bar_3$System)) # %
sum(table(bar_3$System)[c(3,4,7)]) / sum(table(bar_3$System)) #terrestial systems (sum of show caves, caqves and artificial)

bar_3 <- data.frame(table(bar_3$Tested_statistically,bar_3$System))

#rename levels
levels(bar_3$Var2)[1] <- "Not specific"
levels(bar_3$Var2)[2] <- "Anchialine\n& Marine"
levels(bar_3$Var2)[5] <- "Fissural\nsystems"

#reoder
bar_3$Var2 <- factor(bar_3$Var2,levels = c("Anchialine\n& Marine","Groundwater",
                                           "Cave","Show cave","Fissural\nsystems",
                                           "Artificial",
                                           "Not specific"))

(bar_p2 <-  ggplot(bar_3, aes(x=Var2,y=Freq, fill=Var1)) +
    
    geom_bar(stat="identity",position=position_dodge(), color = "grey20")+
    geom_text(aes(label=Freq), vjust=-1, color="black",
              position = position_dodge(0.9), size=2.5)+
    
    ylim(0,750)+
    
    scale_fill_manual("",labels=c("Not tested", "Tested"),values=col_fig2)+
    labs(title=NULL, subtitle = NULL,x=NULL, y = "Frequency")+
    theme_niwot()+
    theme(legend.position =  "none",
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          plot.margin = unit(c(0.2,0.2,0.2,0.2), 'cm')) 
  )

## Tested statistically by taxon
bar_4 <- semi_colon_splitter(input1 = db$Taxon_Group,
                             input2 = db$Tested_statistically, 
                             names  = c("Taxon_Group","Tested_statistically"))

#summary stats
table(bar_4$Taxon_Group) # tot
table(bar_4$Taxon_Group) / sum(table(bar_4$Taxon_Group)) # %

bar_4 <- data.frame(table(bar_4$Tested_statistically, bar_4$Taxon_Group))

bar_4$Var2 <- factor(bar_4$Var2, levels = c("Arthropoda","Other invertebrates",
                                           "Bats","Other vertebrates",
                                           "Plants","Microorganisms",
                                           "Not specific"))

(bar_p3 <-  ggplot(bar_4, aes(x=Var2, y=Freq, fill=Var1)) +
    
    geom_bar(stat="identity",position=position_dodge(), alpha=1, color = "grey20")+
    geom_text(aes(label=Freq), vjust=-1, color="black",
              position = position_dodge(0.9), size=2.5)+
    ylim(0,400)+
    
    scale_fill_manual("",labels=c("Not tested", "Tested"),values=col_fig2)+
    labs(title=NULL, subtitle = NULL,x=NULL, y = NULL)+
    theme_niwot()+
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          plot.margin = unit(c(0.2,0.2,0.2,0.2), 'cm'))
  )

## Tested statistically by conservation action
bar_5 <- data.frame(table(db$Tested_statistically,db$Conservation_Action))

#summary stats
table(db$Conservation_Action) # tot
table(db$Conservation_Action) / sum(table(db$Conservation_Action)) # %

(bar_p4 <-  ggplot(bar_5, aes(x=Var2,y=Freq, fill=Var1)) +
    
    geom_bar(stat="identity",position=position_dodge(), alpha=1,color = "grey20")+
    geom_text(aes(label=Freq), vjust=-1, color="black",
              position = position_dodge(0.9), size=2.5)+
    ylim(0,400)+
    
    scale_fill_manual("",labels=c("Not tested", "Tested"),values=col_fig2)+
    labs(title=NULL, subtitle = NULL,x=NULL, y = NULL)+
    theme_niwot()+
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,
                                     face = c(rep("plain",2),"italic",rep("plain",14))),
          plot.margin = unit(c(0.2,0.2,0.2,0.2), 'cm'))
)

# Arrange on a grid -------------------------------------------------------

pdf(file = "Figure/Figure_2.pdf", width = 10, height =11)

ggpubr::ggarrange(map2,
                  ggpubr::ggarrange(bar_p2, 
                                    bar_p3, 
                  bar_p1, 
                  bar_p4, 
                  ncol = 2, nrow = 2, align = "hv",  labels = c("B","C","D","E")), 
                  common.legend = FALSE,
                  labels = "A",
                  nrow=2)

dev.off()

#####################################
### FIGURE 3:::::::::::::::::::::::::
#####################################

# Chord diagrams ----------------------------------------------------------

## CHORD DIAGRAM 1 - Conservation action vs Impact
## https://www.data-to-viz.com/graph/chord.html

chord_plot1  <- function(){
  require("gridGraphics")
  
  chord <- data.frame(Conservation_Action = db$Conservation_Group, Impact = db$Impact)
  
  ## Adding 'all' values
  
  chord_all <- chord[chord$Impact  == "All",]
  chord     <- chord[!chord$Impact == "All",] ; chord$Impact <- droplevels(chord$Impact)
  
  for(i in nrow(chord_all)) {
    
    chord_i <- chord_all[i,]
    
    chord_i2 <- data.frame(Conservation_Action = rep(chord_i$Conservation_Action, nlevels(chord$Impact)),
                           Impact  = levels(chord$Impact))
    
    chord_i2 <- chord_i2[!chord_i2$Impact == "None",]
    
    chord <- rbind(chord,chord_i2)
    
  }
  
  #renaming levels
  levels(chord$Impact) <- c("Alien species\nPathogens",
                            "Climate\nchange",
                            "None",
                            "Alien species\nPathogens",
                            "Overexploitation\n& Poaching",
                            "Pollution",
                            "    Habitat change\n    (subterranean)",
                            "Habitat change\n(surface)",
                            "Visitors")
  
  mat <- data.frame(table(chord$Conservation_Action, chord$Impact))
  
  colnames(mat) <- c("from", "to", "value")
  
  levels(mat$from)
  levels(mat$to)
  
  nlevels(mat$from)
  nlevels(mat$to)
  
  order_l <- c("Education", "Restoration", "Regulation",
               "Monitoring", "Assessment", "Protection", 
               "Alien species\nPathogens","Pollution", 
               "Climate\nchange", "Habitat change\n(surface)",
               "Visitors" ,"Overexploitation\n& Poaching",
               "    Habitat change\n    (subterranean)","None")
  
  col_l <- c("lightcyan4","brown4","darkgoldenrod2","darkblue","cyan4","blueviolet",rep("grey30",nlevels(mat$to)))
  
  circos.clear()
  circos.par(start.degree =  -8,
             gap.after=c(rep(3,nlevels(mat$from)-1),10,rep(3,nlevels(mat$to)-1),10))
  
  circlize::chordDiagram(mat,
                         grid.border =1,
                         order = order_l,
                         grid.col = col_l
  )
}

chord_plot1()

#Convert as ggplot2 object
chord_plot1  <- cowplot::as_grob(chord_plot1)
chord_plot1  <- ggpubr::as_ggplot(chord_plot1)


pdf(file = "Figure/Figure_3_no_silouhette.pdf", width = 10, height = 10)
par(mar=c(2,2,2,2))
ggpubr::ggarrange(chord_plot1, nrow= 1, ncol = 1)
dev.off()

# Silouhettes where prepared by Irene Frigo with Adobe Illustrator and have been added to Figure 3 outside R, using the program Inkscape.


#####################################
### FIGURE 4:::::::::::::::::::::::::
#####################################

# Modeling temporal trends ---------------------------------------

#Overall temporal trend in testing
bar_1 <- data.frame(table(db_unique$Year_publication,db_unique$Tested_statistically)) ; colnames(bar_1) <- c("yr","Tested","N")
bar_1$yr <- as.numeric(as.character(bar_1$yr))

# Has the frequency of tested paper changed over time?
glm <- data.frame(yr = unique(bar_1$yr),
                  tested = bar_1[bar_1$Tested=="yes",]$N, 
                  untested = bar_1[bar_1$Tested=="no",]$N)

glm <- glm[glm$yr > 1999,] #selecting last 20 years

m1  <- glm(cbind(tested,untested) ~ yr, data = glm, family = "binomial")

performance::check_overdispersion(m1)
performance::check_model(m1)

(pM1 <- parameters::model_parameters(m1))

y <- seq(from = min(glm$yr), to = max(glm$yr), 1)

(inset <- 
    ggplot()+
    xlab(NULL)+ 
    ylab("Proportion tested")+
    geom_line(aes(y = logisticline(y, m1), x = y), colour = "blueviolet",linetype="solid",size=1)+
    
    geom_ribbon(aes(ymax = logisticline_max(y, m1),
                    ymin = logisticline_min(y, m1),x = y),alpha = 0.3,fill="black")+
    
    geom_point(aes(y= ((glm$tested)/(glm$tested + glm$untested)),x=glm$yr), col = "black", alpha = 0.5)+
    
    annotate("text", x = 2011, y = 0.35,
             label =  paste("GLM (n= 22): ",
                            round(pM1$Coefficient[2],2),
                            " ± ",
                            round(pM1$SE[2],2),
                            "; p= ",
                            round(pM1$p[2],3),
                            sep=''),size = 4)+
    
    theme_niwot() )
  
#add missing years for the plot
yr <- seq(from = range(as.numeric(as.character(bar_1$yr)))[1], to = range(as.numeric(as.character(bar_1$yr)))[2], by = 1)
yr <- yr[!yr %in% bar_1$yr]

yr2 <- data.frame(yr = rep(yr,2),
                  Tested = c(rep("no",length(yr)),rep("yes",length(yr))),
                  N = rep(0,2*length(yr)))

bar_1 <- rbind(bar_1,yr2)

#plot
Plot_trend0 <- ggplot(bar_1, aes(x=yr, y=as.numeric(N), fill=Tested))+
  geom_bar(stat="identity",color = "grey20")+
  scale_x_continuous(breaks = seq(from = range(bar_1$yr)[1], 
                                  to = range(bar_1$yr)[2], by = 3),
                     labels = as.character(seq(from = range(bar_1$yr)[1],
                                               to = range(bar_1$yr)[2], by = 3))
                     
  )+ 
  scale_fill_manual("",labels=c("Not tested", "Tested"),values=col_fig2)+
  labs(title=NULL, 
       x=NULL, 
       y = "Frequency",
       subtitle = NULL)+
  theme_niwot()+
  theme(legend.position = "top",
        plot.subtitle = element_text(size = 12, color = "#222222"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), 'cm'))

#add the inset

(Plot_trend0 <- Plot_trend0  + annotation_custom(ggplotGrob(inset),xmin=1965,xmax=1990,ymin=20,ymax=60))

## B - Conservation actions

tot <- data.frame(table(db$Year_publication)) #total number of publication/year

# Creating a database
db_yr2    <- data.frame(table(db$Year_publication,db$Conservation_Group))
db_yr2    <- data.frame(db_yr2, tot = rep(tot$Freq,6)) ; colnames(db_yr2) <- c("yr","Cons","N","tot")
db_yr2$yr <- as.numeric(as.character(db_yr2$yr))

# Selecting papers from 2000 onward
db_yr2 <- db_yr2[db_yr2$yr > 1999,]

# Modelling the temporal trends
model   <- list()
par     <- list()

for (i in levels(factor(db_yr2$Cons))) {
  
  db_i <- db_yr2[db_yr2$Cons==i, ]
  model[[i]]   <- glm(cbind(N,tot) ~ yr, data = db_i, family = "binomial")
  par[[i]] <- parameters::model_parameters(model[[i]])
  
}  

# Model validation
performance::check_model(model[[1]])
performance::check_model(model[[2]])
performance::check_model(model[[3]])
performance::check_model(model[[4]])
performance::check_model(model[[5]])
performance::check_model(model[[6]])

# Model summary 

for (i in 1:nlevels(factor(db_yr2$Cons))) {
  
  message(paste("::::::  ",levels(factor(db_yr2$Cons))[i],"  :::::"))
  
  print(parameters::model_parameters(model[[i]]))   
  
}  #Restoration and regulation significant

#Temporal series of interest
y2 <- seq(from = min(db_yr2$yr), to = max(db_yr2$yr), 1) #temporal series of interest

(Plot_trend1 <- ggplot() +
    ylab(NULL) + xlab(NULL)+ ylim(0,0.3)+
    #trend lines
    geom_line(aes(y = logisticline(y2,model[[1]]), x = y2), colour = COL[1],linetype="solid",size=1.1,alpha=alpha1)+
    geom_line(aes(y = logisticline(y2,model[[2]]), x = y2), colour = COL[2],linetype="solid",size=1.1,alpha=alpha1)+
    geom_line(aes(y = logisticline(y2,model[[3]]), x = y2), colour = COL[3],linetype="solid",size=1.1,alpha=alpha1)+
    geom_line(aes(y = logisticline(y2,model[[4]]), x = y2), colour = COL[4],linetype="solid",size=1.1,alpha=alpha1)+
    geom_line(aes(y = logisticline(y2,model[[5]]), x = y2), colour = COL[5],linetype="solid",size=1.1)+
    geom_line(aes(y = logisticline(y2,model[[6]]), x = y2), colour = COL[6],linetype="solid",size=1.1)+
    #confidence intervals
    geom_ribbon(aes(ymax = logisticline_max(y2, model[[01]]),
                    ymin = logisticline_min(y2, model[[01]]),x = y2),alpha = alpha2,fill=COL[1])+
    geom_ribbon(aes(ymax = logisticline_max(y2, model[[02]]),
                    ymin = logisticline_min(y2, model[[02]]),x = y2),alpha = alpha2,fill=COL[2])+
    geom_ribbon(aes(ymax = logisticline_max(y2, model[[03]]),
                    ymin = logisticline_min(y2, model[[03]]),x = y2),alpha = alpha2,fill=COL[3])+
    geom_ribbon(aes(ymax = logisticline_max(y2, model[[04]]),
                    ymin = logisticline_min(y2, model[[04]]),x = y2),alpha = alpha2,fill=COL[4])+
    geom_ribbon(aes(ymax = logisticline_max(y2, model[[05]]),
                    ymin = logisticline_min(y2, model[[05]]),x = y2),alpha = 0.5, fill=COL[5])+
    geom_ribbon(aes(ymax = logisticline_max(y2, model[[06]]),
                    ymin = logisticline_min(y2, model[[06]]),x = y2),alpha = 0.5, fill=COL[6])+
    
    #Text
    annotate(geom="text", hjust = 0,vjust = 0.3,
             x= 2021.5, y= logisticline_max(y2, model[[01]])[21], 
             label = levels(factor(db_yr2$Cons))[1],
             color="black",alpha=alpha3)+
    
    annotate(geom="text", hjust = 0,vjust = 0,
             x= 2021.5, y= logisticline_max(y2, model[[02]])[21], 
             label = levels(factor(db_yr2$Cons))[2],
             color="black",alpha=alpha3)+
    
    annotate(geom="text", hjust = 0,vjust = 0,
             x= 2021.5, y= logisticline_max(y2, model[[03]])[21], 
             label = levels(factor(db_yr2$Cons))[3],
             color="black",alpha=alpha3)+
    
    annotate(geom="text", hjust = 0,vjust = 0,
             x= 2021.5, y= logisticline_max(y2, model[[04]])[21], 
             label = levels(factor(db_yr2$Cons))[4],
             color="black",alpha=alpha3)+
    
    annotate(geom="text", hjust = 0,vjust = 0,
             x= 2021.5, y= logisticline_max(y2, model[[05]])[21], 
             label = levels(factor(db_yr2$Cons))[5],
             color=COL[5])+
    
    annotate(geom="text", hjust = 0,vjust = 0,
             x= 2021.5, y= logisticline_max(y2, model[[06]])[21], 
             label = levels(factor(db_yr2$Cons))[6],
             color=COL[6])+
    
    #annotation_custom(rasterGrob(img_Assessment), xmin = 2022, ymin = 0.2) +
    
    coord_cartesian(xlim = c(2000, 2021), # This focuses the x-axis on the range of interest
                    clip = 'off') +   # This keeps the labels from disappearing
    
    theme_niwot() + theme(plot.margin = unit(c(0.5,4,0.5,0.5), 'cm'))
)

## C - Impacts

# Creating a database
db_yr3 <- data.frame(table(db$Year_publication,db$Impact2))
db_yr3 <- data.frame(db_yr3, tot = rep(tot$Freq, nlevels(db_yr3$Var2))) ; colnames(db_yr3) <- c("yr","Impact","N","tot")
db_yr3$yr <- as.numeric(as.character(db_yr3$yr))

# Selecting years and removing all
db_yr3 <- db_yr3[db_yr3$yr > 1999,]
db_yr3 <- db_yr3[!db_yr3$Impact == "All",]

# Modelling the temporal trends
model2 <- list()
par2   <- list()

for (i in levels(factor(db_yr3$Impact))) {
  
  db_i <- db_yr3[db_yr3$Impact==i, ]
  model2[[i]] <- glm(cbind(N,tot) ~ yr, data = db_i, family = "binomial")
  par2[[i]] <- parameters::model_parameters(model2[[i]])
  
}  

# Model summary 

for (i in 1:nlevels(factor(db_yr3$Impact))) {
  
  message(paste("::::::  ",levels(factor(db_yr3$Impact))[i],"  :::::"))
  
  print(parameters::model_parameters(model2[[i]]))   
  
}  #Pollution, Pathogens, Visitors significant. Climate change borderline significant.

#Renaming levels
levels(db_yr3$Impact)[1] <- "Alien species\n& Pathogens"
levels(db_yr3$Impact)[3] <- "Climate change"
levels(db_yr3$Impact)[5] <- "Overexploitation\n& Poaching"

(Plot_trend2 <- ggplot() +
    ylab("Proportion of conservation actions") + xlab(NULL) + ylim(0,0.3)+
    #trend lines
    geom_line(aes(y = logisticline(y2,model2[[03]]), x = y2), colour = COL2[3],linetype="solid",size=1.1,alpha=alpha1)+
    geom_line(aes(y = logisticline(y2,model2[[04]]), x = y2), colour = COL2[4],linetype="solid",size=1.1,alpha=alpha1)+
    geom_line(aes(y = logisticline(y2,model2[[06]]), x = y2), colour = COL2[6],linetype="solid",size=1.1,alpha=alpha1)+
    geom_line(aes(y = logisticline(y2,model2[[07]]), x = y2), colour = COL2[7],linetype="solid",size=1.1,alpha=alpha1)+
    
    geom_line(aes(y = logisticline(y2,model2[[08]]), x = y2), colour = COL2[8],linetype="solid",size=1.1)+
    geom_line(aes(y = logisticline(y2,model2[[01]]), x = y2), colour = COL2[1],linetype="solid",size=1.1)+
    geom_line(aes(y = logisticline(y2,model2[[05]]), x = y2), colour = COL2[5],linetype="solid",size=1.1)+
    geom_line(aes(y = logisticline(y2,model2[[02]]), x = y2), colour = COL2[2],linetype="solid",size=1.1)+
    
    #confidence intervals
    geom_ribbon(aes(ymax = logisticline_max(y2, model2[[03]]),
                    ymin = logisticline_min(y2, model2[[03]]),x = y2),alpha = alpha2,fill=COL2[3])+
    geom_ribbon(aes(ymax = logisticline_max(y2, model2[[04]]),
                    ymin = logisticline_min(y2, model2[[04]]),x = y2),alpha = alpha2,fill=COL2[4])+
    geom_ribbon(aes(ymax = logisticline_max(y2, model2[[06]]),
                    ymin = logisticline_min(y2, model2[[06]]),x = y2),alpha = alpha2,fill=COL2[6])+
    geom_ribbon(aes(ymax = logisticline_max(y2, model2[[07]]),
                    ymin = logisticline_min(y2, model2[[07]]),x = y2),alpha = alpha2,fill=COL2[7])+
    
    geom_ribbon(aes(ymax = logisticline_max(y2, model2[[08]]),
                    ymin = logisticline_min(y2, model2[[08]]),x = y2),alpha = 0.4,fill=COL2[8])+
    geom_ribbon(aes(ymax = logisticline_max(y2, model2[[01]]),
                    ymin = logisticline_min(y2, model2[[01]]),x = y2),alpha = 0.4,fill=COL2[1])+
    geom_ribbon(aes(ymax = logisticline_max(y2, model2[[05]]),
                    ymin = logisticline_min(y2, model2[[05]]),x = y2),alpha = 0.4,fill=COL2[5])+
    geom_ribbon(aes(ymax = logisticline_max(y2, model2[[02]]),
                    ymin = logisticline_min(y2, model2[[02]]),x = y2),alpha = 0.4,fill=COL2[2])+
    
    #Text
    annotate(geom="text", hjust = 0,vjust = 0,
             x= 2021.5, y= logisticline_max(y2, model2[[04]])[21]-0.015, 
             label = levels(factor(db_yr3$Impact))[4],
             color=COL2[4],alpha=alpha3)+
    
    annotate(geom="text", hjust = 0,vjust = 0,
             x= 2021.5, y= logisticline_max(y2, model2[[03]])[21]+0.003, 
             label = paste(levels(factor(db_yr3$Impact))[3],"| Habitat change"),
             color=COL2[3],alpha=alpha3)+
    
    annotate(geom="text", hjust = 0,vjust = 0.3,
             x= 2021.5, y= logisticline_max(y2, model2[[01]])[21]-0.044, 
             label = levels(factor(db_yr3$Impact))[1],
             color=COL2[1])+
    
    annotate(geom="text", hjust = 0,vjust = 0.3,
             x= 2021.5, y= logisticline_max(y2, model2[[05]])[21], 
             label = levels(factor(db_yr3$Impact))[5],
             color=COL2[5])+
    
    annotate(geom="text", hjust = 0,vjust = 0.3,
             x= 2021.5, y= logisticline_max(y2, model2[[08]])[21], 
             label = levels(factor(db_yr3$Impact))[8],
             color=COL2[8])+
    
    annotate(geom="text", hjust = 0,vjust = 0.3,
             x= 2021.5, y= logisticline_max(y2, model2[[02]])[21], 
             label = levels(factor(db_yr3$Impact))[2],
             color=COL2[2])+
    
    coord_cartesian(xlim = c(2000, 2021), # This focuses the x-axis on the range of interest
                    clip = 'off') +     # This keeps the labels from disappearing
    
    theme_niwot() + theme(plot.margin = unit(c(0.5,4,0.5,0.5), 'cm'))
)


# Arrange on a grid -------------------------------------------------------

pdf(file = "Figure/Figure_4.pdf", width = 14, height = 10)

ggpubr::ggarrange(Plot_trend0, ggpubr::ggarrange(Plot_trend2, Plot_trend1,
                                                 ncol = 2, nrow = 1, 
                                                 labels  = c("B","C"),align = "hv"),
                  ncol = 1, nrow = 2, 
                  labels  = c("A"), 
                  heights = c(1,1),
                  widths  = c(1,1),
                  align = "v",
                  common.legend = FALSE) 

dev.off()

## End of analysis

## EXTRA
## CHORD DIAGRAM 2 - Conservation action vs Taxonomy

chord_plot2 <- function(){
  require("gridGraphics")
  
  chord2 <- data.frame(Conservation_Action = db$Conservation_Group, Taxonomy = db$Taxon_Group)
  
  ## Duplicating multiple level variables
  chord2 <- semi_colon_splitter(input1 = chord2$Taxonomy,
                                input2 = chord2$Conservation_Action, 
                                names =c("Taxonomy","Conservation_Action"))
  
  # removing none specific
  chord2     <- chord2[!chord2$Taxonomy == "Not specific",] ; chord2$Taxonomy <- droplevels(chord2$Taxonomy)
  chord2$Taxonomy <- as.factor(chord2$Taxonomy)
  
  #renamic levels
  levels(chord2$Taxonomy) <- c("Arthropods","Bats","Microorgansims","Other\ninvertebrates","Other\nvertebrates","Plants")
  
  mat2 <- data.frame(table(chord2$Conservation_Action,chord2$Taxonomy))
  
  colnames(mat2) <- c("from", "to", "value")
  
  levels(mat2$from)
  levels(mat2$to)
  
  nlevels(mat2$from)
  nlevels(mat2$to)
  
  order_l2 <- c("Education", "Restoration", "Regulation",
                "Monitoring", "Assessment", "Protection", 
                levels(mat2$to)[3],
                levels(mat2$to)[1],
                levels(mat2$to)[4],
                levels(mat2$to)[2],
                levels(mat2$to)[5],
                levels(mat2$to)[6])
  
  col_l2 <- c("lightcyan4","brown4","darkgoldenrod2","darkblue","cyan4","blueviolet",rep("grey30",nlevels(mat2$to)))
  
  circos.clear()
  circos.par(gap.after=c(rep(3,nlevels(mat2$from)-1),10,rep(3,nlevels(mat2$to)-1),10))
  
  circlize::chordDiagram(mat2,
                         order = order_l2,
                         grid.col = col_l2)
}

chord_plot2()

#Convert as ggplot2 object
chord_plot2  <- cowplot::as_grob(chord_plot2)
chord_plot2  <- ggpubr::as_ggplot(chord_plot2)

## End of analysis -------------------------------
