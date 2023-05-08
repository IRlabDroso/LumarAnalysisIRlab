

# library(ggplot2)
# library(reshape2)
# library(dplyr)
# library(tseries)
# library(patchwork)
# library(hash)
# library(zoo)
# library(rlang)
# library(ggforce)
# library(tidyverse)
# library(rstatix)
# library(ggpubr)



#' Read the exp_info.csv file
#'
#' @title Read_Exp_info
#'
#' @param file Default = NULL, Path to the exp_info.csv file. If not specified search into the current directory.
#'
#' @examples
#' exp.info = Read_Exp_info
#' exp.info$exp_id
#' @export
Read_Exp_info = function(file=NULL){

  ### find .csv files in the dir ###
  if(is.null(file)){
    doss = list.files(getwd(), pattern = "*.csv")
    file = doss[grep("info",doss,ignore.case = T)]
  }

  ### generate outputs ###
  info.exp<-read.csv(file)
  condNum<-info.exp[1,"Conditions"]
  numOdorant<-info.exp[1,"Odorant.Trials"]
  only.trial<-dplyr::select(info.exp, contains("Trial."))
  odorant<-unlist( only.trial[1:numOdorant])
  conditions = unlist(info.exp[1,c(3:(2+condNum))])
  return(list(exp_info = info.exp,exp_id = info.exp$Exp_id,conditions =conditions ,condNum=condNum,numOdorant=numOdorant,odorant=odorant))
}


####################### Split raw csv by odors ########################


#' add the odorant columns and change names of columns
#'
#' @title Split_CSV
#'
#' @param directory Path to the raw csv file.
#' @param length_exp The duration in seconds of a trial for a single odor.
#'
#' @examples
#' csv = Split_CSV(getwd())
#'
#' @export
Split_CSV = function(directory,length_exp = 184){
  ### Open all csv files in the directory
  setwd(directory)
  doss <- list.files(directory, pattern = "*.csv")

  ### open Exp_info ###
  info.exp = Read_Exp_info()

  ### open Raw data CSV ###
  doss <- doss[!(doss %in% c("Exp_Infos.csv"))]
  csv = doss

  tabl <- read.table(csv, header=TRUE, sep=",")[-1,]
  tabl <- na.omit(tabl)

  ### Keep only the useful columns and rename them ###
  list_of_cols = colnames(tabl)[1]
  new_colnames = c("Timing")
  for(i in 1:((ncol(tabl)-1)/4)){
    list_of_cols = c(list_of_cols,sprintf("Channel1_R%s_IntensityMeanThrs..EGFP_R%s_IntensityMeanThrs..R",i,i))
    new_colnames = c(new_colnames,sprintf("Antenna_%s",i))
  }

  tabl = tabl[,list_of_cols]
  colnames(tabl) = new_colnames


  ### Add new column odorant to specify in which odorant we are for each time point ###
  create_odorant_col = function(num,length_exp=184){
    if(num == info.exp$numOdorant){
      theory_timing_start = (length_exp*num - length_exp) / 60
      theory_timing_stop = tabl$Timing[nrow(tabl)]# take last time point if last odor
    }else{
      theory_timing_start = (length_exp*num - length_exp) / 60
      theory_timing_stop = (length_exp*num) / 60
    }

    y = rep(as.character(info.exp$odorant[num]),nrow(tabl[tabl$Timing>=theory_timing_start & tabl$Timing<=theory_timing_stop,]))
    return(y)
  }

  tabl$odorant = NA
  tabl$odorant = unlist(lapply((1:info.exp$numOdorant),function(x){create_odorant_col(x)}))

  return(tabl)
}

####################### Format Raw data set ########################


#' normalize the dataset produced by Split_CSV
#'
#' @title format_csv
#'
#' @import dplyr
#' @import zoo
#' @import tidyr
#'
#'
#' @param csv_input table output of \code{\link{Split_CSV}}.
#'
#' @examples
#' csv_split = Split_csv(getwd())
#' df = format_csv(csv_input = csv_split)
#'
#' @export
format_csv = function(csv_input){

  ### read csv ###
  exp.info = Read_Exp_info()
  #csv = read.csv(csv)
  csv = csv_input[,-ncol(csv_input)]

  ### Normalize with the last two ROI (forehead of flies) ###
  csv[,-1] = csv[,-1] / rowMeans(csv[,c(ncol(csv)-1,ncol(csv))]) #divide all the values by the mean of last two ROI aka foreheads ROI

  ### delta F ###
  colmean = apply(csv[,-1],2,mean,na.rm=TRUE)
  f = sweep(csv[,-1],2,colmean,"-")
  DF = sweep(f,2,colmean,"/")

  melted = reshape2::melt(csv,id="Timing")
  melted = melted %>% group_by(variable) %>% mutate(ma=rollmean(as.numeric(value),25, align = "right", fill = NA))
  rollwithdf<-mutate(melted, deltaf=((as.numeric(value)-ma)/ma))
  rollwithdf2 <- dplyr::select(rollwithdf, Timing,variable,deltaf)
  rollwithdf2$deltaf[is.na(rollwithdf2$deltaf)] = 0
  new_csv = pivot_wider(rollwithdf2,names_from=variable,values_from = deltaf, names_sort=FALSE  )

  ### Add Timing column that was not included for calculate delta F ###
  csv_DF = cbind(csv$Timing,DF)
  colnames(csv_DF)[1] = "Timing"

  ### set first value as 0 to have comparables values ###
  csv[,-1] = sweep(csv[,-1],2,unlist(colMeans(csv[1:20,-1])),"-")
  csv_DF[,-1] = sweep(csv_DF[,-1],2,unlist(colMeans(csv_DF[1:20,-1])),"-")
  new_csv[,-1] = sweep(new_csv[,-1],2,unlist(colMeans(new_csv[1:20,-1])),"-")

  print(ggplot(csv_DF,aes(x=Timing,y=Antenna_8,colour="Normalized scale"))+
          geom_line()+
          geom_line(data=csv,aes(x=Timing,y=Antenna_8,colour="Original scale"))+
          geom_line(data=new_csv,aes(Timing,Antenna_8,colour="rolling mean"))+
          labs(title = exp.info$exp_id,
               subtitle = "Comparaison original scale vs Normalized scale")+
          ylab("Antenna 8 intensity")+
          scale_color_manual(name = "Legend", values = c("Normalized scale" = "blue","Original scale"="black","rolling mean"="orange")))

  csv_DF = cbind(csv_DF,"odorant" = csv_input$odorant)
  new_csv = cbind(new_csv,"odorant" = csv_input$odorant)

  return(list(normalized = csv_DF,rolling_mean=new_csv))
}

####################### Finding pics ########################

rollingSlope.lm <- function(vector) {

  a <- coef(lm(vector ~ seq(vector)))[2]
  return(a)

}

#' Find the pics in the data and superpose them for 1 odorant
#'
#' @title Finding_pics
#'
#' @import rlang
#'
#' @param csv_input rolling_mean table output of \code{\link{format_csv}}.
#' @param exp_odorant An odorant name present in the csv_input.
#'
#' @examples
#' csv_df = format_csv(csv)
#' Finding_pics(csv_df$rolling_mean, exp_odorant="water 1")
#'
#' @export
Finding_pics =function(csv_input,exp_odorant = NULL){

  ### Read exp_info ###
  exp.info = Read_Exp_info()

  ### Load csv and modify it ###
  csv = csv_input[csv_input$odorant == exp_odorant,]
  csv = csv[,-c(ncol(csv))]
  csv[,-c(1,ncol(csv))] = as.data.frame(scale(csv[,-c(1,ncol(csv))],scale = F))
  csv$Timing = rank(csv$Timing)# to get matching timings between pulses

  ### Create new csv that keep pulse information ###
  pulse_csv = data.frame()

  ### apply each conditions to each antenna ###
  for(i in 1:exp.info$condNum){
    if(i==1){
      colnames(csv)[2:(exp.info$exp_info$Antennas+1)] = paste0(colnames(csv[,-1])[1:exp.info$exp_info$Antennas],sprintf(".%s",exp.info$conditions[i]))
    }else{
      colnames(csv)[((i-1)*exp.info$exp_info$Antennas+2):(i*exp.info$exp_info$Antennas+1)] = paste0(colnames(csv)[((i-1)*exp.info$exp_info$Antennas+2):(i*exp.info$exp_info$Antennas+1)],sprintf(".%s",exp.info$conditions[i]))
    }
  }

  for(i in 2:(ncol(csv)-2)){
    name_antenna = colnames(csv)[i]

    no_response = F# to know if we have a response for this odor or not

    ### Finding the max autocorelation for the first pics ###
    ACF = acf(csv[,i],pl=F,lag = 500)
    ACF = as.data.frame(cbind(ACF$lag,ACF$acf))

    # Find the max value
    max_1 = which(ACF$V2 == max(ACF$V2[ACF$V1>380]))
    if(max_1 < 380 | max(ACF$V2[ACF$V1>380]) < 0.35){
      no_response = T
      # 409 = 184 / 60 / 3 / 0.0025 aka how much frame are in each pulse test
      max_1 = 409 #default value that seems to works properly, This is used only for no response signals so don't care much
    }

    ### Finding the max autocorelation for the second and third pics ###
    ACF = acf(csv[csv$Timing>csv$Timing[max_1],i],pl=F,lag = 500)
    ACF = as.data.frame(cbind(ACF$lag,ACF$acf))
    max_2 = which(ACF$V2 == max(ACF$V2[ACF$V1>100]))

    if(max_2 < 380 | max(ACF$V2[ACF$V1>380]) < 0.35){
      # 409 = 184 / 60 / 3 / 0.0025 aka how much frame are in each pulse test
      no_response = T
      max_2 = 409 #default value that seems to works properly, This is used only for no response signals so don't care much
    }

    First_pic = csv[csv$Timing<=csv$Timing[max_1],]
    Second_pic = csv[csv$Timing>csv$Timing[max_1] & csv$Timing<=csv$Timing[(max_1+max_2)],]
    Third_pic = csv[csv$Timing>csv$Timing[(max_1+max_2)],]

    ### Calculate rolling slope for each pulse ###
    First_pic = First_pic %>%mutate(Slope.lm = First_pic[[name_antenna]])%>% mutate(Slope.lm = rollapply(Slope.lm, width=5, FUN=rollingSlope.lm, fill=NA))
    Second_pic = Second_pic %>%mutate(Slope.lm = Second_pic[[name_antenna]])%>% mutate(Slope.lm = rollapply(Slope.lm, width=5, FUN=rollingSlope.lm, fill=NA))
    Third_pic = Third_pic %>%mutate(Slope.lm = Third_pic[[name_antenna]])%>% mutate(Slope.lm = rollapply(Slope.lm, width=5, FUN=rollingSlope.lm, fill=NA))

    ### Find the max slope and set as 0 to synchronize the pulses ###
    if(no_response){#theory should be each 409
      First_pic$Timing = First_pic$Timing - First_pic$Timing[(nrow(First_pic)/2)]
      Second_pic$Timing = Second_pic$Timing - Second_pic$Timing[(nrow(Second_pic)/2)]
      Third_pic$Timing = Third_pic$Timing - Third_pic$Timing[(nrow(Third_pic)/2)]
    }else{
      First_pic$Timing = First_pic$Timing - First_pic$Timing[which.max(First_pic$Slope.lm)]
      Second_pic$Timing = Second_pic$Timing - Second_pic$Timing[which.max(Second_pic$Slope.lm)]
      Third_pic$Timing = Third_pic$Timing - Third_pic$Timing[which.max(Third_pic$Slope.lm)]
    }


    ### Normalize the values with new 0 as reference ###
    First_0 = which(First_pic$Timing==0)
    Second_0 = which(Second_pic$Timing==0)
    Third_0 = which(Third_pic$Timing==0)
    First_pic[,i] = First_pic[,i] - mean(First_pic[(First_0-15):(First_0-2),i])
    Second_pic[,i] = Second_pic[,i] - mean(Second_pic[(Second_0-15):(Second_0-2),i])
    Third_pic[,i] = Third_pic[,i] - mean(Third_pic[(Third_0-15):(Third_0-2),i])


    ### Merge pics + create a melted DF ###
    dfs = list(First_pic[,c("Timing",name_antenna)],Second_pic[,c("Timing",name_antenna)],Third_pic[,c("Timing",name_antenna)])

    melted_CSV = Reduce(function(x,y) merge(x,y,by="Timing"),dfs)
    colnames(melted_CSV) = c("Timing",paste(name_antenna,c("pic_1","pic_2","pic_3"),sep = "/"))#use "/" to split in next function
    melted_CSV_long = reshape2::melt(melted_CSV,id.vars = "Timing")

    ### return the merged CSV + combined scaled pics csv ###
    if(is_empty(pulse_csv)){#if first antenna
      pulse_csv = melted_CSV

      output_csv = rbind(First_pic[,-which(colnames(First_pic) %in% c("Timing","Slope.lm"))],
                         Second_pic[,-which(colnames(Second_pic) %in% c("Timing","Slope.lm"))],
                         Third_pic[,-which(colnames(Third_pic) %in% c("Timing","Slope.lm"))])
    }else{
      pulse_csv = merge(pulse_csv,melted_CSV,by="Timing")

      output = rbind(First_pic[,-which(colnames(First_pic) %in% c("Timing","Slope.lm"))],
                     Second_pic[,-which(colnames(Second_pic) %in% c("Timing","Slope.lm"))],
                     Third_pic[,-which(colnames(Third_pic) %in% c("Timing","Slope.lm"))])

    }

    output_csv = rbind(First_pic[,-which(colnames(First_pic) %in% c("Timing","Slope.lm"))],
                       Second_pic[,-which(colnames(Second_pic) %in% c("Timing","Slope.lm"))],
                       Third_pic[,-which(colnames(Third_pic) %in% c("Timing","Slope.lm"))])

  }

  return(list(csv=csv,pulse_csv=pulse_csv))
}

####################### Calculating the z_score #######################

#' z score calculation for each odorant
#'
#' @title Z_score_calculation
#'
#' @param pulse_csv pulse_csv table output of \code{\link{Finding_pics}}
#'
#' @examples
#' pics = Finding_pics(csv$rolling_mean, exp_odorant="water 1")
#' z_score_calculation(pics$pulse_csv)
#'
#' @export
Z_score_calculation = function(pulse_csv){
  ### Melt the dataset ###
  melted_csv= reshape2::melt(pulse_csv,id.vars = "Timing")

  ### Formatting csv adding pulse information ###
  pulse = unique(sapply(strsplit(as.character(melted_csv$variable),"/"), `[`, 2))
  melted_csv$pulse = unlist(pulse[match(sapply(strsplit(as.character(melted_csv$variable),"/"), `[`, 2),pulse)])
  melted_csv$pulse = as.factor(melted_csv$pulse)

  ### Formatting csv adding conditions information ###
  melted_csv$variable = sapply(strsplit(as.character(melted_csv$variable),"/"), `[`, 1)
  conditions = unique(sapply(strsplit(as.character(melted_csv$variable),"[.]"), `[`, 2))
  melted_csv$conditions = unlist(conditions[match(sapply(strsplit(as.character(melted_csv$variable),"[.]"), `[`, 2),conditions)])
  melted_csv$conditions = factor(melted_csv$conditions,levels=unique(melted_csv$conditions))

  ### Taking rid of long names to keep antenna id ###
  melted_csv$variable = sapply(strsplit(as.character(melted_csv$variable),"[.]"), `[`, 1)
  melted_csv$variable = factor(melted_csv$variable,levels = unique(melted_csv$variable))
  melted_csv$value = melted_csv$value+1

  ### calculation of standard deviation + mean of each antenna + condition ###
  z_sd = melted_csv[melted_csv$Timing <0 & melted_csv$Timing >-150,] %>% group_by(variable,conditions) %>% summarize(sd = sd(value))
  z_mean = melted_csv[melted_csv$Timing <0 & melted_csv$Timing >-150,] %>% group_by(variable,conditions) %>% summarize(mean = mean(value))
  # z_sd = melted_csv %>% group_by(variable,pulse,conditions) %>% summarize(sd = sd(value))
  # z_mean = melted_csv %>% group_by(variable,pulse,conditions) %>% summarize(mean = mean(value))

  ### compute z_score with usual formula ###
  melted_csv = left_join(melted_csv,z_sd,by = c("variable","conditions")) %>%
    left_join(z_mean,by = c("variable","conditions")) %>%
    mutate(z_score = (value-mean)/sd)

  ### add max value for normal + z_score ###
  melted_csv = melted_csv %>% group_by(variable,pulse,conditions) %>% mutate(max_value = max(value),max_zscore = max(z_score))
  melted_csv = melted_csv %>% group_by(variable,pulse,conditions) %>% mutate(mean_value = mean(value[Timing>-10&Timing<50]),mean_z_score = mean(z_score[Timing>-10&Timing<50]))

  return(melted_csv)
}

####################### Creating the plots #######################

### 1 First plots (plot_trace) for 1 odor ###

#' plot the trace of fluorescence intensity raw or computed
#'
#' @title PlotTrace
#'
#' @param csv table output of \code{\link{Z_score_calculation}}.
#' @param combined Default = True, Plot the trace uncut (False) or cut at each pulse and overlapped (True).
#' @param z_score Default = False, plot the trace on a raw scale (False) or z score scale (True).
#'
#' @examples
#' pics = Finding_pics(csv$rolling_mean, exp_odorant="water 1")
#' pics$pulse_csv = z_score_calculation(pics$pulse_csv)
#' PlotTrace(pics,combined=F,z_score=T)
#'
#' @export
PlotTrace = function(csv,combined = T,z_score = F) {

  if(combined){
    ### load correct data ###
    melted_csv = csv$pulse_csv

    ### create xtitle variable aka conditions ###
    xtitle = unique(melted_csv$conditions)

    ### create ggplot ###
    if(z_score){
      p1 = ggplot(melted_csv,aes(Timing,z_score,col=pulse))
    }else{
      p1 = ggplot(melted_csv,aes(Timing,value,col=pulse))
    }


  }else{
    ### Take off forehead ROIs ###
    csv$csv = csv$csv[,-c(ncol(csv$csv)-1,ncol(csv$csv))]

    ### Melt the dataset ###
    melted_csv = reshape2::melt(csv$csv,id.vars = "Timing")

    ### Taking rid of long names to keep antenna id ###
    xtitle = unique(sapply(strsplit(as.character(melted_csv$variable),"[.]"), `[`, 2))
    melted_csv$variable = sapply(strsplit(as.character(melted_csv$variable),"[.]"), `[`, 1)
    melted_csv$variable = factor(melted_csv$variable,levels = unique(melted_csv$variable))

    ### create ggplot ###
    p1 = ggplot(melted_csv,aes(Timing,value))
  }



  ### finish ggplot ###
  for(p in 1:length(xtitle)){#number of pages
    print(p1+
            geom_line()+
            xlab(xtitle[p])+
            theme(axis.title.y = element_blank(),
                  axis.text.y = element_text(size=13),
                  axis.text.x = element_text(size=12),
                  #axis.ticks.y = element_blank(),
                  strip.text.y = element_text(angle = 0),
                  axis.title.x = element_text(size = 15, face="bold" ))+
            facet_wrap_paginate(~ variable,ncol = 2,nrow = 3,page = p)
    )
  }
}

####################### Creating the plots #######################
### 2 combined all antennas of same conditions ###

#' plot the trace of fluorescence intensity overlapped for each conditions
#'
#' @title PlotTraceCondition
#'
#' @param csv table output of \code{\link{Z_score_calculation}}.
#' @param groupBy Default = "pulse", color by pulse the overlapped traces ("pulse") or by antenna ID ("antenna").
#' @param z_score Default = True, plot the trace on a raw scale (False) or z score scale (True).
#'
#' @examples
#' pics = Finding_pics(csv$rolling_mean, exp_odorant="water 1")
#' pics$pulse_csv = z_score_calculation(pics$pulse_csv)
#' PlotTraceCondition(pics,groupBy="pulse",z_score=T)
#'
#' @export
PlotTraceCondition = function(csv,groupBy = "pulse",z_score = T){
  ### load dataset ###
  melted_csv = csv$pulse_csv

  ### create xtitle variable aka conditions ###
  xtitle = unique(melted_csv$conditions)

  if(groupBy == "pulse"){
    if(z_score){
      p1 = ggplot(melted_csv,aes(Timing,z_score,group=interaction(pulse,variable),col=pulse))
      label = "z score"
    }else{
      p1 = ggplot(melted_csv,aes(Timing,value,group=interaction(pulse,variable),col=pulse))
      label = "Normalized score"
    }
    print(p1+
            geom_line(size=0.3)+
            labs(ylab = label)+
            theme(strip.text = element_text(size=5))+
            facet_wrap(~ conditions)

    )
  }else if(groupBy == "antenna"){
    ### create a list of ggplot for each conditions ###
    if(z_score){
      ggList <- lapply(split(melted_csv, melted_csv$conditions), function(i) {
        ggplot(i, aes(Timing,z_score,group=interaction(pulse,variable),colour=variable)) +
          geom_line(size=0.3)+
          labs(ylab = "z_score",title = unique(i$conditions),colour = "Antenna ID")+
          theme(plot.title = element_text(size =5,face = "bold"))
      })
    }else{
      ggList <- lapply(split(melted_csv, melted_csv$conditions), function(i) {
        ggplot(i, aes(Timing,value,group=interaction(pulse,variable),colour=variable)) +
          geom_line(size=0.3)+
          labs(ylab = "normalized score",title = unique(i$conditions),colour = "Antenna ID")+
          theme(plot.title = element_text(size =5,face = "bold"))
      })
    }

    ### grid them together ###
    print(cowplot::plot_grid(plotlist = ggList,align="h",ncol = 2))
  }
}

####################### Creating the plots #######################

### 3 resume all responses ###

#' Resume plot with boxplots representing max values of each pulse for each antenna or condition
#'
#' @title PlotResume
#'
#' @param csv table output of \code{\link{Z_score_calculation}}.
#' @param groupBy Default = "antenna", separate each antenna wrapped by conditions ("antenna"), If another string than antenna is given the global mean of each condition will be displayed.
#' @param z_score Default = True, plot the trace on a raw scale (False) or z score scale (True).
#'
#' @examples
#' pics = Finding_pics(csv$rolling_mean, exp_odorant="water 1")
#' pics$pulse_csv = z_score_calculation(pics$pulse_csv)
#' PlotResume(pics,groupBy="antenna",z_score=T)
#'
#' @export
PlotResume = function(csv,groupby = "antenna",z_score = T){

  ### load dataset ###
  melted_csv = csv$pulse_csv

  ### compute the maximum by each category ###
  if(z_score){
    max_value_csv = melted_csv%>% group_by(variable,pulse,conditions) %>% summarize(max = max(z_score))
    ylabel = "Z score max"
  }else{
    max_value_csv = melted_csv%>% group_by(variable,pulse,conditions) %>% summarize(max = max(value))
    ylabel = "normalized max"
  }

  ### plot the results based on users choices ###
  if(groupby=="antenna"){
    p1 = ggplot(max_value_csv,aes(variable,max))+
      geom_boxplot(coef=0)+
      geom_jitter(aes(colour=pulse),shape=16, position=position_jitter(0))+
      labs(ylab=ylabel)+
      theme(axis.text.x = element_text(angle = 45,hjust = 1),
            axis.title.x = element_blank())+
      facet_grid(~ conditions,scales="free_x")
    print(p1)

  }else{
    p1 = ggplot(max_value_csv,aes(conditions,max))+
      geom_boxplot()+
      geom_jitter(aes(colour=pulse),shape=16, position=position_jitter(0))+
      labs(ylab=ylabel)+
      theme(axis.text.x = element_text(angle = 45,hjust = 1),
            axis.title.x = element_blank())
    print(p1)

  }

  return(p1)
}

####################### Generating complete pipeline #######################

#' Pipeline to generate a single PDF containing all plots doable with this package for 1 odor.
#'
#' @title create_pdf
#'
#' @param csv_DF table output of \code{\link{format_csv}}.
#' @param odorant Corresponding odorant to generate the plots (Should be present in the exp_info.csv file).
#' @param z_score Default = True, plot the trace on a raw scale (False) or z score scale (True).
#' @param name_file File name of the PDF.
#'
#' @examples
#' exp.info = Read_Exp_info()
#' csv_split = Split_CSV(getwd())
#' csv_DF = format_csv(csv=csv_split)
#' create_pdf(csv_DF,exp.info$odorant[1],z_score = T,name_file="Odorant_1.pdf")
#'
#' @export
create_pdf = function(csv_DF,odorant,z_score = F,name_file){

  dataset = Finding_pics(csv_input = csv_DF$rolling_mean,exp_odorant = as.character(odorant))
  dataset_normalized = Finding_pics(csv_input = csv_DF$normalized,exp_odorant = as.character(odorant))
  dataset$pulse_csv = Z_score_calculation(dataset$pulse_csv)

  ### Create PDF ###
  pdf(name_file)
  PlotTrace(csv = dataset_normalized,combined = F,z_score = z_score)
  PlotTrace(csv = dataset,combined = T,z_score = z_score)
  PlotTraceCondition(csv=dataset,groupBy = "pulse",z_score = z_score)
  PlotTraceCondition(csv=dataset,groupBy = "antenna",z_score = z_score)
  res_antenna = PlotResume(csv=dataset,groupby = "antenna",z_score = z_score)
  res_cond = PlotResume(csv=dataset,groupby = "condition",z_score = z_score)
  dev.off()
  name_antenna_resume = paste0("Antenna_",name_file)
  name_condition_resume = paste0("Condition_",name_file)

  return(df = dataset$pulse_csv)
}


#' Pipeline to generate the complete analysis
#'
#' @title pipeline
#'
#' @param z_score Default = False, plot the trace on a raw scale (False) or z score scale (True).
#'
#' @examples
#' pipeline(z_score=T)
#'
#' @export
pipeline = function(z_score = F){

  ### load all datasets needed ###
  exp.info = Read_Exp_info()
  csv = Split_CSV(getwd())
  csv_DF = format_csv(csv=csv)

  ### create file names ###
  nametags <- LETTERS
  for(i in 1:25){
    nametags <- append(nametags,LETTERS)
  }
  nametags2 <- nametags
  nametags <- sort(nametags)
  for (i in 1:length(nametags)){
    nametags[i] <- paste(nametags[i], nametags2[i], sep="")
  }

  ### Create the directory of the plots if not exist ###
  if(!dir.exists("Plots")){
    dir.create("Plots")
  }

  ### compute graph for each odorant ###
  for_resume_final = c()
  for(i in 1:exp.info$numOdorant){
    ### Create PDF ###
    name_file = paste0("Plots/",paste(nametags[i],exp.info$odorant[i],exp.info$exp_id,sep = "_"),".pdf")
    print(name_file)
    for_resume = create_pdf(csv_DF,exp.info$odorant[i],z_score = z_score,name_file=name_file)

    ### for resume formatting + saved ###
    for_resume = for_resume %>% select(variable,conditions,pulse,max_zscore,mean_z_score,max_value,mean_value)
    for_resume$odorant = exp.info$odorant[i]
    for_resume_final = rbind(for_resume_final,for_resume)

  }

  ### Create short conditions names ###
  cond_names =lapply(strsplit(as.character(levels(for_resume_final$conditions)),"[_]"), `[`, c(1,7,8))
  short_cond_names = sapply(cond_names,function(x) paste(unlist(x, use.names = TRUE), collapse = "_"))
  #levels(for_resume_final$conditions) = short_cond_names
  for_resume_final$short_cond_names = for_resume_final$conditions
  levels(for_resume_final$short_cond_names) = short_cond_names

  ### Convert odorant to factor + change antenna names by their number only ###
  for_resume_final$odorant = factor(for_resume_final$odorant,levels = unique(for_resume_final$odorant))
  levels(for_resume_final$variable) = c(1:length(levels(for_resume_final$variable)))

  ### eliminates duplicates ###
  for_resume_final = unique(for_resume_final)

  ### generate samescale plot ###
  pdf("Plots/samescale_summary.pdf",width = 18)
  print(ggplot(for_resume_final,aes(short_cond_names,max_zscore,fill=short_cond_names))+
          geom_boxplot()+
          theme(axis.text.x = element_text(angle=45,hjust=1),
                axis.title.x = element_blank(),
                legend.title = element_blank(),
                legend.position="bottom",
                legend.direction="horizontal")+
          facet_wrap(~ odorant,scales="free_x",nrow = 1)
  )

  plot_list = c()
  max_value = max(for_resume_final$max_zscore)
  for(i in 1:length(levels(for_resume_final$short_cond_names))){

    ### create new df and add score column aka if double the mean of water 1 value = green color ###
    df =for_resume_final[for_resume_final$short_cond_names == levels(for_resume_final$short_cond_names)[i],]
    color = df %>% group_by(odorant,variable,short_cond_names) %>% mutate(Mean = max(max_zscore))
    water_mean = mean(color$Mean[color$odorant=="water 1"])
    color = color %>% mutate(score = as.factor(ifelse(Mean > water_mean*2,1,0)))

    if(i == 1){
      p = ggplot(color,aes(variable,max_zscore,fill=score))+
        geom_boxplot()+
        ylab(levels(for_resume_final$short_cond_names)[i])+
        scale_y_continuous(limits = c(0,max_value))+
        theme(axis.title.x = element_blank(),
              strip.text = element_text(size=7),
              legend.position = "none")+
        facet_grid(~odorant,scales="free")+
        scale_fill_manual(values = c("white","green"))

    }else{
      p = ggplot(color,aes(variable,max_zscore,fill=score))+
        geom_boxplot()+
        scale_y_continuous(limits = c(0,max_value))+
        ylab(levels(for_resume_final$short_cond_names)[i])+
        theme(axis.title.x = element_blank(),
              strip.text = element_text(size=5),
              strip.background = element_blank(),
              strip.text.x = element_blank(),
              legend.position = "none")+
        facet_grid(~odorant,scales="free")+
        scale_fill_manual(values = c("white","green"))
    }

    plot_list = c(plot_list,list(p))
  }

  ### assemble together each condition plot
  print(cowplot::plot_grid(plotlist = plot_list, nrow =length(levels(for_resume_final$short_cond_names)) ))

  ### create diff from endogenous plot
  endo_cond_name = levels(for_resume_final$short_cond_names)[grep("endo",levels(for_resume_final$short_cond_names),ignore.case = T)]
  if(length(endo_cond_name)==0){

  }else{
    endo_normalized = for_resume_final %>% group_by(odorant) %>% mutate(fold_change = (max_zscore- max_zscore[short_cond_names == endo_cond_name])/max_zscore[short_cond_names == endo_cond_name])

    stat.test = endo_normalized %>%
      group_by(odorant) %>%
      t_test(fold_change ~ short_cond_names, ref.group = endo_cond_name) %>%
      adjust_pvalue(method = "bonferroni") %>%
      add_significance()%>%
      add_xy_position(x="short_cond_names")

    print(ggplot(endo_normalized,aes(short_cond_names,fold_change,color=short_cond_names))+
            geom_boxplot()+
            theme(axis.text.x = element_text(angle=45,hjust=1),
                  axis.title.x = element_blank(),
                  legend.title = element_blank(),
                  legend.position="bottom",
                  legend.direction="horizontal")+
            facet_wrap(~ odorant,scales="free_x",nrow = 1)+
            stat_pvalue_manual(stat.test,hide.ns = T)
    )
  }

  dev.off()
  return(for_resume_final)
}
