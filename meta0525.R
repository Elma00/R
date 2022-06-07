##input data
library(forestplot)
library(readxl)
library(meta)
library(metafor)
library(dplyr)
library(ggplot2)
library(DataCombine)

#update the meta=data
Meta_data <- read_excel("meta_0426.xlsx")
ID_note <- read_excel("ID_note.xlsx")

##############################recalculate data : just correct the direction##############################
#dput(names(Meta_data))
Data_all <- Meta_data
Data_all <- Data_all %>% filter(is.na(Delete)==TRUE,use==1)
ID <- unique(Data_all[,c("ID", "Study")])  #check 96
#Data_all <- Data_all %>% filter(is.na(Population_health)==TRUE)

Data <- Data_all[,1:57]
#Data$disease <- ifelse(is.na(Data$Population_health)==FALSE,1,0)
table(Data$disease)
source("transform.R")
source("plot_meta.R")

forest_data <- Data
sexHormone <- c("TT","BT","SHBG","E2","DHEAS")
sexHormone_name <- c("TT","Bio-T","SHBG","E2","DHEAS")
Disease <- c("CVD","CHD","Stroke","HF")
sex_name <- c("Men","Women")


##over lap
forest_data <- subset(forest_data,!(reference %in% c("m3","m4","m5")))  
##with CHD
forest_data <- subset(forest_data,ID!="74")  #
forest_data <- subset(forest_data,ID!="121")  #
forest_data <- subset(forest_data,ID!="235")  #
forest_data <- subset(forest_data,ID!="246")  #
forest_data <- subset(forest_data,ID!="334")  #
forest_data <- subset(forest_data,ID!="202102")  #


#overlap 
forest_data <- subset(forest_data,ID!="15")  #overlap IID 13,14
forest_data <- subset(forest_data,ID!="139")  #overlap ID138 and with short follow up
forest_data <- subset(forest_data,ID!="140")  #overlap ID138 and with short follow up
forest_data <- subset(forest_data,ID!="355")  #overlap ID138 and with short follow up

forest_data <- subset(forest_data,ID!="147") # Delete due to overlap 148 and without related exposure
##low quality
forest_data <- subset(forest_data,ID!="13")  #low quality
forest_data <- subset(forest_data,ID!="14")  #low quality for DHEAS
forest_data <- subset(forest_data,ID!="15")  #low quality for DHEAS
forest_data <- subset(forest_data,ID!="79")  #low quality ??
forest_data <- subset(forest_data,ID!="179")  #low quality ??

forest_data <- subset(forest_data,ID!="177")  #low quality ??
##
forest_data <- subset(forest_data,ID!="210905")  #UKB double202201
forest_data <- subset(forest_data,ID!="210902")  #overlap ID125
forest_data <- subset(forest_data,ID!="107")  #overlap ID125
forest_data <- subset(forest_data,ID!="310")  #M j shape, Q5 vs Q2
forest_data <- subset(forest_data,ID!="213")  #overlap ID 325 CHD

forest_data <- subset(forest_data,id!="14602")  #M j shape, Q5 vs Q2
forest_data <- subset(forest_data,ID!="208")  #case-conttrol
forest_data <- subset(forest_data,ID!="162")  #referebce hwu Q2-Q4 vs Q1
forest_data <- subset(forest_data,ID!="156") #TT, FT overlap ID354-356, delete due to 354 cover more patient


###all  new
Result_CHD <- data.frame()
for(x in 2){
  for(n in 1:2){
    for(i in 1:length(sexHormone)){
      data_meta <- subset(forest_data,Exposure==sexHormone[i] & sex==sex_name[n] & Outcome==Disease[x])
      #data_meta <- subset(data_meta,ID!="156")    #TT, FT overlap ID354-356, delete due to 354 cover more patient
      ID_note_meta <- subset(ID_note,sex==sex_name[n])
      Data_meta <- merge(data_meta,ID_note_meta,by=c("ID","id","sex","type"),all.x=TRUE)
      disease <- Disease[x]
      sexH <- sexHormone[i]
      gender <- sex_name[n]
      if(dim(data_meta)[1]>0){
        
        m <- metagen(TE=lnhr, seTE, 
                     studlab=paste(Author,"et al, ",year), 
                     n.e=Events, n.c=N, 
                     comb.fixed=T, comb.random=T,
                     label.e=gs('Sample size'), label.c=gs('Cases'), 
                     sm="RR", method.tau = "DL", 
                     hakn=F, prediction = F, backtransf=T, 
                     #byvar=`BMI `, 
                     #byvar=bmi_bigger27, #TT
                     overall=T,
                     data=Data_meta)
        filename_plot <-  paste0("plot/",sex_name[n],"_",Disease[x],"&",sexHormone[i],".pdf")
        pdf(file=filename_plot,height=10,width=10)
        if(n==1){
          if(i==1){meta_plot_g(m,0.3,2,10,10)}
          if(i==2){meta_plot_g(m,0.5,2.5,10,10)}
          if(i==3){meta_plot_g(m,0.5,1.5,10,10)}
          if(i==4){meta_plot_g(m,0.5,1.5,10,10)}
          if(i==5){meta_plot_g(m,0.4,1.25,5,25)}
        }
        if(n==2){
          if(i==1){meta_plot_women(m,0.4,2.5,10,10)}
          if(i==2){meta_plot_women(m,0.2,2,5,15)}
          if(i==3){meta_plot_women(m,0.4,2.5,10,10)}
          if(i==4){meta_plot_women(m,0.4,2.5,10,10)}
          if(i==5){meta_plot_women(m,0.2,2,10,20)}
        }
        dev.off()
        bias_p <- round(metabias(m,k.min=dim(Data_meta)[1])$p.value,3)
        if(bias_p<0.05){
          a <- trimfill(m)
          trim_hr <- round(exp(a$TE.random),2)
          trim_lhr <- round(exp(a$lower.random),2)
          trim_uhr <- round(exp(a$upper.random),2)
          trim <- paste0(trim_hr," (",trim_lhr,",",trim_uhr,")")
        }
        if(bias_p>=0.05){trim <- NA}  
        m <- metagen(TE=lnhr, seTE, 
                     studlab=paste0(Author,", ",year), 
                     n.e=Events, n.c=N, 
                     comb.fixed=T, comb.random=T,
                     label.e=gs('Sample size'), label.c=gs('Cases'), 
                     sm="RR", method.tau = "DL", 
                     hakn=F, prediction = F, backtransf=T, 
                     #byvar=`BMI `, 
                     #byvar=bmi_bigger27, #TT
                     overall=T,
                     data=Data_meta)
        filename_plot <-  paste0("plot/Funnel_",sexvar[n],"_",sexHormone[i], "&", Disease[x],".tiff")
        tiff(file=filename_plot,width=500,height=500)
        funnel(m,studlab = TRUE,label=TRUE,
               #level = 0.95, contour = c(0.9, 0.95, 0.99),
               xlab=paste0("RR for association of ",sexHormone[i]," with ", Disease[x])) 
        dev.off()
        
        if(dim(data_meta)[1]>2){
          
          filename_plot <-  paste0("plot/",sex_name[n],"_",sexHormone[i],"&",Disease[x],"_cum.pdf")
          pdf(file=filename_plot)
          forest(metacum(m,sortvar = year),
                 leftcols =c('studlab'),
                 leftlabs=c('Adding Source (k = No. of Studies)'),
                 rightcols = c( "effect", "ci"),
                 rightlabs = c( "RR", "[95%CI]"),
                 #xlab= paste0("RR (95% CI) for ",sexHormone[i],"  and CHD in ",sex_name[n]),
                 ff.random='bold',
                 ff.lr="bold",
                 ff.xlab = 'bold')
          dev.off()
          ##leave one out
          leave_plot <- paste0("plot/",sex_name[n],"_",sexHormone[i],"&",Disease[x], "_Leaveone.pdf")
          pdf(file=leave_plot)
          forest(metainf(m,sortvar = year),
                 leftcols =c('studlab'),
                 leftlabs=c('Omitting Source'),
                 rightcols = c( "effect", "ci"),
                 rightlabs = c( "RR", "[95%CI]"),
                 #xlab= paste0("RR (95% CI) for ",sexHormone[i],"  and CHD in ",sex_name[n]),
                 ff.random='bold',
                 ff.lr="bold",
                 ff.xlab = 'bold',
          )
          dev.off()
        }
        #meta result
        sm <- summary(m)
        No <- sm$k
        num <- unique(Data_meta[,c("sex","Events","N")] %>%
                        group_by(sex) %>%
                        mutate(all_cases=sum(Events),
                               all_sample=sum(N)) %>% 
                        select(sex,all_cases,all_sample))
        num$total <- paste0(num[2],"/",num[3])
        #random effect
        hr_r <- sprintf("%0.2f",exp(sm$random$TE))
        HR_l_r <- sprintf("%0.2f",exp(sm$random$lower))
        HR_u_r <- sprintf("%0.2f",exp(sm$random$upper))
        HR_r <- paste0(hr_r," (",HR_l_r,", ",HR_u_r,")")
        #fix effect
        hr_f <- sprintf("%0.2f",exp(sm$fixed$TE))
        HR_l_f <- sprintf("%0.2f",exp(sm$fixed$lower))
        HR_u_f <- sprintf("%0.2f",exp(sm$fixed$upper))
        HR_f <- paste0(hr_f," (",HR_l_f,", ",HR_u_f,")")
        I2 <- paste0(sprintf("%0.1f",m$I2*100),"%")
        Q <- round(sm$Q,0)
        
        result <- data.frame(gender,sexH,disease,No,num[4],HR_r,HR_f,Q,I2,bias_p,trim)
        Result_CHD <- rbind(Result_CHD,result)
      }
    }
  }
}
Result_CHD$sexH <- factor(Result_CHD$sexH,levels=c("TT","BT","E2","DHEAS","SHBG"))
Result_CHD <- Result_CHD %>% arrange(gender,sexH)
write.csv(Result_CHD,"result/Result_CHD.csv")

#exlcuding UKB
###all  new
Result_CHD <- data.frame()
for(x in 2){
  for(n in 1:2){
    for(i in 1:length(sexHormone)){
      data_meta <- subset(forest_data,Exposure==sexHormone[i] & sex==sex_name[n] & Outcome==Disease[x])
      data_meta <- subset(data_meta,ID!="202113")    
      data_meta <- subset(data_meta,ID!="202201")    
      data_meta <- subset(data_meta,ID!="202211")    
      ID_note_meta <- subset(ID_note,sex==sex_name[n])
      Data_meta <- merge(data_meta,ID_note_meta,by=c("ID","id","sex","type"),all.x=TRUE)
      disease <- Disease[x]
      sexH <- sexHormone[i]
      gender <- sex_name[n]
      if(dim(data_meta)[1]>0){
        
        m <- metagen(TE=lnhr, seTE, 
                     studlab=paste(Author,"et al, ",year), 
                     n.e=Events, n.c=N, 
                     comb.fixed=T, comb.random=T,
                     label.e=gs('Sample size'), label.c=gs('Cases'), 
                     sm="RR", method.tau = "DL", 
                     hakn=F, prediction = F, backtransf=T, 
                     #byvar=`BMI `, 
                     #byvar=bmi_bigger27, #TT
                     overall=T,
                     data=Data_meta)
        filename_plot <-  paste0("plot/noUKB_",sex_name[n],"_",Disease[x],"&",sexHormone[i],".pdf")
        pdf(file=filename_plot,height=10,width=10)
        if(n==1){
          if(i==1){meta_plot_g(m,0.3,2,10,10)}
          if(i==2){meta_plot_g(m,0.5,2.5,10,10)}
          if(i==3){meta_plot_g(m,0.5,1.5,10,10)}
          if(i==4){meta_plot_g(m,0.5,1.5,10,10)}
          if(i==5){meta_plot_g(m,0.4,1.25,5,25)}
        }
        if(n==2){
          if(i==1){meta_plot_women(m,0.4,2.5,10,10)}
          if(i==2){meta_plot_women(m,0.2,2,5,15)}
          if(i==3){meta_plot_women(m,0.4,2.5,10,10)}
          if(i==4){meta_plot_women(m,0.4,2.5,10,10)}
          if(i==5){meta_plot_women(m,0.2,2,10,20)}
        }
        dev.off()
        bias_p <- round(metabias(m,k.min=dim(Data_meta)[1])$p.value,3)
        if(bias_p<0.05){
          a <- trimfill(m)
          trim_hr <- round(exp(a$TE.random),2)
          trim_lhr <- round(exp(a$lower.random),2)
          trim_uhr <- round(exp(a$upper.random),2)
          trim <- paste0(trim_hr," (",trim_lhr,",",trim_uhr,")")
        }
        if(bias_p>=0.05){trim <- NA}  
        m <- metagen(TE=lnhr, seTE, 
                     studlab=paste0(Author,", ",year), 
                     n.e=Events, n.c=N, 
                     comb.fixed=T, comb.random=T,
                     label.e=gs('Sample size'), label.c=gs('Cases'), 
                     sm="RR", method.tau = "DL", 
                     hakn=F, prediction = F, backtransf=T, 
                     #byvar=`BMI `, 
                     #byvar=bmi_bigger27, #TT
                     overall=T,
                     data=Data_meta)
        filename_plot <-  paste0("plot/noUKB_Funnel_",sexvar[n],"_",sexHormone[i], "&", Disease[x],".tiff")
        tiff(file=filename_plot,width=500,height=500)
        funnel(m,studlab = TRUE,label=TRUE,
               #level = 0.95, contour = c(0.9, 0.95, 0.99),
               xlab=paste0("RR for association of ",sexHormone[i]," with ", Disease[x])) 
        dev.off()
        
        if(dim(data_meta)[1]>2){
          
          filename_plot <-  paste0("plot/noUKB_",sex_name[n],"_",sexHormone[i],"&",Disease[x],"_cum.pdf")
          pdf(file=filename_plot)
          forest(metacum(m,sortvar = year),
                 leftcols =c('studlab'),
                 leftlabs=c('Adding Source (k = No. of Studies)'),
                 rightcols = c( "effect", "ci"),
                 rightlabs = c( "RR", "[95%CI]"),
                 #xlab= paste0("RR (95% CI) for ",sexHormone[i],"  and CHD in ",sex_name[n]),
                 ff.random='bold',
                 ff.lr="bold",
                 ff.xlab = 'bold')
          dev.off()
          ##leave one out
          leave_plot <- paste0("plot/noUKB_",sex_name[n],"_",sexHormone[i],"&",Disease[x], "_Leaveone.pdf")
          pdf(file=leave_plot)
          forest(metainf(m,sortvar = year),
                 leftcols =c('studlab'),
                 leftlabs=c('Omitting Source'),
                 rightcols = c( "effect", "ci"),
                 rightlabs = c( "RR", "[95%CI]"),
                 #xlab= paste0("RR (95% CI) for ",sexHormone[i],"  and CHD in ",sex_name[n]),
                 ff.random='bold',
                 ff.lr="bold",
                 ff.xlab = 'bold',
          )
          dev.off()
        }
        #meta result
        sm <- summary(m)
        No <- sm$k
        num <- unique(Data_meta[,c("sex","Events","N")] %>%
                        group_by(sex) %>%
                        mutate(all_cases=sum(Events),
                               all_sample=sum(N)) %>% 
                        select(sex,all_cases,all_sample))
        num$total <- paste0(num[2],"/",num[3])
        #random effect
        hr_r <- sprintf("%0.2f",exp(sm$random$TE))
        HR_l_r <- sprintf("%0.2f",exp(sm$random$lower))
        HR_u_r <- sprintf("%0.2f",exp(sm$random$upper))
        HR_r <- paste0(hr_r," (",HR_l_r,", ",HR_u_r,")")
        #fix effect
        hr_f <- sprintf("%0.2f",exp(sm$fixed$TE))
        HR_l_f <- sprintf("%0.2f",exp(sm$fixed$lower))
        HR_u_f <- sprintf("%0.2f",exp(sm$fixed$upper))
        HR_f <- paste0(hr_f," (",HR_l_f,", ",HR_u_f,")")
        I2 <- paste0(sprintf("%0.1f",m$I2*100),"%")
        Q <- round(sm$Q,0)
        
        result <- data.frame(gender,sexH,disease,No,num[4],HR_r,HR_f,Q,I2,bias_p,trim)
        Result_CHD <- rbind(Result_CHD,result)
      }
    }
  }
}
Result_CHD$sexH <- factor(Result_CHD$sexH,levels=c("TT","BT","E2","DHEAS","SHBG"))
Result_CHD <- Result_CHD %>% arrange(gender,sexH)
write.csv(Result_CHD,"result/Result_CHD_noUKB.csv")

##meta-regression
Result_r <- data.frame()
for(x in 2){
  for(n in 1:2){
    for(i in 1:length(sexHormone)){
      data_meta <- subset(forest_data,Exposure==sexHormone[i] & sex==sex_name[n] & Outcome==Disease[x])
      data_meta <- subset(data_meta,ID!="156") 
      ID_note_meta <- subset(ID_note,sex==sex_name[n])
      Data_meta <- merge(data_meta,ID_note_meta,by=c("ID","id","sex","type"),all.x=TRUE)
      disease <- Disease[x]
      sexH <- sexHormone[i]
      gender <- sex_name[n]
      Data_meta$bmi_mean[is.na(Data_meta$bmi_mean)==TRUE] <- mean(Data_meta$bmi_mean,na.rm=T)
      Data_meta$age_mean[is.na(Data_meta$age_mean)==TRUE] <- mean(Data_meta$age_mean,na.rm=T)
      Data_meta$age_bigger60 <- ifelse(Data_meta$age_mean<=60,0,1)
      Data_meta$bmi_bigger27 <- ifelse(Data_meta$bmi_mean<27,0,1)
      Data_meta$duration <- as.numeric(Data_meta$duration)
      Data_meta$duration_bigger5 <- ifelse(Data_meta$duration<7,0,1)
      Data_meta$N_biggre1000 <- ifelse(Data_meta$N<1000,0,1) #TT
      var <- names(Data_meta)[130:dim(Data_meta)[2]]
      if(dim(data_meta)[1]>0){
        m <- metagen(TE=lnhr, seTE, 
                     studlab=paste0(Author,", ",year), 
                     n.e=Events, n.c=N, 
                     comb.fixed=T, comb.random=T,
                     label.e=gs('Sample size'), label.c=gs('Cases'), 
                     sm="RR", method.tau = "DL", 
                     hakn=F, prediction = F, backtransf=T, 
                     overall=T,
                     data=Data_meta)
        I2 <- paste0(sprintf("%0.1f",m$I2*100),"%")
        M_reg[[1]] <- metareg(m, ~  age_mean, method.tau = "DL")
        M_reg[[2]] <- metareg(m, ~  bmi_mean, method.tau = "DL")
        M_reg[[3]] <- metareg(m, ~  duration, method.tau = "DL")
        M_reg[[4]] <- metareg(m, ~  N, method.tau = "DL")
        for(b in 1:length(var)){
          m1 <- M_reg[[b]]
          m1_p <- m1$pval[2]
          m1_b <- paste0(sprintf("%0.2f",m1$b[2])," (",sprintf("%0.2f",m1$ci.lb[2]),", ",sprintf("%0.2f",m1$ci.ub[2]),")")
          
          HR_adjust <- paste0(sprintf("%0.2f",exp(m1$beta[1]))," (",
                              sprintf("%0.2f",exp(m1$ci.lb[1])),", ",
                              sprintf("%0.2f",exp(m1$ci.u[1])),")")
          I2_adjust <- paste0(sprintf("%0.1f",m1$I2),"%")
          bubble_plot <- paste0("plot/","bubble_",var[b],"_",sex_name[n],"_",sexHormone[i],"&",Disease[x], ".pdf")
          pdf(file=bubble_plot)
          bubble(m1, studlab = TRUE,ylab=paste0("RR of CHD and ",sexHormone[i]))
          dev.off()
          m1_r <- c(gender,sexH,var[b],HR_adjust,I2,I2_adjust,m1_b,sprintf("%0.2f",m1_p))
          Result_r <- rbind(Result_r,m1_r)
        }
        
      }
    }
  }
}

###subgroup
Data <- merge(forest_data,ID_note,by=c("ID","id","sex","type"),all.x=TRUE)
Data <- Data %>% filter(Outcome=="CHD")
Data <- Data %>% filter(Exposure %in% c("TT","BT","E2","DHEAS","SHBG"))
#dput(table(Data$Country))
table(Data$Country)
Data$region[Data$Country %in% c("China","Korea","Japan","Iran","Australia","India")] <- "Asia-Pacific"
Data$region[Data$Country %in% c("Denmark","European areas","Poland","Spain","Finland","Italy","Netherlands","Norway",
                                "France","Germany","Greece", "Sweden","Switzerland","Turkey","UK")] <- "European"
Data$region[Data$Country %in% c("Canada","US","Brazil" )] <- "Americas"
table(Data$region)
table(is.na(Data$region))

##age
table(Data$age_mean)
table(is.na(Data$age_mean))
Data$age_bigger60 <- ifelse(Data$age_mean<=60,"≤ 60 year","> 60 year")
table(Data$age_bigger60)
###BMI
table(Data$bmi_mean)
table(is.na(Data$bmi_mean))
Data$BMI_bigger27_NA <- ifelse(is.na(Data$bmi_mean)==TRUE,"NA",
                               ifelse(Data$bmi_mean<=27,"≤ 27 kg/m2","> 27 kg/m2"))
table(Data$BMI_bigger27_NA)
Data$BMI_m_all <- Data$bmi_mean
Data$BMI_m_all[is.na(Data$bmi_mean)==TRUE] <- median(Data$bmi_mean,na.rm=TRUE)
Data$BMI_bigger27 <- ifelse(Data$BMI_m_all<=27,"≤ 27 kg/m2","> 27 kg/m2")
table(Data$BMI_bigger27)
table(Data$BMI_bigger27,Data$BMI_bigger27_NA)
###type and duration
table(is.na(Data$type))
##durtaion
table(is.na(Data$duration))
plot(Data$duration,Data$duration_new)
Data$duration <- as.numeric(Data$duration)
mean(Data$duration)
median(Data$duration)
hist(Data$duration)
Data$duration_new <- as.numeric(Data$duration_new)
mean(Data$duration_new)
median(Data$duration_new)
hist(Data$duration_new)
Data$duration_bigger8 <- ifelse(Data$duration_new<=8,"≤ 8 years","> 8 years")
table(Data$duration_bigger8)
##sample size
table(Data$N)
table(is.na(Data$N))
median(Data$N)
Data$N_bigger2000 <- ifelse(Data$N<2000,"< 2000","≥ 2000")
table(Data$N_bigger2000 )
##qulity
quality_CHD <- read_excel("quality_CHD.xlsx")
names(quality_CHD)
quality_CHD <- unique(quality_CHD %>% select("ID","NOS","Quqlity"))

Data <- merge(Data,quality_CHD,by="ID",all.x = T)
table(Data$Quqlity)
table(is.na(Data$Quqlity))

table(Data$sexhormones_used)
table(Data$adjust_sexH)
table(is.na(Data$adjust_sexH))
table(is.na(Data$adjust_disease))
table(Data$adjust_disease)
table(is.na(Data$adjust_BMI))
table(Data$adjust_BMI)

subgroup_var <- c("type","region","age_bigger60","BMI_bigger27_NA","BMI_bigger27","N_bigger2000","duration_bigger8","Quqlity","adjust_sexH","adjust_disease","adjust_BMI")
subgroup_names <- subgroup_var 
Result <- data.frame()
for(j in 2){
  for(k in 1:2){
    for(i in 1:length(sexHormone)){
      data_meta <- subset(Data,Outcome==Disease[j] &
                            Exposure==sexHormone[i]
                          & sex==sexvar[k])
      disease <- Disease[j]
      sexH <- sexHormone[i]
      num_strudy <- dim(data_meta)[1]
      Sex <- sexvar[k]
      if(dim(data_meta)[1]>1){
        data_meta$sub <- NA
        for(n in 1:length(subgroup_var)){
          #meta analysis
          data_meta$sub <- data_meta[,subgroup_var[n]]
          m <- metagen(TE=lnhr, seTE, 
                       studlab=paste(Study,sep=","), 
                       n.e=Events, n.c=N, 
                       comb.fixed=F, comb.random=T,
                       label.e=gs('Sample size'), label.c=gs('Cases'), 
                       sm="HR", method.tau = "DL", 
                       hakn=F, prediction = F, backtransf=T, 
                       byvar=sub, 
                       overall=T,
                       data=data_meta)
          #meta result
          sm <- summary(m)
          No <- sm$k.w
          subs <- sm$bylevs
          hr <- sprintf("%0.2f",exp(sm$within.random$TE))
          HR_l <- sprintf("%0.2f",exp(sm$within.random$lower))
          HR_u <- sprintf("%0.2f",exp(sm$within.random$upper))
          HR <- paste0(hr," (",HR_l,", ",HR_u,")")
          I2 <- paste0(sprintf("%0.1f",sm$I2.w$TE*100),"%")
          p <- sprintf("%0.3f",sm$pval.Q.b.random)
          P <- c(p,NA)
          sub_name <- data.frame(Sex,sexH,disease,subgroup_names[n],NA,NA,NA,p)
          result <- data.frame(Sex,sexH,disease,subs,No,HR,I2,NA)
          names(sub_name) <- names(result)
          Result0 <- rbind(sub_name,result)
          Result0$goup <- subgroup_var[n]
          Result <- rbind(Result,Result0)
        }
      }
    }
  }
}

write.csv(Result,"result/submeta_result.csv")


#manuscript
d <- Data %>% select(ID,Study)
d <- unique(d)
length(unique(Data$ID)) #27
#total participant
participant <- unique(Data[,c("id","sex","N")])
participant <- participant %>% group_by(id,sex) %>%
  mutate(n=max(N))
participant <- unique(participant[,-3])
participant <- participant %>% group_by(sex) %>% mutate(total=sum(n))
#events
event <- unique(Data[,c("id","sex","Outcome","Events")])
event <- event %>% group_by(sex,Outcome) %>% mutate(total=sum(Events))
event <- unique(event[,c(2,3,5)])
#sex
sex_num <- unique(Data[,c("ID","sex")])
table(sex_num$sex)

type_num <-  unique(Data[,c("ID","sex","type","Cohort_short_namename")])
d <- table(type_num$Cohort_short_namename[type_num$sex=="Men"],type_num$type[type_num$sex=="Men"])
d <- table(type_num$Cohort_short_namename[type_num$sex=="Women"],type_num$type[type_num$sex=="Women"])

are_num <-  unique(Data[,c("ID","region")])
table(are_num$region)


men <- unique(Data[,c("ID","sex","type","Cohort_short_namename")])
table(men$sex,men$type)
men <- unique(Data[,c("ID","sex","region")])
table(men$sex,men$region)
#quarlity
qua <- unique(Data[,c("ID","quality")])
table(qua$quality)
#length
time <- as.numeric(Data$duration)
max(time,na.rm=TRUE)
min(time,na.rm=TRUE)

#age
max(Data$Age[Data$sex=="Men"])
age_man <- Data[Data$sex=="Men",]
table(Data$sex)