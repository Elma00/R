##transform 
Data$se <- (log(Data$uphr)-log(Data$lowhr))/(2*1.96)
#dput(table(Data$reference))
Data$HR_0 <- NA
Data$HRlow_0 <- NA
Data$HRup_0 <- NA

ref_h <- c("h2u", "h2","h3", "h4", "h5", "h10","hcc", "hsd")
ref_l <- c("l2u","l3", "l4", "l5","lc", "lcc", "lcu", "lsd","m3", "m4", "m5" )


#change the direction:h4,h2u
Data$HR_0[Data$reference %in% ref_h] <- exp(log(1/Data$hr[Data$reference %in% ref_h]))
Data$HRlow_0[Data$reference %in% ref_h] <- exp(log(1/Data$hr[Data$reference %in% ref_h])-Data$se[Data$reference %in% ref_h]*1.96)
Data$HRup_0[Data$reference %in% ref_h] <- exp(log(1/Data$hr[Data$reference %in% ref_h])+Data$se[Data$reference %in% ref_h]*1.96)

#unchange ther direction
Data$HR_0[Data$reference %in% ref_l] <- Data$hr[Data$reference %in% ref_l]
Data$HRlow_0[Data$reference %in% ref_l] <- Data$lowhr[Data$reference %in% ref_l]
Data$HRup_0[Data$reference %in% ref_l] <- Data$uphr[Data$reference %in% ref_l]

Data$Exposure <- Data$Exposure2
Data$Exposure[Data$Exposure %in% c("FT","TSHBG")] <- "BT"
Data$Exposure[Data$Exposure %in% c("DHEA","DHT")] <- "DHEAS"
Data$Exposure[Data$Exposure %in% c("ESHBG","FE2")] <- "BE2"
table(Data$Exposure)


Data$HR <- Data$HR_0
Data$HRlow <- Data$HRlow_0
Data$HRup <- Data$HRup_0
#Q2: h2,l2
Q2 <- c("h2","l2")
Data$HR[Data$reference %in% Q2] <- exp(-log(Data$hr[Data$reference %in% Q2])*2.54/1.695)
Data$HRlow[Data$reference %in% Q2] <- exp(log(Data$HR)[Data$reference %in% Q2]-Data$se[Data$reference %in% Q2]*2.54/1.695)
Data$HRup[Data$reference %in% Q2] <- exp(log(Data$HR)[Data$reference %in% Q2]+Data$se[Data$reference %in% Q2]*2.54/1.695)
#Q3: h3, l3
Q3 <- c("h3","l3")
Data$HR[Data$reference %in% Q3] <- exp(log(Data$HR_0[Data$reference %in% Q3])*2.54/2.18)
Data$HRlow[Data$reference %in% Q3] <- exp(log(Data$HR)[Data$reference %in% Q3]-Data$se[Data$reference %in% Q3]*2.54/2.18*1.96)
Data$HRup[Data$reference %in% Q3] <- exp(log(Data$HR)[Data$reference %in% Q3]+Data$se[Data$reference %in% Q3]*2.54/2.18*1.96)#Q5
#Q5
Q5 <- c("l5","h5")
Data$HR[Data$reference %in% Q5] <- exp(log(Data$HR_0[Data$reference %in% Q5])*2.54/2.8)
Data$HRlow[Data$reference %in% Q5] <- exp(log(Data$HR)[Data$reference %in% Q5]-Data$se[Data$reference %in% Q5]*2.54/2.8*1.96)
Data$HRup[Data$reference %in% Q5] <- exp(log(Data$HR)[Data$reference %in% Q5]+Data$se[Data$reference %in% Q5]*2.54/2.8*1.96)
#Q10
Q10 <- c("h10","l10")
Data$HR[Data$reference %in% Q10] <- exp(-log(Data$hr[Data$reference %in% Q10])*2.54/4.7)
Data$HRlow[Data$reference %in% Q10] <- exp(log(Data$HR)[Data$reference %in% Q10]-Data$se[Data$reference %in% Q10]*2.54/4.7)
Data$HRup[Data$reference %in% Q10] <- exp(log(Data$HR)[Data$reference %in% Q10]+Data$se[Data$reference %in% Q10]*2.54/4.7)
#per sd
Psd <- c("lsd","hsd")
Data$HR[Data$reference %in% Psd] <- exp(log(Data$HR_0[Data$reference %in% Psd])*2.54)
Data$HRlow[Data$reference %in% Psd] <- exp(log(Data$HR)[Data$reference %in% Psd]-Data$se[Data$reference %in% Psd]*2.54*1.96)
Data$HRup[Data$reference %in% Psd] <- exp(log(Data$HR)[Data$reference %in% Psd]+Data$se[Data$reference %in% Psd]*2.54*1.96)
#per lc or hc
c <- c("lc","hc")
Data$exposure_sd <- as.numeric(Data$exposure_sd)
Data$HR[Data$reference %in% c] <- exp(log(Data$HR_0[Data$reference %in% c])*(Data$exposure_sd[Data$reference %in% c]*2.54))
Data$HRlow[Data$reference %in% c] <- exp(log(Data$HR)[Data$reference %in% c]-Data$se[Data$reference %in% c]*(Data$exposure_sd[Data$reference %in% c]*2.54)*1.96) #1.35
Data$HRup[Data$reference %in% c] <- exp(log(Data$HR)[Data$reference %in% c]+Data$se[Data$reference %in% c]*(Data$exposure_sd[Data$reference %in% c]*2.54)*1.96)
#d <- Data %>% filter(reference %in% c("lc","hc"))
##
cc <- c("lcc","hcc")
d <- Data %>% filter(reference %in% c("lcc","hcc"))
table(Data$comparison)
Data$comparison[Data$comparison=="hc50"] <- 50
Data$comparison[Data$comparison=="lc0.01"] <- 0.01
Data$comparison[Data$comparison=="lc10"] <- 10
Data$comparison[Data$comparison=="lc100"] <- 100
Data$comparison <- as.numeric(Data$comparison)
Data$HR[Data$reference %in% cc] <- exp(log(Data$HR_0[Data$reference %in% cc])/Data$comparison[Data$reference %in% cc]*(Data$exposure_sd[Data$reference %in% cc]*2.54))
Data$HRlow[Data$reference %in% cc] <- exp(log(Data$HR)[Data$reference %in% cc]-Data$se[Data$reference %in% cc]/Data$comparison[Data$reference %in% cc]*(Data$exposure_sd[Data$reference %in% cc]*2.54)*1.96)
Data$HRup[Data$reference %in% cc] <- exp(log(Data$HR)[Data$reference %in% cc]+Data$se[Data$reference %in% cc]/Data$comparison[Data$reference %in% cc]*(Data$exposure_sd[Data$reference %in% cc]*2.54)*1.96)
#d <- Data %>% filter(reference %in% c("lcc","hcc"))
 
##
Data$lnhr <- log(Data$HR)
Data$lnlower <- log(Data$HRlow)
Data$lnupper <- log(Data$HRup)
Data$seTE <- (Data$lnupper-Data$lnlower)/(2*1.96)

