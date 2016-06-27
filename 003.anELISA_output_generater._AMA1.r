rm(list=ls())
setwd("C:\\Users\\HelbD\\Documents\\Magic Briefcase\\Namibia Serology\\anELISA\\AMA1")
options(warn=1)


####################################################################################################
print(load(paste(".\\anELISA Outputs\\", gsub("-", "_", Sys.Date()), ".AMA_output.Raw.r", sep="")))


####################################################################################################  
#Correct concentrations if Flag==Low or High
for(i in 1:nrow(output.data)){
  if(output.data[i,"Flag1"]=="Low"){
    output.data[i,"Conc"] <- 0
  }
  if(output.data[i,"Flag1"]=="High"){
    output.data[i,"Conc"] <- 1.5*(output.data[i,"Dilution"]/output.data[i,"Std.MaxDil"])
  }
}


####################################################################################################  
#Investigate replicates & see if the duplicates diverge by too much
AMA1.output <- output.data[duplicated(output.data$Sample.Name)==F, names(output.data) %in% c("Sample.Name", "Barcode", "Notes")]
for(i in AMA1.output$Sample.Name){
  AMA1.output[AMA1.output$Sample.Name==i, "Conc"] <- mean(output.data[output.data$Sample.Name==i, "Conc"], na.rm=T)
  AMA1.output[AMA1.output$Sample.Name==i, "SD"] <- sd(output.data[output.data$Sample.Name==i, "Conc"], na.rm=T)
}
AMA1.output$coef_of_var <- AMA1.output$SD/AMA1.output$Conc
AMA1.output$coef_of_var[AMA1.output$Conc==0 & AMA1.output$SD==0] <- 0

for(i in AMA1.output$Sample.Name){
  if(AMA1.output[AMA1.output$Sample.Name==i,"coef_of_var"]>=1.5 & !is.na(AMA1.output[AMA1.output$Sample.Name==i,"coef_of_var"])){
    AMA1.output[AMA1.output$Sample.Name==i,"Conc"] <- NA
    AMA1.output[AMA1.output$Sample.Name==i,"SD"] <- NA
  }
}


#################################################################################################### 
save(AMA1.output, file=paste(".\\anELISA Outputs\\", gsub("-", "_", Sys.Date()), ".AMA_output.FINAL.r", sep=""))
write.csv(AMA1.output, file=paste(".\\anELISA Outputs\\", gsub("-", "_", Sys.Date()), ".AMA_output.FINAL.csv", sep=""), row.names=F)
