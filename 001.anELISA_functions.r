####################################################################################################
#CREATE A VECTOR OF WELL LOCATIONS IN THE 96 WELL PLATE
wells <- rep(NA,96) #Lookup for well names to indexes 1-96
n <- 1
for (i in 1:8){
  for (j in 1:12){
    wells[n] <- paste(LETTERS[i],formatC(j,width=2,flag='0'),sep="")
    n <- n+1
  }
}


####################################################################################################
get.std <- function(template,ODs){
  n.std.dup <- length(unique(template[grep("STD",template$Sample,ignore.case=T),"Duplicate"]))
  std.dil <- unique(template[grep("STD",template$Sample,ignore.case=T),"Dilution"])
  std <- matrix(ncol=n.std.dup+2,nrow=length(std.dil))
  colnames(std) <- c("Dil","Conc",LETTERS[1:n.std.dup])
  std[,1] <- sort(std.dil, decreasing=T)
  std[,2] <- sort(1/std.dil)
  for (i in 1:n.std.dup){
    for (j in 1:nrow(std)){
      well <- template[((toupper(template$Sample)=="STD") & 
                          (template$Duplicate==LETTERS[i]) & 
                          (template$Dilution==std[j,1])), "ELISA.Well"]
      std[j,(i+2)] <- ODs[wells==well]
    }
  }
  blank.wells <- template[toupper(template$Sample)=="BLANK", "ELISA.Well"]
  blanks <- NULL
  for (i in blank.wells){
    blanks <- c(blanks,ODs[wells==i])
  }
  return(list(std=std, blanks=blanks))
}


####################################################################################################  
get.conc <- function(OD,MaxOD,MinOD,Xmid,Scale){  #Note: expects Xmid to be on log scale
  conc <- exp(Xmid - (Scale * log((MaxOD-MinOD)/(OD-MinOD)-1)))
  return(conc)
}


####################################################################################################  
AnELISA <- function(){
  
  print(paste("Opening master file:",master.filename))
  master <- read.delim(master.filename, stringsAsFactors=F)
  print("File opened.")
  
  output.data <- NULL
  
  for (r in 1:nrow(master)){   #Step through master file one line at a time
    
    #####################
    # Get ODs from file #
    #####################
    print(paste("Opening ELISA file:",master$ELISA.file.name[r]))
    txtlines <- readLines(paste("SoftMaxPro files/", master$ELISA.file.name[r], sep=""), skipNul=T)
    print("File opened.  Reading ODs.")
    d <- NULL
    
    antigen <- paste(master$Antigen.specific.template[r], master$Sample.file.ID[r], sep=".")
    read <- 1
    ODs <- NULL
    ODs <- c(ODs,as.numeric(strsplit(txtlines[44],"\t")[[1]])[3:98])
    d[[antigen]][[read]] <- ODs
    
    if (!is.null(d)) print(paste("ODs read successfully for",paste(names(d),collapse=", "),collapse=" ")) else print("Error reading file.")
    ag.restrict <- master$Antigen.specific.template[r]
    if (!is.na(ag.restrict) & !(ag.restrict=="")){
      print(paste("Restricting data to",ag.restrict))
      if (length(grep(ag.restrict,names(d))) > 0){
        d <- d[grep(ag.restrict,names(d))]
      } else print("Error, unable to find this antigen name")
    }
    
    ########################
    # Show standard curves #
    ########################
    print(paste("Opening template file:",master$Template.file.name[r]))
    template <- read.delim(master$Template.file.name[r], stringsAsFactors=F)
    print("File opened.")
    
    template$Sample <- sub("SAMPLE", "", template$Sample, ignore.case=TRUE)
    
    cols <- brewer.pal(8,"Set1")
    for (antigen in names(d)){
      reads <- d[[antigen]]
      par(mfrow=c(2,ceiling(length(reads)/2)))
      for (read in 1:length(reads)){
        temp <- get.std(template, reads[[read]])
        std <- temp$std
        if(std[7,3]<max(c(std[6,3], std[6,4]))) std[7,3] <- max(c(std[6,3], std[6,4]))
        if(std[7,4]<max(c(std[6,3], std[6,4]))) std[7,4] <- max(c(std[6,3], std[6,4]))
        blanks <- temp$blanks
        plot(log(std[,"Conc"]),std[,3],ylim=c(0,4),type="b",ylab="OD",
             xlab="Dilution factor",col=cols[1],xaxt="n")
        axis(side=1,at=log(std[,"Conc"]),labels=std[,"Dil"])
        mtext(paste("Read",read),line=-2,font=2)
        for (i in 4:ncol(std)){
          points(log(std[,"Conc"]),std[,i],col=cols[i-2],type="b")
        }
        abline(h=blanks,lty=3)
        if (read==1){
          legend(x=log(std[1,"Conc"]),y=3,legend=c("Blanks",colnames(std)[3:ncol(std)]),
                 col=c("black",cols[1:(ncol(std)-2)]),lty=c(3,rep(1,ncol(std))))
          mtext(paste("STD curves for",antigen),outer=T,line=-2)
        }
      }
      
      ################
      ## analyze it ##
      ################
      temp <- get.std(template, ODs)
      std <- temp$std
      if(std[7,3]<max(c(std[6,3], std[6,4]))) std[7,3] <- max(c(std[6,3], std[6,4]))
      if(std[7,4]<max(c(std[6,3], std[6,4]))) std[7,4] <- max(c(std[6,3], std[6,4]))
      blanks <- temp$blanks
      x <- y <- NULL
      for (i in 3:ncol(std)){
        y <- c(y,std[,i])
        x <- c(x,std[,"Conc"])
      }
      fit.data <- data.frame("OD"=y,"Conc"=x)
      MinOD <- mean(blanks)
      if (MinOD < 0) MinOD <- 0
      fit.data$OD <- fit.data$OD - MinOD   #Taking out the MinOD from the fit so MaxOD is now the upper asymptote
      fit.data$Conc <- log(fit.data$Conc)  #x values now log transformed, so Xmid returned will be log transformed
      par(mfrow=c(1,1))
      revise.plot <- TRUE
      while (revise.plot){
        plot(log(std[,"Conc"]),std[,3],ylim=c(0,4),ylab="OD",
             xlab="Dilution factor",col=cols[1],xaxt="n")
        axis(side=1,at=log(std[,"Conc"]),labels=std[,"Dil"])
        mtext(antigen,line=2,font=2)
        for (i in 4:ncol(std)){
          points(log(std[,"Conc"]),std[,i],col=cols[i-2])
        }
        abline(h=blanks,lty=3)
        legend(x=log(std[1,"Conc"]),y=4,legend=c("Blanks",colnames(std)[3:ncol(std)],"Fit"),
               col=c("black",cols[1:(ncol(std)-2)],"black"),lty=c(3,rep(1,ncol(std)-2),2),
               lwd=c(rep(1,ncol(std)-1),2), cex=0.75)
        
        starting.values <- getInitial(OD ~ SSlogis(Conc,MaxOD,Xmid,Scale), data=fit.data)
        fit <- nls(OD ~ (MaxOD / (1 + exp((Xmid - Conc)/Scale))),data=fit.data, start=starting.values)
        params <- as.list(coef(fit))
        params$MaxOD <- params$MaxOD + MinOD
        x <- seq(min(fit.data$Conc),max(fit.data$Conc),length=100)
        y <- predict(fit, newdata=list(Conc=x)) + MinOD
        lines(x,y,lty=2,lwd=2)
        answer <- readline("\nRemove any outliers? (y/[n]) ")
        if (toupper(answer)=="Y"){
          mtext("Select outlier with the mouse.",col="red",cex=1.5)
          outlier <- locator(n=1)
          conc.toremove <- fit.data[which.min(abs(outlier$x - fit.data$Conc)),"Conc"]
          poss.ODs <- fit.data[fit.data$Conc==conc.toremove,"OD"] + MinOD
          OD.toremove <- poss.ODs[which.min(abs(outlier$y - poss.ODs))]
          points(conc.toremove,OD.toremove,pch=4,cex=2,col="red")
          answer <- readline("\nCorrect outlier? ([y]/n) ")
          if (toupper(answer)!="N"){
            fit.data <- fit.data[!(((fit.data$OD + MinOD) == OD.toremove) & (fit.data$Conc == conc.toremove)),]
            std.row.toremove <- which(log(std[,"Conc"])==conc.toremove)
            std[std.row.toremove,std[std.row.toremove,]==OD.toremove] <- NA
          }
        } else revise.plot <- F
      } # while revise.plot
      if (master$Sample.file.name[r] != ""){
        upperbound <- 100
        mtext("Define an upper bound with the mouse.",col="red",cex=1.5)
        upperbound <- locator(n=1)$y
        abline(h=upperbound, col="red", lty=3)
        print(paste("Opening sample file:",master$Sample.file.name[r]))
        samples <- read.delim(master$Sample.file.name[r], stringsAsFactors=F)
        print("File opened.")
        samples <- subset(samples, (Sample.file.ID==master$Sample.file.ID[r]) & (toupper(Sample.Name)!="EMPTY"))
        samples$Sample.Number <- as.integer(samples$Sample.Number)
        options(warn=-1)
        template.samples <- template[!is.na(as.integer(template$Sample)),]
        options(warn=1)
        template.samples$Sample <- as.integer(template.samples$Sample)
        samples.merged <- merge(template.samples, samples, by.x="Sample", by.y="Sample.Number")
        samples.merged$Flag1 <- ""
        for (i in 1:nrow(samples.merged)){
          OD <- ODs[wells==samples.merged$ELISA.Well[i]]
          samples.merged[i,"OD"] <- OD
          if (OD <= MinOD){
            samples.merged[i,"Flag1"] <- "Low"
          } else if (OD >= params$MaxOD){
            samples.merged[i,"Flag1"] <- "High"
          }
        }
        samples.merged$Flag <- ""
        for (i in 1:nrow(samples.merged)){
          OD <- samples.merged[i,"OD"]
          samples.merged[i,"Conc"] <- samples.merged[i,"Dilution"] * get.conc(OD=samples.merged$OD[i], MaxOD=params$MaxOD, MinOD=MinOD, Xmid=params$Xmid, Scale=params$Scale)
          if (OD >= upperbound) samples.merged[i,"Flag"] <- "AboveUpperbound"
          if (OD >= max(fit.data$OD + MinOD)) samples.merged[i,"Flag"] <- "AboveMaxStd"
          if (OD < min(fit.data$OD + MinOD)) samples.merged[i,"Flag"] <- "BelowMinStd"
        }
        samples.merged <- samples.merged[order(samples.merged$Sample.Name, samples.merged$Duplicate),]
        add.to.output <- samples.merged[c("Sample.Name","Study.ID","Barcode","OD","Conc","Flag1", "Flag","Duplicate","ELISA.Well","Sample.file.ID","Dilution", "Coating", "Notes")]
        add.to.output$Antigen <- gsub("[.].+$", "", antigen)
        add.to.output$Read <- read
        add.to.output$ELISA.date <- as.Date(gsub("[.].+$", "", master$ELISA.file.name[r]), "%y%m%d")
        add.to.output$UpperBound <- upperbound
        add.to.output$MinOD <- MinOD
        add.to.output$MaxOD <- params$MaxOD
        add.to.output$MaxStd <- max(fit.data$OD + MinOD)
        add.to.output$MinStd <- min(fit.data$OD + MinOD)
        add.to.output$Std.MaxDil <- min(data.frame(std)$Dil)
        output.data <- rbind(output.data, add.to.output)
        flagged <- which(samples.merged$Flag != "")
        borders <- rep("black",nrow(samples.merged))
        borders[flagged] <- "red"
        densities <- rep(-1,nrow(samples.merged))
        densities[flagged] <- 20
        barplot(samples.merged$Conc, names.arg=paste(samples.merged$Sample.Name,samples.merged$Sample.description,samples.merged$Duplicate),
                las=2, col=rep(cols,each=2), main=paste(antigen,"Read",read),border=borders,density=densities,cex.names=0.8,ylim=c(0,min(50,max(samples.merged$Conc,na.rm=T))))
        readline("\nHit enter to continue.")
      } #calculating concentrations for samples if Sample.file.name specified in master
    } #for each "antigen" in names(d)
    write.csv(output.data,file=paste(".\\anELISA Outputs\\", gsub("-", "_", Sys.Date()), ".AMA_output.RAW.csv", sep=""))
    save(output.data,file=paste(".\\anELISA Outputs\\", gsub("-", "_", Sys.Date()), ".AMA_output.RAW.r", sep=""))
  } #for each row "r" in master
  dev.off()
} #function wrapper


####################################################################################################  
std_curves <- function(){
  
  print(paste("Opening master file:",master.filename))
  master <- read.delim(master.filename, stringsAsFactors=F)
  print("File opened.")
  
  output.data <- NULL
  
  pdf(file=paste("anELISA Outputs/", gsub("-", "_", Sys.Date()), ".AMA_std_curves.pdf", sep=""),width=8,height=8)
  par(mfrow=c(2,2))
  for (r in 1:nrow(master)){   #Step through master file one line at a time
    
    #####################
    # Get ODs from file #
    #####################
    print(paste("Opening ELISA file:",master$ELISA.file.name[r]))
    txtlines <- readLines(paste("SoftMaxPro files/", master$ELISA.file.name[r], sep=""),skipNul=T)
    print("File opened.  Reading ODs.")
    d <- NULL
    
    antigen <- paste(master$Antigen.specific.template[r], master$Sample.file.ID[r], sep=".")
    read <- 1
    ODs <- NULL
    ODs <- c(ODs,as.numeric(strsplit(txtlines[44],"\t")[[1]])[3:98])
    d[[antigen]][[read]] <- ODs
    
    if (!is.null(d)) print(paste("ODs read successfully for",paste(names(d),collapse=", "),collapse=" ")) else print("Error reading file.")
    ag.restrict <- master$Antigen.specific.template[r]
    if (!is.na(ag.restrict) & !(ag.restrict=="")){
      print(paste("Restricting data to",ag.restrict))
      if (length(grep(ag.restrict,names(d))) > 0){
        d <- d[grep(ag.restrict,names(d))]
      } else print("Error, unable to find this antigen name")
    }
    
    ########################
    # Show standard curves #
    ########################
    print(paste("Opening template file:",master$Template.file.name[r]))
    template <- read.delim(master$Template.file.name[r], stringsAsFactors=F)
    print("File opened.")
    
    template$Sample <- sub("SAMPLE", "", template$Sample, ignore.case=TRUE)
    
    cols <- brewer.pal(8,"Set1")
    for (antigen in names(d)){
      reads <- d[[antigen]]
      for (read in 1:length(reads)){
        temp <- get.std(template, reads[[read]])
        std <- temp$std
        blanks <- temp$blanks
        plot(log(std[,"Conc"]),std[,3],ylim=c(0,4),type="b",ylab="OD",
             xlab="Dilution factor",col=cols[1],xaxt="n")
        axis(side=1,at=log(std[,"Conc"]),labels=std[,"Dil"])
        mtext(paste("STD curves for",antigen),line=-2,font=2)
        for (i in 4:ncol(std)){
          points(log(std[,"Conc"]),std[,i],col=cols[i-2],type="b")
        }
        abline(h=blanks,lty=3)
        if (read==1){
          legend(x=log(std[1,"Conc"]),y=3,legend=c("Blanks",colnames(std)[3:ncol(std)]),
                 col=c("black",cols[1:(ncol(std)-2)]),lty=c(3,rep(1,ncol(std))))
        }
      }
    }
  }
  dev.off()
}