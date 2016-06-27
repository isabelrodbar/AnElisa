rm(list=ls())
setwd("~/Dropbox/Malaria/Exposure_R01/Namibia Serology/anELISA/AMA1/")
options(warn=1)


####################################################################################################
#INSTALL REQUIRED PACKAGES: Hmisc, RColorBrewer
library(Hmisc)
library(RColorBrewer)
source("001.anELISA_functions.r")


####################################################################################################
#ASSIGN FILE NAME FOR THE MASTER DOCUMENT
master.filename <- "master.AMA1.munya.txt" 


##################r##################################################################################  
#generate outputs
std_curves()
AnELISA()

