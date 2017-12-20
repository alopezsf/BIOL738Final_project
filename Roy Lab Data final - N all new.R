#--Liana Bruggemann & Alejandro Lopez, 2017

#--Note as of Dec 11, 2017, this code does not use any functions from the seqinr package

#install.packages ("seqinr")

#--Loading library for today's R session

#library(seqinr)

#set working directory 
setwd("~/desktop/Bioinformatics")

#--First read in the data from the .csv file

read.table("Roy Lab Als data .csv",
           sep = "\t", 
           stringsAsFactors = FALSE,
           header = FALSE)-> TblSplicing

N = 2 #--number of people in each Affected and Unaffected set

NTot = N*2  #--Total number of people

#--This for loop goes through each of the rows in TblSplicing and adds in the average of the two A splices divided by total A reads and  
#--adds in the average of the two U splices divided by total U reads 

for (ir in 1:nrow(TblSplicing)) {
  icIn = 2
  icOut = 1+NTot+1+1
  for (i in (1:2)) {
    for (ic in 1:N) {
      lst = strsplit(TblSplicing[ir,icIn],"/")
      icIn = icIn+1
      
      TotReads = as.integer(unlist(lst)[1])+as.integer(unlist(lst)[2])
      
      if (TotReads == 0) 
        TblSplicing[ir,icOut] = -1
      else 
        TblSplicing[ir,icOut] = as.integer(unlist(lst)[1]) / TotReads #--Ratio of Spliced to Total Reads
      
      icOut = icOut+1
    }
    
    TblSplicing[ir,icOut]  = mean(as.numeric(TblSplicing[ir,(icOut-N):(icOut-1)]))#--Average of all ratios
    icOut = icOut+1
  }
  
  if ((TblSplicing[ir,icOut-1] > 0) & (TblSplicing[ir,icOut-1] < 1) & (TblSplicing[ir,icOut-1-N] > 0) & (TblSplicing[ir,icOut-1-N] < 1))
    TblSplicing[ir,icOut]  = TblSplicing[ir,icOut-1] - TblSplicing[ir,icOut-1-N] #--Difference of the A and U RAverage Ratios
  else
    TblSplicing[ir,icOut]  = -10 #--used to indicate not ratios we want to use

}

cn = "Spliceosome"

for (i in 1:N) {
  cn <- c(cn, paste("RatioA",i,sep=""))
}

for (i in 1:N) {
  cn <- c(cn, paste("RatioU",i,sep=""))
}

cn <- c(cn, "GeneTranscript")

for (i in 1:N) {
  cn <- c(cn, paste("AS",i,sep=""))
}

cn <- c(cn, "RatioAavg")

for (i in 1:N) {
  cn <- c(cn, paste("US",i,sep=""))
}

cn <- c(cn, "RatioUavg")

cn <- c(cn, "RatioDiff")

names(TblSplicing) <- cn

options(max.print = 10*ncol(TblSplicing)) #--display up to 10,000 rows

TblSplicing   #--show our resulting table

SARatios <- TblSplicing[order(TblSplicing$RatioAavg),]$RatioAavg     
#--Create num (vector) of spliced affected ratios and sort smallet to biggest
SARatios <- SARatios [((SARatios > 0) & (SARatios < 1) )]   
#--Remove all of the ratios that are either 0 or 1 or -1 (# of A or U reads was 0)
SARatios

SURatios <- TblSplicing[order(TblSplicing$RatioUavg),]$RatioUavg      
#--Create num (vector) of spliced unaffected ratios and sort smallet to biggest
SURatios <- SURatios [((SURatios > 0) & (SURatios < 1)) ]   
#--Remove all of the ratios that are either 0 or 1 or -1 (# of A or U reads was 0)
SURatios

SADdiff <- TblSplicing[order(TblSplicing$RatioDiff),]$RatioDiff
SADdiff <- SADdiff[(SADdiff != -10)]
SADdiff

#--Make the Plot of Affected and Unafected Ratios verses index 
plot(SARatios, type="b", pch=22, lty=2, col="Black", xlab = "", ylab="Splices / Reads", col.lab="red")
lines(SURatios, type="o", pch=23, lty=2, col="red")
title(main=paste(c("Affected and Unaffected Ratios, started with ",nrow(TblSplicing), "total original data points"),collapse=""), col.main="red", font.main=4)
legend("bottomright", 
       inset=c(0.1,0.1),
       bty="o", pt.cex = 1, cex=1.0, col=c("black","red"), pch=22:23, lty=2,
       legend=c(paste(c("Affected ratio, average (of ",length(SARatios) ," points) = " ,format(mean(SARatios),digits=1,nsmall=4)), collapse=""),
                paste(c("Unaffected ratio, average (of ",length(SURatios) ," points) = " ,format(mean(SURatios),digits=1,nsmall=4)), collapse="")))

plot (SADdiff, type="b", pch=22, lty=2, col="Black", xlab = "", ylab="Ratio Difference", col.lab="black")
title(main="Unaffected Ratios - Affected Ratios", col.main="red", font.main=4)

