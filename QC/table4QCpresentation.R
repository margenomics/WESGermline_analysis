#Júlia Perera
#12/02/2020
#Script to generate QC table from multiQC


project <-"20200227_NGarcia_DOvejero_WES"

if (Sys.info()['sysname'] == "Windows") {
  wd <- paste("Z:",project,"QC",sep="/")
  codeDir <-"Z:/R Code"
  setwd(wd)
}  else if (Sys.info()['user'] == "julia") {
  wd <- paste("/bicoh/MARGenomics/",project,"QC",sep="/")
  codeDir <-"/bicoh/MARGenomics/RCode"
  setwd(wd)
} else {
  wd<-paste("/home/marnal/samprojects/Dades_affymetrix/",project,"QC",sep="/")
  codeDir<-"/home/marnal/samprojects/Dades_affymetrix//R Code"
  setwd(wd)
}



library(openxlsx)


table=readWorkbook(file.path(wd,"multiqc_table.xlsx"))
samples=unique(gsub("_L00.*_R.*_001","",table$Sample.Name))


lanes=4
R=2
t=R*lanes

# Prepare table
df=data.frame(matrix(NA,ncol = R*lanes,nrow=length(samples)))
colnames(df)=c("R1_L001","R2_L001","R1_L002","R2_L002","R1_L003","R2_L003","R1_L004","R2_L004")
rownames(df)=samples

# Fill table
for (i in 1:length(samples)){
  df[i,]=t(table[(i*t-(t-1)):(i*t),"M.Seqs"])
}

#Reorder first all R1, then all R2
df=df[,c(seq(1,t,2),seq(2,t,2))]

# Add total Lib.Size per sample
df$`Total (M.reads)`=rowSums(df)

# Save data to file
wb <- createWorkbook()
## Add worksheets
addWorksheet(wb, "Lib.Size")
writeData(wb, "Lib.Size", df,  rowNames = T)
saveWorkbook(wb, file.path(wd,"table4QCpresentation.xlsx"), overwrite = TRUE)
