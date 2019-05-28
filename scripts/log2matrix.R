log2matrix <- function(folder,file.name){
  data<-as.matrix(read.table(paste(folder,file.name,sep="/"),sep=" ",row.names=1,header=T))
  data<-log2(data+1)
  output<-paste(folder,paste("log_",file.name,sep=""),sep="/")
  write.table(data,output,sep=" ", col.names=T, row.names=T)
  system(paste("sed -i '1s/^/probe\t/'", output, sep=" "))
}