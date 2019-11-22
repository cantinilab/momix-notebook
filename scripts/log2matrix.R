## Support function to apply log2(+1) to a matrix
log2matrix <- function(folder, file.name){
    # Read table
    data <- as.matrix(read.table(paste(folder,file.name,sep="/"),sep=" ",row.names=1,header=TRUE))
    # Apply transformation
    data <- log2(data+1)
    # Output file name
    output <- paste(folder,paste0("log_",file.name), sep="/")
    # Export transformed data
    write.table(data,output,sep=" ", col.names=TRUE, row.names=TRUE)
    # ?
    system(paste("sed -i '1s/^/probe\t/'", output, sep=" "))
}