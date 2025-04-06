rm(list = ls())

args <- commandArgs(trailingOnly = TRUE)

print(args)

reads<-read.csv(file=args[1], sep="", header=FALSE)

myhist <- function(x, ..., breaks="Sturges",
                   main = paste("Histogram of", xname),
                   xlab = xname,
                   ylab = "Frequency") {
  xname = paste(deparse(substitute(x), 500), collapse="\n")
  h = hist(x, breaks=breaks, plot=FALSE)
  plot(h$breaks, c(NA,h$counts), type='S', main=main,
       xlab=xlab, ylab=ylab, axes=FALSE, ...)
  axis(1)
  axis(2)
  lines(h$breaks, c(h$counts,NA), type='s')
  lines(h$breaks, c(NA,h$counts), type='h')
  lines(h$breaks, c(h$counts,NA), type='h')
  lines(h$breaks, rep(0,length(h$breaks)), type='S')
  invisible(h)
}

png(paste0(dirname(args[1]), "/", args[2], "_", args[3], "_filter_reads_size_below_10000.png"))
yo = myhist(reads$V1[which(reads$V1<10000)], log="y", breaks = 200, main = paste0(args[2], " - ", args[3], "\n", "Filter - Read length below 10000"), xlab="Read size")
dev.off()

png(paste0(dirname(args[1]), "/", args[2], "_", args[3], "_no_filter.png"))
yo_2 = myhist(reads$V1, log="y", breaks = 200, main =paste0(args[2], " - ", args[3], "\n", "No Filter"), xlab="Read size")
dev.off()



