library(sangerseqR)

setwd("C:/Users/worms/Desktop")
abi <- read.abif("L107 (NNB)6_81.ab1")
sseq <- sangerseq(abi)
chromatogram(sseq)

trace <- traceMatrix(sseq)

plot(trace[1:1000,])

abi$data