library(limma)
library(ggplot2)

hpos <- -log10(0.05)
vpos <- c(1,-1)

targets01  <- readTargets(file="Targets-01-04-07-10.txt")
RG01       <- read.maimages(targets01$Filename, source="agilent")
w01        <- modifyWeights(array(1,dim(RG01)), RG01$genes$ControlType, c("1","-1"), c(0,0))

bgCorr01   <- backgroundCorrect(RG01, method="normexp")
nWithin01  <- normalizeWithinArrays(bgCorr01, method="loess", weights=w01)
nBetween01 <- normalizeBetweenArrays(nWithin01, method="quantile")

design01   <- modelMatrix(targets01,ref="HK")
biolrep01  <- c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8)
corfit01   <- duplicateCorrelation(nBetween01, ndups = 1, block = biolrep01)

fit01      <- lmFit(nBetween01, design01, block = biolrep01, cor = corfit01$consensus)
fit01      <- eBayes(fit01)

write.fit(fit01, results=NULL, "Targets-01-04-07-10_20110208.txt", digits=3, adjust="BH", sep="\t")

logValues01 <- as.matrix.MAList(nWithin01)

write.table(logValues01, file="Targets01_logValues.txt")

fitData01   <- as.data.frame(fit01)
fitData01S  <- subset(fitData01, select = c(genes.ProbeUID,coefficients.KK,coefficients.KLa,coefficients.HLa,p.value.KK,p.value.KLa,p.value.HLa))
fitData01SP <- subset(fitData01S, select = c(genes.ProbeUID,coefficients.HLa,p.value.HLa))

HLaPlot <- ggplot(fitData01SP,aes(coefficients.HLa,-log10(p.value.HLa),color=coefficients.HLa)) + 
    geom_point() + geom_hline(yintercept=hpos,color="grey") + geom_vline(xintercept=vpos,color="grey") + 
    theme_bw() + scale_color_gradient2(low="green",mid="black",high="red",name="") +
    xlab("log2(Fold-change)") + ylab("-log10(P-value)") +
    theme(legend.justification=c(1,1), legend.position=c(1,1))

ggsave(HLaPlot,filename="hla_plot.png",dpi=600)
