library(limma)

hpos <- -log10(0.05)
vpos <- c(1,-1)

targets02  <- readTargets(file="Targets-02-05-08-10.txt")
RG02       <- read.maimages(targets02$Filename, source="agilent")
w02        <- modifyWeights(array(1,dim(RG02)), RG02$genes$ControlType, c("1","-1"), c(0,0))

bgCorr02   <- backgroundCorrect(RG02, method="normexp")
nWithin02  <- normalizeWithinArrays(bgCorr02, method="loess", weights=w02)
nBetween02 <- normalizeBetweenArrays(nWithin02, method="quantile")

design02   <- modelMatrix(targets02,ref="HK")
biolrep02  <- c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8)
corfit02   <- duplicateCorrelation(nBetween02, ndups = 1, block = biolrep02)

fit02      <- lmFit(nBetween02, design02, block = biolrep02, cor = corfit02$consensus)
fit02      <- eBayes(fit02)

write.fit(fit02, results=NULL, "Targets-02-05-08-10_20110208.txt", digits=3, adjust="BH", sep="\t")

logValues02 <- as.matrix.MAList(nWithin02)

write.table(logValues02, file="Targets02_logValues.txt")

fitData02   <- as.data.frame(fit02)
fitData02S  <- subset(fitData02, select = c(genes.ProbeUID,coefficients.KK,coefficients.KE,coefficients.HE,p.value.KK,p.value.KE,p.value.HE))
fitData02SP <- subset(fitData02S, select = c(genes.ProbeUID,coefficients.HE,p.value.HE))

HEPlot <- ggplot(fitData02SP,aes(coefficients.HE,-log10(p.value.HE),color=coefficients.HE)) + 
    geom_point() + geom_hline(yintercept=hpos,color="grey") + geom_vline(xintercept=vpos,color="grey") + 
    theme_bw() + scale_color_gradient2(low="green",mid="black",high="red",name="") +
    xlab("log2(Fold-change)") + ylab("-log10(P-value)") +
    theme(legend.justification=c(1,1), legend.position=c(1,1))

ggsave(HEPlot,filename="he_plot.png",dpi=600)