library(limma)

hpos <- -log10(0.05)
vpos <- c(1,-1)

targets03  <- readTargets(file="Targets-03-06-09-10.txt")
RG03       <- read.maimages(targets03$Filename, source="agilent")
w03        <- modifyWeights(array(1,dim(RG03)), RG03$genes$ControlType, c("1","-1"), c(0,0))

bgCorr03   <- backgroundCorrect(RG03, method="normexp")
nWithin03  <- normalizeWithinArrays(bgCorr03, method="loess", weights=w03)
nBetween03 <- normalizeBetweenArrays(nWithin03, method="quantile")

design03   <- modelMatrix(targets03,ref="HK")
biolrep03  <- c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8)
corfit03   <- duplicateCorrelation(nBetween03, ndups = 1, block = biolrep03)

fit03      <- lmFit(nBetween03, design03, block = biolrep03, cor = corfit03$consensus)
fit03      <- eBayes(fit03)

write.fit(fit03, results=NULL, "Targets-03-06-09-10_20110208.txt", digits=3, adjust="BH", sep="\t")

logValues03 <- as.matrix.MAList(nWithin03)

write.table(logValues03, file="Targets03_logValues.txt")

fitData03   <- as.data.frame(fit03)
fitData03S  <- subset(fitData03, select = c(genes.ProbeUID,coefficients.KK,coefficients.KM,coefficients.HM,p.value.KK,p.value.KM,p.value.HM))
fitData03SP <- subset(fitData03S, select = c(genes.ProbeUID,coefficients.HM,p.value.HM))

HMPlot <- ggplot(fitData03SP,aes(coefficients.HM,-log10(p.value.HM),color=coefficients.HM)) + 
    geom_point() + geom_hline(yintercept=hpos,color="grey") + geom_vline(xintercept=vpos,color="grey") + 
    theme_bw() + scale_color_gradient2(low="green",mid="black",high="red",name="") +
    xlab("log2(Fold-change)") + ylab("-log10(P-value)") +
    theme(legend.justification=c(1,1), legend.position=c(1,1))

ggsave(HMPlot,filename="hm_plot.png",dpi=600)