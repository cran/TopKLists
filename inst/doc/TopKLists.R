## ----echo=FALSE, results='hide', eval=TRUE, message=FALSE, warning=FALSE----
require(knitr)
options(width=60)

## ----message=FALSE, warning=FALSE-------------------------
library(TopKLists)
data(breast)
head(breast)

## ----eval=TRUE, echo=TRUE, results='hide', message=FALSE, warning=FALSE----
library(TopKLists)

## ----tidy=FALSE-------------------------------------------
k = 30
set.seed = 11
x = c(rep(1,k), rbinom(100, 1, 0.2))

## ----tidy=TRUE--------------------------------------------

v.vect=seq(2,length(x), by=2) #setting up a vector of the nu values

resF=c()

for (v in v.vect)
 {
	res=compute.stream(x, const=0.5, v)
	resF=rbind(resF,c(v,paste(res)))
 }

colnames(resF)=c("v", "j0_est", "k","reason.break", "Js", "v.vector")
head(resF)

table(resF[,2])

## ----truncPlot, out.width='8cm', out.height='8cm', fig.pos="ht!", fig.cap="Estimation of $j_0$ for different values of $\\nu$", fig.align='center'----
plot(resF[,1], resF[,2], pch=19, ylim=c(0, length(x)), 
xlab=substitute(nu), ylab=substitute(paste(hat(j)[0])))
abline(a=31, b=0, col="red")
lines(resF[,1], resF[,2])

## ---------------------------------------------------------
L1=1:30
L2=c(1:10,31:40,11:15)
L3=c(1:10,16:20,11:15)
input=list(L1,L2,L3)
space=list(1:40,1:40,1:40)

## ---------------------------------------------------------
outBorda=Borda(input,space,k=50) 
# "space" is explicitly specified; underlying space-dependent

## ---------------------------------------------------------
outBorda1=Borda(input) 
#"space" is not specified; all lists are assumed to come from the common space (objects 1-40)


## ---------------------------------------------------------
outBorda2=Borda(input,space=input) 
# "space = input" indicates that this is the top-k space

## ---------------------------------------------------------
sum(outBorda$Scores-outBorda1$Scores)
sum(outBorda$Scores-outBorda2$Scores)

## ---------------------------------------------------------
as.list(outBorda$TopK)

## ----borda, fig.cap='Plotting Borda\'s scores', fig.subcap=c('Underlying space-dependent analysis','Top-k space analysis'), out.width='.49\\linewidth', fig.pos='ht!', fig.align='center'----
Borda.plot(outBorda, k=40) # plot scores from underlying space-dependent analysis
Borda.plot(outBorda2) # plot scores from top-k space analysis

## ---------------------------------------------------------
outMC=MC(input,space) 
# "space" is explicitly specified; underlying space-dependent

## ---------------------------------------------------------
outMCa=MC(input,k=30) 
# "space" is not specified, so it is the same as common space (1-40)

## ---------------------------------------------------------
outMCb=MC(input,space=input) 
# "space = input" indicates that this is the top-k space

## ---------------------------------------------------------
sum(outMC$MC2.Prob-outMCa$MC2.Prob)
sum(abs(outMC$MC2.TopK-outMCb$MC2.TopK))

## ---------------------------------------------------------
list(outMC$MC1.TopK, outMC$MC2.TopK, outMC$MC3.TopK)

## ----equil, include=TRUE, fig.cap='Equilibrium probabilities', fig.pos='ht!', out.width='9cm', out.height='9cm', fig.align='center'----
MC.plot(outMC)

## ---------------------------------------------------------
set.seed(12345)
outCEMC=CEMC(input, space, N=8000)  
# "space" is explicitly specified; underlying space-dependent

## ---------------------------------------------------------
list(outCEMC$TopK)
outCEMC$ProbMatrix[1:5,1:5] 

## ---------------------------------------------------------
outCEMC$input.par

## ---------------------------------------------------------
KendallMLists(input,space, outBorda$TopK[,1])
all.aggregates=list(outBorda$TopK[,1],outBorda$TopK[,2],outBorda$TopK[,3],outBorda$TopK[,4],outMC$MC1.TopK,outMC$MC2.TopK,outMC$MC3.TopK,outCEMC$TopK)

## ----plotKendall, include=TRUE, tidy=TRUE, echo=TRUE, out.width='8cm', out.height='8cm', fig.pos='H',fig.cap="Comparison of the modified Kendall distances across several algorithms", fig.align='center'----
Kendall.plot(input,all.aggregates,space,algorithm=c("ARM","MED","GEO","L2N","MC1","MC2","MC3","CEMC"))

## ---------------------------------------------------------
data(breast)

## ---------------------------------------------------------
res = j0.multi(breast, d=6, v=10)
res$maxK

## ---------------------------------------------------------
l1=as.character(breast[1:res$maxK,1])
l2=as.character(breast[1:res$maxK,2])
l3=as.character(breast[1:res$maxK,3])
input=list(l1,l2,l3)
common=unique(unlist(input))
space=list(common,common,common)
outBorda=Borda(input,space)
outMC=MC(input,space)
outCEMC=CEMC(input,space,N=2000)

## ---------------------------------------------------------
agg=list(ARM=outBorda$TopK[,1],MED=outBorda$TopK[,2],GEO=outBorda$TopK[,3],L2N=outBorda$TopK[,4],MC1=outMC$MC1.TopK,MC2=outMC$MC2.TopK,MC3=outMC$MC3.TopK,CEMC=outCEMC$TopK)
head(do.call(cbind, agg))

## ----aggplot, out.width='8cm', out.height='8cm', fig.pos="ht!", fig.cap="Results from different algorithms for the breast cancer example data", fig.align='center'----
Kendall.plot(input,agg,space,algorithm=c("ARM","MED","GEO","L2N","MC1","MC2","MC3","CEMC"))

## ---------------------------------------------------------
deltaplot.dir =  paste0(tempdir(), "/deltaplot")
dir.create(deltaplot.dir, showWarnings = FALSE)
subplot.dir = paste0(tempdir(), "/subplot")
dir.create(subplot.dir, showWarnings = FALSE)

## ----fig.show='hide'--------------------------------------
a=deltaplot(breast, deltas = seq(0,300, by=5), directory=deltaplot.dir) 

## ----fig.show='hide'--------------------------------------
a=deltaplot(breast, deltas = 1:50, subset.lists=200, subplot = TRUE, perc.subplot=50, directory=subplot.dir)

## ---------------------------------------------------------
res = j0.multi(breast, d=6, v=10)
sapply(res, head)

## ----eval=FALSE-------------------------------------------
#  TopKListsGUI(breast)

## ---------------------------------------------------------
sessionInfo()

