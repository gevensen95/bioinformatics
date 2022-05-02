library(phytools)
library(phangorn)
library(geiger)
library(RColorBrewer)
library(ggimage)
library(ggtree)

setwd("~/Desktop/Phylogenetics/16S_test") #change to 16S_test for a file that has DUI separate from 
  #modes of reproduction

b.tree <- read.tree(file = "RAxML_bestTree.16S_tree")
b.tree <- multi2di(b.tree)
b.tree <-root.phylo(b.tree, outgroup = "Argopecten_purpuratus", resolve.root = TRUE)
plot.phylo(b.tree, use.edge.length = FALSE, font = 2)

b.tree2 <- b.tree
b.tree2$edge.length<-NULL
b.tree3 <- b.tree
b.tree3$edge.length[150]<-2.5

bvData<-read.csv("sex_mec.csv")
bvData <- as.data.frame(bvData)

#Plot using ggtree
bvMatches <- match(b.tree$tip.label, bvData$Bivalve)
bvMatches <- bvMatches[!is.na(bvMatches)]
bvData2 <- bvData[bvMatches,]

p <- ggtree(b.tree2) %<+% bvData2 + xlim(-.1, 20) 
p2 <- p + geom_tiplab() +
  geom_tippoint(aes(color = mode, shape = DUI), size = 4) + 
  theme(legend.position = "right", legend.title = element_text(color = "black", size = 20),
        legend.text = element_text(color = "black", size = 14))
p2

p3 <- p + geom_tiplab() + 
  geom_tippoint(aes(color = Subclass), size = 4) + geom_text(aes(label=node), hjust=-.3) +
  theme(legend.position = "right", legend.title = element_text(color = "black", size = 20),
        legend.text = element_text(color = "black", size = 14))
p3

#Plot using phytools
data<-as.matrix(read.csv("sex_mec.csv",row.names=1))[,1]
dotTree(b.tree3, data, ftype="b",legend = TRUE, fsize = 0.6, length = 10, prompt = TRUE)

data2<-read.csv("sex_mec.csv",row.names=1)
data2<-setNames(data2[,1],rownames(data2))
cols<-setNames(palette()[1:2],levels(data2))
bivalve.tree<-make.simmap(drop.tip(b.tree3,
                                    setdiff(b.tree$tip.label,names(data2))),data2,
                           model="SYM")

plot(bivalve.tree,cols,fsize=0.8,lwd=3,ftype="i")
add.simmap.legend(prompt = TRUE, colors = cols)

#Now we will utilize stochastic mapping
bivalve.trees<-make.simmap(drop.tip(b.tree3,
                                   setdiff(b.tree$tip.label,names(data2))),data2,
                          model="SYM",nsim=100)
obj <- summary(bivalve.trees, plot = FALSE)
plot(obj,colors=cols,fsize=0.8,cex=c(0.5,0.3))
add.simmap.legend(colors=cols, prompt = TRUE)

#We will do the same as above, but for DUI
data3<-read.csv("sex_mec.csv",row.names=1)
data3<-setNames(data3[,2],rownames(data3))
cols<-setNames(palette()[1:2],levels(data3))
bivalve.tree3<-make.simmap(drop.tip(b.tree3,
                                   setdiff(b.tree3$tip.label,names(data3))),data3,
                          model="SYM")

plot(bivalve.tree3,cols,fsize=0.8,lwd=3,ftype="i")
add.simmap.legend(prompt = TRUE, colors = cols)

  #Stochastic Mapping
bivalve.trees3<-make.simmap(drop.tip(b.tree3,
                                    setdiff(b.tree3$tip.label,names(data3))),data3,
                           model="SYM",nsim=100)
obj.DUI <- summary(bivalve.trees3, plot = FALSE)
plot(obj.DUI,colors=cols,fsize=0.8,cex=c(0.5,0.3))
add.simmap.legend(colors=cols, prompt = TRUE)

#Plot Discrete Characters
modes <- read.csv("sex_mec.csv",row.names=1)
s.modes<-as.factor(setNames(modes[,2],rownames(modes)))
b.trees <- make.simmap(b.tree3, s.modes, nsim = 100)
obj2<-densityMap(b.trees,s.modes,plot=FALSE)

#Model discrete 4 characters
b.Data<-read.csv("sex_mec.csv", row.names = 1)
mec <-b.Data$sex_mec
chk<-name.check(b.tree, b.Data)
summary(chk)
sex_mec<-setNames(as.factor(b.Data[,"sex_mec"]), rownames(b.Data))
head(sex_mec)

#Equal Rates Model
fitER<-fitDiscrete(b.tree, sex_mec, model="ER")
print(fitER,digits=3)
plot(fitER)

#Symmetric Transition Model
fitSYM<-fitDiscrete(b.tree, sex_mec,model="SYM")
print(fitSYM,digits=3)
plot(fitSYM, show.zeros=FALSE)

#All Rates Different Model
fitARD<-fitDiscrete(b.tree,sex_mec,model="ARD")
print(fitARD,digits=3)
plot(fitARD,show.zeros=FALSE)

#Ordered Model
ordered.model<-matrix(c(
  0,1,4,0,
  2,0,0,4,
  0,0,0,0,
  0,5,0,0),4,4,byrow=TRUE,
  dimnames=list(levels(sex_mec),levels(sex_mec)))
ordered.model
fitOrdered<-fitDiscrete(b.tree,sex_mec,model=ordered.model,surpressWarnings=TRUE)
print(fitOrdered,digits=3)
plot(fitOrdered,show.zeros=FALSE,signif=4)

#Directional Model
directional.model<-matrix(c(
  0,1,2,0,
  0,0,0,1,
  0,0,0,0,
  0,0,0,0),4,4,byrow=TRUE,
  dimnames=list(levels(sex_mec),levels(sex_mec)))
directional.model
fitDirectional<-fitDiscrete(b.tree,sex_mec,model=directional.model,surpressWarnings=TRUE)
plot(fitDirectional,show.zeros=FALSE,signif=4)

#Compare Models
aic<-setNames(c(AIC(fitER),AIC(fitDirectional),AIC(fitOrdered),AIC(fitSYM),AIC(fitARD)),
              c("ER","Directional","Ordered","SYM","ARD"))
aic

round(data.frame(
  k=c(fitER$opt$k,fitDirectional$opt$k,
      fitOrdered$opt$k,fitSYM$opt$k,fitARD$opt$k),
  logL=c(logLik(fitER),logLik(fitDirectional),
         logLik(fitOrdered),logLik(fitSYM),logLik(fitARD)),
  AIC=aic,Akaike.w=as.vector(aic.w(aic))),3)

#Reconstructing Ancestral States - G vs.H
b.Data2<-read.csv("sex_mec.csv",row.names = 1)
s.modes2 <- setNames(b.Data2$sex_mec, rownames(b.Data))
fit.s<-fastAnc(b.tree3,s.modes2,vars=TRUE,CI=TRUE)

continuousTraitMap<-contMap(b.tree3,s.modes2,plot=FALSE,
                            lims=c(0,2.5))
continuousTraitMap<-setMap(continuousTraitMap,
                           c("white","grey","black"))
plot(continuousTraitMap,
     sig=2,fsize=c(0.4,0.9),lwd=c(2,4), ftype = "i",
     leg.txt="log(s.modes2)")
par(fg="darkgrey")
errorbar.contMap(continuousTraitMap,lwd=4)

#Reconstructing Ancestral States - DUI
DUI.Data<-read.csv("sex_mec.csv",row.names = 1)
DUI <- setNames(DUI.Data$Mito, rownames(DUI.Data))
fit.DUI<-fastAnc(b.tree3,DUI,vars=TRUE,CI=TRUE)

continuousTraitMap<-contMap(b.tree3,DUI,plot=FALSE,
                            lims=c(0.5,2.6))
continuousTraitMap<-setMap(continuousTraitMap,
                           c("white","grey","black"))
plot(continuousTraitMap,
     sig=2,fsize=c(0.4,0.9),lwd=c(2,4), ftype = "i",
     leg.txt="log(DUI)")
par(fg="darkgrey")
errorbar.contMap(continuousTraitMap,lwd=4)
