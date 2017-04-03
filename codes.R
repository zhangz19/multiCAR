
require(maptools)
library(maps)
library(fields)
library(lattice)
library(ggplot2) #fortify
library(ape)
library(R.matlab)
library(MBA)

## data process
#rm(list=ls())
#dat <- X <- cds <- W <- list(3)
#datafile <- c('health_ind_County.csv', 'health_ind_Census tract.csv', 'health_ind_Block group.csv')
#for(i in 1:3) {
# dat[[i]] <- read.csv(datafile[i])
# X[[i]] <- read.csv(paste('X',i,'.csv',sep=''),na='#DIV/0!')
# cds[[i]] <- read.csv(paste('coords_drain',i,'.txt',sep=''),sep='\t',h=T)
# W[[i]] <- read.csv(paste('W',i,'.txt',sep=''),sep='\t',h=F)
#}
#
#X0 <- X[[1]]; dat0 <- dat[[1]]; cds0 <- cds[[1]]; W0 <- W[[1]]
#X0 <- X0[-4,] #exclude   26037    Clinton County
#X0 <- X0[which(X0$Geo_FIPS %in% dat0$GEOID10),]
#W0 <- W0[which(W0[,1] %in% dat0$GEOID10),-2]; names(W0) <- c("GEOID10",paste('neighbor',1:(ncol(W0)-1),sep=''))
#a <- merge(cds0,dat0, by='GEOID10')
#a <- merge(a, X0, by.x='GEOID10', by.y='Geo_FIPS')
#a <- a[,-c(5:6,29:32)]; a <- a[,c(1,9,2:3,5:8,4,10:ncol(a))]
#a <- merge(a, W0, by='GEOID10')
#write.csv(a, file='fulldat1.csv',row.names=F)
#
#
#X0 <- X[[2]]; dat0 <- dat[[2]]; cds0 <- cds[[2]]; W0 <- W[[2]]
#X0 <- X0[which(X0$Geo_FIPS %in% dat0$GEOID10),]
#cds0 <- cds0[which(cds0$GEOID10 %in% dat0$GEOID10),]
#ID0 <- W0[,1]; W0[,1] <- ID0[W0[,2]]
#ind0 <- which(ID0 %in% dat0$GEOID10)
#W0 <- W0[which(W0[,1] %in% dat0$GEOID10),-2]; names(W0) <- c("GEOID10",paste('neighbor',1:(ncol(W0)-1),sep=''))
#for(i in 1:nrow(W0)){
# for(j in 2:ncol(W0)){ 
#   if(!(W0[i,j]%in%ind0)) W0[i,j] <- 0
#   else W0[i,j] <- ID0[W0[i,j]]
# }
#} 
## write.csv(W0, file='tmp.csv',row.names=F)
#a <- merge(cds0,dat0, by='GEOID10')
#a <- merge(a, X0, by.x='GEOID10', by.y='Geo_FIPS')
#a <- a[,-c(5:6,29:32)]; a <- a[,c(1,9,2:3,5:8,4,10:ncol(a))]
#a <- merge(a, W0, by='GEOID10')
#write.csv(a, file='fulldat2.csv',row.names=F)
#
#
#X0 <- X[[3]]; dat0 <- dat[[3]]; cds0 <- cds[[3]]; W0 <- W[[3]]
#dat0 <- dat0[-1,]
#X0 <- X0[which(X0$Geo_FIPS %in% dat0$GEOID10),]
#cds0 <- cds0[which(cds0$GEOID10 %in% dat0$GEOID10),]
#ID0 <- W0[,2]; W0[,1] <- ID0[W0[,1]]
#ind0 <- which(ID0 %in% dat0$GEOID10)
#W0 <- W0[which(W0[,1] %in% dat0$GEOID10),-2]; names(W0) <- c("GEOID10",paste('neighbor',1:(ncol(W0)-1),sep=''))
#for(i in 1:nrow(W0)){
# for(j in 2:ncol(W0)){ 
#   if(!(W0[i,j]%in%ind0)) W0[i,j] <- 0
#   else W0[i,j] <- ID0[W0[i,j]]
# }
#} 
#a <- merge(cds0,dat0, by='GEOID10')
#a <- merge(a, X0, by.x='GEOID10', by.y='Geo_FIPS')
#a <- a[,-c(5:6,29:32)]; a <- a[,c(1,9,2:3,5:8,4,10:ncol(a))]
#a <- merge(a, W0, by='GEOID10')
#write.csv(a, file='fulldat3.csv',row.names=F)
#
##for(i in 1:3) print(dim(dat[[i]]))
##for(i in 1:3) print(dim(X[[i]]))
#
##a <- merge(dat0, X0, by.x='GEOID10', by.y='Geo_FIPS')
##for(i in 1:3) {write.csv(X[[i]], file=paste('sX',i,'.csv',sep=''), row.names=F)}



# pdf('BoundaryPlots.pdf')
Allshcoords <- list(3)
for(ids in 1:3){
dirs <- './Census_levels/'
cdnams <- c('County_UTM17','Census_UTM17','Block_UTM17')
shape <- readShapePoly(paste(dirs, cdnams[ids],'.shp', sep=""))
shcoords <- fortify(shape); shcoords <- shcoords[,c(1,2)]
Allshcoords[[ids]] <- shcoords
# cds <- read.csv(paste('coords_drain',ids,'.txt',sep=''),sep='\t',h=T)
plot(shcoords)
# points(cds[,c(3,2)], col='red', pch=16)
}
save(file='Allshcoords',Allshcoords)
# shape <- readShapePoly(paste(dirs, 'saginaw/saginaw.shp', sep=""))
# shcoords <- fortify(shape); shcoords <- shcoords[,c(1,2)]
# plot(shcoords)
# dev.off()




## check the neighborhood
#pdf('CheckNeighbor.pdf',h=8,w=11)
#ids <- 2
#plot(Allshcoords[[ids]], col='gray50',cex=.5) 
#a <- read.csv(file=paste('fulldat',ids,'.csv',sep=''))
#points(a[,c('Longitude','Latitude')],col='gray20')
#for(checkid in c(1,150,220,303,330)){
#allid <- a[checkid,c(29:ncol(a))]; allid <- allid[allid!=0]
#points(a[c(checkid,which(a[,1]%in%allid)),c('Longitude','Latitude')],col=c(2,rep(3,length(allid))),pch=16)
#}
#dev.off()



## -- rerun this to exclude PCT_SE_T025_003, PCT_SE_T011_003, PCT_SE_T011_004
#for(ids in 1:3){
# a <- read.csv(file=paste('fulldat',ids,'.csv',sep=''))
# a <- a[,-which(names(a) %in% c('PCT_SE_T025_003','PCT_SE_T011_003','PCT_SE_T011_004'))]
# write.csv(file=paste('fulldat',ids,'.csv',sep=''),a,row.names=F)
#}

# ----------------------------------------------------------------------------------------------------------

fnams <- c('1-Fulldata_County.csv','2-Fulldata_Census.csv','3-Fulldata_Block.csv')

# check the symmetry of the neighborhood matrix
ids <- 3
# a <- read.csv(file=paste('fulldat',ids,'.csv',sep=''))
a <- read.csv(file=fnams[ids])

W1 <- a[,c(26:ncol(a))]; n <- nrow(a)
W2 <- matrix(0,n,n); for(i in 1:n) W2[i,which(a[,1]%in%W1[i,])] <- 1
all(W2==t(W2))


# investigate the links three levels
dat <- list(3); for(ids in 1:3) dat[[ids]] <- read.csv(file=fnams[ids])
sapply(dat,dim)[1,]           
grp2_1 <- character(nrow(dat[[2]]))
for(i in 1:nrow(dat[[1]])) grp2_1[grepl(pattern=as.character(dat[[1]][i,1]), x=as.character(dat[[2]][,1]))] <- dat[[1]][i,1]

grp3_2 <- character(nrow(dat[[3]]))
for(i in 1:nrow(dat[[2]])) grp3_2[grepl(pattern=as.character(dat[[2]][i,1]), x=as.character(dat[[3]][,1]))] <- dat[[2]][i,1]



# a <- cbind(grp3_2, dat[[3]][,1]); print(a[sample(1:nrow(a),size=20),])
head(cbind(aggregate(dat[[3]][,c("FIBI","EPT")], by=list(grp3_2), mean), dat[[2]][,c("FIBI","EPT")]),20)#,"HBI","IBI"  # almost match, as I expected

head(cbind(aggregate(dat[[2]][,c("FIBI","EPT")], by=list(grp2_1), mean), dat[[1]][,c("FIBI","EPT")]),20)



# check variables
# produce the scatterplot matrix
.panelhist <- function(x, ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.5) )
    h <- hist(x, plot = FALSE, breaks=8)
    breaks <- h$breaks; nB <- length(breaks)
    y <- h$counts; y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y, col="cyan4", ...)
    if(length(unique(x)) == 2){
     text(breaks[c(1,nB-1)]+.05, y[c(1,nB-1)]+.2, c("F","M"))
    }
}


.panelcor <- function(x, y, digits=2, prefix="", cex.cor)
{
    if(length(unique(y)) == 2) return(NULL)
    else{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    # r <- abs(cor(x, y, use="complete.obs", method='spearman'))
    tmp <- cor.test(x,y, exact=F, method="spearman")
    r <- tmp$estimate; pvalue <- tmp$p.value
    txt <- format(c(r, 0.123456789), digits=digits)[1]
    txt <- paste(prefix, txt, sep=""); tcol <- 1
    if(r >= 0.7) tcol <- 'red'  #(tmp$p.value < 0.05)
    else if(r <= -0.7) tcol <- 'blue2'
    text(0.5, 0.5, txt, cex = .8/strwidth(txt), col = tcol)
    }
}  



source('getScatterPlot.R')

# pdf('scatterplotmatrix.pdf',w=13,h=9,pointsize=23)
# for(ids in 1:3){
#  #ids <- 1 #change 1 to 2 and 3
#  a <- read.csv(file=fnams[ids])
#  mat <- a[,c(5:8, 9:25)]
#  #print(summary(lm(cbind(family_IBI,EPT_Taxa,HBI,IBI) ~ ., data=mat)))
#  #X11(w=12,h=9); 
#  par(pch=1,xaxt='n',yaxt='n',mar=c(1,1,0,0)+.1,col='blue4') 
#  pairs(mat, #main='scatterplot matrix with pearson correlation coefficients',
#         diag.panel=.panelhist, 
#         upper.panel=.panelcor, lower.panel=panel.smooth,
#         cex=.5, cex.labels=.21, gap=0.2, font.labels=3)
# }
# dev.off()



# investigate spatial distribution
cexs <- c(1,.5,.3); epss <- c(.3,.012,.01)
pdf('SpatialPlot.pdf',h=8,w=11)

for(ids in 1:3){
eps <- epss[ids] # ids=2: 0.012
a <- read.csv(file=fnams[ids])
load('Allshcoords')
shcoords <- Allshcoords[[ids]]
coords <- a[,c('Longitude','Latitude')]
y <- a[,'FIBI']

mybox <- c(range(shcoords[,1])+c(-1,1)*.0, range(shcoords[,2])+c(-1,1)*.0)
surf <- mba.surf(cbind(coords,y), no.X=400, no.Y=400, extend=T,b.box=mybox)$xyz.est
for(i in 1:length(surf$y)){   # row, lat
  tmp <- range(shcoords[which(abs(shcoords[,2]-surf$y[i])<eps),1])
  if(length(tmp)) surf$z[which(surf$x<tmp[1]|surf$x>tmp[2]),i] <- NA
  if(ids > 1){ if(surf$y[i]>43.6 && surf$y[i]<43.8){
    tmp <- range(shcoords[which((abs(shcoords[,2]-surf$y[i])<eps)&(shcoords[,1]< -83.9)),1])
    if(length(tmp)) surf$z[which((surf$x<tmp[1]|surf$x>tmp[2])&(surf$x< -83.9)),i] <- NA
  }  }
}
for(i in 1:length(surf$x)){   # col, lon
  tmp <- range(shcoords[which(abs(shcoords[,1]-surf$x[i])<eps),2])
  if(length(tmp)) surf$z[i,which(surf$y<tmp[1]|surf$y>tmp[2])] <- NA
  if(ids > 1) {if(surf$x[i]> -84.17 && surf$x[i]< -83.7){
    tmp <- range(shcoords[which((abs(shcoords[,1]-surf$x[i])<eps)&(shcoords[,2]< 43.7)),2])
    if(length(tmp)) surf$z[i,which((surf$y<tmp[1]|surf$y>tmp[2])&(surf$y< 43.7))] <- NA
  } }
}
isNA <- is.na(surf$z)

spatPlot <- function(y, mains=''){ 
surf <- mba.surf(cbind(coords,y), no.X=400, no.Y=400, extend=T,b.box=mybox)$xyz.est
surf$z[which(isNA)] <- NA
par(mar=c(0,0,2,0))
image.plot(surf, xlim=range(c(shcoords[,1],coords[,1]),na.rm=T)+c(-1,1)*.05, 
 main=mains, ylim=range(c(shcoords[,2],coords[,2]),na.rm=T)+c(-1,1)*.05, axes=F)
points(shcoords, pch=16, cex=.24,col='gray70')
points(coords, col='gray40', pch=16, cex=cexs[ids]) 
}
# spatPlot(y)

par(mfrow=c(2,2)); for(j in 5:8) spatPlot(a[,j], names(a)[j])
# summary(a[,5:8])
}

dev.off()




## convert data for MATLAB
for(ids in 1:3){
  dat <- read.csv(file=fnams[ids])
  Y <- dat[,5:8]
  X <- dat[,9:25]    # no intercept
  W0 <- dat[,26:ncol(dat)]; n <- nrow(W0)
  W <- matrix(0,n,n); for(i in 1:n) W[i,which(dat[,1]%in%W0[i,])] <- 1
  # exclude missing cases
  if(ids > 1){ misind <- integer(); for(i in 1:nrow(X)) if(any(is.na(X[i,]))) misind <- c(misind,i)
  X <- X[-misind,];  Y <- Y[-misind,];  W <- W[-misind,-misind];
  }
  writeMat(paste('fulldat',ids,'.mat',sep=''),Y=as.matrix(Y),X=as.matrix(X),W=W)
}



## read output from MATLAB and save it as Excel sheet
out <- readMat('final.mat')

load('data4scatterPlot') #dats, Xnams
Xnams <- Xnams[-c(1:4)]
p <- length(Xnams)
mat <- as.data.frame(out$optpara)
row.names(mat)[1:(p+1)] <- c('Intercept',Xnams) 
write.csv(file='tmp.csv',mat)

for(ids in 1:3){
  dat <- read.csv(file=fnams[ids])
  if(ids > 1){ misind <- integer(); for(i in 1:nrow(X)) if(any(is.na(X[i,]))) misind <- c(misind,i) }
  mat <- dats[[ids]][,1:4]; if(length(misind)) mat <- mat[-misind, ]
  mat <- cbind(mat, out$yhats[[ids]])
  names(mat)[5:8] <- paste('predicted',names(mat)[1:4])
  mat <- mat[,c(1,5,2,6,3,7,4,8)]
  write.csv(file=paste('yMat_level',ids,'.csv', sep=''), mat, row.names=F)
}


out <- readMat('dicMat.mat')
write.csv(file='dicMat.csv', out$dicMat, row.names=F)



## for part-1 analysis
out <- readMat('outAll.mat')

row.names(out$Allmat) <- rep(c('Intercept',Xnams,'tau^2','gamma'),3)
write.csv(file='modeloutput.csv', out$Allmat)

dics <- matrix(0,nrow(out$Alldic), ncol(out$Alldic)*4)
for(i in 1:nrow(dics)) dics[i,] <- rep(out$Alldic[i,], each=4)
write.csv(file='modeldic.csv', dics, row.names=F)  # paste in modeloutput.xlxs

dics[,seq(1,ncol(dics),by=4)]



