
prefix <- './'  # folder you want to save images

# Xnams <- c('FIBI', 'EPT', 'HBI', 'IBI', 'S1', 'S2', 'S3', 'S4', 'S5', 'S6', 'S7', 'S8', 'S9', 'S10', 'S11', 'S12', 'S13', 'S14', 'S15', 'S16', 'S17')

#pdf('scatterplotmatrix.pdf',w=13,h=9,pointsize=23,paper='special')
  

# fnams <- c('1-Fulldata_County.csv','2-Fulldata_Census.csv','3-Fulldata_Block.csv')
# dats <- list(3)
# for(ids in 1:3){
#   a <- read.csv(file=fnams[ids])
#   dats[[ids]] <- a[,c(5:8, 9:25)]
# }
# Xnams <- names(dats[[ids]]) # assume all the same
# save(file='data4scatterPlot',dats)


load('data4scatterPlot')
p <- length(Xnams)

for(ids in 1:3){
 #ids <- 1 #change 1 to 2 and 3
 #a <- read.csv(file=paste('fulldat',ids,'.csv',sep=''))
 #mat <- a[,c(5:8, 9:25)]
 #dats[[ids]] <- mat
          
   
 fac <- 1 # make fac smaller, to rescale figure to be smaller   
 tiff(paste(prefix,'scatterplotmatrix_',ids,'.tif',sep=''),height=9*fac, width=13*fac, units='in', pointsize=23*fac, res=600)
 
 par(pch=1,xaxt='s',yaxt='s',mar=c(0,0,0,0),col='black',mgp=c(1.1,.1,0), tck=-0.1,cex.axis=.61) 
 
 source('pairs2.R')
 pairs2(dats[[ids]], #main='scatterplot matrix with pearson correlation coefficients',
        diag.panel=.panelhist, 
        upper.panel=.panelcor, 
        lower.panel=.panelmain, #panel.lmline, #panel.smooth,
        cex.labels=.1, gap=0, font.labels=3, labels=rep(NA,p))
        
  addAxisLabels()
  
  dev.off()      
}
#dev.off()
