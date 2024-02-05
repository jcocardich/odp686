####
setwd("D:\\Publicaciones\\202X ENSO ODP\\Dataset")
#install.packages("qpcR",dependencies=TRUE)
library(pls);	library(ade4);	library(vegan)
library(gclus);	library(ape);	library(ggrepel)
library(quantmod);	library(pracma);	library(TTR)
source("evplot.R")

tableres1=NULL
tableres2=NULL
tsTi=NULL
tsTiamp=NULL
tsTimax=NULL

xrf <- read.csv("segment_Ti.csv",sep=",", row.names=NULL,stringsAsFactors = FALSE)
seg <- colnames(xrf)[-1]
dfz <- xrf[,1]
dfti <- xrf[,-1]
sumxmean <- c(137448.8, 138747.3, 174377.2, 151095.5, 
	129952.5, 157958.8, 182143.5, 164214.7, 165691.4, 
	169095.2, 176070.1, 237608.9, 284502.5)

for (i in 1:length(seg))	{
df <- cbind(dfz, dfti[i])
df <- df[complete.cases(df),]
z <- df[,1]
ti <- df[,2]					
tic <- ti/sumxmean[i]

###MAXIMA AND MINIMA##################
###https://github.com/stas-g/findPeaks
find_peaks <- function (x, m = 3){
    shape <- diff(sign(diff(x, na.pad = FALSE)))
    pks <- sapply(which(shape < 0), FUN = function(i){
       z <- i - m + 1
       z <- ifelse(z > 0, z, 1)
       w <- i + m + 1
       w <- ifelse(w < length(x), w, length(x))
       if(all(x[c(z : i, (i + 2) : w)] <= x[i + 1])) return(i + 1) else return(numeric(0))
    })
     pks <- unlist(pks)
     pks
}

pe <- find_peaks(ti, m=4);		pe <- as.vector(pe)
allpe <- find_peaks(ti, m=1);		spe <- setdiff(allpe,pe)
						spe <- as.vector(spe)	
va <- find_peaks(-ti, m=4);		va <- as.vector(va)
p <- c(pe, va)	

d <- data.frame(id=ti)						
d$flag <- ifelse(ti %in% d$id[pe], 1, ifelse(ti %in% ti[va], -1, 0))
d$z <- z
dd <- subset(d, flag != 0)
rle <- rle(dd$flag)
myrle <- which(rle$lengths >= 2) 				
rlesum = cumsum(rle$lengths)					
ends = rlesum[myrle]						
newi = ifelse(myrle>1, myrle-1, 0)
starts = rlesum[newi] + 1
if (0 %in% newi) starts = c(1,starts)

aa=NULL
for(j in 1:length(starts)) {
	a <- dd[starts[j]:ends[j],]
	if(a$flag[1] == 1) {
	b <- a[which(a$id == min(a$id)),] } else {
	b <- a[which(a$id == max(a$id)),] }
	aa = rbind(aa,b)						
}

del <- which(as.numeric(row.names(dd)) %in% as.numeric(row.names(aa)))		
ddd1 <- dd[-del,]					
z2 <- ddd1$z					
pez <- ddd1[,c(1,3)][ddd1$flag == 1,]	
vaz <- ddd1[,c(1,3)][ddd1$flag == -1,]	
pez$tic <- tic[as.numeric(row.names(pez))]
vaz$tic <- tic[as.numeric(row.names(vaz))]

###ACCUMULATION RATE########################
sed=NULL
for(f in 1:length(pe)) {
dif<-z[pe[rev(1:length(pe))]][f] - z[pe[rev(1:length(pe))]][f+1]
sed<-c(dif,sed)
}
sed <- sed[!is.na(sed)]
mean(sed)
sd(sed)

###AMPLITUDE ESTIMATION#####################
if(ddd1$flag[1] == -1) 	{				
	ddd1 <- ddd1[-1,]	} else {
	ddd1 <- ddd1
	}
if(rev(ddd1$flag)[1] == 1) {				
	ddd1 <- head(ddd1,-1) }  else {
	ddd1 <- ddd1
	}
dd2 <- ddd1$id						

am=NULL							
for(g in 1:(length(dd2)/2)) {					
	b <- abs(dd2[g*2-1]-dd2[g*2])
	am <- c(am,b)
}
am <- am[!is.infinite(am)]			

###DATA CONVERSION TO MONTHLY RESOLUTION###
if(rev(spe)[1] > rev(pe)[1]){
	spen <- spe[which(spe < rev(pe)[1])] 
} else {
	spen <- spe
}								

yr <- 1:dim(pez)[1]					
iyr <- approx(x=pez$z,y=yr, xout=z)$y 		
ime <- approx(x=iyr, y=ti, xout= seq(1,rev(yr)[1],by=1/12))		
year <- ime$x
val <- ime$y
nime <- cbind(year,val)

nvaz=NULL
for(k in 1:(length(yr)-1)) {
	r <- nime[(k*12):(k*12-11),]
	nmin <- r[which(r[,2] == min(r[,2])), ]
	nvaz <- rbind(nvaz,nmin)
}

tsTi <- qpcR:::cbind.na(tsTi,nime)
timax <- pez$corrti  
tsTimax <- qpcR:::cbind.na(tsTimax,timax)

###TIME SERIES SMOOTHING##################
library(dplyr)
library(zoo)
mmed10 <- rollmedian(val, k=120, fill=NA, align='center')	
tis <- val - mmed10							

###EL NINO frequency######################
y <- length(yr)								
itis <- cbind(year,tis, rep(1:12, y-1))
colnames(itis)[3] <- "mes"
sd(tis, na.rm=T)							
Tismax.sd <- sd(itis[,2][itis[,3]==1], na.rm=T)			

sigma=0.7213344
ENf <- length( which(log(am) > mean(log(am)) + 1*sigma) ) / y	
xENf <- length( which(log(am) > mean(log(am)) + 2*sigma) ) / y	

}






