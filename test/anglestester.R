# house_1_108_cutlines_new.r -- Does Plot of Both House/Senate Cuttline Line
#                           Angles for any Congress From DW-NOMINATE Coordinates
#
#  &&&&&&&&&&&&&&&&&&&&&&&&
ncong <- 83
#  &&&&&&&&&&&&&&&&&&&&&&&&
#
#
#  Dimension Weights for 1 - 108 Scaling:  House 0.3463, Senate 0.375
#
WEIGHT <- 0.375
#WEIGHT <- 1.0
#
#
#  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#  Read Roll Call Coordinates
#  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
rcy.file <- "c:/atime2005/SC01108A1.DAT"
#
# Standard fields and their widths -- DW-NOMINATE Roll Call Parameters
#
rcy.fields <- c("cong","icount","dl1","zml1","dl2","zml2") 
rcy.fieldWidths <- c(3,5,7,7,7,7)
#    
# Input Roll Call Coordinates
#
TTT <- read.fwf(file=rcy.file,widths=rcy.fieldWidths,as.is=TRUE,col.names=rcy.fields)
dim(TTT)
#
#
TRC <- TTT[TTT[,1]==ncong,]
#
nrowrc <- length(TRC[,1])
ncolrc <- length(TRC[1,])
#
TANGLE <- rep(0,nrowrc*8)
dim(TANGLE) <- c(nrowrc,8)
#
i <- 0
j <- 0
while (i < nrowrc) {
  i <- i + 1
  DL1 <- TRC[i,3]
  ZM1 <- TRC[i,4]
  DL2 <- TRC[i,5]
  ZM2 <- TRC[i,6]
  circleconstraint <- ZM1**2 + ZM2**2
  if ((abs(DL1) > 0.0 | abs(DL2) > 0.0) & circleconstraint < .95){
      j <- j + 1
#
      veclength <- sqrt(DL1*DL1+WEIGHT*WEIGHT*DL2*DL2)
#  normal vector in weighted metric
      N1W <- DL1/veclength
      N2W <- (DL2*WEIGHT)/veclength
#  cutting plane vector
      CUTVECTOR1 <- N2W
      CUTVECTOR2 <- -N1W
#
     if (CUTVECTOR2 < 0.0){
          CUTVECTOR1 <- -N2W
          CUTVECTOR2 <- N1W
      }
      TANGLE[j,1] <- i
      TANGLE[j,2] <- j
      TANGLE[j,3] <- N1W
      TANGLE[j,4] <- N2W
      TANGLE[j,5] <- CUTVECTOR1
      TANGLE[j,6] <- CUTVECTOR2
      TANGLE[j,7] <- atan2(CUTVECTOR2,CUTVECTOR1)
      TANGLE[j,8] <- TANGLE[j,7]*(180.0/3.1415926536)
#
  }
#
}
#
# Create Vector for Smoothed Histogram
#
cutting.angles <- TANGLE[TANGLE[,8] > 0.00,8]
cuttingdens <- density(cutting.angles)
#
ymax <- max(cuttingdens$y)
ymax <- 1.1*ymax
#
sucker <- hist(cutting.angles,
br=c(0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180),
freq=FALSE,
#freq=NULL,
border="white",lwd=2,
main="",
xlab="",
ylab="",
ylim=c(0,.02),
axes=FALSE,
col="gray60",font=2)
#
year <- c(0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180)
axis(1,year[1:19],
   labels=c(
  '0',
  '10',
  '20',
  '30',
  '40',
  '50',
  '60',
  '70',
  '80',
  '90',
  '100',
  '110',
  '120',
  '130',
  '140',
  '150',
  '160',
  '170',
  '180'),font=2)
tucker <- round(((sucker$counts)/sum(sucker$counts))*100.0,2)
#axis(2,labels=tucker[1:18],font=2)
#axis(2,labels=c('0','3.88','7.75','11.63','15.50'),font=2)
#axis(2,font=2)
#axis(2,labels=c('0','5','10','15','20'),font=2)
axis(2,at=seq(.000,.020,.005),as.character(seq(0,20,by=5)),font=2)
#
# Main title
mtext("83rd Senate 1953-54\nCutting Line Angles",side=3,line=1.00,cex=1.75,font=2)
# x-axis title
mtext("Cutting Line Angle",side=1,line=2.75,cex=1.2)
# y-axis title
mtext("Percent",side=2,line=2.5,cex=1.2)
#
text(120,.0075,"Conservative\nCoalition",font=2,cex=1.2)
text(55,.014,"Party",font=2,cex=1.2)
lines(c(0,180),c(.005,.005),lty=1)
lines(c(0,180),c(.010,.010),lty=1)
lines(c(0,180),c(.015,.015),lty=1)
lines(c(0,180),c(.020,.020),lty=1)




#Code proving the equivalency of Keith's
#angle plots with the ones used here

rcy.file <- "c:/atime2005/SC01108A1.DAT"
rcy.fields <- c("cong","icount","spread1D","midpoint1D","spread2D","midpoint2D") 
rcy.fieldWidths <- c(3,5,7,7,7,7)
rollcalls <- read.fwf(file=rcy.file,widths=rcy.fieldWidths,as.is=TRUE,col.names=rcy.fields)
nomObject<-list(rollcalls=rollcalls[rollcalls$cong==83,],dimensions=2,weights=c(1,0.375))
class(nomObject)<-c("nomObject")
plot.angles(nomObject)

#The result is the exact same as the one from above
