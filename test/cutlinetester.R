################################################################
# The Original Cutline Function                                #
# Not used for now, but we should discuss this (weights)       #
# Adapted to add cutlines for reuse                            #
# Retained primarily for comparison purposes in the code below #
################################################################

add.cutline2 <- function (cutData, main.title="W-NOMINATE Cutting Lines", d1.title="First Dimension", d2.title="Second Dimension", weight1=.5,weight2=.5) 
{

WEIGHT2<-WEIGHT<-weight1/weight2

  DL1 <- cutData[2]    #spread
  ZM1 <- cutData[1]
  DL2 <- cutData[4]
  ZM2 <- cutData[3]
#
#
YEA1 <- ZM1-DL1
YEA2W <- (ZM2-DL2)*WEIGHT
NAY1 <- ZM1+DL1
NAY2W <- (ZM2+DL2)*WEIGHT
#
#  GET NORMAL VECTOR IN WEIGHTED METRIC -- SET YEA TO ORIGIN
#
A1 <- NAY1 - YEA1
A2 <- NAY2W - YEA2W
#
#  SETTING TO UNIT LENGTH VECTOR
#
ALENGTH <- sqrt(A1*A1+A2*A2)
N1W <- A1/ALENGTH
N2W <- A2/ALENGTH
if (N1W < 0){
  N1W <- -N1W
  N2W <- -N2W
}
#
ws <- N1W*ZM1 + N2W*ZM2*WEIGHT
#
xws <- ws*N1W
yws <- ws*N2W
#
N1WORTHOG <- N2W
N2WORTHOG <- ((-N1W)/WEIGHT)*WEIGHT2
#
BLENGTH <- sqrt(N1WORTHOG*N1WORTHOG + N2WORTHOG*N2WORTHOG)
N1NEW <- N1WORTHOG/BLENGTH
N2NEW <- N2WORTHOG/BLENGTH
#
N1 <- N2NEW
N2 <- -N1NEW
if (N1 < 0){
  N1 <- -N1
  N2 <- -N2
}
#
#
ws <- N1*ZM1 + N2*ZM2*WEIGHT2
#
#  Plot Cutting Line
#
xws <- ws*N1
yws <- ws*N2
segments(xws,yws,xws+N2,yws-N1,lwd=2,col="black")
segments(xws,yws,xws-N2,yws+N1,lwd=2,col="black")
#
#
#
}


###################################################
#  Code proving equivalency of the two methods ####
# Any RC can be selected via 'cutline'         ####
###################################################

par(mfrow=c(1,2))
plot(-1:1,-1:1)
cutline<-50
nomObject$result
add.cutline(c(nomObject$rollcalls[cutline,"midpoint1D"],
                    nomObject$rollcalls[cutline,"spread1D"],
                    nomObject$rollcalls[cutline,"midpoint1D"],
                    nomObject$rollcalls[cutline,"spread1D"]),
		    weight=1)
plot(-1:1,-1:1)
add.cutline2(c(nomObject$rollcalls[cutline,"midpoint1D"],
                    nomObject$rollcalls[cutline,"spread1D"],
                    nomObject$rollcalls[cutline,"midpoint1D"],
                    nomObject$rollcalls[cutline,"spread1D"]))
