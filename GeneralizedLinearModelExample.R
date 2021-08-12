# R - code to bring in the data and fit a GLM
AOI <- rep(c("Low","Medium","High"), times=3)
Terr <- rep(c(1,2,3), each=3)
Exposure <- c(7, 108,179,130,126,129,143,126,40)
LossLAE <- c(210.93,4458.05,10565.98,6206.12,8239.95,12063.68,8441.25,10188.7,4625.34)
Premium <- c(335.99,6479.87,14498.71,10399.79,12599.75,17414.65,14871.7,16379.68,7019.86)
Sampdata <- data.frame(AOI,Terr,Exposure,LossLAE,Premium)
Sampdata
Sampdata$AOI = relevel(Sampdata$AOI,ref="Medium")
Sampdata$Terr = relevel(factor(Sampdata$Terr),ref="2")
#  GLM FIT
gamma.model <- glm(LossLAE~AOI+Terr,offset=log(Exposure), 
                   data=Sampdata,family=Gamma(link="log"))
summary(gamma.model)
# EXTRACT COEFFICIENTS AND RELATIVITES
(tempoutput <- cbind(gamma.model$coefficients,exp(gamma.model$coefficients)) )
#  ALTERNATIVE APPROACH USING WEIGHTED REGRESSION
gamma.model1 <- glm(LossLAE/Exposure~AOI+Terr,weights=1/Exposure,
                    data=Sampdata,family=Gamma(link="log"))
summary(gamma.model1)

