#MTL_Cox for Personalized Prediction
library(dplyr)
data0 <- fread("~/zhangshuaijie/UKB/disease_data_set/20211019experiment/MTL_data/CHD_unionVars_data_std.csv")
data0<-data0[,-1]
beta <- read_xlsx("~/zhangshuaijie/UKB/disease_data_set/20211019experiment/MTL_data/20211204MTL_beta.xlsx",sheet = "BETA")
beta <- beta[which(beta$dis==1),]
var_number<- beta$Num
var_number <- as.vector(var_number+3)

risk_quantile <- as.data.frame(matrix(NA,4,3))
MTL_beta <- as.matrix(beta$beta)
data <- data0[order(data0$survtime),]
dt <- as.numeric(data$survtime)
MTL_HH <- as.data.frame(data$eid)

m <- 1
for (t in c(1,3,5,7)) {
  k <- length(which(dt < as.vector(t %*% 365)))
  MTL_eXB <- exp(data.matrix(data[,-c(1:3)]) %*% MTL_beta)
  h <- rep(0,k)
  for (i in 1:k) {
    h[i] <- sum((data$survtime == dt[i]) & (data$I20to25 == 1))/sum(MTL_eXB[data$survtime >= dt[i]])
  }
  MTL_H0 <- sum(h)
  MTL_S0 <- exp(-MTL_H0)
  MTL_S <- MTL_S0 ^ MTL_eXB
  MTL_H <- 1-MTL_S
  MTL_H <-as.data.frame(MTL_H)
  MTL_HH <- cbind(MTL_HH,MTL_H)
  risk_quantile[m,] <- quantile(MTL_HH[,m+1],na.rm = T,probs = c(0.33,0.5,0.66))
  m <- m+1
}

names(MTL_HH)[1] <- "eid"
names(MTL_HH)[c(2,3,4,5)] <- c("t1_CHD","t3_CHD","t5_CHD","t7_CHD")
names(risk_quantile)[1:3] <- c("Low risk","Average risk", "High risk")

MTL_HH_allDiseases <- inner_join(MTL_HH_allDiseases,MTL_HH,by="eid")
load("~/zhangshuaijie/UKB/disease_data_set/20211019experiment/varsdata.RData")
MTL_allDiseases <- left_join(MTL_HH_allDiseases,data[,c("eid","age")],by="eid")

#Calibration for MTL hypertension
MTL_hypertension <- MTL_allDiseases[,c("eid","t1_hypertension","t3_hypertension","t5_hypertension","t7_hypertension")]
MTL_hypertension$t1_hypertension_calibration <- 0.055*MTL_hypertension$t1_hypertension
MTL_hypertension$t3_hypertension_calibration <- 0.165*MTL_hypertension$t3_hypertension
MTL_hypertension$t5_hypertension_calibration <- 0.28*MTL_hypertension$t5_hypertension
MTL_hypertension$t7_hypertension_calibration <- 0.4*MTL_hypertension$t7_hypertension

risk_quantile <- as.data.frame(matrix(NA,4,3))
risk_quantile[1,] <- quantile(MTL_hypertension$t1_hypertension_calibration,na.rm = T,probs = c(0.33,0.5,0.66))
risk_quantile[2,] <- quantile(MTL_hypertension$t3_hypertension_calibration,na.rm = T,probs = c(0.33,0.5,0.66))
risk_quantile[3,] <- quantile(MTL_hypertension$t5_hypertension_calibration,na.rm = T,probs = c(0.33,0.5,0.66))
risk_quantile[4,] <- quantile(MTL_hypertension$t7_hypertension_calibration,na.rm = T,probs = c(0.33,0.5,0.66))

MTL_allDiseases <- left_join(MTL_allDiseases,MTL_hypertension[,c("eid","t1_hypertension_calibration","t3_hypertension_calibration","t5_hypertension_calibration","t7_hypertension_calibration")],by="eid")
save(MTL_allDiseases,file = "~/zhangshuaijie/UKB/disease_data_set/20211019experiment/MTL_data/MTL_allDiseases_year1357.RData")
#------------
