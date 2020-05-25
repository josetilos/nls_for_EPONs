

#  Non-linear regression of delay percentiles in TDM EPONs 
# Jose Alberto Hernandez
# 22 mayo 2020


library(caret)

rm(list=ls())


# setwd("~/Google Drive/Research/Rexperiments")

set.seed(1234)




# 1G-EPON and traffic parameters, delay values in us

Cpon = 1e9; Nonu = 8
avg_packet_size = 7/12*40 + 4/12*576 +1/12*1500; 
sd_packet_size = sqrt(7/12*40^2 + 4/12*576^2 +1/12*1500^2 - avg_packet_size^2)
dtx = avg_packet_size*8/Cpon*1e6


# Loading delays for 1G-EPON, 4km, 8km and 20 km datasets, and normalization

# 4 km
distance = 4; tau = 5*distance;
Delay4 = t(as.matrix(read.csv("Trimodal_1G_4km.csv", header = F)))
Load = as.matrix(c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.92, 0.94, 0.96,  0.98, 0.99),nrow=1)
Delay4 = Delay4[,1:13]; Load = Load[1:13]
Delay4_norm = (Delay4-dtx)/tau

# 8 km
distance = 8; tau = 5*distance;
Delay8 = t(as.matrix(read.csv("Trimodal_1G_8km.csv", header = F)))
Load = as.matrix(c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.92, 0.94, 0.96,  0.98, 0.99),nrow=1)
Delay8 = Delay8[,1:13]; Load = Load[1:13]
Delay8_norm = (Delay8-dtx)/tau


# 20 km
distance = 20; tau = 5*distance;
Delay20 = t(as.matrix(read.csv("Trimodal_1G_20km.csv", header = F)))
Load = as.matrix(c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.92, 0.94, 0.96,  0.98, 0.99),nrow=1)
Delay20 = Delay20[,1:13]; Load = Load[1:13]
Delay20_norm = (Delay20-dtx)/tau


# Aggregation of all datasets
Delay_norm = rbind(Delay4_norm,Delay8_norm,Delay20_norm)
Load = as.matrix(c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.92, 0.94, 0.96,  0.98, 0.99),nrow=1)
Delay_norm = Delay_norm[,1:12]; Load = Load[1:12]; q = as.numeric(Load)


# Booststrap sampling from normalized dataset
Nsample = 1e5
delay_sample = (Delay_norm[sample(1:dim(Delay_norm)[1],size=Nsample,replace=TRUE),])
delay_sample_mean = apply(delay_sample,2,mean)


# ML models: Simple and Full
model_simp <- nls(delay_sample_mean ~ a*2*(2-q)/(1-q), 
                  start = list(a=Nonu), algorithm = "default")
print(model_simp)

model_full <- nls(delay_sample_mean ~ (a+b*q)/(c-q), 
                  start = list(a=Nonu,b=20,c=1), algorithm = "default")
print(model_full)


# Models' evaluation at 4 km

distance = 4; tau = 5*distance

A = (((sd_packet_size*8)^2)/(avg_packet_size*8)+8*avg_packet_size)/(2*Cpon)*1e6
theo_eq = (4*tau + (A-2*tau)*q)/(1-q)+ dtx


cat("Simple model")
R2(delay_sample_mean*tau+dtx,predict(model_simp)*tau+dtx)
RMSE(delay_sample_mean*tau+dtx,predict(model_simp)*tau+dtx)

cat("Full model")
R2(delay_sample_mean*tau+dtx,predict(model_full)*tau+dtx)
RMSE(delay_sample_mean*tau+dtx,predict(model_full)*tau+dtx)

cat("Theo equation")
R2(delay_sample_mean*tau+dtx,theo_eq)
RMSE(delay_sample_mean*tau+dtx,theo_eq)


# Delay-mean plot


pdf("Rplot_1G_DelayMean.pdf") 
plot(q,theo_eq,type="p", pch=23, col="black",
     ylab = 'E(D) (us)', xlab='Load',main='EPON 1G, 8 ONUs, 4 km distance')
lines(q,theo_eq,lty=2,col='black', lwd=1, pch = 22)
points(q,delay_sample_mean*tau + dtx,col='black',pch=19)

points(q,predict(model_full)*tau + dtx,type="p", pch=22, col="blue")
lines(q,predict(model_full)*tau + dtx,lty=2,col="blue",lwd=1)

points(q,predict(model_simp)*tau+dtx ,type="p", pch=25, col="red")
lines(q,predict(model_simp)*tau+dtx ,lty=2,col="red",lwd=1)

legend("topleft",c("Observ.","NLS Model (Full)","NLS Model (Simp)","Theo Equation"),
       col=c("black","blue","red","black"), lty=c(1,2,2,2), lwd=c(2,1,1,1), pch=c(19, 22, 25, 23))
dev.off() 


# Percentiles

# 25th percentile
delay_sample_perc = apply(delay_sample,2,quantile,0.25)

# ML models: Simple and Full
model_simp <- nls(delay_sample_perc ~ a*2*(2-q)/(1-q), 
                  start = list(a=Nonu), algorithm = "default")
print(model_simp)

model_full <- nls(delay_sample_perc ~ (a+b*q)/(c-q), 
                  start = list(a=Nonu,b=20,c=1), algorithm = "default")
print(model_full)


# Models' evaluation at 4 km

distance = 4; tau = 5*distance

cat("Simple model")
R2(delay_sample_perc*tau+dtx,predict(model_simp)*tau+dtx)
RMSE(delay_sample_perc*tau+dtx,predict(model_simp)*tau+dtx)

cat("Full model")
R2(delay_sample_perc*tau+dtx,predict(model_full)*tau+dtx)
RMSE(delay_sample_perc*tau+dtx,predict(model_full)*tau+dtx)


# Plot

pdf("Rplot_25Perc_1G.pdf")
plot(q,delay_sample_perc*tau + dtx,col='black',pch=19,
     ylab = 'Delay perc (us)', xlab='Load',main='25-th percentile')

points(q,predict(model_full)*tau + dtx,type="p", pch=22, col="blue")
lines(q,predict(model_full)*tau + dtx,lty=2,col="blue",lwd=1)

points(q,predict(model_simp)*tau+dtx ,type="p", pch=25, col="red")
lines(q,predict(model_simp)*tau+dtx ,lty=2,col="red",lwd=1)

legend("topleft",c("Observ.","NLS Model (Full)","NLS Model (Simp)"),
       col=c("black","blue","red"), lty=c(1,2,2), lwd=c(2,1,1), pch=c(19, 22, 25))
dev.off() 


# 50th percentile
delay_sample_perc = apply(delay_sample,2,quantile,0.5)

# ML models: Simple and Full
model_simp <- nls(delay_sample_perc ~ a*2*(2-q)/(1-q), 
                  start = list(a=Nonu), algorithm = "default")
print(model_simp)

model_full <- nls(delay_sample_perc ~ (a+b*q)/(c-q), 
                  start = list(a=Nonu,b=20,c=1), algorithm = "default")
print(model_full)


# Models' evaluation at 4 km

distance = 4; tau = 5*distance

cat("Simple model")
R2(delay_sample_perc*tau+dtx,predict(model_simp)*tau+dtx)
RMSE(delay_sample_perc*tau+dtx,predict(model_simp)*tau+dtx)

cat("Full model")
R2(delay_sample_perc*tau+dtx,predict(model_full)*tau+dtx)
RMSE(delay_sample_perc*tau+dtx,predict(model_full)*tau+dtx)


# Plot

pdf("Rplot_50Perc_1G.pdf")
plot(q,delay_sample_perc*tau + dtx,col='black',pch=19,
     ylab = 'Delay perc (us)', xlab='Load',main='50-th percentile')

points(q,predict(model_full)*tau + dtx,type="p", pch=22, col="blue")
lines(q,predict(model_full)*tau + dtx,lty=2,col="blue",lwd=1)

points(q,predict(model_simp)*tau+dtx ,type="p", pch=25, col="red")
lines(q,predict(model_simp)*tau+dtx ,lty=2,col="red",lwd=1)

legend("topleft",c("Observ.","NLS Model (Full)","NLS Model (Simp)"),
       col=c("black","blue","red"), lty=c(1,2,2), lwd=c(2,1,1), pch=c(19, 22, 25))
dev.off() 


# 75th percentile
delay_sample_perc = apply(delay_sample,2,quantile,0.75)

# ML models: Simple and Full
model_simp <- nls(delay_sample_perc ~ a*2*(2-q)/(1-q), 
                  start = list(a=Nonu), algorithm = "default")
print(model_simp)

model_full <- nls(delay_sample_perc ~ (a+b*q)/(c-q), 
                  start = list(a=Nonu,b=20,c=1), algorithm = "default")
print(model_full)


# Models' evaluation at 4 km

distance = 4; tau = 5*distance

cat("Simple model")
R2(delay_sample_perc*tau+dtx,predict(model_simp)*tau+dtx)
RMSE(delay_sample_perc*tau+dtx,predict(model_simp)*tau+dtx)

cat("Full model")
R2(delay_sample_perc*tau+dtx,predict(model_full)*tau+dtx)
RMSE(delay_sample_perc*tau+dtx,predict(model_full)*tau+dtx)


# Plot

pdf("Rplot_75Perc_1G.pdf")
plot(q,delay_sample_perc*tau + dtx,col='black',pch=19,
     ylab = 'Delay perc (us)', xlab='Load',main='75-th percentile')

points(q,predict(model_full)*tau + dtx,type="p", pch=22, col="blue")
lines(q,predict(model_full)*tau + dtx,lty=2,col="blue",lwd=1)

points(q,predict(model_simp)*tau+dtx ,type="p", pch=25, col="red")
lines(q,predict(model_simp)*tau+dtx ,lty=2,col="red",lwd=1)

legend("topleft",c("Observ.","NLS Model (Full)","NLS Model (Simp)"),
       col=c("black","blue","red"), lty=c(1,2,2), lwd=c(2,1,1), pch=c(19, 22, 25))
dev.off() 


# 90th percentile
delay_sample_perc = apply(delay_sample,2,quantile,0.9)

# ML models: Simple and Full
model_simp <- nls(delay_sample_perc ~ a*2*(2-q)/(1-q), 
                  start = list(a=Nonu), algorithm = "default")
print(model_simp)

model_full <- nls(delay_sample_perc ~ (a+b*q)/(c-q), 
                  start = list(a=Nonu,b=20,c=1), algorithm = "default")
print(model_full)


# Models' evaluation at 4 km

distance = 4; tau = 5*distance

cat("Simple model")
R2(delay_sample_perc*tau+dtx,predict(model_simp)*tau+dtx)
RMSE(delay_sample_perc*tau+dtx,predict(model_simp)*tau+dtx)

cat("Full model")
R2(delay_sample_perc*tau+dtx,predict(model_full)*tau+dtx)
RMSE(delay_sample_perc*tau+dtx,predict(model_full)*tau+dtx)


# Plot

pdf("Rplot_90Perc_1G.pdf")
plot(q,delay_sample_perc*tau + dtx,col='black',pch=19,
     ylab = 'Delay perc (us)', xlab='Load',main='90-th percentile')

points(q,predict(model_full)*tau + dtx,type="p", pch=22, col="blue")
lines(q,predict(model_full)*tau + dtx,lty=2,col="blue",lwd=1)

points(q,predict(model_simp)*tau+dtx ,type="p", pch=25, col="red")
lines(q,predict(model_simp)*tau+dtx ,lty=2,col="red",lwd=1)

legend("topleft",c("Observ.","NLS Model (Full)","NLS Model (Simp)"),
       col=c("black","blue","red"), lty=c(1,2,2), lwd=c(2,1,1), pch=c(19, 22, 25))
dev.off() 


# 95th percentile
delay_sample_perc = apply(delay_sample,2,quantile,0.95)

# ML models: Simple and Full
model_simp <- nls(delay_sample_perc ~ a*2*(2-q)/(1-q), 
                  start = list(a=Nonu), algorithm = "default")
print(model_simp)

model_full <- nls(delay_sample_perc ~ (a+b*q)/(c-q), 
                  start = list(a=Nonu,b=20,c=1), algorithm = "default")
print(model_full)


# Models' evaluation at 4 km

distance = 4; tau = 5*distance

cat("Simple model")
R2(delay_sample_perc*tau+dtx,predict(model_simp)*tau+dtx)
RMSE(delay_sample_perc*tau+dtx,predict(model_simp)*tau+dtx)

cat("Full model")
R2(delay_sample_perc*tau+dtx,predict(model_full)*tau+dtx)
RMSE(delay_sample_perc*tau+dtx,predict(model_full)*tau+dtx)


# Plot

pdf("Rplot_95Perc_1G.pdf")
plot(q,delay_sample_perc*tau + dtx,col='black',pch=19,
     ylab = 'Delay perc (us)', xlab='Load',main='95-th percentile')

points(q,predict(model_full)*tau + dtx,type="p", pch=22, col="blue")
lines(q,predict(model_full)*tau + dtx,lty=2,col="blue",lwd=1)

points(q,predict(model_simp)*tau+dtx ,type="p", pch=25, col="red")
lines(q,predict(model_simp)*tau+dtx ,lty=2,col="red",lwd=1)

legend("topleft",c("Observ.","NLS Model (Full)","NLS Model (Simp)"),
       col=c("black","blue","red"), lty=c(1,2,2), lwd=c(2,1,1), pch=c(19, 22, 25))
dev.off() 



# 99th percentile
delay_sample_perc = apply(delay_sample,2,quantile,0.99)

# ML models: Simple and Full
model_simp <- nls(delay_sample_perc ~ a*2*(2-q)/(1-q), 
                  start = list(a=Nonu), algorithm = "default")
print(model_simp)

model_full <- nls(delay_sample_perc ~ (a+b*q)/(c-q), 
                  start = list(a=Nonu,b=20,c=1), algorithm = "default")
print(model_full)


# Models' evaluation at 4 km

distance = 4; tau = 5*distance

cat("Simple model")
R2(delay_sample_perc*tau+dtx,predict(model_simp)*tau+dtx)
RMSE(delay_sample_perc*tau+dtx,predict(model_simp)*tau+dtx)

cat("Full model")
R2(delay_sample_perc*tau+dtx,predict(model_full)*tau+dtx)
RMSE(delay_sample_perc*tau+dtx,predict(model_full)*tau+dtx)


# Plot

pdf("Rplot_99Perc_1G.pdf")
plot(q,delay_sample_perc*tau + dtx,col='black',pch=19,
     ylab = 'Delay perc (us)', xlab='Load',main='99-th percentile')

points(q,predict(model_full)*tau + dtx,type="p", pch=22, col="blue")
lines(q,predict(model_full)*tau + dtx,lty=2,col="blue",lwd=1)

points(q,predict(model_simp)*tau+dtx ,type="p", pch=25, col="red")
lines(q,predict(model_simp)*tau+dtx ,lty=2,col="red",lwd=1)

legend("topleft",c("Observ.","NLS Model (Full)","NLS Model (Simp)"),
       col=c("black","blue","red"), lty=c(1,2,2), lwd=c(2,1,1), pch=c(19, 22, 25))
dev.off() 







# Simple ML model: Modeling beta parameter

pp = seq(from=0.20,to=0.99, by=0.025)
beta = NA*(1:length(pp));
R2coef = NA*(1:length(pp));
RMSE = NA*(1:length(pp));


for (ii in c(1:length(pp))) {
  daux = apply(delay_sample,2,quantile,pp[ii]);  
  
  model <- nls(daux ~ a*2*(2-q)/(1-q), 
               start = list(a=Nonu), algorithm = "default")
  
  beta[ii] = coef(model)[1]

  R2coef[ii] = R2(daux*tau+dtx,predict(model)*tau+dtx)
  RMSE[ii] = RMSE(daux*tau+dtx,predict(model)*tau+dtx)
  
}


fit_be <- lm(beta ~ poly(pp,3))
print(summary(fit_be))


pdf("SimplML_Beta_1G.pdf")
plot(pp,beta,col='black',ylab = 'beta', xlab='Percentile',main="1G EPON Generalized Simple ML model",pch=21)
lines(pp,predict(fit_be),lty=2,col="red",lwd=2,pch=19)
legend("topleft",c("Beta coef","Regression Model"),col=c("black","red"), lty=c(1,2), lwd=c(2,1), pch=c(21,19))
dev.off() 

pdf("SimplML_R2_1G.pdf")
plot(pp,R2coef,ylab="R2 coef",xlab="Percentile",main="R2 for 1G EPON")
dev.off()






# 10G EPON


Cpon = 10e9;
dtx = avg_packet_size*8/Cpon*1e6

distance = 4; tau = 5*distance;
Delay4 = t(as.matrix(read.csv("Trimodal_10G_4km.csv", header = F)))
Load = as.matrix(c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.92, 0.94, 0.96,  0.98, 0.99),nrow=1)
Delay4 = Delay4[,1:13]; Load = Load[1:13]
Delay4_norm = (Delay4-dtx)/tau

distance = 8; tau = 5*distance;
Delay8 = t(as.matrix(read.csv("Trimodal_10G_8km.csv", header = F)))
Load = as.matrix(c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.92, 0.94, 0.96,  0.98, 0.99),nrow=1)
Delay8 = Delay8[,1:13]; Load = Load[1:13]
Delay8_norm = (Delay8-dtx)/tau

distance = 20; tau = 5*distance;
Delay20 = t(as.matrix(read.csv("Trimodal_10G_20km.csv", header = F)))
Load = as.matrix(c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.92, 0.94, 0.96,  0.98, 0.99),nrow=1)
Delay20 = Delay20[,1:13]; Load = Load[1:13]
Delay20_norm = (Delay20-dtx)/tau


Delay_norm = rbind(Delay4_norm,Delay8_norm,Delay20_norm)
Load = as.matrix(c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.92, 0.94, 0.96,  0.98, 0.99),nrow=1)
Delay_norm = Delay_norm[,1:12]; Load = Load[1:12]


q = as.numeric(Load)


Nsample = 1e5
delay_sample = (Delay_norm[sample(1:dim(Delay_norm)[1],size=Nsample,replace=TRUE),])
delay_sample_mean = apply(delay_sample,2,mean)


# ML models: Simple and Full
model_simp <- nls(delay_sample_mean ~ a*2*(2-q)/(1-q), 
                  start = list(a=Nonu), algorithm = "default")
print(model_simp)

model_full <- nls(delay_sample_mean ~ (a+b*q)/(c-q), 
                  start = list(a=Nonu,b=20,c=1), algorithm = "default")
print(model_full)


# Models' evaluation at 4 km

distance = 4; tau = 5*distance

A = (((sd_packet_size*8)^2)/(avg_packet_size*8)+8*avg_packet_size)/(2*Cpon)*1e6
theo_eq = (4*tau + (A-2*tau)*q)/(1-q)+ dtx


cat("Simple model")
R2(delay_sample_mean*tau+dtx,predict(model_simp)*tau+dtx)
RMSE(delay_sample_mean*tau+dtx,predict(model_simp)*tau+dtx)

cat("Full model")
R2(delay_sample_mean*tau+dtx,predict(model_full)*tau+dtx)
RMSE(delay_sample_mean*tau+dtx,predict(model_full)*tau+dtx)

cat("Theo equation")
R2(delay_sample_mean*tau+dtx,theo_eq)
RMSE(delay_sample_mean*tau+dtx,theo_eq)


# Delay-mean plot

pdf("Rplot_10G_DelayMean.pdf") 
plot(q,theo_eq,type="p", pch=23, col="black",
     ylab = 'E(D) (us)', xlab='Load',main='EPON 10G, 8 ONUs, 4 km distance')
lines(q,theo_eq,lty=2,col='black', lwd=1, pch = 22)
points(q,delay_sample_mean*tau + dtx,col='black',pch=19)

points(q,predict(model_full)*tau + dtx,type="p", pch=22, col="blue")
lines(q,predict(model_full)*tau + dtx,lty=2,col="blue",lwd=1)

points(q,predict(model_simp)*tau+dtx ,type="p", pch=25, col="red")
lines(q,predict(model_simp)*tau+dtx ,lty=2,col="red",lwd=1)

legend("topleft",c("Observ.","NLS Model (Full)","NLS Model (Simp)","Theo Equation"),
       col=c("black","blue","red","black"), lty=c(1,2,2,2), lwd=c(2,1,1,1), pch=c(19, 22, 25, 23))
dev.off() 



# Simple ML model: Modeling beta parameter

pp = seq(from=0.20,to=0.99, by=0.025)
beta = NA*(1:length(pp));
R2coef = NA*(1:length(pp));
RMSE = NA*(1:length(pp));


for (ii in c(1:length(pp))) {
  daux = apply(delay_sample,2,quantile,pp[ii]);  
  
  model <- nls(daux ~ a*2*(2-q)/(1-q), 
               start = list(a=Nonu), algorithm = "default")
  
  beta[ii] = coef(model)[1]
  
  R2coef[ii] = R2(daux*tau+dtx,predict(model)*tau+dtx)
  RMSE[ii] = RMSE(daux*tau+dtx,predict(model)*tau+dtx)
  
}


fit_be <- lm(beta ~ poly(pp,3))
print(summary(fit_be))


pdf("SimplML_Beta_10G.pdf")
plot(pp,beta,col='black',ylab = 'beta', xlab='Percentile',main="Generalized Simple ML model (10G EPON)",pch=21)
lines(pp,predict(fit_be),lty=2,col="red",lwd=2,pch=19)
legend("topleft",c("Beta coef","Regression Model"),col=c("black","red"), lty=c(1,2), lwd=c(2,1), pch=c(21,19))
dev.off() 

pdf("SimplML_R2_10G.pdf")
plot(pp,R2coef,ylab="R2 coef",xlab="Percentile",main="R2 for 10G EPON")
dev.off()






