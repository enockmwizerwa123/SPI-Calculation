library(SPEI)
library(foreign)
library(gplots)
library(plm)
library(psych)
library(lmtest)
library(ggplot2)
library(ggpubr)
library(rcompanion)
library(nlme)
library(lme4)
library(car)
library(tseries)
library(xtable)
library(forecast)
library(TSA)
library(biwavelet)

New_data<-read.csv("Enock.csv") #load the data
View(New_data) # view the dataset
str(New_data) # stusture of dataset
nmiss <- (is.na(New_data)) # check for missing Values
nmiss
dim(New_data) # dimension of dataset
summary(New_data)
xtable(summary(New_data))


# Normalization of dataset using Log
min(New_data$SPI)
min(New_data$SPEI)
New_data$Precip<-log(New_data$Precip)
New_data$PET<-log(New_data$PET)
New_data$SPI<-log(New_data$SPI+3.89)
New_data$SPEI<-log(New_data$SPEI+3.702)
New_data$TB_cases<-log(New_data$TB_cases)
New_data$TB_Inc<-log(New_data$TB_Inc)
fr<-New_data[-c(1,2)]
hist(fr)
#new_data1<-New_data[c("New_data$Precip","New_data$PET","New_data$SPI","New_data$SPEI","New_data$TB_Inc","New_data$")]
View(New_data)
df<-New_data[-c(10,11,12,13,14)]
View(df)
plotNormalHistogram( New_data$TB_cases, prob = FALSE,
                     main = "Normal Distribution overlay on Histogram",
                     length = 1000 )

plotNormalHistogram( New_data$TB_Inc, prob = FALSE,
                     main = "Normal Distribution overlay on Histogram",
                     length = 1000 )


p1 <- ggplot(df, aes(Year, Precip)) +
  geom_point() +
  theme_bw() +
  labs(x = "Year") +
  labs(y = "Precip") +
  geom_smooth()
p2<-ggplot(df, aes(Year, PET)) +
  geom_point() +
  theme_bw() +
  labs(x = "Year") +
  labs(y = "PET") +
  geom_smooth()
p3<-ggplot(df, aes(Year, SPI)) +
  geom_point() +
  theme_bw() +
  labs(x = "Year") +
  labs(y = "SPI") +
  geom_smooth()
p4<-ggplot(df, aes(Year, SPEI)) +
  geom_point() +
  theme_bw() +
  labs(x = "Year") +
  labs(y = "SPEI") +
  geom_smooth()
p5<-ggplot(df, aes(Year, TB_cases)) +
  geom_point() +
  theme_bw() +
  labs(x = "Year") +
  labs(y = "TB_cases") +
  geom_smooth()
p6<-ggplot(df, aes(Year, TB_Inc)) +
  geom_point() +
  theme_bw() +
  labs(x = "Year") +
  labs(y = "TB_Inc") +
  geom_smooth()
# display plot3
ggpubr::ggarrange(p1, p2,p3,p4,p5,p6, ncol = 3, nrow = 2)

hist(New_data)
fr<-New_data(-c(1,2,2))
##### Performing PCA #######################
enock<-read.csv("PCA.csv")
View(enock)
enock<-enock[-1]
View(enock)
enock1<-enock[-6]
enock2<-enock[-5]
View(enock2)
pc1 <- principal(enock1, nfactors =5 , rotate = "none")
pc1 # the full model to see the loadings, Eigenvalues  (SS loadings) and variance accounted for
par(mfrow = c(1,2))
pve <- 100* pc1$values/sum(pc1$values) 
plot(pc1$values, type ="o" , ylab ="Eigenvalues " , xlab ="Principal Component" ,
     col ="blue", main = 'screeplot')

plot (cumsum(pve) , type ="o" , ylab ="Cumulative PVE " , xlab ="Principal Component " , col ="brown3 ",
      main = 'Percentage Variance Explained (PVE) plot')
dev.off()
pc2 <- principal(enock2, nfactors =5 , rotate = "none")
pc2
xtabs(pc2$values)
# the full model to see the loadings, Eigenvalues  (SS loadings) and variance accounted for
par(mfrow = c(1,2))
pve <- 100* pc2$values/sum(pc2$values)
plot(pc2$values, type ="o" , ylab ="Eigenvalues " , xlab ="Principal Component" ,
     col ="blue", main = 'screeplot')

plot (cumsum(pve) , type ="o" , ylab ="Cumulative PVE " , xlab ="Principal Component " , col ="brown3 ",
      main = 'Percentage Variance Explained (PVE) plot')
dev.off()


## Panel Data model for TB_cases as response Variable#################
coplot(TB_cases~Year|Country,type="l",data = New_data) # timeseries plot for TB incidence for each country
plotmeans(TB_Inc~Country,main="heterogeneity accross countries",data=New_data) #plot of 95% confidence interval around the means
# OLS Model
Ols<-lm(TB_cases~Precip+PET+SPEI+SPI,data = df)

summary(Ols$residuals)
xtable(summary(Ols))

par(mfrow=c(2,2))
plot(Ols)



#test for nonconstant variance
bptest(Ols)
#The Breusch-Pagan test indicates that there is some evidence of nonconstant variance is significant enough at the 0.05 level.

# Fixed Effects Model
fixed<-plm(TB_cases~Precip+PET+SPEI+SPI,data =df,index = c("Country","Year"),  model = "within")
summary(fixed)
xtable(summary(fixed))
#Check for model fit
par(mfrow=c(2,2))
plot(fixed)
E1<-fixed$residuals
ggplot(New_data, aes(sample=E1))+
  geom_qq()+
  geom_qq_line()
bptest(fixed)
pcdtest(fixed)
fixef(fixed) #display the fixed effects (constant for each country)
# F test for individuals effect
pFtest(fixed,Ols)# if p-value <5% then fixed model is better than OLS




#Random Effects Model
Random<-plm(TB_cases~Precip+PET+SPEI+SPI,data = df,index =c("Country","Year"),model = "random")
summary(Random)
E2<-Random$residuals
ggplot(New_data, aes(sample=E2))+
  geom_qq()+
  geom_qq_line()
bptest(Random)


# F test for individuals effect
phtest(fixed,Random)# if p-value <5% then fixed model is better than Random



## Panel Data model for TB_Inc as response Variable#################
coplot(TB_Inc~Year|Country,type="l",data = New_data) # timeseries plot for TB incidence for each country
plotmeans(TB_Inc~Country,main="heterogeneity accross countries",data=New_data) #plot of 95% confidence interval around the means
# OLS Model
OlsM<-lm(TB_Inc~Precip+PET+SPEI+SPI,data = df)
summary(OlsM)

xtable(summary(OlsM))
E12<-Ols1$residuals
ggplot(New_data, aes(sample=E12))+
  geom_qq()+
  geom_qq_line()
bptest(Ols1)
par(mfrow=c(2,2))
plot(Ols1)
par(mfrow=c(1,1))


# Fixed Effects Model

fixed_mod1<-lm(TB_Inc~Precip+PET+SPEI+factor(Country)-1,data=New_data)
summary(fixed_mod1)
xtable(summary(fixed_mod1))
fixed1<-plm(TB_Inc~Precip+PET+SPEI+SPI,data = df,index =c("Country","Year"),model = "within")
summary(fixed1)
fixef(fixed) #display the fixed effects (constant for each country)
#par(mfrow=c(2,2))
x<-df[C("Precip","PET","SPI"),"SPEI"]
yhat<-fixed1$fitted
plot(yhat)
plot(fixed1)
par(mfrow=c(1,1))
bptest(fixed1)
# F test for individuals effect
pFtest(fixed,Ols)# if p-value <5% then fixed model is better than OLS




#Random Effects Model
Random1<-plm(TB_Inc~Precip+PET+SPEI+SPI,data = df,index =c("Country","Year"),model = "random")
summary(Random1)
par(mfrow=c(2,2))
plot(Random1)
par(mfrow=c(1,1))

E23<-Random1$residuals
ggplot(New_data, aes(sample=E23))+
  geom_qq()+
  geom_qq_line()
bptest(Random1)
plot(Random1$residuals)
# F test for individuals effect
phtest(fixed1,Random1)# if p-value <5% then fixed model is better than Random


el<-read.csv("precip.csv")
el[is]

#Plot of SPI for each Country
par(mfrow=c(3,2))
New_data1<-read.csv("Burundi12.csv")
View(New_data1)
MB<-spi(New_data1$PRE,12)
plot(spi(ts(New_data1$PRE,freq=12,start=c(2000,1)),scale =12),main = "SPI characteristic In Burundi since 2000-2018")


Newr<-read.csv("Rwanda.csv")
MR<-spi(Newr$PRE,12)
plot(spi(ts(Newr$PRE,freq=12,start=c(2000,1)),scale =12),main = "SPI characteristic In Rwanda since 2000-2018")

plot.spei(MR,main = "SPI characteristic In Rwanda since 2000-2018")

Newk<-read.csv("Kenya13.csv")
MK<-spi(Newk$PRE,12)
plot.spei(MK,main = "SPI characteristic In Kenya since 2000-2018")
plot(spi(ts(Newk$PRE,freq=12,start=c(2000,1)),scale =12),main = "SPI characteristic In Kenya since 2000-2018")


Newc<-read.csv("Congo.csv")
MC<-spi(Newc$PRE,12)
plot.spei(MC,main = "SPI characteristic In DRC since 2000-2018")
plot(spi(ts(Newc$PRE,freq=12,start=c(2000,1)),scale =12),main = "SPI characteristic In DRC since 2000-2018")



Newt<-read.csv("Tanzania12.csv")
MT<-spi(Newt$PRE,12)
plot.spi(MT,main = "SPI characteristic In Tanzania since 2000-2018")

plot(spi(ts(Newt$PRE,freq=12,start=c(2000,1)),scale =12),main = "SPI characteristic In Tanzania since 2000-2018")




Newu<-read.csv("Uganda.csv")
MU<-spi(Newu$PRE,12)
plot.spei(MU,main = "SPI characteristic In Uganda since 2000-2018")
plot(spi(ts(Newu$PRE,freq=12,start=c(2000,1)),scale =12),main = "SPI characteristic In UGANDA since 2000-2018")
dev.off()





Briano<-read.csv("Data12.csv")
View(Briano)




TB_Inc <- ts(Briano$Log_tbinc, start = c(2000,1), end = c(2018,12), frequency = 12)
PRET <- ts(Briano$Log_pet, start = c(2000,1), end = c(2018,12), frequency = 12)
SPEI_s <- ts(Briano$Log_spei, start = c(2000,1), end = c(2018,12), frequency = 12)
SPI_s <- ts(Briano$Log_spi, start = c(2000,1), end = c(2018,12), frequency = 12)
TB_cases <- ts(Briano$Log_tbcases, start = c(2000,1), end = c(2018,12), frequency = 12)
PREsT <- ts(Briano$Log_pre, start = c(2000,1), end = c(2018,12), frequency = 12)
par(mfcol = c(1,1))
forecast::Ccf(TB_Inc,PRET)
dev.off()
par(mfrow=c(2,3))

forecast::Acf(PREsT, lag.max =36)
forecast::Acf(PRET, lag.max = 36)
forecast::Acf(SPEI_s, lag.max = 36)
forecast::Acf(SPI_s, lag.max = 36)
forecast::Acf(TB_cases, lag.max =36)
forecast::Acf(TB_Inc, lag.max =36)
dev.off()








## Panel Data model for TB_cases as response Variable#################
coplot(Log_tbcases~Year|Country,type="l",data = Briano) # timeseries plot for TB incidence for each country
plotmeans(Log_tbinc~Country,main="heterogeneity accross countries",data=New_data) #plot of 95% confidence interval around the means
# OLS Model
Olsmodel<-lm(Log_TB_Case~Log_PRE+Log_PET+Log_SPEI+Log_SPI,data = Briano)
summary(Olsmodel)
xtable(summary(Olsmodel))
E2b<-Olsb$residuals

ggplot(Briano, aes(sample=E2b))+
  geom_qq()+
  geom_qq_line()
plot(Olsb)
lmtest::bptest(Olsmodel)  # Breusch-Pagan test

# Fixed Effects Model
fixedmodel<-plm(Log_TB_Case~Log_PRE+Log_PET+Log_SPEI+Log_SPI,data = Briano,index =c("Country","Year"),model = "within")
summary(fixedmodel)
xtable(summary(fixedmodel))
lmtest::bptest(fixedb)  # Breusch-Pagan test
fixef(fixed) #display the fixed effects (constant for each country)
# F test for individuals effect
pFtest(fixedmodel,Olsmodel)# if p-value <5% then fixed model is better than OLS
E2c<-fixedb$residuals
View(E2c)
plot(fixedb)


#Random Effects Model
Randommmodel<-plm(Log_TB_Case~Log_PRE+Log_PET+Log_SPEI+Log_SPI,data = Briano,index =c("Country","Year"),model = "random")
summary(Randommmodel)
xtable(summary(Random1))
# F test for individuals effect
phtest(fixedmodel,Randommmodel)# if p-value <5% then fixed model is better than Random



## Panel Data model for TB_Inc as response Variable#################
coplot(Log_tbinc~Year|Country,type="l",data = Briano) 
coplot(Log_tbcases~Year|Country,type="l",data = Briano) 

# timeseries plot for TB incidence for each country
plotmeans(TB_Inc~Country,main="heterogeneity accross countries",data=New_data) #plot of 95% confidence interval around the means
# OLS Model
Olscmodel<-lm(Log_TB_Inc~Log_PRE+Log_PET+Log_SPEI+Log_SPI,data = Briano)
summary(Olscmodel)
xtable(summary(Olscmodel))


lmtest::bptest(Olscmodel)  # Breusch-Pagan test
dev.off()
# Fixed Effects Model
fixedmodel1<-plm(Log_TB_Inc~Log_PRE+Log_PET+Log_SPEI+Log_SPI,data = Briano,index =c("Country","Year"),model = "within")
summary(fixedmodel1)
xtable(summary(fixedmodel1))
fixef(fixed) #display the fixed effects (constant for each country)
# F test for individuals effect
pFtest(fixedmodel1,Olsmodel)# if p-value <5% then fixed model is better than OLS




#Random Effects Model
Randommodel2<-plm(Log_TB_Inc~Log_PRE+Log_PET+Log_SPEI+Log_SPI,data = Briano,index =c("Country","Year"),model = "random")
summary(Randommodel2)
# F test for individuals effect
phtest(fixedmodel1,Randommodel2)# if p-value <5% then fixed model is better than Random

cor(New_data[, c('Precip', 'PET', 'SPI','SPEI')])


lmtest::bptest(lmMod)  # Breusch-Pagan test
#studentized Breusch-Pagan test

ggplot(Briano,                            # Draw ggplot2 time series plot
       aes(x = year,
           y = X,
           col = Country_id)) +
woe<-read.csv("Woe.csv")


X<-cbind(Briano$Log_pre,Briano$Log_pet,Briano$Log_speicat,Briano$Log_spicat)



plot(Briano$Log_tbinc, Briano$Log_pet,                                   # Create plot with groups
     main = "This is my Plot",
     xlab = "X-Values",
     ylab = "Y-Values",
     col = Briano$Country_id,
     pch = Briano$Year)




ggplot2::ggplot(New_data1) +
  geom_bar(aes(x = Year, y =PRE, col = sign, fill = sign),
           show.legend = F, stat = "identity") +
  scale_color_manual(values = c("pos" = "darkblue", "neg" = "red")) +
  scale_fill_manual(values = c("pos"  = "darkblue", "neg" = "red")) +
  scale_y_continuous(limits = c(-3, 3), 
                     breaks = -3:3) +
  ylab("SPEI") + ggtitle("12-Month SPEI") +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5))



ggplot(data = New_data, aes(x = Year, y = Precip)) +
  geom_col(data = New_data["PRE" <= 0], fill = "red") +
  geom_col(data = New_data["PRE" >= 0], fill = "blue") +
  theme_bw()

plot(spi(ts(New_data$Precip,freq=12,start=c(2000,1)),scale = 12))




















































































































































































