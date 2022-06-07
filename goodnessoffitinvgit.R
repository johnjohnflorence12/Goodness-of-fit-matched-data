load("S:/grups/EPIC_CERVIX/2_Bases de dades/epic_cx_cas_control_inv.Rdata")


library(rms)
library(survival)
library(Epi)
library(epiDisplay)

rownames(epic_cx_cas_control_inv)
rownames(epic_cx_cas_control_inv) <- seq(1,dim(epic_cx_cas_control_inv)[1])



## Tenemos:

## ICC - 513 observaciones, 171 estratos, 3 sujetos por estrato


modeloinv2 <- clogit(cc_inv~ Bmi_C_new + Pa_Index_new + mar.stat.new + pill.dur.new + dur.tmq.smk.new + l.school.new2 + n.ftp.new1 + hpv.1.L1+resultat.chlam.tracho.2+resultat.HHV2.t1t2.1 + strata(caseset),epic_cx_cas_control_inv)
print(summary(clogit(cc_inv~ Bmi_C_new + Pa_Index_new + mar.stat.new + pill.dur.new + dur.tmq.smk.new + l.school.new2 + n.ftp.new1 + hpv.1.L1+resultat.chlam.tracho.2+resultat.HHV2.t1t2.1 + strata(caseset),epic_cx_cas_control_inv)))
anova(modeloinv2)
# Para no tener problemas con los NAs

NAs2 <- which(is.na(modeloinv2$coefficients)) # identificamos cuales de los coeficientes son valores perdidos, para posteriormente quitar las
# variables a las que corresponden
## precisamente son las que habiamos intentado quitar de la base de datos
cfs2 <- na.omit(modeloinv2$coefficients) # quitamos los NAs de los coeficientes
X2 <- model.matrix(modeloinv2)[, -NAs2] # quitamos las columnas correspondientes a las variables con NAs
X2 %*% cfs2 # nos da el vector que queriamos, con valores lógicos
(vecexp2 <- (exp(X2 %*% cfs2)))
length(vecexp2)
row.names(vecexp2)
sujetos2 <- as.numeric(names(vecexp2[,]))
length(sujetos2)


#Variable respuesta

epic_cx_cas_control_inv <- epic_cx_cas_control_inv[with(epic_cx_cas_control_inv, order(epic_cx_cas_control_inv$caseset, decreasing =T)), ] # Orden directo 

rownames(epic_cx_cas_control_inv)
length(rownames(epic_cx_cas_control_inv))

epic_cx_cas_control_inv <- epic_cx_cas_control_inv[as.integer(rownames(epic_cx_cas_control_inv)) %in% sujetos2,]
rownames(epic_cx_cas_control_inv)
length(rownames(epic_cx_cas_control_inv))

yvec2 <- epic_cx_cas_control_inv$cc_inv

epic_cx_cas_control_inv[1:20,c("idepic","caseset")]
head(vecexp2)


# vector de sumas exponenciales

n2 <- 513

coefestratos2 <- sort(unique(epic_cx_cas_control_inv$caseset))
length(coefestratos2)



sumas2 <- c()

for (j in coefestratos2){
  sumatorioinv=0
  for (i in 1:length(sujetos2)){
    individuo2 <- sujetos2[i]
    if (epic_cx_cas_control_inv[rownames(epic_cx_cas_control_inv)==individuo2,]$caseset==j){
      sumatorioinv=sumatorioinv+vecexp2[rownames(vecexp2)==individuo2]
    }
  }
  sumas2 <- c(sumas2, sumatorioinv)
}

length(sumas2)


#construcción del vector tita

titaaux2 <- matrix(NA, nrow=n2, ncol = 3) #numeric(n)


tita2 <- numeric(n2) #numeric(n)

for (j in coefestratos2){
  for (k in 1:length(sujetos2)){
    individuo2 <- sujetos2[k]
    if (epic_cx_cas_control_inv[rownames(epic_cx_cas_control_inv)==individuo2,]$caseset==j){
      tita2[k] <- vecexp2[k]/sumas2[match(j,coefestratos2)] # tita[k]#porque no queremos la j (el estrato), sino la posición. Ya que en el vector suma
      # la posición 1 será para el estrato 1, pero la posición 3 será para el estrato 7
      titaaux2[k,1] <- tita2[k]
      titaaux2[k,2] <- individuo2
      titaaux2[k,3] <- j
    }
  }
}

tita2 <- t(t(tita2))

kkkk2 <- data.frame(titaaux2)

#construcción del vector xtilda

matsumas2 = matrix(0, nrow = 171, ncol = 35)
k2<-0
for (j in coefestratos2){
  k2 <- k2+1 #??? esto es el length para que recorra hasta 248(650)
  sumatorioinv2 <- numeric(30)
  for (i in 1:length(sujetos2)){ #recorremos todas las filas de la matriz de diseño
    individuo2 <- sujetos2[i]
    if (epic_cx_cas_control_inv[rownames(epic_cx_cas_control_inv)==individuo2,]$caseset==j){
      sumatorioinv2=sumatorioinv2+(X2[i,]*tita2[i]) #ahora el sumatorio no es un valor sino una fila entera
    }
  }
  matsumas2[k2, ] <- sumatorioinv2
}

# Ahora hacemos la resta

matXtilda2 = matrix(0, nrow = n2, ncol = 35)

for (j in coefestratos2){
  for (i in 1:n2){
    individuo2 <- sujetos2[i]
    if (epic_cx_cas_control_inv[rownames(epic_cx_cas_control_inv)==individuo2,]$caseset==j){
      matXtilda2[i,] <- X2[i,]-matsumas2[match(j,coefestratos2),]
    }
  }
}




X2[1,] # fila 1 todas las columnas



## y la matriz diagonal U

U2 <- diag(as.numeric(tita2), n2)
dim(U2)



## Los valores de apalancamiento

# apal: nx1, tita: nx1, xtilda:px1, mat1:nxp, U: pxp

apal2 <- numeric(n2)

transpuesta2 <- t(matXtilda2)

M2 <- (transpuesta2 %*% U2 %*% matXtilda2)
M2 <- solve(M2)

apalaux2 <- matrix(NA, nrow=n2, ncol = 2) #numeric(n)


for (i in 1:n2){
  apal2[i] <- tita2[i]*(t(matXtilda2[i,]) %*% M2 %*% (matXtilda2[i,]))
  apalaux2[i,1] <- tita2[i]
  apalaux2[i,2] <- apal2[i]
}

apal2
length(apal2)


#Gráfico1

windows()
ggplot(epic_cx_cas_control_inv, aes(x = tita2, y = apal2)) +
  geom_point(
    aes(color = cc_inv, shape =cc_inv),
    size = 1.5, 
    alpha = 0.8 # It's nice to add some transparency because there may be overlap.
  )+
  # Use custom colors
  scale_color_manual(
    values = c("blue", "red")) +
  # Use custom colors
  
  labs(x = "Probabilidad" , y = "Apalancamiento") 



epic_cx_cas_control_inv$cc_inv <- as.factor(epic_cx_cas_control_inv$cc_inv) #por si el gráfico no funciona hacer esto


windows(width = 8)
par(las = 1, font = 2, font.lab = 4, font.axis = 2)
plot(apal2 ~ tita2, xlab = "Probabilidad", ylab = "Apalancamiento")
#plot(apal ~ tita, xlab = "Probabilidad", ylab = "Apalancamiento", ylim = c(0, 0.25))

identify(tita2,apal2)

# 508 169 393

#Bmi_C_new + Pa_Index_new + mar.stat.new + pill.dur.new + dur.tmq.smk.new + l.school.new2 + n.ftp.new1 + hpv.1.L1+resultat.chlam.tracho.2+resultat.HHV2.t1t2.1


# 508 tita 9.659895e-09 y apal 1 estrato 608

epic_cx_cas_control_inv[508,]$Pa_Index_new
epic_cx_cas_control_inv[508,]$mar.stat.new 

#169 tita 1.639634e-08 y apal 9.021989e-01 estrato 191

epic_cx_cas_control_inv[169,]$Pa_Index_new 
epic_cx_cas_control_inv[169,]$n.ftp.new1 


# 393 tita 6.267593e-09 y apal 6.860819e-01 estrato 481

epic_cx_cas_control_inv[393,]$pill.dur.new 
epic_cx_cas_control_inv[393,]$n.ftp.new1 




kapal2 <- data.frame(apalaux2)



epic_cx_cas_control_inv[508,]$cc_inv
epic_cx_cas_control_inv[508,]$caseset
epic_cx_cas_control_inv[508,]$country.x
epic_cx_cas_control_inv[508,]$Bmi_C_new
epic_cx_cas_control_inv[508,]$Pa_Index_new #NA
epic_cx_cas_control_inv[508,]$Bmi_Adj #NA

epic_cx_cas_control_inv[169,]$cc_inv
epic_cx_cas_control_inv[169,]$caseset
epic_cx_cas_control_inv[169,]$country.x
epic_cx_cas_control_inv[169,]$Bmi_C_new
epic_cx_cas_control_inv[169,]$Pa_Index_new
epic_cx_cas_control_inv[169,]$Bmi_Adj

epic_cx_cas_control_inv[393,]$cc_inv
epic_cx_cas_control_inv[393,]$caseset
epic_cx_cas_control_inv[393,]$country.x
epic_cx_cas_control_inv[393,]$Bmi_C_new
epic_cx_cas_control_inv[393,]$Pa_Index_new
epic_cx_cas_control_inv[393,]$Bmi_Adj


# Standardized Pearson residual

res2 <- numeric(n2)

for (i in 1:n2){
  res2[i] <- (yvec2[i]-tita2[i])/((tita2[i]*(1-apal2[i]))^(1/2))
}

res2

# Diagnósticos de falta de ajuste y influyentes

falta2 <- res2^2

# ones <- rep(1,n)
# influ1 <- falta1*(apal/(ones-apal))

diagnos2 <- matrix(NA, nrow=n2, ncol = 4) #numeric(n)


influ2 <- numeric(n2)
for (i in 1:n2){
  influ2[i] <- falta2[i]*(apal2[i]/(1-apal2[i]))
  diagnos2[i,1] <- tita2[i]
  diagnos2[i,2] <- apal2[i]
  diagnos2[i,3] <- falta2[i]
  diagnos2[i,4] <- influ2[i]
  
}

diagnosdf2 <- data.frame(diagnos2)


diagnosf2 <- matrix(NA, nrow=171, ncol = 2) #numeric(n)

faltatot2 <- c()
p2 <-0
for (j in coefestratos2){
  p2 <- p2+1
  sumatorioinv3=0
  for (i in 1:length(sujetos2)){
    individuo2 <- sujetos2[i]
    if (epic_cx_cas_control_inv[rownames(epic_cx_cas_control_inv)==individuo2,]$caseset==j){
      sumatorioinv3=sumatorioinv3+((falta2[i]))
    }
  }
  faltatot2 <- c(faltatot2, sumatorioinv3)
  diagnosf2[p2,1] <- faltatot2[p2]
  diagnosf2[p2,2] <- j
}

diagnosfdt2 <- data.frame(diagnosf2)



windows()
par( font = 2, font.lab = 4, font.axis = 2, las = 1,
     oma = c(0, 0, 1, 0), mar = c(5, 5, 4, 3))
#plot(coefestratos2,faltatot2,pch=19,col="black", xlab = "Nº de estrato", ylab = "Cambio en Pearson chi-square")
#identify(coefestratos2,faltatot2) # 57 24 58
plot(coefestratos2,faltatot2,pch=19,col="black", ylim=c(0,70), xlab = "Nº de estrato", ylab = "Cambio en Pearson chi-square")
#identify(coefestratos2,faltatot2) #89 73 20
#plot(coefestratos2,faltatot2,pch=19,col="black", ylim=c(0,20), xlab = "Nº de estrato", ylab = "Cambio en Pearson chi-square")
title("Diagnóstico de falta de ajuste", outer = T, cex.main = 1.5)


####falta de daignostico: 57 24 58 89 73 20

# 57 fa:9.378857e+09 estrato 190 , 5.567969e+07
###individuo 170 9.033761e-01, apal: 3.216427e-03, falta: 9.062911e-01, influ: 2.924425e-03 caso
###individuo 177 9.662393e-02, apal: 3.007167e-02, falta: 9.961966e-02, influ: 3.088609e-03 control
###individuo 181 1.072558e-10, apal: 5.901688e-03, falta: 9.378857e+09, influ: 5.567969e+07 control

epic_cx_cas_control_inv$caseset==190
epic_cx_cas_control_inv[169,]$cc_inv 
epic_cx_cas_control_inv[170,]$cc_inv
epic_cx_cas_control_inv[171,]$cc_inv# lo quitamos por tremenda falta de ajuste y influencia

# 24 fa:2.810488e+09 estrato 63 , 1.189249e+08
###individuo 69 3.708661e-10, apal: 4.059686e-02, falta: 2.810488e+09, influ: 1.189249e+08 coontrol
###individuo 78 1, apal: 9.854267e-11, falta: 9.002209e-19, influ: 8.871017e-29 caso
###individuo 82 5.779336e-10, apal: 6.326349e-02, falta: 6.169650e-10, influ: 4.166738e-11 control

epic_cx_cas_control_inv$caseset==63
epic_cx_cas_control_inv[70,]$cc_inv # lo quitamos por tremenda falta de ajuste
epic_cx_cas_control_inv[71,]$cc_inv
epic_cx_cas_control_inv[72,]$cc_inv


# 58 fa:6.236049e+08 estrato 191 ,5.752653e+09
###individuo 169 1.639634e-08, apal: 9.021989e-01, falta: 6.236049e+08, influ: 5.752653e+09 coontrol
###individuo 171 1, apal: 1.795987e-08, falta: 1, influ: 1.795987e-08 caso
###individuo 180 1.670156e-09, apal: 9.189936e-02, falta: 1.839175e-09, influ: 1.861237e-10 control

epic_cx_cas_control_inv$caseset==191
epic_cx_cas_control_inv[172,]$cc_inv # lo quitamos por tremenda falta de ajuste y influencia
epic_cx_cas_control_inv[173,]$cc_inv
epic_cx_cas_control_inv[174,]$cc_inv

# 89 fa:5.340366e+01 estrato 285 ,1.168618e+00
###individuo 268 0.95021770, apal: 1.872004e-03, falta: 2.613007e-03, influ: 4.900732e-06 caso
###individuo 271 0.03133473, apal: 2.221223e-02, falta: 3.204655e-02, influ: 7.279958e-04 control
###individuo 279 0.01844757, apal: 2.141459e-02, falta: 5.336900e+01, influ: 1.167885e+00 control

epic_cx_cas_control_inv$caseset==285
epic_cx_cas_control_inv[265,]$cc_inv 
epic_cx_cas_control_inv[266,]$cc_inv
epic_cx_cas_control_inv[267,]$cc_inv # lo quitamos por tremenda falta de ajuste y influ

# 73 fa:4.009703e+01 estrato 230 ,1.040509e+00
###individuo 217 0.18079376, apal: 9.607916e-02, falta: 4.106508e+00, influ: 4.364872e-01 caso
###individuo 221 0.02733932, apal: 1.626708e-02, falta: 3.517691e+01, influ: 5.816879e-01 control
###individuo 222 0.79186691, apal: 2.671759e-02, falta: 8.136045e-01, influ: 2.233427e-02 control

epic_cx_cas_control_inv$caseset==230
epic_cx_cas_control_inv[217,]$cc_inv 
epic_cx_cas_control_inv[218,]$cc_inv # lo quitamos por tremenda falta de ajuste y influ
epic_cx_cas_control_inv[219,]$cc_inv 


# 20 fa:3.417269e+01 estrato 56 ,8.674786e-01
###individuo 59 0.19527489, apal: 8.011939e-02, falta: 2.122829e-01, influ: 1.848933e-02 control
###individuo 60 0.02913125, apal: 2.442157e-02, falta: 3.316651e+01, influ: 8.302546e-01 control
###individuo 61 0.77559386, apal: 2.305438e-02, falta: 7.938967e-01, influ: 1.873471e-02 caso

epic_cx_cas_control_inv$caseset==56
epic_cx_cas_control_inv[58,]$cc_inv 
epic_cx_cas_control_inv[59,]$cc_inv # lo quitamos por tremenda falta de ajuste y influ
epic_cx_cas_control_inv[60,]$cc_inv 





diagnosinf2 <- matrix(NA, nrow=171, ncol = 2) #numeric(n)



influtot2 <- c()
kte2 <- 0
for (j in coefestratos2){
  sumatorioinv4=0
  kte2 <- kte2+1
  for (i in 1:length(sujetos2)){
    individuo2 <- sujetos2[i]
    if (epic_cx_cas_control_inv[rownames(epic_cx_cas_control_inv)==individuo2,]$caseset==j){
      sumatorioinv4=sumatorioinv4+influ2[i]
    }
  }
  influtot2 <- c(influtot2, sumatorioinv4)
  diagnosinf2[kte2,1] <- influtot2[kte2]
  diagnosinf2[kte2,2] <- j
}

diagnosinfdt2 <- data.frame(diagnosinf2)


windows()
par( font = 2, font.lab = 4, font.axis = 2, las = 1,
     oma = c(0, 0, 1, 0), mar= c(5,5,1,5))
#plot(coefestratos2,influtot2,pch=19,col="black", xlab = "Nº de estrato", ylab = "Distancia de Cook")
#identify(coefestratos2, influtot2) # 58
plot(coefestratos2,influtot2,pch=19,col="black", ylim=c(0,5), xlab = "Nº de estrato", ylab = "Distancia de Cook")
identify(coefestratos2, influtot2) #50 106 133
title("Diagnóstico influyente", outer = T, cex.main = 1.5)

# 50 fa:1.803058e+01 estrato 167 ,3.043293e+00
###individuo 151 0.05822277, apal: 1.452366e-01, falta: 1.782205e+01, influ: 3.028222e+00 control
###individuo 154 0.14356757, apal: 7.950150e-02, falta: 1.559672e-01, influ: 1.347056e-02 control
###individuo 155 0.79820966, apal: 2.953406e-02, falta: 5.256582e-02, influ: 1.599729e-03 caso

epic_cx_cas_control_inv$caseset==167
epic_cx_cas_control_inv[148,]$cc_inv # lo quitamos por tremenda falta de ajuste y influ
epic_cx_cas_control_inv[149,]$cc_inv 
epic_cx_cas_control_inv[150,]$cc_inv 


# 106 fa:1.488877e+01 estrato 357 ,2.956411e+00
###individuo 286 0.03162888, apal: 4.347169e-02, falta: 3.306633e-02, influ: 1.502778e-03 control
###individuo 287 0.89399006, apal: 1.706903e-02, falta: 9.095146e-01, influ: 1.579412e-02 caso
###individuo 288 0.07438105, apal: 1.740635e-01, falta: 1.394619e+01, influ: 2.939114e+00 control

epic_cx_cas_control_inv$caseset==357
epic_cx_cas_control_inv[316,]$cc_inv 
epic_cx_cas_control_inv[317,]$cc_inv 
epic_cx_cas_control_inv[318,]$cc_inv # lo quitamos por tremenda falta de ajuste y influ

# 133 fa:6.201952e+00 estrato 474 ,1.594679e+00
###individuo 389 0.3023789, apal: 1.347040e-01, falta: 1.860043e+00, influ: 2.895601e-01 caso
###individuo 406 0.4825394, apal: 1.203997e-01, falta: 5.485894e-01, influ: 7.509090e-02 control
###individuo 408 0.2150817, apal: 2.448622e-01, falta: 3.793319e+00, influ: 1.230028e+00 control

epic_cx_cas_control_inv$caseset==474
epic_cx_cas_control_inv[397,]$cc_inv 
epic_cx_cas_control_inv[398,]$cc_inv 
epic_cx_cas_control_inv[399,]$cc_inv # lo quitamos por influ


listfa2 <- c(508, 169, 393,
             171, 70, 172,267, 218, 59,
             148, 318,399)

modelocisnew2 <- clogit(as.numeric(cc_inv)~ Bmi_C_new + Pa_Index_new + mar.stat.new + pill.dur.new + dur.tmq.smk.new + l.school.new2 + n.ftp.new1 + hpv.1.L1+resultat.chlam.tracho.2+resultat.HHV2.t1t2.1 + strata(caseset),epic_cx_cas_control_inv[-c(listfa2),])
print(summary(clogit(as.numeric(cc_inv)~ Bmi_C_new + Pa_Index_new + mar.stat.new + pill.dur.new + dur.tmq.smk.new + l.school.new2 + n.ftp.new1 + hpv.1.L1+resultat.chlam.tracho.2+resultat.HHV2.t1t2.1 + strata(caseset),epic_cx_cas_control_inv[-c(listfa),])))



