load("S:/grups/EPIC_CERVIX/2_Bases de dades/epic_cx_cas_control_cis.Rdata")



library(rms)
library(survival)
library(Epi)
library(epiDisplay)

rownames(epic_cx_cas_control_cis)
rownames(epic_cx_cas_control_cis) <- seq(1,dim(epic_cx_cas_control_cis)[1])



## Tenemos:

## CIS - 1263 observaciones, 421 estratos, 3 sujetos por estrato, 907 variables

class(epic_cx_cas_control_cis$cc_cis)
epic_cx_cas_control_cis$cc_cis <- as.numeric(epic_cx_cas_control_cis$cc_cis)


modelocis2 <- clogit(cc_cis~ Bmi_C_new + Pa_Index_new + mar.stat.new + pill.dur.new + dur.tmq.smk.new + l.school.new2 + n.ftp.new1 + hpv.1.L1+resultat.chlam.tracho.2+resultat.HHV2.t1t2.1 + strata(caseset),epic_cx_cas_control_cis)
print(summary(clogit(cc_cis~ Bmi_C_new + Pa_Index_new + mar.stat.new + pill.dur.new + dur.tmq.smk.new + l.school.new2 + n.ftp.new1 + hpv.1.L1+resultat.chlam.tracho.2+resultat.HHV2.t1t2.1 + strata(caseset),epic_cx_cas_control_cis)))

# Para no tener problemas con los NAs

NAs <- which(is.na(modelocis2$coefficients)) # identificamos cuales de los coeficientes son valores perdidos, para posteriormente quitar las
# variables a las que corresponden
## precisamente son las que habiamos intentado quitar de la base de datos
cfs <- na.omit(modelocis2$coefficients) # quitamos los NAs de los coeficientes
X <- model.matrix(modelocis2)[, -NAs] # quitamos las columnas correspondientes a las variables con NAs
X %*% cfs # nos da el vector que queriamos, con valores lógicos
(vecexp <- (exp(X %*% cfs)))
length(vecexp)
row.names(vecexp)
sujetos <- as.numeric(names(vecexp[,]))
length(sujetos)


#Variable respuesta

epic_cx_cas_control_cis <- epic_cx_cas_control_cis[with(epic_cx_cas_control_cis, order(epic_cx_cas_control_cis$caseset)), ] # Orden directo 

rownames(epic_cx_cas_control_cis)
length(rownames(epic_cx_cas_control_cis))

epic_cx_cas_control_cis <- epic_cx_cas_control_cis[as.integer(rownames(epic_cx_cas_control_cis)) %in% sujetos,]
rownames(epic_cx_cas_control_cis)
length(rownames(epic_cx_cas_control_cis))

yvec <- epic_cx_cas_control_cis$cc_cis

epic_cx_cas_control_cis[1:20,c("idepic","caseset")]
head(vecexp)


# vector de sumas exponenciales

n <- 1263

coefestratos <- sort(unique(epic_cx_cas_control_cis$caseset))
length(coefestratos)

sumas <- c()

for (j in coefestratos){
  sumatorio=0
  for (i in 1:length(sujetos)){
    individuo <- sujetos[i]
    if (epic_cx_cas_control_cis[rownames(epic_cx_cas_control_cis)==individuo,]$caseset==j){
      sumatorio=sumatorio+vecexp[rownames(vecexp)==individuo]
    }
  }
  sumas <- c(sumas, sumatorio)
}

length(sumas)


#construcción del vector tita


for (j in coefestratos){
  for (k in 1:length(sujetos)){
    individuo <- sujetos[k]
    if (epic_cx_cas_control_cis[rownames(epic_cx_cas_control_cis)==individuo,]$caseset==j){
      tita[k] <- vecexp[k]/sumas[match(j,coefestratos)] # tita[k]#porque no queremos la j (el estrato), sino la posición. Ya que en el vector suma
      # la posición 1 será para el estrato 1, pero la posición 3 será para el estrato 7
    }
  }
}


#identifica sujetos

#matsumas = matrix(0, nrow = 248, ncol = 28)

tita <- numeric(n) #numeric(n)

titaaux <- matrix(NA, nrow=n, ncol = 3) #numeric(n)


for (j in coefestratos){
  for (k in 1:length(sujetos)){
    individuo <- sujetos[k]
    if (epic_cx_cas_control_cis[rownames(epic_cx_cas_control_cis)==individuo,]$caseset==j){
      tita[k] <- vecexp[k]/sumas[match(j,coefestratos)] # tita[k]#porque no queremos la j (el estrato), sino la posición. Ya que en el vector suma
      # la posición 1 será para el estrato 1, pero la posición 3 será para el estrato 7
      titaaux[k,1] <- tita[k]
      titaaux[k,2] <- individuo
      titaaux[k,3] <- j
    }
  }
}

tita <- t(t(tita))
titaaux


sujetos
coefestratos

#construcción del vector xtilda

matsumas = matrix(0, nrow = 421, ncol = 37)
k<-0
for (j in coefestratos){
  k <- k+1 #??? esto es el length para que recorra hasta 248(650)
  sumatorio2 <- numeric(28)
  for (i in 1:length(sujetos)){ #recorremos todas las filas de la matriz de diseño
    individuo <- sujetos[i]
    if (epic_cx_cas_control_cis[rownames(epic_cx_cas_control_cis)==individuo,]$caseset==j){
      sumatorio2=sumatorio2+(X[i,]*tita[i]) #ahora el sumatorio no es un valor sino una fila entera
    }
  }
  matsumas[k, ] <- sumatorio2
}

# Ahora hacemos la resta

matXtilda = matrix(0, nrow = n, ncol = 37)

for (j in coefestratos){
  for (i in 1:n){
    individuo <- sujetos[i]
    if (epic_cx_cas_control_cis[rownames(epic_cx_cas_control_cis)==individuo,]$caseset==j){
      matXtilda[i,] <- X[i,]-matsumas[match(j,coefestratos),]
    }
  }
}




X[1,] # fila 1 todas las columnas



## y la matriz diagonal U

U <- diag(as.numeric(tita), n)
dim(U)



## Los valores de apalancamiento

# apal: nx1, tita: nx1, xtilda:px1, mat1:nxp, U: pxp

apal <- numeric(n)

transpuesta <- t(matXtilda)

M <- (transpuesta %*% U %*% matXtilda)
M <- solve(M)

apalaux <- matrix(NA, nrow=n, ncol = 2) #numeric(n)


for (i in 1:n){
  apal[i] <- tita[i]*(t(matXtilda[i,]) %*% M %*% (matXtilda[i,]))
  apalaux[i,1] <- tita[i]
  apalaux[i,2] <- apal[i]
}

apal
length(apal)

apalaux

kapal <- data.frame(apalaux)





windows()
ggplot(epic_cx_cas_control_cis, aes(x = tita, y = apal)) +
  geom_point(
    aes(color = cc_cis, shape =cc_cis),
    size = 1.5, 
    alpha = 0.8 # It's nice to add some transparency because there may be overlap.
  )+
  # Use custom colors
  scale_color_manual(
    values = c("blue", "red")) +
  # Use custom colors
  
  labs(x = "Probabilidad" , y = "Apalancamiento") 



epic_cx_cas_control_cis$cc_cis <- as.factor(epic_cx_cas_control_cis$cc_cis)


windows(width = 8)
par(las = 1, font = 2, font.lab = 4, font.axis = 2)
plot(apal ~ tita, xlab = "Probabilidad", ylab = "Apalancamiento")
#plot(apal ~ tita, xlab = "Probabilidad", ylab = "Apalancamiento", ylim = c(0, 0.25))

identify(tita,apal)

# observaciones 17, 1045

# el 17 tiene tita: 3.909492e-08 y apal: 1.0000000061, individuo 17 estrato 12

epic_cx_cas_control_cis[17,]$pill.dur.new
epic_cx_cas_control_cis[17,]$n.ftp.new1 
epic_cx_cas_control_cis[17,]$resultat.chlam.tracho.2 



#Bmi_C_new + Pa_Index_new + mar.stat.new + pill.dur.new + dur.tmq.smk.new + l.school.new2 + n.ftp.new1 + hpv.1.L1+resultat.chlam.tracho.2+resultat.HHV2.t1t2.1

# el 1045 tiene tita: 1.303446e-07 y apal: 0.781988326, individuo 1045 estrato 537

epic_cx_cas_control_cis[1045,]$mar.stat.new



kkkk <- data.frame(titaaux)










# Standardized Pearson residual

res <- numeric(n)

for (i in 1:n){
  res[i] <- (yvec[i]-tita[i])/((tita[i]*(1-apal[i]))^(1/2))
}

res

# Diagnósticos de falta de ajuste y influyentes

falta1 <- res^2


diagnos <- matrix(NA, nrow=n, ncol = 4) #numeric(n)


influ1 <- numeric(n)
for (i in 1:n){
  influ1[i] <- falta1[i]*(apal[i]/(1-apal[i]))
  diagnos[i,1] <- tita[i]
  diagnos[i,2] <- apal[i]
  diagnos[i,3] <- falta1[i]
  diagnos[i,4] <- influ1[i]
  
}

diagnosdf <- data.frame(diagnos)


diagnosf <- matrix(NA, nrow=421, ncol = 2) #numeric(n)


faltatot <- c()
p <-0
for (j in coefestratos){
  p <- p+1
  sumatorio3=0
  for (i in 1:length(sujetos)){
    individuo <- sujetos[i]
    if (epic_cx_cas_control_cis[rownames(epic_cx_cas_control_cis)==individuo,]$caseset==j){
      sumatorio3=sumatorio3+((res[i])^2)
    }
  }
  faltatot <- c(faltatot, sumatorio3)
  diagnosf[p,1] <- faltatot[p]
  diagnosf[p,2] <- j
}

diagnosfdt <- data.frame(diagnosf)


diagnosinf <- matrix(NA, nrow=421, ncol = 2) #numeric(n)



influtot <- c()
kte <- 0
for (j in coefestratos){
  sumatorio4=0
  kte <- kte+1
  for (i in 1:length(sujetos)){
    individuo <- sujetos[i]
    if (epic_cx_cas_control_cis[rownames(epic_cx_cas_control_cis)==individuo,]$caseset==j){
      sumatorio4=sumatorio4+(influ1[i])
    }
  }
  influtot <- c(influtot, sumatorio4)
  diagnosinf[kte,1] <- influtot[kte]
  diagnosinf[kte,2] <- j
}


diagnosinfdt <- data.frame(diagnosinf)



library(Lock5Data)

max(faltatot)

windows()
par( font = 2, font.lab = 4, font.axis = 2, las = 1,
     oma = c(0, 0, 2, 0), mar = c(5, 5, 4, 3))
plot(coefestratos,faltatot,pch=19,col="black", xlab = "Nº de estrato", ylab = "Cambio en Pearson chi-square")
identify(coefestratos, faltatot)
#plot(coefestratos,faltatot,pch=19,col="black", ylim=c(0,20), xlab = "Nº de estrato", ylab = "Cambio en Pearson chi-square")
#identify(coefestratos, faltatot)
#plot(coefestratos,faltatot,pch=19,col="black", ylim=c(0,25), xlab = "Nº de estrato", ylab = "Cambio en Pearson chi-square")
#identify(coefestratos, faltatot)

title("Diagnóstico de falta de ajuste", outer = T, cex.main = 1.5)


####falta de daignostico: 112 88 140 362 371 203 336

# 112 fa:48.2113577 estrato 168 , 0.5419437
###individuo 329 0.02020316, apal: 0.011011388, falta: 48.04647018, influ: 5.349489e-01 control
###individuo 332 0.13119859, apal: 0.046785764, falta: 0.13763809, influ: 6.755568e-03 caso
###individuo 338 0.84859825, apal: 0.00870492, falta: 0.02724938, influ: 2.392867e-04 control

epic_cx_cas_control_cis$caseset==168
epic_cx_cas_control_cis[334,]$cc_cis # lo quitamos por tremenda falta de ajuste
epic_cx_cas_control_cis[335,]$cc_cis
epic_cx_cas_control_cis[336,]$cc_cis



#88 fa: 39.1503595 estrato 140, 0.9995782
###individuo 262 0.89638654, apal: 0.004877023, falta: 0.01203539, influ: 5.898456e-05 caso
###individuo 284 0.03695832, apal: 0.024405099, falta: 25.72221939, influ: 6.755568e-03 control
###individuo 302 0.06665514, apal: 0.025853761, falta: 13.41610475, influ: 3.560623e-01 control

epic_cx_cas_control_cis$caseset==140
epic_cx_cas_control_cis[262,]$cc_cis
epic_cx_cas_control_cis[263,]$cc_cis # lo quuitamos por falta de ajuste mayo que 20
epic_cx_cas_control_cis[264,]$cc_cis

# 140 fa:26.6089416 estrato 237, 0.6061444
###individuo 417 0.8386814, apal: 0.008275575, falta: 0.84567987, influ: 0.0070568869 caso
###individuo 418 0.1114296, apal: 0.037346190, falta: 7.36059470, influ: 0.2855545417 control
###individuo 434 0.0498890, apal: 0.016751957, falta: 18.40266699, influ: 0.3135329715 control

epic_cx_cas_control_cis$caseset==237
epic_cx_cas_control_cis[418,]$cc_cis
epic_cx_cas_control_cis[419,]$cc_cis
epic_cx_cas_control_cis[420,]$cc_cis # podriamos quitar

# 362 fa: 26.9915166 estrato 554 0.5546913
###individuo 1085 0.85130349, apal: 0.006721530, falta: 8.570643e-01, influ: 5.799766e-03 control
###individuo 1153 0.09621727, apal: 0.023860492, falta: 8.696873e+00, influ: 2.125840e-01 caso
###individuo 1187 0.05247924, apal: 0.018921439, falta: 1.743758e+01, influ: 3.363075e-01 control

epic_cx_cas_control_cis$caseset==554
epic_cx_cas_control_cis[1083,]$cc_cis
epic_cx_cas_control_cis[1084,]$cc_cis 
epic_cx_cas_control_cis[1085,]$cc_cis#quitamos control

#371 fa: 26.3694536 estrato 563 0.5367061

###individuo 1099 0.08199965, apal: 0.023814433, falta: 1.052789e+01, influ: 2.568320e-01 caso
###individuo 1120 0.86081456, apal: 0.005694324, falta: 2.263383e-02, influ: 1.296225e-04 control
###individuo 1157 0.05718579, apal: 0.017376865, falta: 1.581893e+01, influ: 2.797445e-01 control

epic_cx_cas_control_cis[1111,]$caseset
epic_cx_cas_control_cis[1111,]$cc_cis #es un caso quitamos todo el estrato
epic_cx_cas_control_cis[1112,]$cc_cis
epic_cx_cas_control_cis[1113,]$cc_cis

#203 fa: 23.9433982 estrato 328 0.223943
###individuo 619 0.17293395, apal: 0.049936960, falta: 0.18202366, influ: 9.567479e-03 control
###individuo 620 0.78661603, apal: 0.013190929, falta: 0.79713092, influ: 1.065545e-02 caso
###individuo 623 0.04045002, apal: 0.008793180, falta: 22.96424367, influ: 2.037201e-01 control

epic_cx_cas_control_cis[609,]$caseset
epic_cx_cas_control_cis[607,]$cc_cis
epic_cx_cas_control_cis[608,]$cc_cis
epic_cx_cas_control_cis[609,]$cc_cis #quitamos control

#336 fa: 20.6837880 estrato 496 0.4444609

###individuo 800 0.52807025, apal: 0.039774715, falta: 0.54994412, influ: 0.0227799358 control
###individuo 893 0.04786686, apal: 0.019448072, falta: 19.31478153, influ: 0.3830855469-02 control
###individuo 1008 0.42406289, apal: 0.045001010, falta: 8.190623e-01, influ: 3.859547e-02 caso

epic_cx_cas_control_cis[1006,]$caseset
epic_cx_cas_control_cis[1006,]$cc_cis
epic_cx_cas_control_cis[1007,]$cc_cis #quitamos control
epic_cx_cas_control_cis[1008,]$cc_cis


############# diagnostico influyente: 60 88 41 45 40 140 

#6 di: NA estrato: 12 fa:NA
###individuo 14 9.605768e-01, apal: 0.000698745, falta: 0.961248437, influ: 6.721380e-04 caso
###individuo 17 3.909492e-08, apal: 1.0000000061, falta: NA, influ: NA control
###individuo 18 3.942319e-02, apal: 0.0170245722, falta: 0.040105979, influ: 6.946126e-04 control
epic_cx_cas_control_cis$caseset==12
epic_cx_cas_control_cis[16,]$cc_cis
epic_cx_cas_control_cis[17,]$cc_cis #ya hemos quitado este individuo
epic_cx_cas_control_cis[18,]$cc_cis

#60 di: 1.29236621 estrato: 110 fa 6.728351
###individuo 177 0.4729737, apal: 0.053391709, falta: 0.499650949, influ: 2.818190e-02 caso
###individuo 179 0.3767978, apal: 0.055363093, falta: 0.398881142, influ: 2.337755e-02 control
###individuo 190 0.1502284, apal: 0.175487543, falta: 5.829819158, influ: 1.240807e+00 control
epic_cx_cas_control_cis$caseset==110
epic_cx_cas_control_cis[178,]$cc_cis
epic_cx_cas_control_cis[179,]$cc_cis
epic_cx_cas_control_cis[180,]$cc_cis #quitamos control

#88 di: 0.253160329 estrato: 132 fa 8.545823
###individuo 235 0.1075267, apal: 0.029705752, falta: 7.63432180, influ: 2.337263e-01 control
###individuo 237 0.5212119, apal: 0.020654964, falta: 0.53220460, influ: 1.122451e-02 control
###individuo 238 0.3712613, apal: 0.021185553, falta: 0.37929694, influ: 8.209539e-03 caso
epic_cx_cas_control_cis$caseset==132
epic_cx_cas_control_cis[238,]$cc_cis
epic_cx_cas_control_cis[239,]$cc_cis #dejamos
epic_cx_cas_control_cis[240,]$cc_cis


#41 di: 0.785433914 estrato: 78 fa 11.8221
###individuo 127 0.05362331, apal: 0.038099390, falta: 0.055747251, influ: 2.208062e-03 control
###individuo 128 0.86364869, apal: 0.014746355, falta: 0.876574979, influ: 1.311975e-02 caso
###individuo 129 0.08272799, apal: 0.066047495, falta: 10.889777422, influ: 7.701061e-01 control
epic_cx_cas_control_cis$caseset==78
epic_cx_cas_control_cis[121,]$cc_cis
epic_cx_cas_control_cis[122,]$cc_cis
epic_cx_cas_control_cis[123,]$cc_cis #quitamos control

#45 di: 0.678566475 estrato: 87 fa 7.667558
###individuo 138 0.2895729, apal: 0.028259150, falta: 0.297993996, influ: 8.665949e-03 control
###individuo 139 0.5861586, apal: 0.028283284, falta: 0.603219647, influ: 1.755762e-02 caso
###individuo 141 0.1242684, apal: 0.087932391, falta: 6.766344323, influ: 6.523429e-01 control
epic_cx_cas_control_cis$caseset==87
epic_cx_cas_control_cis[133,]$cc_cis
epic_cx_cas_control_cis[134,]$cc_cis
epic_cx_cas_control_cis[135,]$cc_cis #quitamos control

#40 di: 0.604378800 estrato: 76 fa 7.245873


###individuo 124 0.2584399, apal: 0.029126115, falta: 0.266193088, influ: 7.985765e-03 control
###individuo 125 0.1300296, apal: 0.083534871, falta: 6.351129484, influ: 5.788990e-01 caso
###individuo 126 0.6115305, apal: 0.027078614, falta: 0.628550810, influ: 1.749400e-02 control
epic_cx_cas_control_cis[120,]$caseset
epic_cx_cas_control_cis[118,]$cc_cis
epic_cx_cas_control_cis[119,]$cc_cis #quitamos todo el estrato
epic_cx_cas_control_cis[120,]$cc_cis

#140 di: 0.606144400 estrato: 237 fa 26.60894
###individuo 417 0.8386814, apal: 0.008275575, falta: 0.84567987, influ: 0.0070568869 caso
###individuo 418 0.1114296, apal: 0.037346190, falta: 7.36059470, influ: 0.2855545417 control
###individuo 434 0.0498890, apal: 0.016751957, falta: 18.40266699, influ: 0.3135329715 control
epic_cx_cas_control_cis$caseset==237
epic_cx_cas_control_cis[418,]$cc_cis
epic_cx_cas_control_cis[419,]$cc_cis
epic_cx_cas_control_cis[420,]$cc_cis #se repite




listfa <- c(17, 1045, 
            334, 263, 420,1085, 1111, 1112, 1113, 609, 1007,
            180, 123,135, 118, 119, 120)

modelocisnew <- clogit(as.numeric(cc_cis)~ Bmi_C_new + Pa_Index_new + mar.stat.new + pill.dur.new + dur.tmq.smk.new + l.school.new2 + n.ftp.new1 + hpv.1.L1+resultat.chlam.tracho.2+resultat.HHV2.t1t2.1 + strata(caseset),epic_cx_cas_control_cis[-c(listfa),])
print(summary(clogit(as.numeric(cc_cis)~ Bmi_C_new + Pa_Index_new + mar.stat.new + pill.dur.new + dur.tmq.smk.new + l.school.new2 + n.ftp.new1 + hpv.1.L1+resultat.chlam.tracho.2+resultat.HHV2.t1t2.1 + strata(caseset),epic_cx_cas_control_cis[-c(listfa),])))

# Al quitar las observaciones con falta de ajuste o influyentes se vuelve a hacer el modelo de regresión logística quitando estas observaciones



windows(width = 12, height = 10)
par(font.lab = 2, font.axis = 2, font = 2)
plot(faltatot ~ coefestratos, xlab = "Probabilidad", ylab = "Apalancamiento", pch = 15, cex = 1, xlim = c(0, 1))

windows()
ggplot(epic_cx_cas_control_cis, aes(x = coefestratos, y = faltatot)) +
  geom_point(
    aes(color = as.factor(cc_cis), shape = as.factor(cc_cis)),
    size = 1.5, 
    alpha = 0.8 # It's nice to add some transparency because there may be overlap.
  ) +
  # Use custom colors
  scale_color_manual(
    values = c("red", "blue")
  )





windows()
par(font = 2, font.lab = 4, font.axis = 2, las = 1,
    oma = c(0, 0, 2, 0), mar= c(5,5,1,5))
plot(coefestratos,influtot,pch=19,col="black", xlab = "Nº de estrato", ylab = "Distancia de Cook")
identify(coefestratos,influtot)
title("Diagnóstico influyente", outer = T, cex.main = 1.5)





