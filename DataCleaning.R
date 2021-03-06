workingDir <- "C:/Users/usuario/Desktop/Curso13-14/02Semestre/DataMining"
dataDir <- file.path(workingDir, "Datos")
resultsDir <- file.path(workingDir, "Resultados")
setwd(workingDir)
misDatos <- read.csv("DatosClinicos.csv")

#importamos el fichero de marcadores (csv2)
markers <- read.csv2(file=paste("Markers.csv"))
#unimos las columnas de markers que nos interesan a los datos clinicos anteriores
misDatos <- cbind(misDatos, markers[,2:5])

#para no tener que especificar m�s en qu� tabla estamos trabajando podemos usar:
#attach(misDatos)

# INFORMACI�N DE LOS CAMPOS DEL FICHERO:

#grupo A: no recidiva
#grupo B: recidva temprana
#grupo C: recidiva tardia

#estadog: estado global -> SG(supervivencia global) -> vivo o muerto
#sgd: tiempo de seguimiento (d�as)
#sgm: tiempo de seguimiento (meses) -> de �ltimo control (de muerte en caso de estadog=1)

#recid: tiempo de recidiva -> SLE(supervivencia libre de enfermedad) ->sano o recidiva
#iled: tiempo de recidiva (d�as)
#ilem: tiempo de recidiva (meses)

#DATA CLEANING:

#analizamos como se han guardado los datos 
#(en este caso la variable tama�o -tam- se ha guardado como Factor cuando deber�a 
#ser numerico continuo)
class(misDatos)
str(misDatos)

#para plotear todas las variables cruzadas de la tabla (visi�n general):
#pairs(misDatos)

#el valor DESCONOCIDO en el campo tam es el que provoca que se hayan guardado 
#los datos como Factor
#para evitar esto, los valores perdidos en R se etiquetan como NA
#una de las opciones que hacer con los valores perdidos es la imputaci�n (p.ej. poner la media)
summary(misDatos)

#buscamos los registros que contengan el valor DESCONOCIDO en el campo tam
misDatos$tam
#es lo mismo que:
# misDatos[,12]
# misDatos[,"tam"]

#con la funci�n levels se recuperan todos los valores posibles de una variable tipo Factor
levels(misDatos$tam)

#utilizamos una variable auxiliar para calcular el valor medio de tama�o
tamAux<-misDatos$tam
tamAux<-as.numeric(as.character(tamAux))
meanTam<-summary(tamAux)[4]
#accedemos a las posiciones en las que se encuentra DESCONOCIDO y las sustituimos por NA
misDatos$tam[misDatos$tam=="DESCONOCIDO"] <- NA 
#cambiamos el tipo de datos del tama�o de factor a numerico
misDatos$tam <- as.numeric(as.character(misDatos$tam))
#sustituimos los NA's por la el valor medio del tama�o
misDatos$tam[is.na(misDatos$tam)] <- meanTam

#str(misDatos)

#ya que solo tenemos un caso de fenotipo Luminal-HER2 y �ste es un subtipo de Luminal B
#lo sustituimos:
misDatos$fenotipo[misDatos$fenotipo=="Luminal-HER2"] <- "Luminal B"
misDatos$fenotipo

#como soluci�n provisional a los casos sin fenotipo, los incluimos en el tipo Luminal B
misDatos$fenotipo[misDatos$fenotipo==""] <- "Luminal B"
summary(misDatos$fenotipo)


misDatos$fenotipo<-factor(misDatos$fenotipo, levels=c("Basal like", "HER2 enriched", "Luminal A", "Luminal B"))

# QUEREMOS SABER DE QU� GRUPO SON LAS PACIENTES QUE HAN RECIDIVADO
# MOSTRAMOS UNA SUBMATRIZ CON GRUPO, ESTADOG, ILEM 
#(importate el uso del comando 'c' -concatenaci�n-)

misDatos[misDatos$recid==1,c("grupo","estadog","ilem")]

#para renombrar las columas utilizamos:
#colnames(misDatos) <- c("nombre1","nombre2","nombre3",...)

# Modelo de regresi�n log�stica para estimar un modelo predictivo de recidiva de pacientes 
#en base a diferentes factores pron�stico de la enfermedad: 
#tama�o del tumor, n�mero de ganglios, edad del paciente y fenotipo

regLog <- glm(recid~tam+ngang+edad+fenotipo, data=misDatos, family=binomial)
summary(regLog)
#plot(regLog)
pred <- predict(regLog, misDatos, type="response")

# utilizamos el umbral 0.5 para asignar los valores 0 o 1 a la prediccion
pred[pred<0.5]<-0
pred[pred>=0.5]<-1

pred<-factor(pred, levels=c("0", "1"))
#obtenemos la matriz de confusion y calculamos la accuracy (acc)
confusionMatrix <- table(pred,misDatos$recid)
hits <- confusionMatrix[1,1]+confusionMatrix[2,2]
total <- sum(confusionMatrix)
acc <- hits/total

acc

#obtenemos ACC=0.56

#Esquema de validaci�n hold-out (60% training, 40% test) para calcular el 
#porcentaje de clasificaci�n correcta (ACC, accuracy) promedio en training y test 
#(generalizaci�n). Estudiar c�mo var�a el ACC para test en funci�n del n�mero de 
#repeticiones del hold-out.


#creamos una funci�n para el hold-out
#argumentos: porcentage de training, n�mero de repeticiones

holdOut <- function(prc.train, repNum) {
  holdOut.acc<-NULL
  for(i in 1:repNum) {
  
    indt <- sample(total, total*prc.train/100)
    dat.train <- misDatos[indt,] #conjunto de training
    dat.test <- misDatos[-indt,] #conjunto de test
    
    rl.train <- glm(recid~tam+ngang+edad+fenotipo, data=dat.train, family=binomial, control=glm.control(epsilon=1e-8, maxit=100, trace=FALSE))
    rl.train$xlevels[["fenotipo"]] <- union(rl.train$xlevels[["fenotipo"]], levels(dat.test$fenotipo))
    pred.test <- predict(rl.train, dat.test, type="response")
    
    # utilizamos el umbral 0.5 para asignar los valores 0 o 1 a la prediccion
    pred.test[pred.test<0.5]<-0
    pred.test[pred.test>=0.5]<-1
    pred.test<-factor(pred.test, levels= c("0","1"))
    
    #calculamos la matriz de confusion y la accuracy
    confmatrix<-table(pred.test, dat.test$recid)
    hits<-confmatrix[1,1]+confmatrix[2,2]
    total<-sum(confmatrix)
    accuracy<-hits/total
    holdOut.acc<-c(holdOut.acc,accuracy)
  }
  
  #calculamos la accuracy promedio (generalizaci�n)
  accG<-mean(holdOut.acc,na.rm=TRUE)
  return (accG)
}

#probamos a llamar a la funci�n con 60% y 10 repeticiones
holdOut(60,10)

#para estudiar la variaci�n del acc en funci�n del n�mero de repeticiones
#llamamos a la funci�n holdOut en bucle (de 10 a 100 repetiones aumentando de 10 en 10)
var.accG<-NULL
for(i in seq(5,100,by=5)){
var.accG <- c(var.accG,holdOut(60,i))
}
var.accG
hist(var.accG)




#Calcular el ACC promedio siguiendo un esquema de validaci�n cruzada 
#(10-fold cross validation) y comparar los resultados de forma gr�fica 
#con los obtenidos utilizando el esquema de hold-out.

#para el kfold: reordenar primero todos los datos con un sample de todos los datos
#despu�s dividir en 10 cajas e ir utilizando cada una como test en una iteraci�n
Kfold <- function(k) {
  fold.acc<-NULL
  dat.merge <- sample(total,total)
  fold.indt <- dat.merge[1:(total/k)]
  fold.train <- misDatos[fold.indt,]
  fold.test <- misDatos[-fold.indt,]
  fold.rl <- glm(recid~tam+ngang+edad+fenotipo, data=fold.train, family=binomial)
  fold.pred <- predict(fold.rl, fold.test, type="response")
  fold.pred[fold.pred<0.5]<-0
  fold.pred[fold.pred>=0.5]<-1
  confusionMatrix <- confusionMatrix(fold.test, fold.test$recid)
  accuracy <- confusionMatrix$overall[1]
  fold.acc<-c(fold.acc,accuracy)
  for(i in 1:(k-1)) {
    fold.indt <- dat.merge[((total/k)*i)+1:(total/k)*(i+1)]
    fold.train <- misDatos[fold.indt,]
    fold.test <- misDatos[-fold.indt,]
    fold.rl <- glm(recid~tam+ngang+edad+fenotipo, data=fold.train, family=binomial)
    fold.pred <- predict(fold.rl, fold.test, type="response")
    fold.pred[fold.pred<0.5]<-0
    fold.pred[fold.pred>=0.5]<-1
    tabla<-table(fold.pred, fold.test$recid)
    hits<-(tabla[1,1]+tabla[2,2])
    total<-sum(tabla)
    accuracy<-hits/total
    fold.acc<-c(fold.acc,accuracy)    
  }
  fold.accG<-mean(fold.acc,na.rm=TRUE)  
  return (fold.accG)  
}
Kfold(10)


fold.acc<-NULL
dat.merge <- sample(total,total)
fold.indt <- dat.merge[1:(total/10)]
fold.train <- misDatos[fold.indt,]
fold.test <- misDatos[-fold.indt,]
fold.rl <- glm(recid~tam+ngang+edad+fenotipo, data=fold.train, family=binomial)
fold.rl$xlevels[["fenotipo"]] <- union(fold.rl$xlevels[["fenotipo"]], levels(fold.test$fenotipo))
fold.pred <- predict(fold.rl, fold.test, type="response")
fold.pred[fold.pred<0.5]<-0
fold.pred[fold.pred>=0.5]<-1
confusionMatrix <- confusionMatrix(fold.test, fold.test$recid)
accuracy <- confusionMatrix$overall[1]
fold.acc<-c(fold.acc,accuracy)
for(i in 1:(k-1)) {
  fold.indt <- dat.merge[((total/k)*i)+1:(total/k)*(i+1)]
  fold.train <- misDatos[fold.indt,]
  fold.test <- misDatos[-fold.indt,]
  fold.rl <- glm(recid~tam+ngang+edad+fenotipo, data=fold.train, family=binomial)
  fold.pred <- predict(fold.rl, fold.test, type="response")
  fold.pred[fold.pred<0.5]<-0
  fold.pred[fold.pred>=0.5]<-1
  tabla<-table(fold.pred, fold.test$recid)
  hits<-(tabla[1,1]+tabla[2,2])
  total<-sum(tabla)
  accuracy<-hits/total
  fold.acc<-c(fold.acc,accuracy)    
}
fold.accG<-mean(fold.acc,na.rm=TRUE)  

#�rea bajo la curva ROC:
#en lugar de aplicar un umbral para la predicci�n de 0.5, 
#ahora cogemos todos los valores obtenidos al hacer la predcci�n, los ordenamos y los
#vamos utilizando como umbral cada uno de ellos para:
#1. calcular la sensibilidad y la especificidad
#2. graficar cada punto en una gr�fica de sensibilidad frente a (1-especificidad)
#3. unir los puntos de la gr�fica en una curva
#El valor del �rea bajo la curva (entre 0 y 1) indica la eficacia del modelo. 
#R
#library(pROC)
#utilizar: obj.roc <- roc(response, predictor, smooth=FALSE, auc=TRUE,...)
#response=vector de recidiva real
#predictor=vector de predicci�n de recidiva
#para calcular el �rea bajo la curva (AUC):
#auc(obj.roc,...) � auc(response, predictor,...)

# M�S ADELANTE...:
#se calcula la media de las AUC de todas las predicciones
#se calcula la curva ROC media --> paquete Daim


obj.roc <- roc(misDatos$recid, pred, smooth=FALSE, auc=TRUE)
auc(obj.roc)
plot(obj.roc)

# REDES NEURONALES: (MLP)
#Se refuerzan o se debilitan las conexiones entre neuronas dependiendo de la experiencia
#LTD: depresi�n sin�ptica
#LTP: potenciaci�n sin�ptica
#la funci�n de cada neurona ser�a el sumatorio de los impulsos que le llegan a trav�s de
#las dendritas multiplicada cada una por el grado de importancia de esa conexi�n
#la salida final ser�a un logaritmo (entre 0 y 1)
#las funciones de las neuronas intermedias son tangentes hiperb�licas (entre -1 y 1)
#Algoritmo de retropropagaci�n de errores: si la salida no es la deseada se vuelve hacia
#atr�s modificando las conexiones para que la pr�xima vez, esa misma entrada genere una
#m�s parecida a la deseada.
#Es importante ajustar bien la constante de aprendizaje para que no oscile demasiado pero
#tampoco tarde demasiado en aprender.
# R: 
#library(nnet)
#funci�n de entrenamiento:
#nn.fit <- nnet(y~, data=x, size=5, entropy= TRUE, maxit=1000, decay=5e-4)
#param size: n�mero de neuronas en la capa oculta
#param entropy: a TRUE si la salida es binaria
#param maxit: n�mero de iteraciones de cada patr�n
#param decay: (lo dejamos como est� por ahora) para no llegar a sobreentrenamiento
#funci�n de generalizaci�n.
#nn.pred <- predict(nn.fit,x,type="raw")

nn.fit<-nnet(recid~tam+ngang+edad+fenotipo,data=misDatos,size=5,entropy= TRUE,maxit=1000,decay=5e-4)
nn.pred <- predict(nn.fit, misDatos, type="raw")
nn.roc <- roc(misDatos$recid, nn.pred, smooth=FALSE, auc=TRUE)
plot(nn.roc)
auc(nn.roc)

# M�QUINAS DE SOPORTE VECTORIAL: (SVM)
#se trata de encontrar los vectores de soporte, es decir, los patrones de cada tipo 
#m�s cercanos al plano de separaci�n, y maximizar la distancia de �stos al plano
#cuando la separaci�n no puede ser lineal se utiliza una funci�n kernel no lineal
# R:
#library(e1071)
#funci�n de aprendizaje:
#svm.fit <- svm(...)
#funci�n de predicci�n:
#svm.pred <- predict(svm.fit,x,probability=TRUE)

svm.fit <- svm(recid~tam+ngang+edad+fenotipo,data=misDatos,cost=1000,gamma=1,probability=TRUE)
svm.pred <- predict(svm.fit, misDatos, probability=TRUE)
svm.roc <- roc(misDatos$recid, svm.pred, smooth=FALSE, auc=TRUE)
plot(svm.roc)
auc(svm.roc)

# �RBOLES DE DECISI�N: (DT)
# R:
library(rpart)
#funci�n de construcci�n del �rbol:
#dt.fit <- rpart(...)
#param control: controla en tema de la entrop�a
#funci�n de predicci�n:
#dt.pred <- predict(dt.fit, x, type="prob")

dt.fit <- rpart(recid~tam+ngang+edad+fenotipo,data=misDatos,y=TRUE,control=rpart.control(cp=.05))
dt.pred <- predict(dt.fit, misDatos, type="prob")

#crear: funci�n cr(modelo, formula, datos)
#De esta forma se puede utilizar esta funci�n para aplicar la validaci�n cruzada con
#cualquier modelo y conjunto de datos. Al tener la f�rmula (y~) como par�metro, tambi�n se
#pueden elegir las variables que se utilizan.

#Implementar una funci�n holdOut generalizada:
#holdOut(formula, datos, nSim)
#debe devolver un objeto de tipo list con: 
#AUCg(vector) y AUCg promedio
#AUCt(vector) y AUCt promedio
#ACCg(vector) y ACCg promedio
#ACCt(vector) y ACCt promedio
# g se refiere al conjunto de test (generalizaci�n)
# t se refiere al conjunto de training

holdOut <- function(formula, datos, nSim, prc.train) {
  AUCg <- NULL
  #AUCt <- NULL
  ACCg <- NULL
  #ACCt <- NULL
  for(i in 1:nSim) {
    
    indt <- sample(total, total*prc.train/100)
    dat.train <- datos[indt,] #conjunto de training
    dat.test <- datos[-indt,] #conjunto de test
    
    rl.train <- glm(formula, data=dat.train, family=binomial)
    pred.test <- predict(rl.train, dat.test, type="response")
    
    #calculamos la AUC (area bajo la curva ROC)
    obj.auc <- auc(dat.test$recid, pred.test)
    AUCg <- c(AUCg,obj.auc)
    
    # utilizamos el umbral 0.5 para asignar los valores 0 o 1 a la prediccion
    pred.test[pred.test<0.5]<-0
    pred.test[pred.test>=0.5]<-1 
    
    pred.test<-factor(pred.test, levels=c("0","1"))
    
    #calculamos la matriz de confusion y la accuracy
    confMatrix<-confusionMatrix(pred.test, dat.test$recid)
    accuracy <- confMatrix$overall[1]
    ACCg <- c(ACCg,accuracy)
   
  }
  #calculamos la auc promedio (generalizaci�n)
  aucG <- mean(AUCg)
  #calculamos la accuracy promedio (generalizaci�n)
  accG<-mean(ACCg,na.rm=TRUE)
  
  obj.list <- list(AUCg=AUCg, aucG=aucG, ACCg=ACCg, accG=accG)
  return (obj.list)
}

lr.edad <- holdOut(recid~edad, misDatos, 10, 60)
var.edad <- var(lr.edad[[1]])

#si utilizamos la variable "grupo" podemos comprobar claramente que es 100% discriminativa,
#realiza la predicci�n con un 100% de precisi�n, por lo que AUC = 1:
lr.grupo <- holdOut(recid~grupo, misDatos, 10, 60)
var.grupo <- var(lr.grupo[[1]])


#Construir una tabla de contingencia entre la variable "grupo" y la recidiva para luego
#aplicarle un test de fisher y ver si las variables son dependientes:
tablaGrupo <- table (misDatos$grupo, misDatos$recid)
chi2.tablaGrupo <- fisher.test(tablaGrupo)
#hip�tesis nula: las variables son independientes
#resultado: p-valor muy bajo -> se rechaza la hip�tesis nula

var.auc <- NULL
varG <- NULL
for(i in 1:10) {
  for(j in 1:10) {
    lr.edad <- holdOut(recid~edad, misDatos, 5*i, 60)
    var.auc <- c(var.auc, lr.edad[[2]])
  }  
  varG <- c(varG,var(var.auc))
}
varG
plot(varG, type="l")

#Observando la gr�fica podemos ver que cu�nto mayor es el n�mero de repeticiones del holdOut,
#menor es la varianza o desviaci�n est�ndar de la AUC promedio. Es decir, m�s se ajusta la 
#distribuci�n de AUC a una NORMAL.
#La desviaci�n est�ndar de la media muestral es el error est�ndar de la media:
#sd(AUCg-promedio) = var(AUCg)/sqrt(n)
#donde n es el tama�o de la muestra, en este caso el n�mero de repts del holdOut.

# �! Leer en detalle: Libro Seefeld - cap�tulo 13.
#(de alguna forma habr� que sustituir el rexp por la funci�n holdOut)

#Comparar la distribuci�n de las medias para n=5, n=20, n=200, n=2000 
#(un histograma para cada una).
#Seg�n el teorema del l�mite la distribuci�n de las medias a partir de n=30, se ajustar�
#siempre cada vez m�s a una normal sin importar la distribuci�n inicial de los datos.

# DISCRETIZACI�N DE VARIABLES:
#uso del paquete discretization
#funciones CAIM y CCAC para discretizar las variables de edad y tama�o, por ejemplo.
#Discretizar la variable edad en: < 40, 40-55, > 55.
#Discretizar la variable tama�o en: <=2, 2-5, >=5
#Analizar la diferencia entre los resultados obtenidos con las variables con y sin discretizar


