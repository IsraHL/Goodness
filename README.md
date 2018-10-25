# Goodness

## ¿Qué hace?

El paquete Goodness incluye la función GOF que nos brinda una solución a las pruebas de bondad de ajuste de manera compuesta. 
Practicamente estamos asumiendo que no conocemos nada acerca de nuestra distribución.

## Instalación


```[R project]
if (!require(devtools)) {
    install.packages("devtools")
}
devtools::install_github('IsraHL/Goodness')
```

## Ejemplos

Los siguientes datos fueron obtenido de Mood(1974) y representan 37 nacimientos consecutivos observados:

```
X=c(1142,84,1125,906,308,922,745,679,1388,505,546,386,
349,969,833,991,236,842,957,1004,392,1186,645,492,
1426,460,26,760,148,1425,520,607,182,857,810,606,775)
```
Ocupamos la GOF de la siguiene manera:

```
GOF(X,"continuo",10000)
```
Finalmente los resultados del programa GOF ocupando 10000 muestras simuladas para el cálculo del p-valor,
cuyo tiempo estimado de ejecución es de 3 minutos con 53 segundos nos índica que la muestra puede ser representada por una distribución Uniforme
con parámetros a = 26 y b = 1426 de haber elegido alguna de las estadísticas diferentes a la estadística Anderson-Darling. 
Los datos podrían ser representados como una distribución Log-Normal con parámetros mu = 6.3428 y  sigma = 0.81486, 
de haber elegido las estadíssticas Kolmogorov-Smirnov, Cramer-von Mises y Anderson-Darling. 
De haber asignado a cualquiera de las estadísticas como nuestro indicador los datos podrían ser representados mediante una distribución Normal con parámetos mu = 709;027 y sigmna = 365;270, 
como una distribución Cauchy con parámetros a = 745 y b = 248,5 o como una distribución Gamma con parámetros a = 2;415 y b = 0;0034.


Graficando en R tenemos

```
hist(X,xlab="Minutos del dia",ylab="Densidad", freq=F, main=""
,col="gray95", ylim=c(0,.0013))
lines(0:1500, dnorm(0:1500,709.027,365.270), lwd=2.5, col= "cyan3", lty=3)
lines(0:1500, dcauchy(0:1500,745,248.5), lwd=2.5, col= "plum2", lty=4)
lines(0:1500, dlnorm(0:1500,6.3428,0.81486), lwd=2.5, col= "khaki3" , lty=5)
lines(0:1500, dgamma(0:1500,2.4151,0.0034062), lwd=2.5, col= "red4" , lty=6)
legend(1100,0.0012,c("Normal","Cauchy", "Log-Normal","Gamma"),
col=c("cyan3", "plum2", "khaki3","red4"), lty= c(3,4,5,6), lwd=2.5, cex = .8,
bty="o",bg="grey96", box.lty=.1)

```
![Nacimientos](URL https://github.com/IsraHL/Goodness/blob/master/NNN3.jpeg)

