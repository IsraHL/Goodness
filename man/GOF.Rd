\name{GOF}
\alias{GOF}
\title{Goodness-of-fit-test}
\usage{
GOF(x , tipo, a )
}
\description{
GOF tiene como objetivo allar la función de distribución (caso continuo) o función masa de probabilidad (caso discreto) que más se asemeje a nuestra muestra. Dando como resultado los diferentes parámetros de las funciones y brindando el
calculo de las estadísticas Kolmogorov-Smirnov, Kuiper, Cramer-von Mises, Anderson-Darling y Watson en el caso continuo o las estadísticas Cramer-von Mises, Anderson-Darling o Watson para el caso discreto.
}
\arguments{
    \item{x}{numerico, muestra representativa.}
    \item{tipo}{character, "discreto", "continuo"}
    \item{a}{numeric, número de simulaciones para el p-valor}
}
\author{Israel Hernández Luna}
\references{Hernandez, Israel (2018). Solución en R para las pruebas de bondad de ajuste. UNAM.}
\examples{
x=c(rep(0,109),rep(1,65),rep(2,22),rep(3,3),4) #muestra
GOF(x,"discreto")
GOF(x,"discreto",10000)
}
