GOF=function(x, tipo, a){

  #' @title Goodness-of-fit-test
  #' @description This function estimates the parameters and performs goodness-of-fit-test for the common density functions and probability mass functions.
  #' @param x Data.
  #' @param tipo Data Type: discreto (discrete) or continuo (continuous).
  #' @param a Number of simulations for the p-value, 4000 for default.
  #' @return Results of parameter estimates and the Goodness of fit tests
  #' @export GOF
  #' @author Israel Hernandez Luna
  #' @references Hernandez, Israel (2018). Solución en R para las pruebas de bondad de ajuste. UNAM.
  #' @examples
  #' x=c(rep(0,109),rep(1,65),rep(2,22),rep(3,3),4) #muestra
  #' GOF(x,"discreto")
  #' GOF(x,"discreto",10000)
  #'


  if (missing(x) || mode(x) != "numeric")
    stop("'x' deberia ser un vector numerico no vacio")
  if (any(!is.finite(x)))
    stop("'x' contiene valores perdidos o infinitos")
  if (!is.character(tipo) || missing(tipo))
    stop("Falta el tipo de dato (continuo o discreto)")
  if (missing(a))
    a=4000
  if(a < 1000)
    stop("El numero de simulaciones para el p-valor debe de ser mayor a 1000")

  a=floor(a)
  fun=tolower(tipo)


  if (fun == "discreto" || fun =="d")
  {
    if (any(x < 0))
      stop("existen valores menores a 0, las distribuciones Binomial,Poisson, Bin Negativa, Geometrica tienen valores positivos")
    else
    {

      #Estadistica W^2 Cramer-von Mises

      W=function(x,pj)
      {
        n=length(x)
        a=sort(unique(x))
        k=length(a)
        ej=pj*n
        oj=numeric(k) #cuantos hay de cada xi
        for( i in 1:k)
        {
          oj[i]=sum(x==a[i])
        }
        zj=cumsum(oj-ej)
        w=sum((zj^2)*pj)/n
        return(w)
      }


      #Estadistica U^2 Watson

      U=function(x,pj)
      {
        n=length(x)
        a=sort(unique(x))
        k=length(a)
        ej=pj*n
        oj=numeric(k)
        for( i in 1:k)
        {
          oj[i]=sum(x==a[i])
        }
        zj=cumsum(oj-ej)
        z=sum(zj*pj)
        u=sum(((zj-z)^2)*pj)/n
        return(u)
      }


      #Estadistica A^2 Anderson-Darling

      A=function(x,pj)
      {
        n=length(x)
        a=sort(unique(x))
        k=length(a)
        ej=pj*n
        hj=cumsum(pj)
        oj=numeric(k)
        for( i in 1:k)
        {
          oj[i]=sum(x==a[i])
        }
        zj=cumsum(oj-ej)
        aa=(zj^2)*pj/(hj*(1-hj))
        aa[aa == 1/0] = 0
        Am=sum(aa)/n
        if(Am == 1/0) return(1) else return(Am)
      }


      parametros=rep(0,7)
      names(parametros)=c("n.Uniforme" ,"size.Binom", "p.Binom" ,"l.Poisson" , "p.Geom", "size.BN", "p.BN")
      Wm=rep(0,5)
      Um=rep(0,5)
      Am=rep(0,5)
      names(Wm)=c("Uniforme", "Binomial" ,"Poisson", "Geometrica", "Binomial.Neg")
      Wp.valor=rep(0,5)
      Up.valor=rep(0,5)
      Ap.valor=rep(0,5)

      n=length(x)
      m=mean(x)
      v=var(x)


      #Distribucion Uniforme

      dis.unif=function(x,N) c(rep(1,length(x)))/N
      N.u=max(x)
      parametros[1]=N.u
      prob0=dis.unif(sort(unique(x)),N.u)
      Wm[1]=W(x,prob0)
      Um[1]=U(x,prob0)
      Am[1]=A(x,prob0)
      xx=unique(x)

      #p.valor
      W0=numeric(length(a))
      U0=numeric(length(a))
      A0=numeric(length(a))
      for(i in 1:a)
      {
        ru=sample(xx,n,replace=TRUE)
        W0[i]=W(ru,dis.unif(sort(unique(ru)),N.u))
        U0[i]=U(ru,dis.unif(sort(unique(ru)),N.u))
        A0[i]=A(ru,dis.unif(sort(unique(ru)),N.u))
      }

      Wp.valor[1]=sum(W0>Wm[1])/a
      Up.valor[1]=sum(U0>Um[1])/a
      Ap.valor[1]=sum(A0>Am[1])/a

      #Distribucion Binomial

      nb=0
      nbb=floor((max(x)^1.95)*(v^.95)/((m^.95)*(max(x)-m)^.95)+.5)
      if (max(x) > nbb) nb=max(x) else nb=nbb
      pb=m/nb
      parametros[2]=nb
      parametros[3]=pb
      prob1=dbinom(sort(unique(x)),nb,pb)
      Wm[2]=W(x,prob1)
      Um[2]=U(x,prob1)
      Am[2]=A(x,prob1)

      #p.valor
      W0=numeric(length(a))
      U0=numeric(length(a))
      A0=numeric(length(a))
      for(i in 1:a)
      {
        rb=rbinom(n,nb,pb)
        W0[i]=W(rb,dbinom(sort(unique(rb)),nb,pb))
        U0[i]=W(rb,dbinom(sort(unique(rb)),nb,pb))
        A0[i]=W(rb,dbinom(sort(unique(rb)),nb,pb))
      }

      Wp.valor[2]=sum(W0>Wm[2])/a
      Up.valor[2]=sum(U0>Um[2])/a
      Ap.valor[2]=sum(A0>Am[2])/a


      #Distribucion Poisson

      l=m
      prob2=dpois(sort(unique(x)),l)
      Wm[3]=W(x,prob2)
      Um[3]=U(x,prob2)
      Am[3]=A(x,prob2)
      parametros[4]=l

      #p.valor
      W0=numeric(length(a))
      U0=numeric(length(a))
      A0=numeric(length(a))
      for(i in 1:a)
      {
        rp=rpois(n,l)
        W0[i]=W(rp,dpois(sort(unique(rp)),l))
        U0[i]=U(rp,dpois(sort(unique(rp)),l))
        A0[i]=A(rp,dpois(sort(unique(rp)),l))
      }

      Wp.valor[3]=sum(W0>Wm[3])/a
      Up.valor[3]=sum(U0>Um[3])/a
      Ap.valor[3]=sum(A0>Am[3])/a


      #Distribucion Geometrica

      pg= 1/(1+mean(x))
      prob3=dgeom(sort(unique(x)),pg)
      Wm[4]=W(x,prob3)
      Um[4]=U(x,prob3)
      Am[4]=A(x,prob3)
      parametros[5]=pg

      #p.valor
      W0=numeric(length(a))
      U0=numeric(length(a))
      A0=numeric(length(a))
      for(i in 1:a)
      {
        rg=rgeom(n,pg)
        W0[i]=W(rg,dgeom(sort(unique(rg)),pg))
        U0[i]=U(rg,dgeom(sort(unique(rg)),pg))
        A0[i]=A(rg,dgeom(sort(unique(rg)),pg))
      }

      Wp.valor[4]=sum(W0>Wm[4])/a
      Up.valor[4]=sum(U0>Um[4])/a
      Ap.valor[4]=sum(A0>Am[4])/a


      #Distribucion Binomial Negativa

      if( v > m)
      {
        size=m^2/(v-m)
        pbn=size/(size+m)
        prob4=dnbinom(sort(unique(x)),size,pbn)
        Wm[5]=W(x,prob4)
        Um[5]=U(x,prob4)
        Am[5]=A(x,prob4)
        parametros[6]=size
        parametros[7]=pbn

        #p.valor
        W0=numeric(length(a))
        U0=numeric(length(a))
        A0=numeric(length(a))
        for(i in 1:a)
        {
          rbn=rnbinom(n,size,pbn)
          W0[i]=W(rbn,dnbinom(sort(unique(rbn)),size,pbn))
          U0[i]=U(rbn,dnbinom(sort(unique(rbn)),size,pbn))
          A0[i]=A(rbn,dnbinom(sort(unique(rbn)),size,pbn))
        }

        Wp.valor[5]=sum(W0>Wm[5])/a
        Up.valor[5]=sum(U0>Um[5])/a
        Ap.valor[5]=sum(A0>Am[5])/a
      }
      else
      {
        Wm[4]=0
        Um[4]=0
        Am[4]=0
        Wp.valor[4]=0
        Up.valor[4]=0
        Ap.valor[4]=0
      }

    }
    return(list(rbind(Wm,Wp.valor,Um,Up.valor,Am,Ap.valor),parametros))
  }

  if(fun == "continuo" || fun == "c")
  {

    #Estadistica Dn de Kolmogorov-Smirnov

    Dn=function(pr)
    {
      pr=sort(pr)
      n=length(pr)
      a=numeric(n)
      aa=numeric(n)
      b=seq(1:n)
      for( i in 1:n)
      {
        a[i]=(b[i]/n)-pr[i]
        aa[i]=pr[i]-((b[i]-1)/n)
      }
      d=max(a)
      dd=max(aa)
      return(max(c(d,dd)))
    }

    #Estadistica Vn de Kuiper

    Vn=function(pr)
    {
      pr=sort(pr)
      n=length(pr)
      a=numeric(n)
      aa=numeric(n)
      b=seq(1:n)
      for( i in 1:n)
      {
        a[i]=(b[i]/n)-pr[i]
        aa[i]=pr[i]-((b[i]-1)/n)
      }
      return(max(a)+max(aa))
    }


    #Estadistica Wn^2 de Cramer-von Mises

    Wn=function(pr)
    {
      pr=sort(pr)
      n=length(pr)
      a=numeric(n)
      b=seq(1:n)
      for( i in 1:n)
      {
        a[i]=(pr[i]-(2*b[i]-1)/(2*n))^2
      }
      return(sum(a)+1/(12*n))
    }


    #Estadistica An^2 Anderson-Darling

    An=function(pr)
    {
      pr=sort(pr)
      n=length(pr)
      a=numeric(n)
      b=seq(1:n)
      for(i in 1:n){
        a[i]=(1/n)*((2*b[i]-1)*log(pr[i])+(2*n+1-2*b[i])*log(1-pr[i]))
      }
      a=a[a!=-Inf]
      n=length(a)
      return(-n-sum(a))
    }


    #Estadistica Un^2 de Watson

    Un=function(pr)
    {
      pr=sort(pr)
      n=length(pr)
      a=numeric(n)
      b=seq(1:n)
      for( i in 1:n)
      {
        a[i]=(pr[i]-(2*b[i]-1)/(2*n))^2
      }
      c=(sum(a)+1/(12*n))
      return(c-n*(mean(pr)-.5)^2)
    }

    parametros=rep(0,17)
    names(parametros)=c("min.U", "max.U", "l.Exp", "mu.Normal", "sd.Normal" , "mu.LogN", "sd.LogN", "a.Gamma" , "b.Gamma", "a.Cauchy" , "b.Cauchy" , "a.Beta" , "b.Beta" , "a.Weibull" , "b.Weibull" , "k.Ji^2" , "k.T" )
    Dm=rep(0,10)
    Vm=rep(0,10)
    Wm=rep(0,10)
    Am=rep(0,10)
    Um=rep(0,10)
    names(Dm)=c("Uniforme", "Exponencial", "Normal" , "Log-Normal" , "Gamma" , "Cauchy" , "Beta", "Weibull" , "Ji^2" , "T" )
    Dp.valor=rep(0,10)
    Vp.valor=rep(0,10)
    Wp.valor=rep(0,10)
    Ap.valor=rep(0,10)
    Up.valor=rep(0,10)

    n=length(x)
    m=mean(x)
    v=var(x)


    #Distribucion Uniforme

    mx=ceiling(max(x))
    mn=floor(min(x))
    prob1=punif(sort(x),mn,mx)
    Dm[1]=Dn(prob1)
    Vm[1]=Vn(prob1)
    Wm[1]=Wn(prob1)
    Am[1]=An(prob1)
    Um[1]=Un(prob1)
    parametros[1]=mn
    parametros[2]=mx

    #p.valor
    D0=numeric(length(a))
    V0=numeric(length(a))
    W0=numeric(length(a))
    A0=numeric(length(a))
    U0=numeric(length(a))
    for(i in 1:a)
    {
      ru=runif(n,mn,mx)
      p=punif(sort(ru),mn,mx)
      D0[i]=Dn(p)
      V0[i]=Vn(p)
      W0[i]=Wn(p)
      A0[i]=An(p)
      U0[i]=Un(p)
    }

    Dp.valor[1]=sum(D0>Dm[1])/a
    Vp.valor[1]=sum(V0>Vm[1])/a
    Wp.valor[1]=sum(W0>Wm[1])/a
    Ap.valor[1]=sum(A0>Am[1])/a
    Up.valor[1]=sum(U0>Um[1])/a

    #Distribucion Exponencial

    if (any(x <= 0))
    {
      Dm[2]=0
      Vm[2]=0
      Wm[2]=0
      Am[2]=0
      Um[2]=0
      parametros[3]=0
      Dp.valor[2]=0
      Vp.valor[2]=0
      Wp.valor[2]=0
      Ap.valor[2]=0
      Up.valor[2]=0
    }
    else
    {
      lambda.e= 1/m
      prob2=pexp(sort(x),lambda.e)
      Dm[2]=Dn(prob2)
      Vm[2]=Vn(prob2)
      Wm[2]=Wn(prob2)
      Am[2]=An(prob2)
      Um[2]=Un(prob2)
      parametros[3]=lambda.e

      #p.valor
      D0=numeric(length(a))
      V0=numeric(length(a))
      W0=numeric(length(a))
      A0=numeric(length(a))
      U0=numeric(length(a))
      for(i in 1:a)
      {
        rex=rexp(n,lambda.e)
        p=pexp(sort(rex),lambda.e)
        D0[i]=Dn(p)
        V0[i]=Vn(p)
        W0[i]=Wn(p)
        A0[i]=An(p)
        U0[i]=Un(p)
      }

      Dp.valor[2]=sum(D0>Dm[2])/a
      Vp.valor[2]=sum(V0>Vm[2])/a
      Wp.valor[2]=sum(W0>Wm[2])/a
      Ap.valor[2]=sum(A0>Am[2])/a
      Up.valor[2]=sum(U0>Um[2])/a
    }

    #Distribucion Normal

    xx=(x-m)^2
    sd0=sqrt((1/n)*sum(xx))
    prob3=pnorm(sort(x),m,sd0)
    Dm[3]=Dn(prob3)
    Vm[3]=Vn(prob3)
    Wm[3]=Wn(prob3)
    Am[3]=An(prob3)
    Um[3]=Un(prob3)
    parametros[4]=m
    parametros[5]=sd0

    #p.valor
    D0=numeric(length(a))
    V0=numeric(length(a))
    W0=numeric(length(a))
    A0=numeric(length(a))
    U0=numeric(length(a))
    for(i in 1:a)
    {
      rn=rnorm(n,m,sd0)
      p=pnorm(sort(rn),m,sd0)
      D0[i]=Dn(p)
      V0[i]=Vn(p)
      W0[i]=Wn(p)
      A0[i]=An(p)
      U0[i]=Un(p)
    }

    Dp.valor[3]=sum(D0>Dm[3])/a
    Vp.valor[3]=sum(V0>Vm[3])/a
    Wp.valor[3]=sum(W0>Wm[3])/a
    Ap.valor[3]=sum(A0>Am[3])/a
    Up.valor[3]=sum(U0>Um[3])/a

    #Distribucion Lognormal

    if (any(x <= 0))
    {
      Dm[4]=0
      Vm[4]=0
      Wm[4]=0
      Am[4]=0
      Um[4]=0
      parametros[6]=0
      parametros[7]=0
      Dp.valor[4]=0
      Vp.valor[4]=0
      Wp.valor[4]=0
      Ap.valor[4]=0
      Up.valor[4]=0
    }
    else
    {
      mu.l=sum(log(x))/n
      lxx=(log(x)-mu.l)^2
      sd.l=sqrt((1/n)*sum(lxx))
      prob4=plnorm(sort(x),mu.l,sd.l)
      Dm[4]=Dn(prob4)
      Vm[4]=Vn(prob4)
      Wm[4]=Wn(prob4)
      Am[4]=An(prob4)
      Um[4]=Un(prob4)
      parametros[6]=mu.l
      parametros[7]=sd.l

      #p.valor
      D0=numeric(length(a))
      V0=numeric(length(a))
      W0=numeric(length(a))
      A0=numeric(length(a))
      U0=numeric(length(a))
      for(i in 1:a)
      {
        rl=rlnorm(n,mu.l,sd.l)
        p=plnorm(sort(rl),mu.l,sd.l)
        D0[i]=Dn(p)
        V0[i]=Vn(p)
        W0[i]=Wn(p)
        A0[i]=An(p)
        U0[i]=Un(p)
      }

      Dp.valor[4]=sum(D0>Dm[4])/a
      Vp.valor[4]=sum(V0>Vm[4])/a
      Wp.valor[4]=sum(W0>Wm[4])/a
      Ap.valor[4]=sum(A0>Am[4])/a
      Up.valor[4]=sum(U0>Um[4])/a
    }

    #Distribucion Gamma

    if (any(x <= 0))
    {
      Dm[5]=0
      Vm[5]=0
      Wm[5]=0
      Am[5]=0
      Um[5]=0
      parametros[8]=0
      parametros[9]=0
      Dp.valor[5]=0
      Vp.valor[5]=0
      Wp.valor[5]=0
      Ap.valor[5]=0
      Up.valor[5]=0
    }
    else
    {

      Newton=function(x,a0){
        m=mean(x)
        ml=mean(log(x))
        a1=a0*(1+(digamma(a0)+log(m/a0)-ml)/(1-a0*trigamma(a0)))
      }

      alpha=function(x){
        a00=(mean(x)^2)/var(x)
        error=1
        while(error>0.0001){
          a1=Newton(x,a00)
          error=abs(a1-a00)
          a00=a1
        }
        a1
      }

      a.g=alpha(x)
      b.g=a.g/m
      prob5=pgamma(sort(x),a.g,b.g)
      Dm[5]=Dn(prob5)
      Vm[5]=Vn(prob5)
      Wm[5]=Wn(prob5)
      Am[5]=An(prob5)
      Um[5]=Un(prob5)
      parametros[8]=a.g
      parametros[9]=b.g

      #p.valor
      D0=numeric(length(a))
      V0=numeric(length(a))
      W0=numeric(length(a))
      A0=numeric(length(a))
      U0=numeric(length(a))
      for(i in 1:a)
      {
        rg=rgamma(n,a.g,b.g)
        p=pgamma(sort(rg),a.g,b.g)
        D0[i]=Dn(p)
        V0[i]=Vn(p)
        W0[i]=Wn(p)
        A0[i]=An(p)
        U0[i]=Un(p)
      }

      Dp.valor[5]=sum(D0>Dm[5])/a
      Vp.valor[5]=sum(V0>Vm[5])/a
      Wp.valor[5]=sum(W0>Wm[5])/a
      Ap.valor[5]=sum(A0>Am[5])/a
      Up.valor[5]=sum(U0>Um[5])/a
    }

    #Distribucion Cauchy

    a.c=median(x)
    b.c=IQR(x)/2
    prob6=pcauchy(sort(x),a.c,b.c)
    Dm[6]=Dn(prob6)
    Vm[6]=Vn(prob6)
    Wm[6]=Wn(prob6)
    Am[6]=An(prob6)
    Um[6]=Un(prob6)
    parametros[10]=a.c
    parametros[11]=b.c

    #p.valor
    D0=numeric(length(a))
    V0=numeric(length(a))
    W0=numeric(length(a))
    A0=numeric(length(a))
    U0=numeric(length(a))
    for(i in 1:5000)
    {
      rc=rcauchy(n,a.c,b.c)
      p=pcauchy(sort(rc),a.c,b.c)
      D0[i]=Dn(p)
      V0[i]=Vn(p)
      W0[i]=Wn(p)
      A0[i]=An(p)
      U0[i]=Un(p)
    }

    Dp.valor[6]=sum(D0>Dm[6])/a
    Vp.valor[6]=sum(V0>Vm[6])/a
    Wp.valor[6]=sum(W0>Wm[6])/a
    Ap.valor[6]=sum(A0>Am[6])/a
    Up.valor[6]=sum(U0>Um[6])/a

    #Distribucion Beta

    #Beta

    if (any(x > 1) || any(x < 0))
    {
      Dm[7]=0
      Vm[7]=0
      Wm[7]=0
      Am[7]=0
      Um[7]=0
      parametros[12]=0
      parametros[13]=0
      Dp.valor[7]=0
      Vp.valor[7]=0
      Wp.valor[7]=0
      Ap.valor[7]=0
      Up.valor[7]=0
    }
    else
    {
      a.beta=-(m*v + m^3 - m^2)/v
      b.beta=((m-1)*v + m^3 - 2*m^2 + m)/v
      prob7=pbeta(sort(x),a.beta,b.beta)
      Dm[7]=Dn(prob7)
      Vm[7]=Vn(prob7)
      Wm[7]=Wn(prob7)
      Am[7]=An(prob7)
      Um[7]=Un(prob7)
      parametros[12]=a.beta
      parametros[13]=b.beta

      #p.valor
      D0=numeric(length(a))
      V0=numeric(length(a))
      W0=numeric(length(a))
      A0=numeric(length(a))
      U0=numeric(length(a))
      for(i in 1:a)
      {
        rb=rbeta(n,a.beta,b.beta)
        p=pbeta(sort(rb),a.beta,b.beta)
        D0[i]=Dn(p)
        V0[i]=Vn(p)
        W0[i]=Wn(p)
        A0[i]=An(p)
        U0[i]=Un(p)
      }

      Dp.valor[7]=sum(D0>Dm[7])/a
      Vp.valor[7]=sum(V0>Vm[7])/a
      Wp.valor[7]=sum(W0>Wm[7])/a
      Ap.valor[7]=sum(A0>Am[7])/a
      Up.valor[7]=sum(U0>Um[7])/a
    }

    #Distribucion Weibull

    Newton.w=function(x,a0){
      n=length(x)
      a1= a0 - ( sum((x^a0)*log(x))/sum(x^a0) - 1/a0 - sum(log(x))/n )/
        ( sum((x^a0)*(log(x)^2))/sum(x^a0) - ((sum((x^a0)*log(x)))^2)/((sum(x^a0))^2) +1/a0^2 )
    }

    a.w=function(x){
      a00=median(x)
      error=1
      while(error>0.0001){
        a1=Newton.w(x,a00)
        error=abs(a1-a00)
        a00=a1
      }
      a1
    }


    if ( any(x < 0) || is.na(Newton.w(x,median(x)))==TRUE)
    {
      Dm[8]=0
      Vm[8]=0
      Wm[8]=0
      Am[8]=0
      Um[8]=0
      parametros[14]=0
      parametros[15]=0
      Dp.valor[8]=0
      Vp.valor[8]=0
      Wp.valor[8]=0
      Ap.valor[8]=0
      Up.valor[8]=0
    }
    else
    {
      a.we=a.w(x)
      b.we=(sum(x^a.we)/n)^(1/a.we)
      prob8=pweibull(sort(x),a.we,b.we)
      Dm[8]=Dn(prob8)
      Vm[8]=Vn(prob8)
      Wm[8]=Wn(prob8)
      Am[8]=An(prob8)
      Um[8]=Un(prob8)
      parametros[14]=a.we
      parametros[15]=b.we

      #p.valor
      D0=numeric(length(a))
      V0=numeric(length(a))
      W0=numeric(length(a))
      A0=numeric(length(a))
      U0=numeric(length(a))
      for(i in 1:a)
      {
        rw=rweibull(n,a.we,b.we)
        p=pweibull(sort(rw),a.we,b.we)
        D0[i]=Dn(p)
        V0[i]=Vn(p)
        W0[i]=Wn(p)
        A0[i]=An(p)
        U0[i]=Un(p)
      }

      Dp.valor[8]=sum(D0>Dm[8])/a
      Vp.valor[8]=sum(V0>Vm[8])/a
      Wp.valor[8]=sum(W0>Wm[8])/a
      Ap.valor[8]=sum(A0>Am[8])/a
      Up.valor[8]=sum(U0>Um[8])/a
    }


    #Distribucion Ji-cuadrada

    if ( any(x < 0) )
    {
      Dm[9]=0
      Vm[9]=0
      Wm[9]=0
      Am[9]=0
      Um[9]=0
      parametros[16]=0
      Dp.valor[9]=0
      Vp.valor[9]=0
      Wp.valor[9]=0
      Ap.valor[9]=0
      Up.valor[9]=0
    }
    else
    {

      k=floor(mean(x)+.5)
      prob9=pchisq(sort(x),k)
      Dm[9]=Dn(prob9)
      Vm[9]=Vn(prob9)
      Wm[9]=Wn(prob9)
      Am[9]=An(prob9)
      Um[9]=Un(prob9)
      parametros[16]=k

      #p.valor
      D0=numeric(length(a))
      V0=numeric(length(a))
      W0=numeric(length(a))
      A0=numeric(length(a))
      U0=numeric(length(a))
      for(i in 1:a)
      {
        rq=rchisq(n,k)
        p=pchisq(sort(rq),k)
        D0[i]=Dn(p)
        V0[i]=Vn(p)
        W0[i]=Wn(p)
        A0[i]=An(p)
        U0[i]=Un(p)
      }

      Dp.valor[9]=sum(D0>Dm[9])/a
      Vp.valor[9]=sum(V0>Vm[9])/a
      Wp.valor[9]=sum(W0>Wm[9])/a
      Ap.valor[9]=sum(A0>Am[9])/a
      Up.valor[9]=sum(U0>Um[9])/a
    }


    #Distribucion t

    k.t=floor(2*v/(v-1)+.5)
    if (k.t<=.5)
    {
      Dm[9]=0
      Vm[9]=0
      Wm[9]=0
      Am[9]=0
      Um[9]=0
      parametros[16]=0
      Dp.valor[9]=0
      Vp.valor[9]=0
      Wp.valor[9]=0
      Ap.valor[9]=0
      Up.valor[9]=0
    }
    else
    {
      prob10=pt(sort(x),k.t)
      Dm[10]=Dn(prob10)
      Vm[10]=Vn(prob10)
      Wm[10]=Wn(prob10)
      Am[10]=An(prob10)
      Um[10]=Un(prob10)
      parametros[17]=k.t

      #p.valor
      D0=numeric(length(a))
      V0=numeric(length(a))
      W0=numeric(length(a))
      A0=numeric(length(a))
      U0=numeric(length(a))
      for(i in 1:a)
      {
        r.t=rt(n,k.t)
        p=pt(sort(r.t),k.t)
        D0[i]=Dn(p)
        V0[i]=Vn(p)
        W0[i]=Wn(p)
        A0[i]=An(p)
        U0[i]=Un(p)
      }

      Dp.valor[10]=sum(D0>Dm[10])/a
      Vp.valor[10]=sum(V0>Vm[10])/a
      Wp.valor[10]=sum(W0>Wm[10])/a
      Ap.valor[10]=sum(A0>Am[10])/a
      Up.valor[10]=sum(U0>Um[10])/a
    }


    return(list(rbind(Dm,Dp.valor, Vm,Vp.valor, Wm,Wp.valor, Am,Ap.valor, Um,Up.valor ),parametros))

  }

  else {stop("en tipo debe de ir 'continuo' o 'discreto'")}

}
