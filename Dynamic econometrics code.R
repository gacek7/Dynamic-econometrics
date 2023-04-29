#usuwa znaczenie zmiennych i rozpoczyna nowa sesje
rm(list = ls())
.rs.restartR()

############ DODATKOWE PAKIETY ###################

if(!require(tseries)){install.packages('tseries')}
if(!require(forecast)){install.packages('forecast')}
if(!require(plotrix)){install.packages('plotrix')}
if(!require(readxl)){install.packages('readxl')}
if(!require(portes)){install.packages('portes')}
if(!require(HDtest)){install.packages('HDtest')}
if(!require(lmtest)){install.packages('lmtest')}
if(!require(MTS)){install.packages('MTS')}
if(!require(MVN)){install.packages('MVN')}
if(!require(MSQC)){install.packages('MSQC')}
if(!require(matrixcalc)){install.packages('matrixcalc')}


library(forecast)#wykorzystane funkcje - seasonaldummy
library(plotrix)
library(readxl)
library(portes)
library(HDtest)
library(lmtest)
library(MTS)
library(MSQC)
library(matrixcalc)

#### WLASNE FUNKCJE ######

estymacjaVAR <- function(dane,k=2,const=0,trend=0,lsez=0,X={}) {
  
  dane = as.matrix(dane)
  y = dane[-c(1:k),]#usuwamy wiersze, od 1 do k
  n = ncol(y)#liczba analizowanych zmiennych
  T_mod = nrow(y)#liczba modelowanych obserwacji
  Z = matrix(0,T_mod,k*n)#'szkielet' do macierzy opoznien
  
  #wypelniamy macierz opoznien
  for (i in 1:k) {
    Z[,((i-1)*n+1):(i*n)] = dane[(k+1-i):(nrow(dane)-i),]
  }
  
  #dolaczamy stala
  if (const == 1) {
    Z = cbind(Z, matrix(1,T_mod,1))
  }
  
  #dolaczamy zmienna 'czas'
  if (trend == 1) {
    Z = cbind(Z, seq(1,T_mod,1))
  }
  
  #sztuczne zmienne sezonowe
  if (lsez > 1) {
    sezonowe = seasonaldummy(ts(1:nrow(dane),freq=lsez,start=1))
    sezonowe = sezonowe[,1:(sez-1)]-1/lsez*matrix(1,nrow(dane),lsez-1)#ortogonalne (wysrodkowane, centered)
    Z = cbind(Z,sezonowe[(k+1):nrow(dane),])
  }
  
  #pozostale zmienne egzogeniczne
  Z = cbind(Z,X)# X - macierz zmiennych egzogenicznych, ktora nalezy wczesniej utworzyc
  
  #estymacja
  hat_A = solve(crossprod(Z,Z))%*%crossprod(Z,y)
  y_teor = Z%*%hat_A
  hat_E = y-y_teor
  tilde_Sigma = crossprod(hat_E,hat_E)/T_mod
  hat_Gamma = crossprod(Z,Z)/T_mod
  
  #logarytm maksymalnej wartosci funkcji wiarygodnosci
  ln_wiary = -T_mod/2*log(det(tilde_Sigma))-n*T_mod/2*(log(2*pi)+1)
  
  out <- list(hat_A = hat_A,
              tilde_Sigma = tilde_Sigma,
              hat_Gamma = hat_Gamma,
              ln_wiary = ln_wiary,
              T_mod = T_mod,
              n = n,
              k = k,
              const = const,
              trend = trend,
              lsez = lsez,
              y = y,
              dane = dane,
              y_teor = y_teor,
              hat_E = hat_E,
              Z = Z,
              X = X)
  return(out)
}

KI <- function(dane,k_max = 2,const = 0,trend = 0,lsez = 0,X = {}){
  kryteria = matrix(0,4,k_max)
  rownames(kryteria) <- c("k","AIC","HQ", "SC")
  for (k in 1:k_max){
    model = estymacjaVAR(dane[-c(1:(k_max-k)),],k,const,trend,lsez,X)
    ln_det_Sig = log(det(model$tilde_Sigma))
    n = model$n
    T_mod = model$T_mod
    kara = 2*(n^2)*k/T_mod
    kryteria[1,k] = k
    kryteria[2,k] = ln_det_Sig+kara #AIC
    kryteria[3,k] = ln_det_Sig+kara*log(log(T_mod)) #HQ
    kryteria[4,k] = ln_det_Sig+kara*log(T_mod) #SC
  }
  wybor = matrix(0,1,3)
  colnames(wybor) <- c("AIC","HQ", "SC")
  wybor[1,1] = which.min(kryteria[2,])
  wybor[1,2] = which.min(kryteria[3,])
  wybor[1,3] = which.min(kryteria[4,])
  out <- list(kryteria = kryteria,
              wybor = wybor)
  return(out)
}


w_wlasne <- function(model,modul = T){
  hat_A = model$hat_A
  n = model$n
  k = model$k
  tilde_A = t(hat_A)
  tilde_A = rbind(tilde_A[,1:(n*k)],cbind(diag(n*k-n),matrix(0,n*k-n,n)))
  ww = eigen(tilde_A,symmetric = F)$value #wartosci wlasne macirzy stowarzyszonej
  if (modul == T) ww = abs(ww) #moduly wartosci wlasnych
  return(ww)
}

LJBtest <- function(model,multivariate = T){
  E = model$hat_E
  E = E - apply(E,2,mean)
  TT = model$T_mod
  Sigma = model$tilde_Sigma
  P = t(chol(Sigma))
  omega = E%*%t(solve(P))
  b1 = 1/TT*apply(omega^3,2,sum)
  b2 = 1/TT*apply(omega^4,2,sum)
  LJBmulti_skos = TT/6*crossprod(b1,b1)
  LJBmulti_kurtoza = TT/24*crossprod(b2-3,b2-3)
  LJBmulti = LJBmulti_skos+LJBmulti_kurtoza
  if (multivariate == T){
    out <- list(LJBmulti = c(LJBmulti,2*model$n,1 - pchisq(LJBmulti,2*model$n)),
                skosnosc = c(LJBmulti_skos,model$n,1 - pchisq(LJBmulti_skos,model$n)),
                kurtoza = c(LJBmulti_kurtoza,model$n,1 - pchisq(LJBmulti_kurtoza,model$n))
    )
  }
  if (multivariate == F){
    wynik = matrix(0,model$n,6)
    rownames(wynik) <- colnames(model$y)
    colnames(wynik) <- c("LJB","p_value","skosnosc","p_value","kurtoza","p_value")
    for (i in 1:model$n){
      omega = as.matrix(E[,i])/(Sigma[i,i])^0.5
      b1 = 1/TT*apply(omega^3,2,sum)
      b2 = 1/TT*apply(omega^4,2,sum)
      ls = TT/6*crossprod(b1,b1)
      lk = TT/24*crossprod(b2-3,b2-3)
      LJB = ls+lk
      wynik[i,] = c(LJB,1-pchisq(LJB,2),ls,1-pchisq(ls,1),lk,1-pchisq(lk,1))
    }
    out <- list(LJBmulti = c(LJBmulti,2*model$n,1 - pchisq(LJBmulti,2*model$n)),
                skosnosc = c(LJBmulti_skos,model$n,1 - pchisq(LJBmulti_skos,model$n)),
                kurtoza = c(LJBmulti_kurtoza,model$n,1 - pchisq(LJBmulti_kurtoza,model$n)),
                jednowymiarowe = wynik)
  }
  return(out)
} 

testWalda <- function(model,mC, c){
  E = model$hat_E
  A = model$hat_A
  Z = model$Z
  Sigma = 1/(model$T_mod-nrow(A))*crossprod(E,E)
  q = nrow(c)
  #funkcja 'vec' pochodzi z pakietu 'matrixcalc'
  lambda_w = t(mC%*%vec(t(A))-c)%*%solve(mC%*%kronecker(solve(crossprod(Z,Z)),Sigma)%*%t(mC))%*%(mC%*%vec(t(A))-c)
  lambda_F = lambda_w/q
  out <- list(chi2_test = paste("lambda_W = ",lambda_w, ", df = ",q, ", p-value = ", pchisq(lambda_w,q,lower.tail = F)),
              F_test = paste("lambda_F = ",lambda_F, ", df1 = ",q, ", df2 = ",model$T_mod*model$n-nrow(vec(A)), ", p-value = ", pf(lambda_F,q,model$T_mod*model$n-nrow(vec(A)),lower.tail = F)))
  return(out)
}

mojeIRF <-function(model,h = 10,CIRF = F){
  n = model$n
  k = model$k
  irf = matrix(0,h+1,n^2)
  rownames(irf) = seq(0,h,1)
  zmienne = colnames(model$hat_A)
  nazwy = {}
  for (i in 1:n){
    for (j in 1:n)
      nazwy = c(nazwy,paste("w", i, "->", zmienne[j]))
  }
  colnames(irf) = nazwy
  tildeA = t(model$hat_A[1:(n*k),])
  tildeA = rbind(tildeA,cbind(diag(n*(k-1)),matrix(0,n*(k-1),n)))
  P = t(chol(model$tilde_Sigma))
  phi = diag(n*k)
  irf[1,] = t(vec((phi[1:n,1:n])%*%P))
  for (i in 1:h){
    phi = phi%*%tildeA
    irf[i+1,] = t(vec((phi[1:n,1:n])%*%P))
  }
  if (CIRF==T) irf = apply(irf,2,cumsum)
  return(irf)
}

mojeFEVD <-function(model,h = 10){
  n = model$n
  zmienne = colnames(model$hat_A)
  Phi2 = apply(mojeIRF(model,h-1,F)^2,2,cumsum)
  FEVD = matrix(0,h,n^2)
  rownames(FEVD) = seq(1,h,1)
  nazwy = {}
  for (i in 1:n){
    pomoc = {}
    for (j in 1:n){
      pomoc = cbind(pomoc,Phi2[,(i+(j-1)*n):(i+(j-1)*n)])
      nazwy = c(nazwy,paste("w", j, "->", zmienne[i]))
    }
    mianownik = as.matrix(apply(pomoc,1,sum))
    FEVD[,((i-1)*n+1):(i*n)] = apply(pomoc, 2, '/',mianownik)
  }
  colnames(FEVD) = nazwy
  return(FEVD)
}

procJohansena <-function(dane, k = 2, ect_const = 0, ect_trend = 0, const = 0, trend = 0, lsez = 0, X = {}){
  dane = as.matrix(dane)
  n = ncol(dane)#liczba analizowanych zmiennych
  przyrosty = diff(dane)
  
  Z0 = przyrosty[-c(1:(k-1)),]#usuwamy wiersze, od 1 do k-1
  T_mod = nrow(Z0)#liczba modelowanych obserwacji
  #opoznione poziomy
  Z1 = dane[-c(1:(k-1)),]#usuwamy wiersze, od 1 do k-1
  Z1 = Z1[-nrow(Z1),]#usuwamy ostatni wiersz
  #skladowe determnistyczne w relacji kointegrujacej
  #stala
  if (ect_const == 1) {
    Z1 = cbind(Z1, matrix(1,T_mod,1))
  }
  #trend
  if (ect_trend == 1) {
    Z1 = cbind(Z1, seq(0,(T_mod-1),1))
  }
  
  m = ncol(Z1)
  
  #'szkielet' do macierzy opoznionych przyrostow
  Z2 = matrix(0,T_mod,(k-1)*n)
  
  #wypelniamy macierz opoznien
  for (i in 1:(k-1)) {
    Z2[,((i-1)*n+1):(i*n)] = przyrosty[(k-i):(nrow(przyrosty)-i),]
  }
  
  #dolaczamy nieograniczona stala
  if (const == 1) {
    Z2 = cbind(Z2, matrix(1,T_mod,1))
  }
  
  #dolaczamy niograniczony trend
  if (trend == 1) {
    Z2 = cbind(Z2, seq(1,T_mod,1))
  }
  
  #sztuczne zmienne sezonowe
  if (lsez > 1) {
    sezonowe = seasonaldummy(ts(1:nrow(dane),freq=lsez,start=1))
    sezonowe = sezonowe[,1:(lsez-1)]-1/lsez*matrix(1,nrow(dane),lsez-1)#ortogonalne (wysrodkowane, centered)
    Z2 = cbind(Z2,sezonowe[(k+1):nrow(dane),])
  }
  
  #pozostale zmienne egzogeniczne
  Z2 = cbind(Z2,X)
  
  #macierze momentow z proby (M_ij dla i, j = 0,1,2)
  for (i in 0:2){
    for (j in 0:2){
      nam <- paste("M", i, j, sep = "")
      assign(nam, 1/T_mod*crossprod(get(paste0 ("Z", i)),get(paste0 ("Z", j))))
    }
  }
  
  #macierze S_ij dla i, j = 0,1)
  for (i in 0:1){
    for (j in 0:1){
      nam <- paste("S", i, j, sep = "")
      assign(nam, get(paste0 ("M", i,j))-get(paste0 ("M", i,2))%*%solve(get(paste0 ("M", 2,2)))%*%get(paste0 ("M", 2,j)))
    }
  }
  
  #pierwiastek macierzy S11
  wektory = eigen(S11, symmetric = T)$vectors
  wartosci = diag(m)
  diag(wartosci) = eigen(S11, symmetric = T)$values
  inv_pierwiastek = solve(wektory%*%wartosci^0.5%*%solve(wektory))
  
  macierz = inv_pierwiastek%*%S10%*%solve(S00)%*%S01%*%inv_pierwiastek
  lambda = (eigen(macierz, symmetric = T)$values)[1:n]
  #potencjane wektory kointegrujace
  beta = (inv_pierwiastek%*%eigen(macierz, symmetric = T)$vectors)
  #potencjane wspolczynniki dostosowan
  alpha = (S01%*%beta)[,1:n]
  
  dotestu = as.matrix(-T_mod*log(1-lambda))
  test = cbind(as.matrix(seq(0,n-1,1)),dotestu,as.matrix(rev(cumsum(rev(dotestu)))),n-as.matrix(seq(0,n-1,1)))
  colnames(test) = c('H0', 'test w. wlasnej', 'test sladu', 'n-r')
  rownames(test) = rep('r = ',length.out=n)
  
  out <- list(dane = dane,
              Z0 = Z0,
              Z1 = Z1,
              Z2 = Z2,
              X = X,
              T_mod = T_mod,
              n = n,
              k = k,
              ect_const = ect_const,
              ect_trend = ect_trend,
              const = const,
              trend = trend,
              lsez = lsez,
              beta = beta[,1:n],
              alpha = alpha,
              test = test
  )
  return(out)
  
}

estymacjaVEC <- function(wynik_testu, r = 1){
  Z0 = wynik_testu$Z0
  Z1 = wynik_testu$Z1
  Z2 = wynik_testu$Z2
  k = wynik_testu$k
  n = wynik_testu$n
  if (r > 0) {
    hat_alpha = (wynik_testu$alpha)[,1:r]
    hat_beta = (wynik_testu$beta)[,1:r]
  } else {
    hat_alpha = matrix(0,n,1)
    hat_beta = matrix(0,ncol(Z1),1)
  }
  hat_mPi = hat_alpha%*%t(hat_beta)
  hat_Psi_t = solve(crossprod(Z2,Z2))%*%crossprod(Z2,Z0-Z1%*%hat_beta%*%t(hat_alpha))
  hat_E = Z0-Z1%*%hat_beta%*%t(hat_alpha)-Z2%*%hat_Psi_t
  hat_Sigma = 1/wynik_testu$T_mod*crossprod(hat_E,hat_E)
  
  out <- list(hat_alpha = hat_alpha,
              hat_beta = hat_beta,
              hat_mPi = hat_mPi,
              hat_Psi = t(hat_Psi_t),
              hat_E = hat_E,
              hat_Sigma = hat_Sigma,
              k = k,
              n = n,
              ect_const = wynik_testu$ect_const,
              ect_trend = wynik_testu$ect_trend,
              const = wynik_testu$const,
              trend = wynik_testu$trend,
              lsez = wynik_testu$lsez,
              dane = wynik_testu$dane,
              X = wynik_testu$X
  )
  return(out)
}

VECdoVAR <- function(model){
  k = model$k
  n = model$n
  dane = model$dane
  Gammy = (model$hat_Psi)[,1:(n*(k-1))]
  hat_Phi = (model$hat_Psi)[,(n*(k-1)+1):ncol(model$hat_Psi)]
  mPi = model$hat_mPi
  A = mPi[,1:n]+diag(n)+Gammy[,1:n]
  for (i in 2:(k-1)){
    A = cbind(A,Gammy[,((i-1)*n+1):(i*n)]-Gammy[,((i-2)*n+1):((i-1)*n)])
  }
  A = cbind(A,-Gammy[,(n*(k-2)+1):(n*(k-1))])
  hat_A = rbind(t(A),t(mPi[,(n+1):ncol(mPi)]),t(hat_Phi))
  colnames(hat_A) = rownames(model$hat_Psi)
  out <- list(hat_A = hat_A,
              tilde_Sigma = model$hat_Sigma,
              k = k,
              n = n,
              ect_const = model$ect_const,
              ect_trend = model$ect_trend,
              const = model$const,
              trend = model$trend,
              lsez = model$lsez,
              dane = model$dane,
              X = model$X
  )
  return(out)
}

prognoza <- function(model, H = 4, pu = 0.95, Xprog = {}){
  X = model$X
  TT = nrow(model$X)
  n = model$n
  k = model$k
  const = model$const
  trend = model$trend
  lsez = model$lsez
  if (is.null(model$ect_const)) ect_const = 0 else ect_const = model$ect_const
  if (is.null(model$ect_trend)) ect_trend = 0 else ect_trend = model$ect_trend
  dane = model$dane
  tilde_Sigma = model$tilde_Sigma
  
  y = dane[-c(1:k),]#usuwamy wiersze, od 1 do k
  Z = matrix(0,TT,k*n)#'szkielet' do macierzy opoznien
  
  #wypelniamy macierz opoznien
  for (i in 1:k) {
    Z[,((i-1)*n+1):(i*n)] = dane[(k+1-i):(nrow(dane)-i),]
  }
  
  #dolaczamy stala
  if (const == 1) {
    Z = cbind(Z, matrix(1,TT,1))
  }
  
  #dolaczamy zmienna 'czas'
  if (trend == 1) {
    Z = cbind(Z, seq(1,TT,1))
  }
  
  #sztuczne zmienne sezonowe
  if (lsez > 1) {
    sezonowe = seasonaldummy(ts(1:(nrow(dane)+H),freq=lsez,start=1))
    sezonowe = sezonowe[,1:(lsez-1)]-1/lsez*matrix(1,nrow(dane)+H,lsez-1)#ortogonalne (wysrodkowane, centered)
    Z = cbind(Z,sezonowe[(k+1):nrow(dane),])
  }
  
  #pozostale zmienne egzogeniczne
  Z = cbind(Z,X)# X - macierz zmiennych egzogenicznych, ktora nalezy wczesniej utworzyc
  
  egzo = {}
  if (ect_const==1) egzo = cbind(egzo,matrix(1,H,1))
  if (ect_trend==1) egzo = cbind(egzo,seq(TT,TT-1+H,1))
  if (const==1) egzo = cbind(egzo,matrix(1,H,1))
  if (trend==1) egzo = cbind(egzo,seq(TT+1,TT+H,1))
  if (lsez > 1) egzo = cbind(egzo,sezonowe[(nrow(dane)+1):(nrow(dane)+H),])
  egzo = cbind(egzo,Xprog)
  A = model$hat_A
  yop = cbind(t(as.matrix(y[TT,])),t(as.matrix(Z[TT,1:(n*(k-1))])))
  p_punktowa = matrix(0,H,n)
  niepewnosc = matrix(0,H,n)
  colnames(p_punktowa) <- colnames(dane)
  colnames(niepewnosc) <- colnames(dane)
  tildeA = t(A[1:(n*k),])
  tildeA = rbind(tildeA,cbind(diag(n*(k-1)),matrix(0,n*(k-1),n)))
  MSE = matrix(0,n,n)
  phi = diag(n*k)
  for (i in 1:H){
    Zprog = cbind(yop,t(as.matrix(egzo[i,])))
    yp = Zprog%*%A
    p_punktowa[i,] = yp
    MSE = MSE + phi[1:n,1:n]%*%tilde_Sigma%*%t(phi[1:n,1:n])
    niepewnosc[i,] = diag(MSE)^0.5
    yop = cbind(yp,t(as.matrix(yop[,1:(n*(k-1))])))
    phi = phi%*%tildeA
  }
  #prognoza przedzialowa
  kwantyl = qnorm((pu+1)/2)
  lewy = p_punktowa - kwantyl*niepewnosc
  prawy = p_punktowa + kwantyl*niepewnosc
  
  p_przedzialowa = matrix(0,H,2*n)
  nazwy = {}
  for (i in 1:n){
    p_przedzialowa[,2*i-1] = lewy[,i]
    p_przedzialowa[,2*i] = prawy[,i]
    nazwy = c(nazwy,paste('dolny',colnames(dane)[i]),paste('gorny',colnames(dane)[i]))
  } 
  colnames(p_przedzialowa) = nazwy
  
  out <- list(p_punktowa = p_punktowa,
              niepewnosc = niepewnosc,
              p_przedzialowa = p_przedzialowa)
  return(out)
}
