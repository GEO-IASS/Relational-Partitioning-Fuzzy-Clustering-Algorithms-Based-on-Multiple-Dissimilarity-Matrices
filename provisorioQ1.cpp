#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
mat UvectorC(List Gk, mat lambda, mat d1, mat d2, float m=1.6){
  
  mat U(d1.n_rows, lambda.n_rows);
  U.zeros();
  float somad1;
  float somad2;
  float numerador;
  float somaGh1 = 0;
  float somaGh2 = 0;
  
  for(int i=0; i<U.n_rows; i++){
    for(int k=0; k<U.n_cols; k++){
      
      somad1 = 0;
      somad2 = 0;
      rowvec clusterK = Gk[k];
      numerador = 0;
      
      for(int e=0; e<clusterK.n_elem; e++){
        //cout << clusterK(e) << endl;
        //cout << i << endl;
        //cout << d1(i,clusterK(e)) << endl;
        somad1 += d1(i,clusterK(e)-1);
        somad2 += d2(i,clusterK(e)-1);
      }
      
      numerador = lambda(k,0)*somad1 + lambda(k,1)*somad2;
      //cout << numerador << endl;
      
      rowvec denominador(lambda.n_rows);
      denominador.zeros();
      
      for(int h=0; h<lambda.n_rows; h++){
        
        somaGh1 = 0;
        somaGh2 = 0;
        rowvec clusterH = Gk[h];
        
        for(int e1=0; e1<clusterH.n_elem; e1++){
          somaGh1 += d1(i,clusterH(e1)-1);
          somaGh2 += d2(i,clusterH(e1)-1);
        }
        
        denominador(h) = lambda(h,0) * somaGh1 + lambda(h,1)* somaGh2;
        
      }
      
      U(i,k) = pow(sum( pow( (numerador/denominador), 1/(m-1) ) ),-1);
    }
  }
  
  
  return U;
}

// [[Rcpp::export]]
List bestPrototypeMC(mat U, mat lambda, mat d1, mat d2, int q, int K, float m=1.6, int s=1){
  
  List listaPrototipos(K);
  uvec Gmin;
  vec l(d1.n_rows);
  double J;
  
  for(int k=0; k<K; k++){
    for(int h=0; h<d1.n_rows; h++){
      J=0;
      for(int i=0; i<d1.n_rows; i++){
        J = J + pow(U(i,k),m) * (pow(lambda(k,0),s) * d1(i,h) + pow(lambda(k,1),s) * d2(i,h));
      }
      l(h) = J;
    }
    
    //cout << l << endl;
    Gmin = sort_index(l);
    Gmin = Gmin.subvec(0,q-1);
    //Gmin.reshape(1,3);
    listaPrototipos(k) = Gmin+1; 
  }
  
  return listaPrototipos;
}

// [[Rcpp::export]]
double JcalcMC(mat U, mat d1, mat d2, mat lambda, List Gk, float m=1.6){
  
  double soma4 = 0;
  
  for(int k=0; k<lambda.n_rows; k++){
    
    double soma3 = 0;
    for(int i=0; i<U.n_rows; i++){

      double soma1 = 0;
      double soma2 = 0; 
      vec clusterK = Gk(k);
      
      for(int e=0; e<clusterK.n_elem; e++){
        soma1 += d1(i,clusterK(e)-1);
        soma2 += d2(i,clusterK(e)-1);
        //cout << clusterK << endl;
      }
      
      soma3 += pow(U(i,k),m) * (lambda(k,0)*soma1 + lambda(k,1)*soma2);
  
    }
    
    soma4 += soma3;
  }
  return soma4;
}

// [[Rcpp::export]]
mat relevantWvectorC(mat U, List Gk, mat d1, mat d2, double m=1.6){
  
  mat lambda(U.n_cols,2);
  double numerador; 
  double denominador;
  
  for(int k=0; k<lambda.n_rows; k++){
    for(int j=0; j<lambda.n_cols; j++){
      
      double soma1=0;
      double soma2=0;
      
      for(int i=0; i<d1.n_rows; i++){
        
        double soma11=0;
        double soma21=0;
        vec clusterK = Gk(k);
        
        for(int e=0; e<clusterK.n_elem;e++){
          soma11 += d1(i,clusterK(e)-1);
          soma21 += d2(i,clusterK(e)-1);
        }
        
        soma1 += (pow(U(i,k),m) * soma11);
        soma2 += (pow(U(i,k),m) * soma21);
      }
      
      numerador = pow(soma1*soma2, 0.5);
      //cout << numerador << endl;
      
      if(j==0){
        denominador = soma1;
      } 
      if(j==1){
        denominador = soma2;
      }
      
      lambda(k,j) = (numerador/denominador);
    }
  }
  return lambda;
}

// [[Rcpp::export]]
List algoritmo(List Gk, mat d1, mat d2, int K=7, float m=1.6, int s=1, int Maxit = 1000, 
                double epsilon=0.001, int q = 3, int t = 0){
  
  mat U(d1.n_rows,K);
  mat lambda(K,2);
  lambda.ones();
  double J0, J1;
  Function UvectorC("UvectorC");
  Function bestPrototypeMC("bestPrototypeMC");
  Function JcalcMC("JcalcMC");
  Function relevantWvectorC("relevantWvectorC");
  
  U = as<arma::mat>(UvectorC(Gk, lambda, d1, d2, m));
  J0 = as<double>(JcalcMC(U,d1,d2,lambda,Gk));

      
  for(int e=0; e<Maxit; e++){
    
    t++;
    Gk = bestPrototypeMC(U, lambda, d1, d2, q, K, m, s=1);
      
    // PASSO 2: COMPUTATION OF BEST RELEVANCE WEIGHT VECTOR
    lambda = as<arma::mat>(relevantWvectorC(U, Gk, d1, d2, m=1.6));
        
    // PASSO 3: DEFINICAO DA MELHOR PARTICAO FUZZY
    U = as<arma::mat>(UvectorC(Gk, lambda, d1, d2, m=1.6));
          
    // CRITERIO DE PARADA: compute J
    J1 = as<double>(JcalcMC(U,d1,d2,lambda,Gk));
    
    if(std::abs(J1 - J0) < epsilon){
      break;
    } else if(J1 > J0){
      cout << "Aumentou: ";
      cout << std::setprecision(12) << J1 << endl; 
      J0 = J1;
    } else if(J1 <= J0){
      cout << "Diminuiu: ";
      cout << std::setprecision(12) << J1 << endl;
      J0 = J1;
    }
    
  }
  
  cout << "------------- Próxima iteração -------------------" << endl;
  
  List output; 
  output["criterio"] = J1; output["U"] = U; output["lambda"] = lambda;
  output["Gk"] = Gk;
  
  return output;
}
  

  
  
/*** R
library(mclust)
# library(microbenchmark)
# 
# microbenchmark(j1 = Uvector(Gk, lambda, d1, d2),  
#                j2 = UvectorC(Gk, lambda, d1, d2), times = 1L)
# 
# microbenchmark(
# prot1 = bestPrototypeMC(U, lambda, d1, d2, q, K, m, s=1),
# prot2 = bestPrototypeM(U, lambda, d1, d2, q, K, m, s=1), times = 1L)
# 
# 
# microbenchmark(a1 = JcalcM(U,d1,d2,lambda,Gk),
# a2 = JcalcMC(U,d1,d2,lambda,Gk), times=1L)
# 
# microbenchmark(
# e1 = relevantWvector(U,Gk,d1,d2,m),
# e2 = relevantWvectorC(U,Gk,d1,d2,m), times = 1L)

library(plyr)
library(dplyr)
library(tidyr)
library(MVA)   # para a matriz de distancias euclidianas 
library(caret) # para criar os folds
  
# Carregando o conjunto de dados:
dados <- read.table(file.choose(), sep=",")
dim(dados)
  
  
# converter as categorias em números:
  
# BRICKFACE   -- 1
# SKY         -- 6
# FOLIAGE     -- 3
# CEMENT      -- 2
# WINDOW      -- 7
# PATH        -- 5
# GRASS       -- 4
  
dados[ ,1] = dados[ ,1] %>% factor(.) %>% as.numeric(.)
# Dividir em dois views:
d = sample(1:2100, 2100, replace=FALSE)
#d = 1:2100
view1 = dados[d, c(1,2:10)]   #shape  
view2 = dados[d, c(1,11:20)]  #rgb 
labels = dados[d,1]

d1 = as.matrix(dist(view1[,-1])); dim(d1)  # matriz de dissimilaridade 1 
d2 = as.matrix(dist(view2[,-1])); dim(d2)  # matriz de dissimilaridade 2
    
# selecionar aleatoriamente K protótipos distintos Gk
options(max.print=999999) 

#set.seed(21)
K=7
index = sample(1:nrow(cbind(view1, view2)), K*q, replace = FALSE)
Gk = list()

cont = 1
for(i in 1:K){
  Gk[[i]] <- index[cont:(cont+q-1)]
  cont = cont+q
}

Gk 

# para cada objeto e_i, calcule o grau de pertencimento no cluster k

lambda = matrix(1, nrow = K, ncol = 2)  #  Matriz de pesos de relevancia inicial

U = UvectorC(Gk,lambda,d1,d2,m=1.6)
rowSums(U)
# compute J:

J0 = JcalcMC(U,d1,d2,lambda,Gk)

Jfinal = algoritmo(Gk, d1, d2, Maxit = 1000, epsilon = 10^-3, q=3)$criterio

for(i in 1:100){
  
  index = sample(1:nrow(cbind(view1, view2)), K*q, replace = FALSE)
  Gk = list()
  
  cont = 1
  for(i in 1:K){
    Gk[[i]] <- index[cont:(cont+q-1)]
    cont = cont+q
  }
  
  algo1 = algoritmo(Gk, d1, d2, Maxit = 1000, epsilon = 10^-3, q=3)
  #algo1$lambda[,1] * algo1$lambda[,2]
  #rowSums(algo1$U)
  #algo1$Gk
  if(algo1$criterio < Jfinal){
    Jfinal <- algo1$criterio
    algofinal = algo1
    cat("Novo mínimo:\n", Jfinal, "\n")
  }  
}

algofinal

# partição hard

g = algofinal$U
  
  
for(i in 1:nrow(algofinal$U)){
  for(j in 1:ncol(algofinal$U)){
      
    if(algofinal$U[i,j] == max(algofinal$U[i, ])){
      algofinal$U[i,j] = 1;
    } else{ 
      algofinal$U[i,j] = 0;
    }
      
  }
    
}
  
  
colSums(algofinal$U)
objetos = list("BRICKFACE"  = NULL,
               "CEMENT"= NULL,
               "FOLIAGE"= NULL,
               "GRASS" = NULL,
               "PATH"= NULL,
               "SKY"= NULL,
               "WINDOW"= NULL
)

randVec = NULL;  # vetor das classificações dadas pelo algoritmo;
realVec = labels

for(i in 1:nrow(algofinal$U))
  for(j in 1:ncol(algofinal$U)){
    if(algofinal$U[i,j] == 1){
      objetos[[j]] = c(objetos[[j]],i)
      randVec[i] = j
    }
  }
      
mean(realVec == randVec)      

algofinal[["objetos"]] = objetos
algofinal[["AdjustedRandIndex"]] = mclust::adjustedRandIndex(realVec,randVec)

setwd("/home/kassio/Desktop")

# salvando num arquivo txt:
options(max.print=999999)
capture.output(algofinal, file = "ResultadoFinal.txt")

*/



