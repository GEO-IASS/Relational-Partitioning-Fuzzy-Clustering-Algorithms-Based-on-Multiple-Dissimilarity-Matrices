library(mvtnorm)
library(plyr)
library(dplyr)
library(class)
library(tidyr)
library(MVA)   # para a matriz de distancias euclidianas 
library(caret) # para criar os folds
library(stats) # friedman test
library(ggplot2)
library(PMCMR) # nemenyi test

# Carregando o conjunto de dados:
dados <- read.table(file.choose(), sep=",")
dim(dados)
head(dados,2)
tail(dados,2)
# converter as categorias em números:

# BRICKFACE   -- 1
# SKY         -- 6
# FOLIAGE     -- 3
# CEMENT      -- 2
# WINDOW      -- 7
# PATH        -- 5
# GRASS       -- 4

dados[ ,1] = dados[ ,1] %>% factor(.) %>% as.numeric(.)
V1 = dados[,1]
head(dados, 5)

# classificador bayesiano:
# INPUTS:  treino - amostra de treinamento
#          teste - amostra de teste
#          V1 e V2 - rótulos do treino e do teste
# OUTPUTS: taxas de acerto do treino e do teste
#          e vetores com a classificação. 


classificadorBayesiano <- function(treino, teste, V1, V2, view){
  
  Priori = NULL
  medias = list()
  covariancias = list()
  
  # estimando a probabilidade a priori:
  for(i in 1:7){
    Priori[i] = filter(treino, treino[,1]==i) %>% nrow(.)/nrow(treino)
  }
  
  for(i in 1:7){
    # estimando os parâmetros da normal multivariada
    data = filter(treino, treino[,1]==i)
    
    #media
    medias[[i]] = colMeans(data[,-1])
    #matriz de covariancia
    covariancias[[i]] = diag(apply(data[,-1], 2,var))
  }
  
  # densidade normal multivariada TEM QUE TIRAR A TERCEIRA VARIÁVEL
  #classificacaoTreino = NULL;
  
  # 
  # for(i in 1:nrow(treino)){
  #   
  #   posteriorisTreino = NULL;
  #   for(j in 1:length(medias)){
  #     posteriorisTreino[j] = dmvnorm(x=treino[i,c(-1,-3)],medias[[j]][-3],covariancias[[j]][-3,-3]) * Priori[j]
  #     #posteriorisTreino[j] = dmvnorm(x=treino[i,-1],medias[[j]],covariancias[[j]]) * Priori[j]
  #   }
  #   posteriorisTreino = posteriorisTreino/sum(posteriorisTreino)
  #   classificacaoTreino[i] = posteriorisTreino %>% which.max(.)
  # }
  
  classificacaoTeste = NULL;
  for(i in 1:nrow(teste)){
    
    posteriorisTeste = NULL;
    if (view == 1){
      for(j in 1:length(medias)){
        posteriorisTeste[j] = dmvnorm(x=teste[i,c(-1,-4)],medias[[j]][-3],covariancias[[j]][-3,-3]) * Priori[j]
        #posteriorisTeste[j] = dmvnorm(x=teste[i,-1],medias[[j]],covariancias[[j]]) * Priori[j]
      }}
    else if (view == 2){
      for(j in 1:length(medias)){
        posteriorisTeste[j] = dmvnorm(x=teste[i,-1],medias[[j]],covariancias[[j]]) * Priori[j]
        #posteriorisTeste[j] = dmvnorm(x=teste[i,-1],medias[[j]],covariancias[[j]]) * Priori[j]
      }}
    
    posteriorisTeste = posteriorisTeste/sum(posteriorisTeste) 
    classificacaoTeste[i] = which.max(posteriorisTeste)
  }
  
  return(list(#classificacaoTreino = classificacaoTreino, 
    #taxaAcertoTreino = mean(treino[,1] == classificacaoTreino),
    classificacaoTeste = classificacaoTeste,
    taxaAcertoTeste = mean(teste[,1] == classificacaoTeste)))
  
}

knn1 <- function(treino, teste, k){
  matrizDist = as.matrix(pdist(treino[,-1], teste[,-1])) # Calcula a matriz de distancia dos elementos do teste com os do treino
  classTeste = NULL
  i = 1
  for(i in 1:nrow(teste)){
    kn = treino[order(matrizDist[,i], decreasing = F)[1:k],1] # pega os k primeiro elementos com as menores distancias
    knunique = unique(kn) # classes encontradas
    classTeste[i] = knunique[which.max(tabulate(match(kn,knunique)))] # selecao da classe mais predominante atribuindo no elemento do teste#
  }
  return(list(ClassificaoTeste = classTeste, taxaAcerto = mean(classTeste == teste[,1])))
}





# validação cruzada:
# k-folds:

# funcao para recuperar um fold:
getFold <- function(view = 1, fold){
  
  if(view == 1)
    return(as.matrix(view1[foldsList1[[fold]], ]))
  if(view == 2)
    return(as.matrix(view2[foldsList2[[fold]], ]))
}  

# calculo da moda de um vetor
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}


# Dividir em dois views:
#dados = dados[sample(1:nrow(dados), nrow(dados), replace = FALSE), ]
view1 = dados[ ,c(1:10)]   #shape 
view2 = dados[ ,c(1,11:20)]  #rgb



#### VALIDAÇÃO CRUZADA ##### VIEW1

times = 30 # número de validações cruzadas

erroBayesGeralView1 = erroBayesGeralView2 = NULL;
erroKnnGeralView1 = erroKnnGeralView2 =NULL;
erroVMGeral = NULL;

for(j in 1:times){
  print(j)
  # criando 10 folds
  foldsList1 <- createFolds(view1[,1], k=10, list=TRUE)  
  
  # validação cruzada:
  erroBayes1 = erroBayes2 = NULL;
  erroKnn1 = erroKnn2 = erroVM = NULL;
  for(i in 1:10){
    teste1=view1[foldsList1[[i]], ] # este fold será usado para teste
    rotulosTreino1 = view1[-foldsList1[[i]],1]; length(rotulosTreino1)
    rotulosTeste1 = teste1[,1]; length(rotulosTeste1)
    treino1 = view1[-foldsList1[[i]],]  
    
    teste2=view2[foldsList1[[i]], ] # este fold será usado para teste
    rotulosTreino2 = view2[-foldsList1[[i]],1]; length(rotulosTreino2)
    rotulosTeste2 = teste2[,1]; length(rotulosTeste2)
    treino2 = view2[-foldsList1[[i]],]
    
    # classificação bayesiana
    bayes1 = classificadorBayesiano(treino1,teste1,V1 = rotulosTreino1,V2 = rotulosTeste1, 1)
    erroBayes1[i] = 1 - bayes1$taxaAcertoTeste 
    
    bayes2 = classificadorBayesiano(treino2,teste2,V1 = rotulosTreino2,V2 = rotulosTeste2, 2)
    erroBayes2[i] = 1 - bayes2$taxaAcertoTeste 
    
    # knn
    KNN1 = knn1(treino1, teste1, k=1)
    KNN1 = KNN1$ClassificaoTeste
    erroKnn1[i] = 1-mean(KNN1==rotulosTeste1)
    
    KNN2 = knn1(treino2, teste2, k=1)
    KNN2 = KNN2$ClassificaoTeste
    erroKnn2[i] = 1-mean(KNN2==rotulosTeste2)
    
    # voto majoritario
    deltaKb1 = bayes1$classificacaoTeste
    deltaKb2 = bayes2$classificacaoTeste
    
    erroVM[i] = 1-mean(apply(cbind(deltaKb1, deltaKb2, as.numeric(KNN1), as.numeric(KNN2)), 1, Mode) == rotulosTeste1)
    
    
  }

  erroBayesGeralView1[j] = mean(erroBayes1)
  erroBayesGeralView2[j] = mean(erroBayes2)
  erroKnnGeralView1[j] = mean(erroKnn1)
  erroKnnGeralView2[j] = mean(erroKnn2)
  erroVMGeral[j] = mean(erroVM)
}

# ESTIMATIVAS PONTUAIS PARA TAXA DE ERRO DOS CLASSIFICADORES
alpha = 0.05
a = qt(alpha/2, df = (times-1), lower.tail = F)

estimMatrix = 
  cbind(
    rbind(
mean(erroBayesGeralView1),  # erro do classificador de Bayes para view1
mean(erroKnnGeralView1),    # erro do classificador Knn para view1
mean(erroBayesGeralView2),  # erro do classificador de Bayes para view2
mean(erroKnnGeralView2),    # erro do classificador Knn para view2
mean(erroVMGeral)               # taxa de erro do voto majoritáirio

# INTERVALOS DE CONFIANÇA PARA TAXA DE ERRO DOS CLASSIFICADORES
    ),
rbind(
mean(erroBayesGeralView1)  + c(-a,a)*sd(erroBayesGeralView1)/sqrt(times), 
mean(erroKnnGeralView1)    + c(-a,a)*sd(erroKnnGeralView1)/sqrt(times),
mean(erroBayesGeralView2)  + c(-a,a)*sd(erroBayesGeralView2)/sqrt(times),
mean(erroKnnGeralView2)    + c(-a,a)*sd(erroKnnGeralView2)/sqrt(times),
mean(erroVM)               + c(-a,a)*sd(erroVM)/sqrt(times)
)
  )

library(xtable)
xtable(estimMatrix, digits=4)

sds = c(sd(erroBayesGeralView1), sd(erroKnnGeralView1), sd(erroBayesGeralView2), sd(erroKnnGeralView2), sd(erroVM))
sds
# TESTE DE FRIEDMAN PARA COMPARAR OS CLASSIFICADORES


m1 = as.matrix(cbind(erroBayesGeralView1, erroBayesGeralView2, erroKnnGeralView1, erroKnnGeralView2, erroVMGeral))
friedman.test(m1)

ntest = posthoc.friedman.nemenyi.test(m1)$p.value
xtable(ntest)


m = data.frame(Times = 1:30, erroBayesGeralView1, erroBayesGeralView2, erroKnnGeralView1, erroKnnGeralView2, erroVM)
m = gather(m, value = "Erro", key = "Tipo", erroBayesGeralView1, erroBayesGeralView2, erroKnnGeralView1, erroKnnGeralView2, erroVM)
ggplot(m,aes(Times,y=Erro,color=Tipo)) + geom_line()







