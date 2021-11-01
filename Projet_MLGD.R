### Generation du signal

x <- seq(-4*pi,3*pi, 1/10)
n <- length(x)
f = cos(x) + atan(x)  
Y = f + rnorm(n)
X=ComputeCosMat(n)  # X est donc orthonormale
Z= t(X)%*%Y
beta = 0.25

# Calcul de l'estimateur de Gibbs

fhat_gibbs = 0
for (i in (1:n)){
  fhat_gibbs = fhat_gibbs + (exp(beta*(Z[i]^2))/(exp(2*beta + log(n)) + exp(beta*Z[i]^2)))*Z[i]*X[,i]
}

# Calcul de l'estimateur par sélection de modèles
K <- 1
lambda = K * (1 + sqrt(2*log(n)))^2
#print(lambda)

thre <- Z*Z
Xm = X[,thre > lambda]  # on trouve le meilleur sous-espace en seuillant 
betasm = solve(t(Xm)%*%Xm)%*%t(Xm)%*%Y  # EMV avec sous-espace réduit
betatilde = rep(0,n)
betatilde[thre > lambda] <- betasm # On reprojette sur un vecteur de dimension originale 

fhat_model = X%*%betatilde