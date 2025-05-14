library(rTensor)

X = rand_tensor(c(2,3,4))
U = pracma::randortho(2)
V = pracma::randortho(3)
W = pracma::randortho(4)

Y = ttl(X, list(U, diag(3), diag(4)), 1:3)

Y[,2,3]
