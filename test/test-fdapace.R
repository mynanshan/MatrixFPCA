
set.seed(1)
n <- 100
M <- 51
pts <- seq(0, 1, length.out=M)
mu <- rep(0, length(pts))
sampDense <- fdapace::MakeGPFunctionalData(n, M, mu, K=1, basisType='sin', sigma=0.01)
samp4 <- fdapace::MakeFPCAInputs(tVec=sampDense$pts, yVec=sampDense$Yn)
res4E <- fdapace::FPCA(samp4$Ly, samp4$Lt, list(error=TRUE))
sampN <- fdapace::GetNormalisedSample(res4E, errorSigma=TRUE)

fdapace::CreatePathPlot(subset=1:20, inputData=samp4, obsOnly=TRUE, showObs=FALSE)
fdapace::CreatePathPlot(subset=1:20, inputData=sampN, obsOnly=TRUE, showObs=FALSE)
