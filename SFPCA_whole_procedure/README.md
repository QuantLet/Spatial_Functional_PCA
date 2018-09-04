[<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/banner.png" width="888" alt="Visit QuantNet">](http://quantlet.de/)

## [<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/qloqo.png" alt="Visit QuantNet">](http://quantlet.de/) **SFPCA_whole_procedure** [<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/QN2.png" width="60" alt="Visit QuantNet 2.0">](http://quantlet.de/)

```yaml

﻿
Name of Quantlet: SFPCA_whole_procedure

Published in: Spatial_Functional_PCA

Description: 'Conducts the whole procedure, including preprocessing data, estimation (getting loadings and factors, smoothing factors, updating loadings) and risk attitude analysis, related to the paper Spatial Functional Principal Component Analysis with Applications to Brain Image Data by Yingxing Li, Chen Huang and Wolfgang Härdle.'

Keywords : 
- fmri
- pca
- smoothing
- basis
- factor analysis

Author:          Chen Huang, Yingxing Li

Submitted:        Wednesday, June 20 2018 by Chen Huang

```

### R Code
```r

## install and load packages
libraries = c("R.matlab", "NMOF", "splines", "tictoc", "matrixStats")
lapply(libraries, function(x) if (!(x %in% installed.packages())) {
    install.packages(x)
})
lapply(libraries, library, quietly = TRUE, character.only = TRUE)
## source functions
cube_mean2 = function(Y, m, n, l) {
    M = dim(Y)[1]
    N = dim(Y)[2]
    L = dim(Y)[3]
    Y.new = array(0, dim = c(m, n, l))
    Y.new2 = array(0, dim = c(m, n, l))
    b1 = floor(M/m)
    b2 = floor(N/n)
    b3 = floor(L/l)
    for (mm in 1:m) {
        for (nn in 1:n) {
            for (ll in 1:l) {
                Y.new[mm, nn, ll] = mean(Y[((mm - 1) * b1 + 1):(mm * b1), ((nn - 
                  1) * b2 + 1):(nn * b2), ((ll - 1) * b3 + 1):(ll * b3)])
                Y.new2[mm, nn, ll] = var(Y[((mm - 1) * b1 + 1):(mm * b1), ((nn - 
                  1) * b2 + 1):(nn * b2), ((ll - 1) * b3 + 1):(ll * b3)])
            }
        }
    }
    result = list(Y.new = Y.new, va = mean(Y.new2))
}

genddd = function(N, m, n, l) {
    mn = m * n
    d3 = floor(N/mn) * (N%%mn == 0) + (floor(N/mn) + 1) * (N%%mn != 0)
    temp = N - (d3 - 1) * mn
    d2 = floor(temp/m) * (temp%%mn == 0) + floor(temp/m + 1) * (temp%%mn != 0)
    d1 = temp - (d2 - 1) * m
    d = c(d1, d2, d3)
    return(d)
}

genB = function(n, K, p) {
    B = bs(c((1:n)/n - 1/2/n), df = K + p, degree = p)
    return(B)
}

genP = function(m, K, p) {
    ms = choose(m, 0:(m)) * (-1)^c(0:(m))
    D_m = matrix(0, K + p - m, K + p)
    for (i in 1:(K + p - m)) {
        D_m[i, i:(i + m)] = ms
    }
    P = t(D_m) %*% D_m
    return(P)
}

gentD = function(B, P) {
    eB = eigen(t(B) %*% B)
    # pK=sum(eB$values<=1e-05) if(pK>0)
    # eBvp=c(1/eB$values[which(eB$values>1e-05)],rep(0,pK))
    # Bhalf=eB$vectors%*%diag(sqrt(eBvp))%*%t(eB$vectors)
    Bhalf = eB$vectors %*% diag(1/sqrt(eB$values)) %*% t(eB$vectors)
    eBP = eigen(Bhalf %*% P %*% t(Bhalf))
    D = eBP$values
    D = D * (D > 1e-05)
    V = eBP$vectors
    temp = t(V) %*% Bhalf %*% t(B)
    result = list(D = D, temp = temp)
    return(result)
}

genlamD = function(lambda, D) {
    lamD = 1/(1 + lambda * D)
    return(lamD)
}

GCV = function(lambda, sY2, dim, D1, D2, D3, tildeY2) {
    lamD1 = genlamD(lambda[1], D1)
    lamD2 = genlamD(lambda[2], D2)
    lamD3 = genlamD(lambda[3], D3)
    
    W = outer(outer(lamD1, lamD2), lamD3)
    WtildeY2 = W * tildeY2
    W2tildeY2 = W^2 * tildeY2
    
    GCVnum = sY2 - 2 * sum(WtildeY2) + sum(W2tildeY2)
    GCVnum/(dim[1] * dim[2] * dim[3] - sum(lamD1) * sum(lamD2) * sum(lamD3))
    # GCVnum/((dim[1]-sum(lamD1))*(dim[2]-sum(lamD2))*(dim[3]-sum(lamD3)))
}

################# stage 1 in R
setwd("D:/fmri 1st")
interval = 3
stim = readMat("decartFin2.mat")$decart
stim.scan = ceiling(round(stim[, 2])/2.5) + 1
scan.idx = matrix(0, nrow = 81, ncol = interval)
for (q in 1:81) {
    scan.idx[q, ] = (stim.scan[q] + 1):(stim.scan[q] + interval)
}
writeMat(paste("scanidx.mat", sep = ""), scan.idx = scan.idx)
############### stage 2 in MATLAB 'divide_data.m') stage 3 in R
subjects = c(1:22)[-c(7, 13, 14, 20, 22)]
I = 17 * 81
m = 30
n = 35
l = 30
M = 91
N = 109
L = 91
Yall = array(0, c(I, M, N, L))
aYall = array(0, c(I, m, n, l))
aYva = rep(0, I)
i = 1

tic("total")
tic("preprocess")
for (s in 1:17) {
    subject = subjects[s]
    for (q in 1:81) {
        filename = paste("voxels_", subject, "_", q, ".mat", sep = "")
        voxels = readMat(filename)$temp.q
        print(paste(s, ",", q, ",", i))
        for (d1 in 1:M) {
            for (d2 in 1:N) {
                for (d3 in 1:L) {
                  Yall[i, d1, d2, d3] = mean(voxels[d1, d2, d3, 2:3] - voxels[d1, 
                    d2, d3, 1])
                }
            }
        }
        result = cube_mean2(Yall[i, , , ], m, n, l)
        aYall[i, , , ] = result$Y.new
        aYva[i] = result$va
        i = i + 1
    }
}
toc()

tic("estimation")
tic("getloading")
LAY = matrix(0, I, I)
for (i in 1:I) {
    for (j in i:I) {
        LAY[i, j] = sum(aYall[i, , , ] * aYall[j, , , ])
        # print(paste(i,',',j))
    }
}
LAY = LAY + t(LAY)
diag(LAY) = diag(LAY)/2

eLAY = eigen(LAY/m/n/l - diag(aYva))
r = 19
eA = t(eLAY$vectors[, 1:r])  #get loadings approach 1
# eA=t(eLAY$vectors[,1:r] %*% diag(sqrt(eLAY$values[1:r]))) #get loadings
# approach 2 ???in this case, normalization on eF is not needed)
load.old = eA
toc()

tic("getfactor")
eF.raw = array(0, c(r, M, N, L))
eF = array(0, c(r, M, N, L))
for (k in 1:r) {
    for (i in 1:I) {
        eF.raw[k, , , ] = eF.raw[k, , , ] + eA[k, i] * Yall[i, , , ]
    }
    eF[k, , , ] = eF.raw[k, , , ]/sqrt(sum(eF.raw[k, , , ]^2)) * sqrt(M * N * L)  #normalize eigenfunctions
}
toc()

tic("smoothfactor")
#### smooth eF
K = c(30, 35, 30)
dim = c(M, N, L)
p = 3
mm = 2

B1 = genB(M, K[1], p)
B2 = genB(N, K[2], p)
B3 = genB(L, K[3], p)
P1 = genP(mm, K[1], p)
P2 = genP(mm, K[2], p)
P3 = genP(mm, K[3], p)

res1 = gentD(B1, P1)
D1 = res1$D
temp1 = res1$temp

res2 = gentD(B2, P2)
D2 = res2$D
temp2 = res2$temp

res3 = gentD(B3, P3)
D3 = res3$D
temp3 = res3$temp

eF.smooth = array(0, dim = c(r, M, N, L))
eF.trim = array(0, dim = c(r, M, N, L))
for (k in 1:r) {
    Y2 = eF[k, , , ]^2
    sY2 = sum(Y2)
    
    oldY2 = array(0, c(K[1] + p, K[2] + p, L))
    for (d3 in 1:L) oldY2[, , d3] = temp1 %*% eF[k, , , d3] %*% t(temp2)
    
    tildeY2 = array(0, c(K[1] + p, K[2] + p, K[3] + p))
    for (d1 in 1:(K[1] + p)) for (d2 in 1:(K[2] + p)) tildeY2[d1, d2, ] = temp3 %*% 
        oldY2[d1, d2, ]
    
    tildeY2 = tildeY2^2
    
    lambda0 = gridSearch(GCV, levels = list(exp(seq(-15, 15, length = 15)), exp(seq(-15, 
        15, length = 15)), exp(seq(-15, 15, length = 15))), sY2 = sY2, dim = dim, 
        D1 = D1, D2 = D2, D3 = D3, tildeY2 = tildeY2, printDetail = FALSE)$minlevels
    lambda = nlm(GCV, lambda0, sY2 = sY2, dim = dim, D1 = D1, D2 = D2, D3 = D3, tildeY2 = tildeY2)$estimate
    
    H1 = B1 %*% solve(t(B1) %*% B1 + lambda[1] * P1) %*% t(B1)
    H2 = B2 %*% solve(t(B2) %*% B2 + lambda[2] * P2) %*% t(B2)
    H3 = B3 %*% solve(t(B3) %*% B3 + lambda[3] * P3) %*% t(B3)
    
    oldfac = array(0, dim = dim)
    for (d3 in 1:L) oldfac[, , d3] = H1 %*% eF[k, , , d3] %*% H2
    
    hatfac = array(0, dim = dim)
    for (d1 in 1:M) for (d2 in 1:N) hatfac[d1, d2, ] = H3 %*% oldfac[d1, d2, ]
    
    eF.smooth[k, , , ] = hatfac
    
    #### trim the factors with the quantile bands for(k in 1:r){
    quantiles = quantile(eF.smooth[k, , , ], probs = c(0.001, 0.999))
    for (d1 in 1:M) {
        for (d2 in 1:N) {
            for (d3 in 1:L) {
                if (quantiles[1] < eF.smooth[k, d1, d2, d3] && eF.smooth[k, d1, d2, 
                  d3] < quantiles[2]) {
                  eF.trim[k, d1, d2, d3] = 0
                } else {
                  eF.trim[k, d1, d2, d3] = 1
                }
            }
        }
    }
    
}
toc()

tic("updateloading")
#### update loadings
F = matrix(0, r, r)
for (i in 1:r) {
    for (j in 1:r) {
        F[i, j] = sum(eF.smooth[i, , , ] * eF.smooth[j, , , ])
    }
}
Y.tilde = matrix(0, r, I)
for (i in 1:r) {
    for (j in 1:I) {
        Y.tilde[i, j] = sum(Yall[j, , , ] * eF.smooth[i, , , ])
    }
}
load.new = solve(F) %*% Y.tilde
toc()
toc()
toc()

#### RA analysis
subjects = c(1:22)[-c(7, 13, 14, 20, 22)]
RA = cbind(subjects, c(5.6, 4.8, 5.6, 11.3, 4.1, 3.7, 5, 4.7, 6.3, 3.8, 1.3, 12.6, 
    8.6, 1.8, 16.6, 18.3, 0.6))

load.new.t = t(load.new)
load.mean = matrix(0, 17, r)
load.sd = matrix(0, 17, r)
for (s in 1:17) {
    load.mean[s, ] = colMeans(load.new.t[((s - 1) * 81 + 1):(s * 81), ])
    load.sd[s, ] = colSds(load.new.t[((s - 1) * 81 + 1):(s * 81), ])
}

## model selection




## in-sample fitting m1 (3.3)
summary.lm(lm(RA[, 2] ~ 0 + load.mean[, c(1, 5)]))
AIC(lm(RA[, 2] ~ 0 + load.mean[, c(1, 5)]))
# m2 (3.5)
summary.lm(lm(RA[, 2] ~ 0 + load.mean[, c(1, 5)] + load.sd[, c(1, 4, 5, 9)]))
AIC(lm(RA[, 2] ~ 0 + load.mean[, c(1, 5)] + load.sd[, c(1, 4, 5, 9)]))
# m3 (3.7)
summary.lm(lm(RA[, 2] ~ 0 + load.mean[, c(1, 5)] + load.sd[, c(1, 4, 5, 8)]))
AIC(lm(RA[, 2] ~ 0 + load.mean[, c(1, 5)] + load.sd[, c(1, 4, 5, 8)]))
# final model (3.8)
summary.lm(lm(RA[, 2] ~ 0 + load.mean[, c(1, 4, 5, 7, 9)] + load.sd[, c(1, 2, 4, 
    5, 9)]))
AIC(lm(RA[, 2] ~ 0 + load.mean[, c(1, 4, 5, 7, 9)] + load.sd[, c(1, 2, 4, 5, 9)]))
# m4 (3.9)
summary.lm(lm(RA[, 2] ~ 0 + load.mean[, c(1, 4, 5, 6)] + load.sd[, c(1, 4, 5, 8, 
    9)]))
AIC(lm(RA[, 2] ~ 0 + load.mean[, c(1, 4, 5, 6)] + load.sd[, c(1, 4, 5, 8, 9)]))
## out-of-sample forecasting
RA.hat.os = rep(0, 17)
# m1 (3.3)
for (s in 1:17) {
    RA.hat.os[s] = lm(RA[-s, 2] ~ 0 + load.mean[-s, c(1, 5)])$coefficients %*% c(load.mean[s, 
        c(1, 5)])
}
cor = cbind(RA[, 2], RA.hat.os)
cor(cor, method = "spearman")[1, 2]
cor(cor, method = "kendall")[1, 2]
# m2 (3.5)
for (s in 1:17) {
    RA.hat.os[s] = lm(RA[-s, 2] ~ 0 + load.mean[-s, c(1, 5)] + load.sd[-s, c(1, 4, 
        5, 9)])$coefficients %*% c(load.mean[s, c(1, 5)], load.sd[s, c(1, 4, 5, 9)])
}
cor = cbind(RA[, 2], RA.hat.os)
cor(cor, method = "spearman")[1, 2]
cor(cor, method = "kendall")[1, 2]
# m3 (3.7)
for (s in 1:17) {
    RA.hat.os[s] = lm(RA[-s, 2] ~ 0 + load.mean[-s, c(1, 5)] + load.sd[-s, c(1, 4, 
        5, 8)])$coefficients %*% c(load.mean[s, c(1, 5)], load.sd[s, c(1, 4, 5, 8)])
}
cor = cbind(RA[, 2], RA.hat.os)
cor(cor, method = "spearman")[1, 2]
cor(cor, method = "kendall")[1, 2]
# final model (3.8)
for (s in 1:17) {
    RA.hat.os[s] = lm(RA[-s, 2] ~ 0 + load.mean[-s, c(1, 4, 5, 7, 9)] + load.sd[-s, 
        c(1, 2, 4, 5, 9)])$coefficients %*% c(load.mean[s, c(1, 4, 5, 7, 9)], load.sd[s, 
        c(1, 2, 4, 5, 9)])
}
cor = cbind(RA[, 2], RA.hat.os)
cor(cor, method = "spearman")[1, 2]
cor(cor, method = "kendall")[1, 2]
# m4 (3.9)
for (s in 1:17) {
    RA.hat.os[s] = lm(RA[-s, 2] ~ 0 + load.mean[-s, c(1, 4, 5, 6)] + load.sd[-s, 
        c(1, 4, 5, 8, 9)])$coefficients %*% c(load.mean[s, c(1, 4, 5, 6)], load.sd[s, 
        c(1, 4, 5, 8, 9)])
}
cor = cbind(RA[, 2], RA.hat.os)
cor(cor, method = "spearman")[1, 2]
cor(cor, method = "kendall")[1, 2]

```

automatically created on 2018-09-04