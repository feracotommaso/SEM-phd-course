library(lavaan)
# manic data: build the variance covariance matrix from basic paper info
manic.means <- c(10.09, 12.07, 10.25, 9.96, 10.90, 11.24, 10.30, 10.44) #Vector of means
manic.sd    <- c(3.06,   3.53,  3.18, 2.85,  2.49,  3.95,  3.35,  3.13) #Vector of sds
manic.cor   <- matrix (nrow = 8, ncol = 8) #correlation matrix ...
manic.cor[lower.tri(manic.cor, diag = TRUE )] <- c (1.00 , 0.72 , 0.66 , 0.65 , 0.40 , 0.29 , 0.36 , 0.38 , 
                                                    1.00 , 0.78 , 0.74 , 0.56 , 0.35 , 0.46 , 0.42 , 
                                                    1.00 , 0.75 , 0.57 , 0.39 , 0.49 , 0.49 , 
                                                    1.00 , 0.58 , 0.46 , 0.33 , 0.40 , 
                                                    1.00 , 0.52 , 0.43 , 0.49 , 
                                                    1.00 , 0.47 , 0.47 , 
                                                    1.00 , 0.62 , 
                                                    1.00) #write the lower triangle of the matrix
manic.cor[upper.tri(manic.cor, diag = TRUE )]<- t(manic.cor)[upper.tri(manic.cor, diag = TRUE)] #copy it to the upper triangle
manic.cov <- cor2cov(manic.cor, manic.sd) #from correlations do covariances

# norming data: build the variance covariance matrix from basic paper info
norming.means <- c(10.10, 10.30, 9.80, 10.10, 10.10, 10.10, 9.90, 10.20)
norming.sd    <- c( 3.10,  2.90, 3.00,  2.90,  3.20,  3.30, 3.40,  3.30)
norming.cor   <- matrix(nrow =8 , ncol =8)
norming.cor[lower.tri(norming.cor, diag = TRUE)]<-c (1.00 , 0.65 , 0.68 , 0.49 , 0.45 , 0.34 , 0.50 , 0.42 , 
                                                     1.00 , 0.72 , 0.53 , 0.49 , 0.31 , 0.50 , 0.48 , 
                                                     1.00 , 0.58 , 0.47 , 0.30 , 0.40 , 0.44 , 
                                                     1.00 , 0.40 , 0.30 , 0.33 , 0.33 , 
                                                     1.00 , 0.36 , 0.48 , 0.47 , 
                                                     1.00 , 0.32 , 0.33 , 
                                                     1.00 , 0.59 , 
                                                     1.00)
norming.cor[upper.tri( norming.cor, diag = TRUE )]<-t(norming.cor)[upper.tri(norming.cor, diag = TRUE)]
norming.cov <- cor2cov (norming.cor, norming.sd)

# WE NOW HAVE TWO VARIANCE-COVARIANCE MATRICES.
# WHAT IS MISSING?




# ... the data
# Simulate data with N participants per group
library(MASS) # we will use the MASS package to simulate multivariate normals
set.seed(3)   # set the seed for reproducibility
d.man=mvrnorm(n=150,manic.means,manic.cov,empirical=T) # simulate N manic data?
d.nor=mvrnorm(n=150,norming.means,norming.cov,empirical=T) # simulate N norming data?
dati=rbind(d.man,d.nor) #merge the data
sum(dati<0) # Negative values are impossible. Set them to 0
dati[dati < 0] <- 0
# Create the dmg data frame with id, diagnosis, and colnames
dmg=data.frame(id=c(1:nrow(dati)),diagnosis=factor(rep(c("manic","norming"),c(150,150))),dati)
names(dmg)[3:10]=c("Info","Sim","Vocab","Comp","PicComp","PicArr","BlkDsgn","ObjAsmb")

#Save the data
# save(dmg,file="05.Invariance/data/dmg.RData")
