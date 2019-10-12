# predictionBands
R package to compute distribution-free prediction bands using density estimators according to

Izbicki, R., Shimizu, G. Y., Stern, R. B. (2019). Distribution-free conditional predictive bands using density estimators.

The package estimates conditional densities using [FlexCode](https://github.com/rizbicki/FlexCoDE/). (More on FlexCoDE: Izbicki, R.; Lee, A.B. [Converting High-Dimensional Regression to High-Dimensional Conditional Density Estimation](https://projecteuclid.org/euclid.ejs/1499133755). Electronic Journal of Statistics, 2017)

Two types of bands are available: 'dist-split' returns intervals (ideal for unimodal response distributions); 'cd-split' returns unions of intervals (ideal for multimodal response distributions)

To install the package, run

```R
# install.packages("devtools")
devtools::install_github("rizbicki/FlexCoDE")
devtools::install_github("rizbicki/predictionBands")
```
(this package requires FlexCoDE).

A simple example:

```R
#generate data
n<-1000
n_new<-50
d<-10
data<-matrix(NA,n,d+1)
x<-matrix(rnorm(n*d),n,d)
y<-x[,1]+rnorm(n,0,0.1)
fit<-fit_predictionBands(x,y,0.5,0.4,0.1)

xnew<-matrix(rnorm(n_new*d),n_new,d)
ynew<-xnew[,1]+rnorm(n_new,0,0.1)

#Dist-split
bands<-predict(fit,xnew,type="dist")
bands[[1]]
bands[[2]]
plot(bands)
plot(bands,ynew) # if ynew is available

#CD-split
bands<-predict(fit,xnew,type="cd")
bands[[1]]
bands[[2]]
plot(bands)
plot(bands,ynew) # if ynew is available
```

