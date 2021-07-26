
# finds k nearest neighbors in xTrain of each xTest
which_neighbors <-  function(xTrain,xTest,k){
  return(FNN::get.knnx(data=xTrain,query = xTest,k=k)$nn.index)
}

profile_density <- function(t_grid,y_grid,cde_estimate)
{
  v2 <- cde_estimate[order(cde_estimate)]
  v2s <- rev(cumsum(rev(v2)))*(y_grid[2]-y_grid[1])
  v2s <- v2s[findInterval(t_grid, v2) + 1]
  v2s[which(is.na(v2s))] <- 0
  return(v2s)
}

#' Fit conformal prediction bands based on density estimation for regression
#'
#' @param x Matrix with covariates of training set
#' @param y Vector with the (continuous) responses of training set
#' @param per_train # percentage of samples used for traning density estimator (defaults to 40\%)
#' @param per_val # percentage of samples used for tuning density estimator (defaults to 10\%)
#' @param per_ths  # percentage of samples used for computeing thresholds for the conformal method (defaults to 50\%)
#' @param k # Number of clusters for cd-split. Default to round(per_ths*nrow(as.matrix(x))/100) so that each cluster has on average 100 samples
#' @param regressionFunction # regression function to be used for FlexCode. Defaults to Random Forests. See FlexCode documentation for additional regression methods.
#' @param ... Additional arguments to FlexCoDE::fitFlexCoDE
#'
#' @return Returns an object of the class predictionBands with the following components:
#' \item{density_fit}{Object of the class FlexCoDE with the estimated density}
#' \item{cum_dist_evaluated_train}{Cumulative conditional distribution functions on the training set (for dist-split)}
#' \item{conformity_score_train}{Conformal scores on the training set (for cd-split)}
#' \item{conformity_score_train_hpd}{Conformal scores on the training set (for hpd-split)}
#' \item{t_grid}{Level sets of the densities}
#' \item{g_train}{Profiles of the training sample}
#' \item{center_kmeans}{The center of the clusters found by kmeans (in the profile space)}
#' @export
#'
#' @examples
#'
#' # generate data
#' n <- 1000
#' n_new <- 50
#' d <- 10
#' data <- matrix(NA,n,d+1)
#' x <- matrix(rnorm(n*d),n,d)
#' y <- x[,1]+rnorm(n,0,0.1)
#' fit <- fit_predictionBands(x,y,0.5,0.4,0.1)
#'
#' xnew <- matrix(rnorm(n_new*d),n_new,d)
#' ynew <- xnew[,1]+rnorm(n_new,0,0.1)
#'
#'  # Dist-split
#'  bands <- predict(fit,xnew,type="dist")
#'  bands[[1]]
#'  bands[[2]]
#'  plot(bands)
#'  plot(bands,ynew)
#'
#'  # CD-split
#'  bands <- predict(fit,xnew,type="cd")
#'  bands[[1]]
#'  bands[[2]]
#'  plot(bands)
#'  plot(bands,ynew)
fit_predictionBands <- function(x,y,
                                per_train=0.4,
                                per_val=0.1,
                                per_ths=1-per_train-per_val,
                                k=max(round(per_ths*nrow(as.matrix(x))/100),1),
                                regressionFunction=FlexCoDE::regressionFunction.Forest,
                                ...)
{
  x <- as.matrix(x)
  n_levels <- 1000
  splits <- sample(c("Train","Validation","Threshold"),size = nrow(x),
                   prob=c(per_train,per_val,per_ths),
                   replace=TRUE)

  fit <- FlexCoDE::fitFlexCoDE(xTrain=x[splits=="Train",],
                               zTrain=y[splits=="Train"],
                               xValidation=x[splits=="Validation",],
                               zValidation=y[splits=="Validation"],
                               regressionFunction = regressionFunction,
                               ...)

  pred_train_cde <- FlexCoDE::predict.FlexCoDE(fit,x[splits!="Threshold",])
  t_grid <- seq(0,max(pred_train_cde$CDE),length.out = n_levels)
  g_train_cde <- matrix(NA,nrow(pred_train_cde$CDE),
                        length(t_grid))
  for(ii in 1:nrow(pred_train_cde$CDE))
  {
    g_train_cde[ii,] <- profile_density(t_grid,pred_train_cde$z,
                                        pred_train_cde$CDE[ii,])
  }

  kmeans_result <- try(kmeanspp(g_train_cde,k=k),
                       silent = TRUE)
  if(class(kmeans_result)=="try-error")
  {
    kmeans_result <- kmeans(g_train_cde,centers = k)
  }
  centers_kmeans <- kmeans_result$centers
  rm(g_train_cde)
  rm(pred_train_cde)


  pred_train <- FlexCoDE::predict.FlexCoDE(fit,x[splits=="Threshold",])  ## done as FlexCoDE didn't use .S3method

  # hpd-split
  which_select <- cbind(1:length(y[splits=="Threshold"]),
                        which_neighbors(as.matrix(pred_train$z),
                                        as.matrix(y[splits=="Threshold"]),
                                        k=1))
  which_smaller <- apply(pred_train$CDE<=pred_train$CDE[which_select],1,which)
  conformity_score_train_hpd <- rep(NA,nrow(pred_train$CDE))
  for(ii in 1:nrow(pred_train$CDE))
  {
    conformity_score_train_hpd[ii] <- sum(pred_train$CDE[ii,which_smaller[[ii]]])
  }
  band <- diff(pred_train$z)[1]
  conformity_score_train_hpd <- conformity_score_train_hpd*band


  t_grid <-seq(0,max(pred_train$CDE),length.out = n_levels)
  # observed densities:
  which_select <- cbind(1:length(y[splits=="Threshold"]),
                        which_neighbors(as.matrix(pred_train$z),
                                        as.matrix(y[splits=="Threshold"]),
                                        k=1))
  conformity_score_train <- pred_train$CDE[which_select]

  # Profiles
  g_train <- matrix(NA,length(conformity_score_train),
                    length(t_grid))
  for(ii in 1:length(conformity_score_train))
  {
    g_train[ii,] <- profile_density(t_grid,pred_train$z,
                                    pred_train$CDE[ii,])
  }


  cum_dist_evaluated_train <- cum_dist(pred_train$z,pred_train$CDE,
                                       y[splits=="Threshold"])

  return_value <- NULL
  return_value$density_fit <- fit
  return_value$cum_dist_evaluated_train <- cum_dist_evaluated_train
  return_value$conformity_score_train <- conformity_score_train
  return_value$conformity_score_train_hpd <- conformity_score_train_hpd
  return_value$t_grid <- t_grid
  return_value$band <- band
  return_value$g_train <- g_train
  return_value$centers_kmeans <- centers_kmeans
  class(return_value) <-"predictionBands"
  return(return_value)
}

#' Compute conformal prediction bands based on density estimation on new samples
#'
#' @param cd_split_fit # Object of the class predictionBands, fitted using function \code{\link{fit_predictionBands}}
#' @param xnew # new covariates (one per row) where prediction bands are to be computed
#' @param type # type of prediction bands. dist for 'dist-split' and cd for 'cd-split'.
#' @param alpha # Miscoverage level (defaults to 10\%)
#'
#' @return Object of the class 'bands' with the following components:
#' \item{y_grid}{Grid of values for y}
#' \item{densities}{Matrix where each row contains the estimated density of the corresponding row of xnew on each value of y_grid}
#' \item{ths}{thresholds agains which to compare conformal scores.}
#' \item{prediction_bands_which_belong}{List with logical values with the same size as nrow(xnew); each element contains shows which elements of y_grid are contained on the prediction bands of the corresponding row of xnew}
#' \item{intervals}{List with the same size as nrow(xnew) that contains the prediction bands for each element of xnew. For cd-split, results are shown as a union of disjoint intervals.}
#' \item{type}{Type of prediction band ("cd", "dist" or "hpd")}
#' \item{alpha}{Miscoverage level}
#' @export
#'
#' @examples See \code{\link{fit_predictionBands}}
predict.predictionBands <- function(cd_split_fit,xnew,type="dist",alpha=0.1)
{
  pred_test <- FlexCoDE::predict.FlexCoDE(cd_split_fit$density_fit,xnew)

  if(type=="cd")
  {
    prediction_bands_which_belong <- list()
    intervals <- list()

    ths <- rep(NA,length(cd_split_fit$conformity_score_train))
    g_test <- matrix(NA,nrow(xnew),length(cd_split_fit$t_grid))
    for(ii in 1:nrow(xnew))
    {
      g_test[ii,] <- profile_density(cd_split_fit$t_grid,
                                     pred_test$z,
                                     pred_test$CDE[ii,])
    }
    which_partition_test <- which_neighbors(cd_split_fit$centers_kmeans,
                                            g_test,1)
    which_partition_train <- which_neighbors(cd_split_fit$centers_kmeans,
                                             cd_split_fit$g_train,1)

    ths_partition <- rep(NA,nrow(cd_split_fit$centers_kmeans))
    for(ii in 1:nrow(cd_split_fit$centers_kmeans))
    {
      ths_partition[ii] <-  quantile(cd_split_fit$conformity_score_train[which_partition_train==ii],
                                     probs=alpha)
    }
    ths_partition[is.na(ths_partition)] <- quantile(cd_split_fit$conformity_score_train,probs=alpha)
    ths <- ths_partition[which_partition_test]
    for(ii in 1:nrow(xnew))
    {
      prediction_bands_which_belong[[ii]] <- pred_test$CDE[ii,]>=ths[ii]
      intervals[[ii]] <- compute_intervals(prediction_bands_which_belong[[ii]],
                                           pred_test$z)
    }

    return_value <- list(y_grid=pred_test$z,ths=ths,densities=pred_test$CDE,
                         prediction_bands_which_belong=prediction_bands_which_belong,
                         intervals=intervals,ths=ths,type="cd",alpha=alpha)

  } else if(type=="dist"){
    ths <-  quantile(cd_split_fit$cum_dist_evaluated_train,
                     probs = c(alpha/2,1-alpha/2))
    prediction_bands_which_belong <- list()
    intervals <- list()
    FTest <- matrix(NA,nrow(xnew),length(pred_test$z))
    for (ii in 1:nrow(xnew)){
      FTest[ii,] <- cumsum(pred_test$CDE[ii,])*diff(pred_test$z)[1]
      prediction_bands_which_belong[[ii]] <- FTest[ii,]>=ths[1]&FTest[ii,]<=ths[2]
      intervals[[ii]] <- paste0("(",min(pred_test$z[prediction_bands_which_belong[[ii]]]),
                                ",",max(pred_test$z[prediction_bands_which_belong[[ii]]]),")")
    }
    return_value <-list(y_grid=pred_test$z,ths=ths,densities=pred_test$CDE,
                        prediction_bands_which_belong=prediction_bands_which_belong,
                        intervals=intervals,type="dist",alpha=alpha)

  }   else if(type=="hpd") {
    th <- quantile(cd_split_fit$conformity_score_train_hpd,probs=alpha)
    prediction_bands_which_belong <- list()
    intervals <- list()
    th_hpd <- rep(NA,nrow(xnew))
    for(ii in 1:nrow(xnew))
    {
      th_hpd[ii] <- findThresholdHPD(cd_split_fit$band,pred_test$CDE[ii,],1-th)
      prediction_bands_which_belong[[ii]] <- pred_test$CDE[ii,]>=th_hpd[ii]
      intervals[[ii]] <- compute_intervals(prediction_bands_which_belong[[ii]],
                                           pred_test$z)
    }
    return_value <-list(y_grid=pred_test$z,th_hpd=th_hpd,densities=pred_test$CDE,intervals=intervals,
                        prediction_bands_which_belong=prediction_bands_which_belong,type="hpd",alpha=alpha)
  }  else {
    stop(paste0("Type of distribution not implemented", " (",type,")"))
  }
  class(return_value) <- "bands"
  return(return_value)


}


findThresholdHPD=function(binSize,estimates,confidence)
{
  estimates=as.vector(estimates)
  maxDensity=max(estimates)
  minDensity=min(estimates)
  newCut=(maxDensity+minDensity)/2
  eps=1
  ii=1
  while(ii<=1000)
  {
    prob=sum(binSize*estimates*(estimates>newCut))
    eps=abs(confidence-prob)
    if(eps<0.0000001) break; # level found
    if(confidence>prob) maxDensity=newCut
    if(confidence<prob) minDensity=newCut
    newCut=(maxDensity+minDensity)/2
    ii=ii+1
  }
  return(newCut)
}


compute_intervals <- function(sequence,t_grid) {
  aux <- rle(sequence)
  cum_aux <- cumsum(aux$lengths)
  if(aux$values[1]==TRUE)
  {
    start <- cum_aux[aux$values==FALSE]
    start <-c(1,start[1:sum(aux$values==TRUE)]+1)
    start <- start[-length(start)]
    end <- cum_aux[aux$values==TRUE]
  } else {
    start <- cum_aux[aux$values==FALSE]
    start <-start[1:sum(aux$values==TRUE)]+1
    end <- cum_aux[aux$values==TRUE]
  }
  return(paste(paste0("(",paste(t_grid[start],
                                t_grid[end],sep=","),")"),
               collapse =  "U"))
}


#' Plots prediction bands for new observations
#'
#' @param bands # Predicted bands fitted using \code{\link{predict.predictionBands}}
#' @param ynew # (optional) response values of new samples
#'
#' @return # ggplot object with a plot of the predicition bands (together with ynew if these are available)
#' @export
#'
#' @examples  See \code{\link{fit_predictionBands}}
plot.bands <- function(bands,ynew=NULL)
{
  data_plot_list <- NULL
  for(ii in 1:length(bands$prediction_bands_which_belong))
  {
    pred_region <- bands$y_grid[bands$prediction_bands_which_belong[[ii]]]
    if(length(pred_region)==0)
    {
      data_plot_list[[ii]] <- data.frame(x=ii,y=median(bands$y_grid))
      next;
    }
    data_plot_list[[ii]] <- data.frame(x=ii,y=pred_region)
  }
  data_plot_regions <- do.call("rbind",data_plot_list)

  if(bands$type=="dist")
  {
    type <- "Dist-split"
  } else if(bands$type=="cd") {
    type <- "CD-split"
  } else {
    type <- "HPD-split"
  }
  title <- paste0(type, "; alpha=",bands$alpha)
  g <- ggplot()+
    geom_point(data=data_plot_regions,
               aes(x=x,y=y),alpha=0.02,color="red")+
    theme_bw()+
    xlab("Sample id")+
    theme(axis.text=element_text(size=20),
          axis.title=element_text(size=24,face="bold"),
          legend.title = element_blank(),
          legend.position = c(0.01, 0.3),
          legend.justification = c(0, 1),
          plot.title = element_text(size = 28, face = "bold"),
          legend.text=element_text(size=22))+
    ggtitle(title)

  if(!is.null(ynew))
  {
    data_plot_y <-data.frame(x=1:length(bands$prediction_bands_which_belong),y=ynew)
    g <- g+geom_point(data=data_plot_y,
                      aes(x=x,y=y))
  }


  plot(g)
  invisible(g)
}


cum_dist <- function(y_grid,cde_estimates,y_values)
{
  which_closest <- FNN::get.knnx(data = y_grid,
                                 query=y_values,k=1)$nn.index
  apply(as.matrix(1:nrow(cde_estimates)),1,function(xx){
    return(sum(cde_estimates[xx,1:which_closest[xx]])*diff(y_grid)[1])
  })
}

