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
#' @param regressionFunction # regression function to be used for FlexCode. Defaults to Random Forests. See FlexCode documentation for additional regression methods.
#' @param ... Additional arguments to FlexCoDE::fitFlexCoDE
#'
#' @return Returns an object of the class predictionBands with the following components:
#' \item{density_fit}{Object of the class FlexCoDE with the estimated density}
#' \item{cum_dist_evaluated_train}{Cumulative conditional distribution functions on the training set (for dist-split)}
#' \item{conformity_score_train}{Conformal scores on the training set (for cd-split)}
#' \item{t_grid}{Level sets of the densities}
#' \item{g_train}{Profiles of the training sample}
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
#'  plot(bands)
#'  plot(bands,ynew)
#'
#'  # CD-split
#'  bands <- predict(fit,xnew,type="cd")
#'  plot(bands)
#'  plot(bands,ynew)
fit_predictionBands <- function(x,y,
                                per_train=0.4,
                                per_val=0.1,
                                per_ths=1-per_train-per_val,
                                regressionFunction=FlexCoDE::regressionFunction.Forest,
                                ...)
{
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


  pred_train <- predict(fit,x[splits=="Threshold",])
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


  cum_dist_evaluated_train <- cum_dist(pred_train$z,pred_train$CDE,y[splits=="Threshold"])

  return_value <- NULL
  return_value$density_fit <- fit
  return_value$cum_dist_evaluated_train <- cum_dist_evaluated_train
  return_value$conformity_score_train <- conformity_score_train
  return_value$t_grid <- t_grid
  return_value$g_train <- g_train
  class(return_value) <-"predictionBands"
  return(return_value)
}

#' Compute conformal prediction bands based on density estimation on new samples
#'
#' @param cd_split_fit # Object of the class predictionBands, fitted using function \code{\link{fit_predictionBands}}
#' @param xnew # new covariates (one per row) where prediction bands are to be computed
#' @param type # type of prediction bands. dist for 'dist-split' and cd for 'cd-split'.
#' @param alpha # Miscoverage level (defaults to 10\%)
#' @param k # Number of nearest neighbors for cd-split. Default to 100
#'
#' @return Object of the class 'bands' with the following components:
#' \item{y_grid}{Grid of values for y}
#' \item{densities}{Matrix where each row contains the estimated density of the corresponding row of xnew on each value of y_grid}
#' \item{ths}{thresholds agains which to compare conformal scores.}
#' \item{prediction_bands_which_belong}{List with logical values with the same size as nrow(xnew); each element contains shows which elements of y_grid are contained on the prediction bands of the corresponding row of xnew}
#' \item{intervals}{List with the same size as nrow(xnew) that contains the prediction bands for each element of xnew. For cd-split, results are shown as a union of disjoint intervals.}
#' \item{type}{Type of prediction band}
#' \item{alpha}{Miscoverage level}
#' @export
#'
#' @examples See \code{\link{fit_predictionBands}}
predict.predictionBands <- function(cd_split_fit,xnew,type="dist",alpha=0.1,
                                    k=min(100,length(cd_split_fit$conformity_score_train)))
{
  pred_test <- predict(cd_split_fit$density_fit,xnew)


  if(type=="cd")
  {
    prediction_bands_which_belong <- list()
    intervals <- list()
    if(k<length(cd_split_fit$conformity_score_train))
    {
      ths <- rep(NA,length(cd_split_fit$conformity_score_train))

      g_test <- matrix(NA,nrow(xnew),length(cd_split_fit$t_grid))
      for(ii in 1:nrow(xnew))
      {
        g_test[ii,] <- profile_density(cd_split_fit$t_grid,
                                       pred_test$z,
                                       pred_test$CDE[ii,])
      }
      neighbors <- which_neighbors(cd_split_fit$g_train,
                                   g_test,k=k)
      for(ii in 1:nrow(xnew))
      {
        ths[ii] <- quantile(cd_split_fit$conformity_score_train[neighbors[ii,]],probs=alpha)
        prediction_bands_which_belong[[ii]] <- pred_test$CDE[ii,]>=ths[ii]
        intervals[[ii]] <- compute_intervals(prediction_bands_which_belong[[ii]],
                                             pred_test$z)
      }

    } else {
      ths <- quantile(cd_split_fit$conformity_score_train,
                      probs=alpha)
      for(ii in 1:nrow(xnew))
      {
        prediction_bands_which_belong[[ii]] <- pred_test$CDE[ii,]>=ths
        intervals[[ii]] <- compute_intervals(prediction_bands_which_belong[[ii]],
                                             pred_test$z)
      }
      ths <- rep(ths,nrow(xnew))
    }

    return_value <-list(y_grid=pred_test$z,ths=ths,densities=pred_test$CDE,
                        prediction_bands_which_belong=prediction_bands_which_belong,
                        intervals=intervals,ths=ths,type="cd",alpha=alpha,k=k)
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

  } else{
    stop(paste0("Type of distribution not implemented", " (",type,")"))
  }
  class(return_value) <- "bands"
  return(return_value)


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
    data_plot_list[[ii]] <- data.frame(x=ii,y=bands$y_grid[bands$prediction_bands_which_belong[[ii]]])
  }
  data_plot_regions <- do.call("rbind",data_plot_list)

  type <- ifelse(bands$type=="dist","Dist-split","CD-split")
  title <- paste0(type, "; alpha=",bands$alpha)
  if(bands$type=="cd")
  {
    title <- paste0(title," (k=",bands$k,")")
  }
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
  return(g)
}

cum_dist <- function(y_grid,cde_estimates,y_values)
{
  which_closest <- FNN::get.knnx(data = y_grid,
                                 query=y_values,k=1)$nn.index
  apply(as.matrix(1:nrow(cde_estimates)),1,function(xx){
    return(sum(cde_estimates[xx,1:which_closest[xx]])*diff(y_grid)[1])
  })
}


