#!/usr/bin/env Rscript
# NEXO / spatial-temporal statistical analysis of REAL imported measurements.
# Never manufactures observations; the entrypoint fails if the CSV is absent.
# Base R only. Status not executed unless Rscript is installed.
args <- commandArgs(trailingOnly=TRUE)
input <- if(length(args)>0) args[1] else 'data/processed/sima_normalized.csv'
output <- if(length(args)>1) args[2] else 'results'
if (!file.exists(input)) stop('No normalized REAL CSV found. Run acquisition and preprocessing first.')
set.seed(20261018)
dat <- read.csv(input, stringsAsFactors=FALSE)
required <- c('station_id','timestamp_source','variable','value','unit','latitude','longitude')
if(!all(required %in% names(dat))) stop('Canonical schema mismatch')
dat$time <- as.POSIXct(dat$timestamp_source, tz='UTC')  # IMPORTANT source TZ unresolved; use for ordering ONLY
if(anyNA(dat$time)) stop('Cannot interpret timestamps')
if(!all(is.finite(dat$latitude)&is.finite(dat$longitude))) stop('Bad station coordinates')
dat <- dat[!is.na(dat$value) & is.finite(dat$value),]
if(nrow(dat)<30) stop('Insufficient valid records')
dir.create(file.path(output,'tables'),recursive=TRUE,showWarnings=FALSE)
dir.create(file.path(output,'figures'),recursive=TRUE,showWarnings=FALSE)
# Source timezone must be resolved before inferring diurnal cycle or local meteorological matching.
summary <- aggregate(value ~ station_id+variable+unit,dat,function(x) c(n=length(x),mean=mean(x),sd=sd(x),p95=unname(quantile(x,.95))))
write.csv(summary,file.path(output,'tables','real_descriptives.csv'),row.names=FALSE)
pdf(file.path(output,'figures','real_distributions.pdf'),width=9,height=5)
boxplot(value~station_id,data=dat,las=2,main='Valores ambientales por estación (datos reales; si la fuente fue validada)',ylab='Unidades en tabla')
dev.off()
# A station-blocked split: all observations of each held-out station withheld from training.
stations <- unique(dat$station_id)
if(length(stations)>=6){
    held <- sample(stations,max(1,round(length(stations)*.25)))
    train <- dat[!(dat$station_id %in% held),]
    test <- dat[(dat$station_id %in% held),]
    # Baseline: mean by variable from TRAIN ONLY, not spatially extrapolated.
    m <- aggregate(value~variable,train,mean)
    names(m)[2] <- 'prediction'
    prediction <- merge(test,m,by='variable')
    if(nrow(prediction)>0){
        err <- prediction$value-prediction$prediction
        metrics <- data.frame(n=nrow(prediction),rmse=sqrt(mean(err^2)),mae=mean(abs(err)),
                              heldout_stations=paste(held,collapse=';'),model='train-variable-mean-baseline')
        write.csv(metrics,file.path(output,'tables','real_spatial_holdout_baseline.csv'),row.names=FALSE)
    }
}
# Multivariate correlation/PCA: complete-case station-time aligned values, with adequate common observations only.
wide <- reshape(dat[c('timestamp_source','station_id','variable','value')],
                idvar=c('timestamp_source','station_id'),timevar='variable',direction='wide')
measures <- names(wide)[grepl('^value\\.',names(wide))]
if(length(measures)>=2){
  X <- wide[measures]
  complete <- complete.cases(X)
  if(sum(complete)>=30){
    stable <- vapply(X[complete,,drop=FALSE],sd,0.)>0
    if(sum(stable)>=2){
      pc <- prcomp(X[complete,stable,drop=FALSE],center=TRUE,scale.=TRUE)
      write.csv(pc$rotation,file.path(output,'tables','real_pca_loadings.csv'))
    }
  }
}
cat('R processing completed; station-blocked baseline is distinct from geostatistical model.\n')

# Calendar-date aggregation can be performed on source dates without asserting a UTC offset.
dat$source_date <- as.Date(substr(dat$timestamp_source,1,10))
if(!all(is.na(dat$source_date))){
  daily <- aggregate(value~source_date+variable,dat,mean)
  write.csv(daily,file.path(output,'tables','real_daily_mean_source_dates.csv'),row.names=FALSE)
  pdf(file.path(output,'figures','real_daily_tendencies.pdf'),width=10,height=5)
  groups <- split(daily,daily$variable)
  plot(range(daily$source_date,na.rm=TRUE),range(daily$value,na.rm=TRUE),type='n',xlab='Fecha de la fuente (sin zona horaria resuelta)',ylab='Concentración (unidades del proveedor)',main='Medias diarias observadas (no hay inferencia causal)')
  colors <- seq_along(groups)
  for(i in seq_along(groups))lines(groups[[i]]$source_date,groups[[i]]$value,col=colors[i])
  legend('topright',legend=names(groups),col=colors,lty=1,cex=.7)
  dev.off()
}

# Cluster bootstrap across station averages; interpret as sampling variability only.
if(length(unique(dat$station_id))>=6){
  by_st <- aggregate(value~station_id+variable,dat,mean)
  by_par <- split(by_st,by_st$variable)
  ci <- do.call(rbind,lapply(names(by_par),function(v){
    a <- by_par[[v]]$value
    B <- 300L
    boot <- replicate(B,mean(sample(a,length(a),replace=TRUE)))
    data.frame(variable=v,stations=length(a),mean_of_station_means=mean(a),
      bootstrap_sd=sd(boot),percentile_025=unname(quantile(boot,.025)),
      percentile_975=unname(quantile(boot,.975)),
      warning='Uncalibrated station bootstrap; not predictive interval')
  }))
  write.csv(ci,file.path(output,'tables','real_station_bootstrap.csv'),row.names=FALSE)
}

# Chronological holdout, within each pollutant; train dates preceding test dates.
if(length(unique(dat$source_date[!is.na(dat$source_date)]))>15){
  dates <- sort(unique(dat$source_date[!is.na(dat$source_date)]))
  border <- dates[floor(.75*length(dates))]
  tr <- dat[!is.na(dat$source_date)&dat$source_date<=border,]
  te <- dat[!is.na(dat$source_date)&dat$source_date>border,]
  mu <- aggregate(value~variable,tr,mean);names(mu)[2]<-'prediction'
  tt <- merge(te,mu,by='variable')
  if(nrow(tt)>0){
    err <- tt$value-tt$prediction
    write.csv(data.frame(train_end_as_reported=as.character(border),
       model='variable-mean-train-period',n=nrow(tt),rmse=sqrt(mean(err^2)),mae=mean(abs(err))),
       file.path(output,'tables','real_time_holdout_baseline.csv'),row.names=FALSE)
  }
}
