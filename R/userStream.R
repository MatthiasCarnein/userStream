#' Reference Class userStream
#'
#' Reference class mostly used to expose C class object
#'
#' @field C exposed C class
#'
#' @author Matthias Carnein \email{matthias.carnein@@uni-muenster.de}
#'
#' @export userStream
userStream <- setRefClass("userStream_R", fields = list(
  C ="ANY",
  n = "integer",
  upToDate = "logical",
  macro = "matrix",
  macroWeights = "numeric",
  macroAssignment = "integer",
  k = "integer",
  weighted = "logical"
))

userStream$methods(
  initialize = function(r, lambda=0, tgap=1000, k, weighted=TRUE) {
    C <<- new(UserStream, r, lambda, tgap) ## Exposed C class
    n <<- 0L
    upToDate <<- FALSE
    k <<- as.integer(k)
    weighted <<- weighted
    macro <<- as.matrix(data.frame())
    macroWeights <<- numeric()
    .self
  }
)

userStream$methods(
  serialize = function(){
    list(C=.self$C$serialize(), n = .self$n, upToDate=.self$upToDate, k=.self$k, macro=.self$macro, macroWeights=.self$macroWeights, weighted=.self$weighted)
  }
)


userStream$methods(
  deserialize = function(serialized) {
    C <<- new(UserStream, serialized[["C"]])
    n <<- serialized[["n"]]
    upToDate <<- serialized[["upToDate"]]
    k <<- serialized[["k"]]
    weighted <<- serialized[["weighted"]]
    macro <<- serialized[["macro"]]
    macroWeights <<- serialized[["macroWeights"]]
  }
)


userStream$methods(
  cluster = function(dataframe, userIds, times=NULL){
    dataframe=as.matrix(dataframe)
    userIds = as.matrix(userIds)

    if(!is.null(times)) times = as.matrix(times)
    for(i in seq_len(nrow(dataframe))){
      .self$update(dataframe[i,], userIds[i], times[i])
    }
  }
)


userStream$methods(
  update = function(newUser, userId, time=NULL){

    upToDate <<- FALSE
    n <<- .self$n + 1L

    ## allow all kinds of time formats
    if(is.character(time)){ ## as string
      time = as.integer(as.POSIXct(strptime(time, "%Y-%m-%d %H:%M:%S")))
    } else if (is(time, "POSIXlt")){ ## POSIXlt object
      time = as.integer(as.POSIXct(time))
    } else if (is(time, "POSIXct")){ ## POSIXct object
      time = as.integer(time)
    } else if(is.null(time)){ ## NULL for internal count
      time = .self$n
    } else if(is.numeric(time)){ ## directly as numeric time
      time = as.integer(time)
    }

    .self$C$update(newUser, userId, time)
  }
)


userStream$methods(
  get_microclusters = function(unscale=TRUE) {
    as.data.frame(.self$C$get_microclusters(unscale))
  }
)

userStream$methods(
  get_microweights = function() {
    .self$C$get_microweights()
  }
)

userStream$methods(
  plotCentres = function(i, data=NULL, assignment, last){

    mc = .self$C$get_microclusters(unscale=T)
    weights = .self$C$get_microweights()
    if(length(unique(weights))==1){
      weights <- rep(2.5, length(weights))
    } else{
      weights = (weights - min(weights))/(max(weights)-min(weights)) * (5 - 1) + 1
    }

    ## only works for two dimensions
    col="gray" ## assignment
    plot(as.matrix(last), pch=19, cex=1, col=col, main=i)
    col="red" # 1:nrow(mc)
    points(mc, pch=1, cex=weights, col=col)


    if(!is.null(.self$k)){
      macroC = .self$get_macroclusters()
      macroW = .self$get_macroweights()

      if(length(unique(macroW))==1){
        macroW <- rep(2.5, length(macroW))
      } else{
        macroW = (macroW - min(macroW))/(max(macroW)-min(macroW)) * (5 - 1) + 1
      }

      col="blue" # 1:nrow(macroC)
      points(macroC, pch=3, cex=macroW, col=col)
    }
  }
)



userStream$methods(
  SSQ = function(assignment, last, centres){

    sum(sapply(seq_along(assignment), function(i){
      sqrt(sum((last[i,] - centres[assignment[i],])^2))
    }))
  }
)

userStream$methods(
  silhouette = function(assignment, last){

    if(length(unique(assignment))==1) return(1)
    mean(cluster::silhouette(assignment, dist(last))[,"sil_width"])
  }
)



userStream$methods(
  batch_iterate = function(data, userIds, dates, clients=NULL, horizon=1000, plot=F, sleep=0, evaluate=T, prequential=F){

    n_ = nrow(data)
    do.call(rbind,lapply(seq(1,n_,by=horizon), function(i){

      allRange = 1:min((i+horizon-1),n_)
      range = i:min((i+horizon-1),n_)

      allLast = aggregate(data[allRange,], by = list(userIds[allRange]), FUN = tail, n = 1) ## get most recent observation of each customer
      last = aggregate(data[range,], by = list(userIds[range]), FUN = tail, n = 1) ## get most recent observation for each customer of last batch

      if(prequential){

        means = apply(last[,-1], 2, mean)
        dev = apply(last[,-1], 2, sd)
        dev[dev==0]=1
        scaledLast = matrix(apply(last[,-1], 1, function(x){
          (x - means) / dev
        }), nrow = nrow(last), byrow = T)

        ## before insertion
        assignment = cbind(last[,1], .self$C$getClosest(as.matrix(scaledLast)))
        if(evaluate){

          if(i!=1){
            evalResult = c(
              "points" = min((i+horizon-1),n_),
              "SSQ_micro" = .self$SSQ(assignment[,2], scaledLast, .self$get_microclusters()),
              "SSQ_macro" = .self$SSQ(.self$microToMacro()[assignment[,2]], scaledLast, .self$get_macroclusters()),
              "silhouette_micro" = .self$silhouette(assignment[,2], scaledLast),
              "silhouette_macro" = .self$silhouette(.self$microToMacro()[assignment[,2]], scaledLast),
              "numMicroClusters" = nrow(.self$get_microclusters()),
              "numMacroClusters" = nrow(.self$get_macroclusters())
            )
          } else{
            evalResult = NULL
          }}


      }

      .self$cluster(data[range,], userIds[range], dates[range])

      allAssignment = .self$C$getAssignment() ## after insertion

      if(!prequential){

        ## filter customers that are no longer assigned (because the cluster was removed) Alternaitvely: assign customer to closest cluster
        assignment = allAssignment[allAssignment[,2] >=0 ,]
        relevant = intersect(last[,1], assignment[,1])
        last = last[last[,1] %in% relevant,]
        assignment = assignment[assignment[,1] %in% relevant,]


        means = apply(last[,-1], 2, mean)
        dev = apply(last[,-1], 2, sd)
        dev[dev==0]=1
        scaledLast = matrix(apply(last[,-1], 1, function(x){
          (x - means) / dev
        }), nrow = nrow(last), byrow = T)

        if(evaluate){
          evalResult = c(
            "points" = min((i+horizon-1),n_),
            "SSQ_micro" = .self$SSQ(assignment[,2], scaledLast, .self$get_microclusters()),
            "SSQ_macro" = .self$SSQ(.self$microToMacro(), scaledLast, .self$get_macroclusters()),
            "silhouette_micro" = .self$silhouette(assignment[,2], scaledLast),
            "silhouette_macro" = .self$silhouette(.self$microToMacro()[assignment[,2]], scaledLast),
            "numMicroClusters" = nrow(.self$get_microclusters()),
            "numMacroClusters" = nrow(.self$get_macroclusters())
          )
        }
      }


      if(!is.null(dates)){
        lastDates = aggregate(dates[allRange], by = list(userIds[allRange]), FUN = tail, n = 1)$x ## get dates in same order as last
      } else{
        lastDates=NULL
      }

      if(!is.null(clients)){
        clients = aggregate(clients[allRange,], by = list(userIds[allRange]), FUN = tail, n = 1)[,-1] ## get dates in same order as last
      }

      if(plot==TRUE) .self$plotCentres(min((i+horizon-1),n_), data, allAssignment[,2], allLast[,-1])

      if(sleep!=0) Sys.sleep(sleep)


      if(evaluate){
        return(evalResult)
      }

    }))
  }
)


userStream$methods(
  recluster = function(k){
    mc = .self$C$get_microclusters(unscale=T)

    means = apply(mc, 2, mean)
    dev = apply(mc, 2, sd)
    dev[dev==0]=1

    mc = matrix(apply(mc, 1, function(x){
      (x - means) / dev
    }), nrow = nrow(mc), byrow = T)

    weights = .self$C$get_microweights()

    k_actual = min(c(nrow(mc), .self$k))

    if(.self$weighted){

      clustmem <- sample.int(n = k_actual, size = nrow(mc), replace = TRUE) ## random assignment
      km = lsbclust:::KMeansW(start = clustmem, data = mc, weight = weights, nclust = k_actual) ## weighted kmeans

      macroWeights <<- sapply(seq_len(k_actual), function(x){
        sum(weights[km$cluster == x])
      })

      macro <<- km$centers
      macro <<- matrix(apply(macro, 1, function(x){
        (x * dev) + means
      }), nrow = nrow(macro), byrow = T)
      macroAssignment <<- km$cluster
    } else{

      km = kmeans(mc, .self$k)

      macro <<- km$centers
      macro <<- matrix(apply(macro, 1, function(x){
        (x * dev) + means
      }), nrow = nrow(macro), byrow = T)
      macroWeights <<- sapply(seq_len(.self$k), function(x){
        sum(km$cluster == x)
      })
      macroAssignment <<- km$cluster
    }


  }
)


userStream$methods(
  get_macroclusters = function(){
    if(!.self$upToDate){
      upToDate <<- TRUE
      .self$recluster(.self$k)
    }
    .self$macro
  }
)

userStream$methods(
  get_macroweights = function(){
    if(!.self$upToDate){
      upToDate <<- TRUE
      .self$recluster(.self$k)
    }
    .self$macroWeights
  }
)



userStream$methods(
  microToMacro = function() {
    if(!.self$upToDate){
      upToDate <<- TRUE
      .self$recluster(.self$k)
    }
    return(.self$macroAssignment)
  }
)
