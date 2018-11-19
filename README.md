# userStream - Customer Segmentation using Stream Clustering

userStream is a new clustering algorithm that is applicable to customer segmentation.
Customer segmentation is one of the most important tools in marketing.
It divides a market of potential customers into distinct groups of customers that share similar behaviour or preference.
This allows to create suitable value propositions for each customer segment.
A major drawback of customer segmentation is that it is usually a snapshot analysis where segments are identified at a specific point in time.
To levitate this problem we propose a new stream clustering algorithm which processes a stream of transactions in order to identify and track customer segments over time.

## Installation

The easiest way to install the package is by using devtools and pointing it to the corresponding GitLab:

```R
devtools::install_git("https://wiwi-gitlab.uni-muenster.de/m_carn01/userStream")
```


## Usage

After loading the library, a new userstream object needs to created. The library makes use of R's Reference Classes which allows to initialize the object using the function `new()`. Using the corresponding parameters, the radius threshold (`r`), decay factor (`lambda`), decay time (`tgap`) and number of macro clusters (`k`) can be configured. In addition `weighted` determines whether to perform unweighted or weighted reclustering.

The created object then provides two main methods:
1. The function `cluster()` allows to cluster an entire data frame in a single run and is the fastes method to calculate the end-result of the stream clustering algorithm.
2. The function `batch_iterate()` allows to cluster the data frame in batches of fixed size. This allows to produce intermediate results, calculate quality measures and plot the result.

Both methods require 3 main parameters: the usage-related features, the corresponding customer ID as well as the Date of the transaction. `batch_iterate` further allows to configure the size of the batches, whether to plot and evaluate the intermediate results. In addition, the parameter `prequential` determines whether the quality measures are calculated before or after the insertion of the new points (the former is usually referred to as prequential evaluation or train-then-test evaluation).

Therefore the main usage scenario is as follows:
```R
library(userStream)

userstream = userStream$new(r=1, lambda=0.01, tgap=30, k=5, weighted=F)

userstream$cluster(data[,c("orderSize", "numPurchases")], data[,"Customer"], data[,"Date"])

userstream$batch_iterate(data[,c("orderSize", "numPurchases")], data[,"Customer"], data[,"Date"], clients, export=F, horizon=1000, plot=T, evaluate=T, prequential=F, sleep=0)
```




## Full Example

The following demonstrates a full example how to segments customers based on a transaction log using the retail data from: https://archive.ics.uci.edu/ml/datasets/online+retail

First, the data is loaded and a few usage-related features are computed over time (Return Rate, Number of Purchases, Mean Order Size and Recency of Last Purchase). Then, two of these possible features (Order Size and Number of Purchases) are used to run the clustering algorithm.

```R
# Load Data ---------------------------------------------------------------

library(openxlsx)

## download to temporary directory
tmp = tempfile(fileext = ".xlsx")
download.file(url = "https://archive.ics.uci.edu/ml/machine-learning-databases/00352/Online Retail.xlsx", destfile = tmp, mode="wb")
data = read.xlsx(tmp) ## read file
data$InvoiceDate = trunc(as.integer(convertToDateTime(data$InvoiceDate))/60/60/24) ## convert to date time in days
data = data[!is.na(data$CustomerID),] ## only use customers with ID

## sort chronologically
data = data[order(data$InvoiceDate),]

## select first entries for demonstration purposes
data = data[1:10000,]




# Create Usage-Related Features over Time ---------------------------------

pb <- txtProgressBar(min=1, max=nrow(data))
features = do.call(rbind,lapply(seq_len(nrow(data)), function(i){
  setTxtProgressBar(pb, i)
  previous = data[seq_len(i),] # previous purchases
  previous = previous[previous$CustomerID == data$CustomerID[i],] # from customer

  returnRate = (-1*sum(previous$Quantity[previous$Quantity<=0])) / sum(previous$Quantity[previous$Quantity > 0]) ## returned percentage of products
  if(returnRate==Inf) return(NULL)
  numPurchases = nrow(previous) ## frequency
  orderSize = mean(previous$Quantity) ## mean order size
  recency = data$InvoiceDate[i]-data$InvoiceDate[1] ## last purchase as days since start

  Date = data$InvoiceDate[i]
  Customer = data$CustomerID[i]
  data.frame(Customer, Date, returnRate, numPurchases, orderSize, recency, stringsAsFactors = F) #, orderSize, recency
}))
close(pb)



# Cluster -----------------------------------------------------------------

library(userStream)

## either cluster all observations at once
userstream = userStream$new(r=1, lambda=0.01, tgap=30, k=5, weighted=F)
userstream$cluster(features[,c("orderSize", "numPurchases")], features[,"Customer"], features[,"Date"])


## or cluster in batches to calculate quality measures and plot (with a delay of 1 second between batches for better visualization)
userstream = userStream$new(r=1, lambda=0.01, tgap=30, k=5, weighted=F)
result = userstream$batch_iterate(features[,c("orderSize", "numPurchases")], features[,"Customer"], features[,"Date"], horizon=1000, plot=T, evaluate=T, prequential=F, sleep=1)
```
