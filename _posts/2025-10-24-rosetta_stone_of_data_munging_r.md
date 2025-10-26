---
layout: post
title:  "Rosetta Stone of Data Munging in R"
date:   2025-10-24 10:13:50 -0500
description: "A translation guide to data manipulation using R packages"
tags: [R]
---

## Introduction
In my journey through R, I initially leaned heavily on `Base-R`, not by choice, but by circumstance. At the time, more advanced data manipulation packages like [dplyr](https://dplyr.tidyverse.org) from the [tidyverse](https://tidyverse.org) and [data.table](https://r-datatable.com) were either nonexistent or still in their infancy and not yet widely adopted.

This post isn’t meant to spark (or settle) a debate about which syntax is best. Each paradigm, `Base-R`, `tidyverse`, and `data.table` has its own strengths and weaknesses. While I’ll admit a personal preference and deeper knowledge for data.table, I will still reach for the others when a task calls for it.

That said, I do think there are diminishing returns to mastering all three fluently. This guide isn’t about becoming fluent in every flavor of R, it’s about building a practical bridge between them. Whether you’re transitioning from one to another because you need to read a colleagues code or are just trying to understand an answer in stack-overflow, this post aims to serve as a Rosetta Stone of sorts between these dialects.

Where possible I will try and elaborate on how the syntax is working and what it's doing, however I highly recommend checking out the many detailed guides and tutorials specific to each paradigm of data-munging, I've linked a few at the end I think are particularly interesting or easy to follow.

## Prepare data
To follow along with the examples presented here you'll need to install and load the relevant packages. `Base-R` obviously is already available, `tidyverse` includes both `dplyr` and `tidyr` for data manipulation, and we will need the `data.table` library as well. We will be using flight data from the [nycflights13](https://nycflights13.tidyverse.org) package in order to have a standard dataset to manipulate. please refer to the the "Session Info" section of this post for specific package versions.

I also want to note that I would not normally recommend loading both `data.table` and `dplyr` libraries in a single R session. Many function names conflict between the two packages. As such you might see me refer explicitly to the function in one package or another with the `::` operator. For example both `dplyr` and `data.table` have a function `first()`, in order to explicitly use this function from one package or another we can do `data.table::first()` or `dplyr::first()`. This would not be necessary if loading only one library.
```R
# Install packages
install.packages(c('tidyverse', 'data.table', 'nycflights13'))

# load libraries
library(dplyr) # part of tidyverse
library(tidyr) # part of tidyverse
library(data.table)
library(nycflights13) # data to manipulate
```
### Data Structure
Before diving in, I think it is usefull to get a basic idea of what the data we will be working with looks like. We will primarily be workign with the flights data object, it should be available if you've run the `library(nycflights13)` function call above. This object is currently in a tibble, which is a native format of the `tidyverse`, we'll get to that in the next section, but for now it is worth knowing that all of these basic functions we'll be using, `nrow()`, `ncol()`, `summary()`, etc. are basic, almost universal methods and will function regardless of the data object type, at least for what we are discussing today.

<hr>
To begin we will look at the following:
- `nrow()` will give is the number of rows
- `ncol()` will give us the number of columns
- `dim()` will give us both the rows and columns in a two element vector
```r
#######################
### Basic structure ###
nrow(flights) # get number of rows
ncol(flights) # gen number of columns
dim(flights) # get dimensions of object
```
<div class="clearfix"></div>

<figure class="glow-figure-left">
  <img src="{{ site.baseurl }}/assets/posts/Rosetta_Stone_of_Data_Munging_in_R/basic_structure_1.png"
       style="max-width: 150px;"
       alt="output of structure queries" class="glow-frame-left">
  <figcaption>Output of Basic Structure Commands</figcaption>
</figure>

<hr>

- `str()` will output the basic structure of the object which will tell us
	- the rows and columns
	- columns names
	- column class, integer, character, numeric, etc.
	- the first few values in the colum

- `summary()` will output basic summary statistics for each column
```r
#######################
# Detailed Structure ##
str(flights) # get structure of object
summary(flights) # get summary statistics for object
```
<div class="clearfix"></div>

<figure class="glow-figure-left">
  <img src="{{ site.baseurl }}/assets/posts/Rosetta_Stone_of_Data_Munging_in_R/basic_structure_2.png"
       style="max-width: 1400px;"
       alt="output of detailed queries - 1" class="glow-frame-left">
  <figcaption>Output of str(flights)</figcaption>
</figure>

<div class="clearfix"></div>

<figure class="glow-figure-left">
  <img src="{{ site.baseurl }}/assets/posts/Rosetta_Stone_of_Data_Munging_in_R/basic_structure_3.png"
       style="max-width: 1400px;"
       alt="output of detailed queries - 2" class="glow-frame-left">
  <figcaption>Output of summary(flights)</figcaption>
</figure>

<hr>

- `head()` will output the first few rows of the object
- `tail()` will output the last few rows of the object
- `View()` will pull up the object in a interactive view-pane, assuming you're in rstudio

```r
#######################
#### View Object ######
head(flights) # get first few rows
tail(flights) # get last few rows
view(flights) # view object manually 
```
<div class="clearfix"></div>

<figure class="glow-figure-left">
  <img src="{{ site.baseurl }}/assets/posts/Rosetta_Stone_of_Data_Munging_in_R/basic_structure_4.png"
       style="max-width: 1400px;"
       alt="output of view queries - 1" class="glow-frame-left">
  <figcaption>Output of head(flights)</figcaption>
</figure>

<hr>

As we can see `flights` holds basic flight information for the year 2013 for new yourk airports. The data itself should be fairly self-explanatory but don't worry if you're not sure what each column is, while usefull just to more intuitively understand why we are doing something for the upcoming sections, it's not strictly necessary.
### Formats

```r
#######################
##### Conversions #####
flightsDF <- as.data.frame(flights) # base R
flightsTB <- as_tibble(flights) # dplyr
flightsDT <- as.data.table(flights) # data.table
```
### Basic Data Manipulation

```r
#######################
###### Selection ######

# select all rows for DL (Delta)
flightsDF[flightsDF$carrier == 'DL',]
flightsTB %>% filter(carrier == 'DL')
flightsDT[carrier == 'DL']
```
```r
# select a single column called dest
flightsDF[,'dest']
flightsTB %>% select(dest)
flightsDT[,.(dest)]
```
```r
# select multiple columns by name
flightsDF[,c('origin', 'dest')]
flightsTB %>% select(origin, dest)
flightsDT[,.(origin, dest)]
```
```r
# select multiple columns by column position (arrival and dest)
flightsDF[,c(13,14)]
flightsTB %>% select(c(13, 14))
flightsDT[,c(13, 14)]
```
```r
#############################
##### Creation/Deletion #####

# add a new column called "newCol"
flightsDF$newCol <- paste0(flightsDF$origin, ':', flightsDF$dest)
flightsTB <- flightsTB %>% mutate(newCol = paste0(origin, ':', dest))
flightsDT[,newCol := paste0(origin, ':', dest)]
```
```r
# remove the created column called "newCall"
flightsDF$newCol <- NULL
flightsTB <- flightsTB %>% select(-c(newCol))
flightsDT[,newCol := NULL]
```

### Basic Data Cleaning

```r
# replace missing values
flightsDF[is.na(flightsDF$air_time),'air_time'] <- 0
flightsTB <- flightsTB %>% mutate(air_time = if_else(is.na(air_time), 0, air_time))
flightsDT[is.na(air_time), air_time := 0]
```
```r
# rename columns
names(flightsDF)[names(flightsDF) == 'dest'] <- 'destination'
flightsTB <- flightsTB %>% rename(destination = dest)
setnames(flightsDT, c('dest'), c('destination'))
```
```r
# sort by longest distance flown
flightsDF <- flightsDF[order(-flightsDF$distance),]
flightsTB <- flightsTB %>% arrange(-distance)
flightsDT <- flightsDT[order(-distance)]
```
```r
# re-arrange columns
flightsDF <- flightsDF[,c('carrier', names(flightsDF)[names(flightsDF) != 'carrier'])]
flightsTB %>% relocate(carrier, .before=year)
setcolorder(flightsDT, 'carrier', before='year')
```

### Grouping and Aggregating

```r
# determine the number of flights for each airport
aggregate(flightsDF$origin, by=list(flightsDF$origin), FUN=length) 
flightsTB %>% group_by(origin) %>% summarise(count=n())
flightsDT[,.N, by=.(origin)]
```
```r
# determine the most frequent route flown
aggregate(flightsDF[,c('origin', 'destination')], by=list(flightsDF$origin, flightsDF$destination), FUN=length) |> (\(df) df[order(-df$origin),])() |> (\(df) df[1,])()
flightsTB %>% group_by(origin, destination) %>% summarise(count=n()) %>% ungroup() %>% arrange(-count) %>% slice_head(n=1)
flightsDT[,.N,by=.(origin, destination)][order(-N)][1]
```
```r
# determine the average air-time for the most frequent route for each carrier
aggregate(air_time~carrier,
          data = flightsDF[flightsDF$origin == 'JFK' & flightsDF$destination == 'LAX', ],
          FUN = function(x) c(mean = mean(x, na.rm = TRUE)))
flightsTB %>% filter(origin == 'JFK' & destination == 'LAX') %>% group_by(carrier) %>% summarise(avg=mean(air_time))
flightsDT[origin == 'JFK' & destination == 'LAX', mean(air_time), by=.(carrier)]
```
```r
# determine both the average and median air-time for the most frequent route for each carrier
aggregate(air_time~carrier,
          data = flightsDF[flightsDF$origin == 'JFK' & flightsDF$destination == 'LAX', ],
          FUN = function(x) c(mean = mean(x, na.rm = TRUE), median = median(x, na.rm = TRUE)))
flightsTB %>% filter(origin == 'JFK' & destination == 'LAX') %>% group_by(carrier) %>% summarise(avg=mean(air_time), med=median(air_time))
flightsDT[origin == 'JFK' & destination == 'LAX', .(avg=mean(air_time), med=median(air_time)), by=.(carrier)]
```

### Reshaping

```r
# base R
depatureDeviationDF <- aggregate(dep_delay~carrier + month,
                                 data=flightsDF,
                                 FUN=function(x) c(avgDelay = mean(x, na.rm = TRUE)))

depatureDeviationDF_wide <- reshape(depatureDeviationDF,
                                    timevar = "month",
                                    idvar = "carrier",
                                    direction = "wide")
colnames(depatureDeviationDF_wide) <- gsub('dep_delay.', '', colnames(depatureDeviationDF_wide))

depatureDeviationDF_long <- reshape(depatureDeviationDF_wide, 
                                    direction = "long", 
                                    varying = 2:5, 
                                    v.names = "avgDelay", 
                                    timevar = "month", 
                                    times = names(depatureDeviationDF_wide)[2:5], 
                                    idvar = "carrier")
depatureDeviationDF_long <- depatureDeviationDF_long[,c('carrier', 'month', 'avgDelay')]

```
```r
# dplyr + tidyr
depatureDeviationTB <- flightsTB %>% group_by(month, carrier) %>% summarise(avgDelay=mean(dep_delay, na.rm=TRUE))
departureDeviationTB_wide <- depatureDeviationTB %>% 
  select(carrier, month, avgDelay) %>%
  pivot_wider(names_from = month, values_from = avgDelay)

departureDeviationTB_long <- departureDeviationTB_wide %>%
  pivot_longer(
    cols = 2:5,
    names_to = "month",
    values_to = "avgDelay"
  ) %>%
  select(carrier, month, avgDelay)
```
```r
# determine the most accurate airline based on departure delay for each month
departureDeviationDT <- flightsDT[,.(avgDelay=mean(dep_delay, na.rm=TRUE)),by=.(month, carrier)]
departureDeviationDT_wide <- dcast(departureDeviationDT, carrier~month, value.var='avgDelay')
departureDeviationDT_long <- melt(departureDeviationDT_wide, id.vars=c('carrier'), measure.vars=c(2:5), variable.name='month', value.name='avgDelay')
```


## Session Info
This guide was develped on R version 4.4.1, please find specific package versions

<div class="clearfix"></div>

<figure class="glow-figure">
  <img src="{{ site.baseurl }}/assets/posts/Rosetta_Stone_of_Data_Munging_in_R/sessioninfo.png"
       style="max-width: 1400px;"
       alt="Gel bead-in-EMulsion" class="glow-frame">
  <figcaption>Output of sessionInfo()</figcaption>
</figure>

