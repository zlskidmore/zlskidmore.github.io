################################
########## SETUP ###############

# Install packages
install.packages(c('tidyverse', 'data.table', 'nycflights13'))

# load libraries
library(dplyr) # part of tidyverse
library(tidyr) # part of tidyverse
library(data.table)
library(nycflights13) # data to manipulate

#############################
##### Format Conversions ####

flightsDF <- as.data.frame(flights) # base R
flightsTB <- as_tibble(flights) # dplyr
flightsDT <- as.data.table(flights) # data.table





################################################################################
########################## Basic Data Manipulation #############################

# select rows
flightsDF[flightsDF$carrier == 'DL',]
flightsTB %>% filter(carrier == 'DL')
flightsDT[carrier == 'DL']

# select a single column
flightsDF[,'dest']
flightsTB %>% select(dest)
flightsDT[,.(dest)]

# select multiple columns by name
flightsDF[,c('origin', 'dest')]
flightsTB %>% select(origin, dest)
flightsDT[,.(origin, dest)]

# select multiple columns by position
flightsDF[,c(13,14)]
flightsTB %>% select(c(13, 14))
flightsDT[,c(13, 14)]

# add columns
flightsDF$newCol <- paste0(flightsDF$origin, ':', flightsDF$dest)
flightsTB <- flightsTB %>% mutate(newCol = paste0(origin, ':', dest))
flightsDT[,newCol := paste0(origin, ':', dest)]

# remove columns
flightsDF$newCol <- NULL
flightsTB <- flightsTB %>% select(-c(newCol))
flightsDT[,newCol := NULL]

################################################################################
#################### Basic data cleaning #######################################

# replace missing values
flightsDF[is.na(flightsDF$air_time),'air_time'] <- 0
flightsTB <- flightsTB %>% mutate(air_time = if_else(is.na(air_time), 0, air_time))
flightsDT[is.na(air_time), air_time := 0]

# rename columns
names(flightsDF)[names(flightsDF) == 'dest'] <- 'destination'
flightsTB <- flightsTB %>% rename(destination = dest)
setnames(flightsDT, c('dest'), c('destination'))

# sort by longest distance flown
flightsDF <- flightsDF[order(-flightsDF$distance),]
flightsTB <- flightsTB %>% arrange(-distance)
flightsDT <- flightsDT[order(-distance)]

# re-arrange columns
flightsDF <- flightsDF[,c('carrier', names(flightsDF)[names(flightsDF) != 'carrier'])]
flightsTB %>% relocate(carrier, .before=year)
setcolorder(flightsDT, 'carrier', before='year')

################################################################################
################### group by functions #########################################

# determine the number of flights for each airport
aggregate(flightsDF$origin, by=list(flightsDF$origin), FUN=length) 
flightsTB %>% group_by(origin) %>% summarise(count=n())
flightsDT[,.N, by=.(origin)]

# determine the most frequent route flown
aggregate(flightsDF[,c('origin', 'destination')], by=list(flightsDF$origin, flightsDF$destination), FUN=length) |> (\(df) df[order(-df$origin),])() |> (\(df) df[1,])()
flightsTB %>% group_by(origin, destination) %>% summarise(count=n()) %>% ungroup() %>% arrange(-count) %>% slice_head(n=1)
flightsDT[,.N,by=.(origin, destination)][order(-N)][1]

# determine the average air-time for the most frequent route for each carrier
aggregate(air_time~carrier,
          data = flightsDF[flightsDF$origin == 'JFK' & flightsDF$destination == 'LAX', ],
          FUN = function(x) c(mean = mean(x, na.rm = TRUE)))
flightsTB %>% filter(origin == 'JFK' & destination == 'LAX') %>% group_by(carrier) %>% summarise(avg=mean(air_time))
flightsDT[origin == 'JFK' & destination == 'LAX', mean(air_time), by=.(carrier)]

# determine both the average and median air-time for the most frequent route for each carrier
aggregate(air_time~carrier,
          data = flightsDF[flightsDF$origin == 'JFK' & flightsDF$destination == 'LAX', ],
          FUN = function(x) c(mean = mean(x, na.rm = TRUE), median = median(x, na.rm = TRUE)))
flightsTB %>% filter(origin == 'JFK' & destination == 'LAX') %>% group_by(carrier) %>% summarise(avg=mean(air_time), med=median(air_time))
flightsDT[origin == 'JFK' & destination == 'LAX', .(avg=mean(air_time), med=median(air_time)), by=.(carrier)]

################################################################################
##################### reshapeing ###############################################

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

# determine the most accurate airline based on departure delay for each month
departureDeviationDT <- flightsDT[,.(avgDelay=mean(dep_delay, na.rm=TRUE)),by=.(month, carrier)]
departureDeviationDT_wide <- dcast(departureDeviationDT, carrier~month, value.var='avgDelay')
departureDeviationDT_long <- melt(departureDeviationDT_wide, id.vars=c('carrier'), measure.vars=c(2:5), variable.name='month', value.name='avgDelay')


