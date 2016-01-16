
library(reshape2)
library(colorspace)



# set URL
URL <- "http://www.esapubs.org/archive/ecol/E094/244/"

# get benthic algae/invert data
dat.bi <- read.csv(paste0(URL,"Benthic%20density%20raw%20data.csv"))
colnames(dat.bi) <- tolower(colnames(dat.bi))

# get benthic fish data
dat.bf <- read.csv(paste0(URL,"Benthic%20fish%20density%20raw%20data.csv"))
colnames(dat.bf) <- tolower(colnames(dat.bf))

# get midwater fish data
dat.mf <- read.csv(paste0(URL,"Midwater%20fish%20density%20raw%20data.csv"))
colnames(dat.mf) <- tolower(colnames(dat.mf))


URL <- "https://raw.githubusercontent.com/eric-ward/TVVARSS/master/species_sampled_lookup.csv"
URL <- "https://raw.github.com/eric-ward/TVVARSS/master/species_sampled_lookup.csv"

URL <- "https://raw.githubusercontent.com/eric-ward/TVVARSS/master/species_sampled_lookup.csv?token=AE0I0DvPOlzn38ULW4v9IDcQRnvg7WlOks5Woq-VwA%3D%3D"

# set URL
URL <- "https://raw.githubusercontent.com/eric-ward/TVVARSS/master/"
fileName <- "species_sampled_lookup.csv?token=AE0I0DvPOlzn38ULW4v9IDcQRnvg7WlOks5Woq-VwA%3D%3D"
# get LUT of guild names
guilds <- read.csv(paste0(URL,fileName))
colnames(guilds) <- tolower(colnames(guilds))







# spp to drop/eliminate
spp.out <- c("Piscivorous fishes - pelagic","Large piscivorous fishes")
# drop spp/guilds of no interest
dat <- dat[!(dat$guild %in% spp.out),]


# get means by region (W,N,S)
dat$region <- NA
#dat$region[dat$station==2 | dat$station==3 | dat$station==7] <- "W"
#dat$region[dat$station==1] <- "N"
#dat$region[dat$station==4 | dat$station==5 | dat$station==6] <- "S"
dat$region[dat$station==2 | dat$station==3 | dat$station==7] <- 3
dat$region[dat$station==1] <- 2
dat$region[dat$station==4 | dat$station==5 | dat$station==6] <- 1

# add small random deviate to all densities to address true 0's
#dat$density <- dat$density + runif(dim(dat)[1],0.9,1)/2500






# get mean density by guild & station
dat.m <- aggregate(density ~ guild + station + period, data=dat, "sum")			  

guild.names <- sort(unique(dat.m$guild))

# log of density
dat.m$dens <- log(dat.m$density)
			  
# number of guilds
n.g <- length(guild.names)			  

# number of sites
n.s <- length(unique(dat.m$station))			  

# transform data to wide form
dat.m2 <- dcast(dat.m, period ~ guild + station, value.var="dens")

# getting time periods with no samples (ie, NAs)
per.miss <- seq(max(dat.m2[,"period"]))[!(seq(max(dat.m2[,"period"])) %in% dat.m2[,"period"])]

# insert NAs for missing dates
dat.miss <- cbind(per.miss,matrix(NA,length(per.miss),ncol(dat.m2)-1))
colnames(dat.miss) <- colnames(dat.m2)
dat.m2 <- rbind(dat.m2,dat.miss)

# total length of ts
TT <- dim(dat.m2)[1]

dat.m2 <- dat.m2[order(dat.m2$period),]

