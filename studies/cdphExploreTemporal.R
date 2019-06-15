
library(R.utils)
sourceDirectory('Documents/research/dataInt/R/')

src <- "/Users/brianconroy/Documents/research/cdph/data/"
rodents <- read.csv(paste(src, "CDPH_scurid_updated_full.csv", sep=""), header=T, sep=",")


# view yearly counts and rates
table(rodents$Year)
barplot(table(rodents$Year))

rodents_pos <- rodents[rodents$Res == 'POS',]
rodents_neg <- rodents[rodents$Res == 'NEG',]

years <- sort(unique(rodents$Year))
year_pos <- c()
year_neg <- c()
for (y in years){
  year_pos <- c(year_pos, nrow(rodents_pos[rodents_pos$Year == y,]))
  year_neg <- c(year_neg, nrow(rodents_neg[rodents_neg$Year == y,]))
}

names(year_pos) <- years
names(year_neg) <- years
barplot(year_pos)
barplot(year_neg)

year_rates <- year_pos/(year_pos + year_neg)
plot(year_rates, type='l')


# view binned counts and rates
radius <- 2
start <- years[1]
end <- tail(years, 1)
centers <- seq(start, end, by=(2*radius + 1))

bin_pos <- c()
bin_neg <- c()
for (cen in centers){
  cen_yrs <- seq(cen-radius, cen+radius)
  cen_pos <- 0
  cen_neg <- 0
  for (y in cen_yrs){
    cen_pos <- cen_pos + nrow(rodents_pos[rodents_pos$Year == y,])
    cen_neg <- cen_neg + nrow(rodents_neg[rodents_neg$Year == y,])
  }
  bin_pos <- c(bin_pos, cen_pos)
  bin_neg <- c(bin_neg, cen_neg)
}

names(bin_pos) <- centers
names(bin_neg) <- centers
barplot(bin_pos)
barplot(bin_neg)

bin_rates <- bin_pos/(bin_pos + bin_neg)
plot(bin_rates, type='l')
bin_rates

