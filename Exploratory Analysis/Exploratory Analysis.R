
library(data.table)
library(dplyr)
library(magrittr)
library(ggplot2)
library(lubridate)

data_raw = fread("Glider data/cabot_20200717_114_delayed_corrected_v4.csv")

data_renamed = data_raw %>% rename(
    temp = sea_water_temperature,
    salt = sea_water_practical_salinity,
    density = sea_water_density,
    oxygen = micromoles_of_oxygen_per_unit_mass_in_sea_water
)

# Split time into date and time
data_date_time = data_renamed %>% mutate(
    date = ymd(substr(time, 1, 10)),
    time = hms(substr(time, 12, 19))
)

test = data_renamed[1:10000,]

q = test$time[1] %>% stringr::str_sub( 1, -1)

q = test$time
# Remove the trailing " UTC" from q
w = stringr::str_sub(q, 1, -4)

e = as_datetime(w, format = )

ymd(q)
w = as_datetime(q)

q = q[1]
w = w[1]

w = parse_date_time(q, "%Y-%m-%d %H:%M:%S")

# ---------------------------------------------------------------------------- #
#                              Plot subset of data                             #
# ---------------------------------------------------------------------------- #

test = data[1:2000,]


lat = data$latitude

diffs = diff(lat)

plot(lat[1:10000])

plot(test$latitude)
plot(data$longitude)
plot(data$depth)

plot(test$temp)
plot(test$salt)
plot(test$density)
plot(test$oxygen)


check = test %>% select(-date, -time, -latitude, -longitude)
GGally::ggpairs(check)



# ----------------------- Plot depth colored by oxygen ----------------------- #

test_data_plot = test %>% mutate(ind = row_number())

# Sequential plot of depth
ggplot(test_data_plot, aes(x = ind, y = depth, color = oxygen)) + geom_point()



# ---------------------------------------------------------------------------- #
#                             Find source of jumps                             #
# ---------------------------------------------------------------------------- #

# test = data[100001:110000,]
plot(test$depth)

surface_data = test %>% select(latitude, longitude, depth) %>%
    filter(depth < 2) 

plot(-surface_data$depth)
abline(v=20)
plot(-diff(surface_data$depth))
abline(h=0)
abline(v=20)

filter(test, depth <1)

which(test$depth <1)

plot(test$latitude)




which(test$depth < 2)

windows = c(864, 884, 3874, 3921, 6991, 7047)
plot(test$latitude)
abline(v=windows)

plot(test$latitude[windows[1]:windows[2]])

