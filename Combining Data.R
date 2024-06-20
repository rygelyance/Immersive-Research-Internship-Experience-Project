library("dplyr")
library("data.table")
path<-"C:/Users/rygel/Documents/team-twin-cities/TC_Meteorology_Data/"
days<-dir(path) #makes a vector of folder names

setwd("C:/Users/rygel/Documents/team-twin-cities/TC_Meteorology_Data/")

combined_files <- bind_rows(lapply(days, fread))

setwd("C:/Users/rygel/Documents/team-twin-cities/")

write.csv(combined_files, "Full_Meteorology_Data.csv")

stations <- read.csv("Station Names and IDs.csv")

path<-"C:/Users/rygel/Documents/team-twin-cities/Twin_Cities_PM25/"
months<-dir(path) #makes a vector of folder names

setwd("C:/Users/rygel/Documents/team-twin-cities/Twin_Cities_PM25/")


combined_files <- bind_rows(lapply(months, fread))

write.csv(combined_files, "Full_PM25_Data.csv", row.names = F)

setwd("C:/Users/rygel/Documents/team-twin-cities/")

PM25 <- read.csv("Full_PM25_Data.csv")
Holidays <- read.csv("major_holidays_2000_2025.csv")
Holidays = subset(Holidays, select = -year)
Weather <- read.csv("Full_Meteorology_Data.csv")
Weather = subset(Weather, select = -X)
Weather = subset(Weather, select = -V1)
stations = read.csv("Station Names and IDs.csv")

weather_fixed_date <- Weather %>%
  mutate(date = stringr::str_extract(date, "[0-9]{4}[0-9]{2}[0-9]{2}")) %>%
  mutate(date = paste(substr(date, 1, 4), "-", substr(date, 5, 6), "-", substr(date, 7, nchar(date)), sep = ""))

pm25_fixed_date <- PM25 %>%
  mutate(date = stringr::str_extract(date, "[0-9]{4}[0-9]{2}[0-9]{2}")) %>%
  mutate(date = paste(substr(date, 1, 4), "-", substr(date, 5, 6), "-", substr(date, 7, nchar(date)), sep = ""))

pm25_stations = merge(pm25_fixed_date, stations, by = "station_num")

add_weather = merge(pm25_stations, weather_fixed_date, by = "date", all.x = TRUE, all.y = FALSE)

add_holidays = merge(add_weather, Holidays, by = "date", all.x = TRUE, all.y = FALSE)

full_fixed_cols <- add_holidays %>%
  mutate(holiday = ifelse(is.na(holiday), FALSE, TRUE)) %>%
  mutate(day_of_week = weekdays(as.Date(date))) %>%
  mutate(month = months(as.Date(date)))

date_stations_sorted <- full_fixed_cols[order(full_fixed_cols$date, full_fixed_cols$station_num),]

write.csv(date_stations_sorted, "Mega_Dataframe.csv", row.names = F)
