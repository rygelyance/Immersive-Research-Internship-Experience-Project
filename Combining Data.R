library("dplyr")
library("data.table")
path<-"C:/Users/rygel/Documents/Files/VSCode Workspace/Local R Workspace/team-twin-cities/TC_Meteorology_Data/"
days<-dir(path) #makes a vector of folder names


combined_files <- bind_rows(lapply(days, fread))

write.csv(combined_files, "Full_Meteorology_Data.csv")

stations <- read.csv("stations_with_locations.csv")

path<-"C:/Users/rygel/Documents/Files/VSCode Workspace/Local R Workspace/team-twin-cities/Twin_Cities_PM25/"
months<-dir(path) #makes a vector of folder names

setwd("C:/Users/rygel/Documents/Files/VSCode Workspace/Local R Workspace/team-twin-cities/Twin_Cities_PM25/")


combined_files <- bind_rows(lapply(months, fread))

write.csv(combined_files, "Full_PM25_Data.csv", row.names = F)

setwd("C:/Users/rygel/Documents/Files/VSCode Workspace/Local R Workspace/team-twin-cities/")

PM25 <- read.csv("Full_PM25_Data.csv")
Holidays <- read.csv("major_holidays_2000_2025.csv")
Weather <- read.csv("TC_Meteorology_Data/Full_Meteorology_Data.csv")

weather_fixed_date <- Weather %>%
  mutate(date = stringr::str_extract(date, "[0-9]{4}[0-9]{2}[0-9]{2}")) %>%
  mutate(date = paste(substr(date, 1, 4), "-", substr(date, 5, 6), "-", substr(date, 7, nchar(date)), sep = ""))

pm25_fixed_date <- PM25 %>%
  mutate(date = stringr::str_extract(date, "[0-9]{4}[0-9]{2}[0-9]{2}")) %>%
  mutate(date = paste(substr(date, 1, 4), "-", substr(date, 5, 6), "-", substr(date, 7, nchar(date)), sep = ""))

Holidays_weather = merge(weather_fixed_date, Holidays, by = "date", all.x = TRUE, all.y = FALSE)
stations = read.csv("Station Names and IDs.csv")

part_merge = merge(Holidays_weather, stations, by = "date", all.x = TRUE)