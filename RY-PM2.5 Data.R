library("tidyverse")
library("ggmap")

geo <- read.csv("stations_with_locations.csv")

library("terra")

sample_coords <- cbind(geo$lon2, geo$lat2)
lr_stations <- vect(sample_coords)
geom(lr_stations)
crdref <- "+proj=longlat +datum=WGS84"
pts <- vect(sample_coords, crs=crdref)
pts_buffer <- buffer(pts, width = 500) #TEMPORARY, DO NOT USE UNTIL FINALIZED!


path<-"G:/Shared drives/2024 FIRE Light Rail/DATA/PM25/"
months<-dir(path) #makes a vector of folder names

for (m in 1:length(months)) {
  print(months[m])
  days<-dir(paste0(path,months[m])) #makes a vector of filenames within each folder


  days_output<-c()
  for (d in 1:length(days)) {
    print(days[d])

    #read tif file
    r<-rast(paste0(path, months[m], "/", days[1]))

    #changes the crs system
    buffer_project<-terra::project(pts_buffer,  crs(r))

    #pts_buffer is the buffer around stations
    #crops raster to contain only buffers around stations
    int<-crop(r, buffer_project,
            snap="in",
            mask=TRUE)

    #convert cropped raster into dataframe and fine average value
    cntrl_df<-terra::extract(int, buffer_project, fun="mean", na.rm=TRUE)

    #rename columns
    names(cntrl_df)<-c("station_num","pm25")

    #create a dataframe date, shape index, and pm25
    output <- as.data.frame(c("date"=days[1], cntrl_df))

    #combine output with previous loop
    days_output<-rbind(days_output, output)
  }

  write.csv(days_output, paste0("Twin_Cities_PM25/", months[m],".csv"), row.names = F)

}
