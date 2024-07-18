library("tidyverse")
library("knitr")

df<-read.csv("Mega_Dataframe.csv") 

df2<-df %>%
  mutate(date=as.Date(date, format='%Y-%m-%d'))

#period of analysis
startdate<-as.Date("2000-06-01", format='%Y-%m-%d')
enddate<-as.Date("2008-06-01", format='%Y-%m-%d')

#opening data of metro
opendate<-as.Date("2004-06-14", format='%Y-%m-%d')

#date when groundbreak starts
conststart <- as.Date("2001-01-17", format = "%Y-%m-%d")

#Heavy-Duty Engine and Vehicle Standards
HD_Engine <- as.Date("2007-01-01", format = "%Y-%m-%d")
#Nonroad Diesel Rule
NR_Diesel <- as.Date("2004-06-29", format = "%Y-%m-%d")
#Minnesota Mercury Reduction Act
MC_Reduction <- as.Date("2006-05-11", format = "%Y-%m-%d")
#Next Gen Energy Act
NG_Energy <- as.Date("2007-05-25", format = "%Y-%m-%d")
#2030 Regional Development Framework
RDF <- as.Date("2004-01-01", format = "%Y-%m-%d")
#Clean Air Minnesota Initiative
CA_Init <- as.Date("2003-01-01", format = "%Y-%m-%d")
#Minnesota Renewable Energy Standard
RE_Standard <- as.Date("2001-01-01", format = "%Y-%m-%d")
#Vehicle Emissions Inspection and Maintenance Programs
VE_Programs <- as.Date("2001-04-05", format = "%Y-%m-%d")

df3 <- df2 %>%
  filter(date>=startdate & date <= enddate) %>%
  mutate(MetroOpen = ifelse(date >= opendate, 1, 0)) %>%
  mutate(dow = wday(date)) %>%
  mutate(construction = ifelse(date > conststart & date < opendate, 1, 0)) %>%
  group_by(station_num) %>%
  arrange(station_num, date) %>%
  mutate(t = as.numeric(date-startdate)) %>%
  mutate(t2 = t^2, t3 = t^3, t4 = t^4) %>%
  mutate(l_tair = lag(Tair_f_tavg)) %>%
  mutate(l_tair_2 = l_tair^2) %>%
  mutate(l_tair_3 = l_tair^3) %>%
  mutate(l_tair_4 = l_tair^4) %>%
  mutate(l_qair = lag(Qair_f_tavg)) %>%
  mutate(l_qair_2 = l_qair^2) %>%
  mutate(l_qair_3 = l_qair^3) %>%
  mutate(l_qair_4 = l_qair^4) %>%
  mutate(l_wind = lag(Wind_f_tavg)) %>%
  mutate(l_wind_2 = l_wind^2) %>%
  mutate(l_wind_3 = l_wind^3) %>%
  mutate(l_wind_4 = l_wind^4) %>%
  mutate(CA_Init = ifelse(date >= CA_Init, 1, 0)) %>%
  mutate(HD_Engine = ifelse(date >= HD_Engine, 1, 0)) %>%
  mutate(MC_Reduction = ifelse(date >= MC_Reduction, 1, 0)) %>%
  mutate(NG_Energy = ifelse(date >= NG_Energy, 1, 0)) %>%
  mutate(NR_Diesel = ifelse(date >= NR_Diesel, 1, 0)) %>%
  mutate(RDF = ifelse(date >= RDF, 1, 0)) %>%
  mutate(RE_Standard = ifelse(date >= RE_Standard, 1, 0)) %>%
  mutate(VE_Programs = ifelse(date >= VE_Programs, 1, 0))

summary(m1 <- lm(log(pm25) ~ MetroOpen:as.factor(station_num) +
                   construction +
                   Tair_f_tavg +
                   l_tair +
                   l_tair_2 +
                   l_tair_3 +
                   l_tair_4 +
                   Qair_f_tavg +
                   l_qair +
                   l_qair_2 +
                   l_qair_3 +
                   l_qair_4 +
                   Wind_f_tavg +
                   l_wind +
                   l_wind_2 +
                   l_wind_3 +
                   l_wind_4 +
                   holiday +
                   t +
                   t2 +
                   t3 +
                   t4 +
                   as.factor(month) +
                   as.factor(dow)
                 , data = df3))

c <- coef(m1)
len_coef<-length(coef(m1))

#get coefficients of the station-level effect
coef<-coef(m1)[(len_coef-36): len_coef]
station_num<-c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17,18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37)
coefdf<-as.data.frame(cbind(coef, station_num))


output4<-coefdf |>
  mutate(FID=station_num-1)

output5<-output4 |>
  mutate(coef=as.numeric(coef)*100)

buff<-vect("Station_Buffers.shp")
shape<-tigris::blocks(state="MN", county="Ramsey", class="sp", year=2000)
shape2<-tigris::blocks(state="MN", county="Hennepin", class="sp", year=2000)
shapevect1<-vect(shape)
shapevect2<-vect(shape2)
shapevect<-rbind(shapevect1, shapevect2)

buff2<-merge(buff, output5, by="FID")
buffdf = as.data.frame(buff2)

library("maptiles")

osmpos <- create_provider(name = "CARTO.POSITRON",
                          url = "https://{s}.basemaps.cartocdn.com/light_all/{z}/{x}/{y}.png",
                          sub = c("a", "b", "c", "d"),
                          citation = "© OpenStreetMap contributors © CARTO ")

buff3<-buffer(buff2, width=2000)
bg <- get_tiles(ext(buff3),provider = osmpos, crop = TRUE)

plot(bg, alpha=0.05)
plot(shapevect, add=TRUE)
plot(buff2, 
     "coef",
     type="interval",
     breaks=c(-30, -25, -20, -15, -10, -5),
     col=map.pal("water"),
     cex.main=2.125,
     plg=list( # parameters for drawing legend
       title = "Change in PM2.5 \n (in Percents)",
       title.cex = 1.5, # Legend title size
       cex = 2),
     add=TRUE) #legend text size

library(ggplot2)
library(tidyverse)

age_groups <- read.csv("age_group_changes.csv")
age_groups_fixed <- pivot_longer(age_groups, cols = 2:3, names_to = "gender", values_to = "coef")
ggplot(data = age_groups_fixed, aes(x = X, y = coef, fill = gender)) +
  geom_bar(position = "dodge", stat = "identity") +
  scale_x_discrete(limits=age_groups_fixed$X) +
  labs(y= "Reduction of PM2.5 \n (in Percent)", x = "Age Group") +
  scale_fill_manual(values=c("#b84279","#3751e5"))
