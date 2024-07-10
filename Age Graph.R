ba<-blockagedf |>
  mutate(aveage= (
    (total_male_under_5_years + total_female_under_5_years) * 2.5 +
      (total_male_5_to_9_years + total_female_5_to_9_years) * 7 +
      (total_male_10_to_14_years + total_female_10_to_14_years) * 12 +
      (total_male_15_to_17_years + total_female_15_to_17_years) * 16 +
      (total_male_18_and_19_years + total_female_18_and_19_years) * 18.5 +
      (total_male_20_years + total_female_20_years) * 20 +
      (total_male_21_years + total_female_21_years) * 21 +
      (total_male_22_to_24_years + total_female_22_to_24_years) * 23 +
      (total_male_25_to_29_years + total_female_25_to_29_years) * 27 +
      (total_male_30_to_34_years + total_female_30_to_34_years) * 32 +
      (total_male_35_to_39_years + total_female_35_to_39_years) * 37 +
      (total_male_40_to_44_years + total_female_40_to_44_years) * 42 +
      (total_male_45_to_49_years + total_female_45_to_49_years) * 47 +
      (total_male_50_to_54_years + total_female_50_to_54_years) * 52 +
      (total_male_55_to_59_years + total_female_55_to_59_years) * 57 +
      (total_male_60_and_61_years + total_female_60_and_61_years) * 60.5 +
      (total_male_62_to_64_years + total_female_62_to_64_years) * 63 +
      (total_male_65_and_66_years + total_female_65_and_66_years) * 65.5 +
      (total_male_67_to_69_years + total_female_67_to_69_years) * 68 +
      (total_male_70_to_74_years + total_female_70_to_74_years) * 72 +
      (total_male_75_to_79_years + total_female_75_to_79_years) * 77 +
      (total_male_80_to_84_years + total_female_80_to_84_years) * 82 +
      (total_male_85_years_and_over + total_female_85_years_and_over) * 85
  )
  /(total_male+total_female)
         )

shape<-tigris::blocks(state="MN", county="Ramsey", class="sp", year=2000)
shape2<-tigris::blocks(state="MN", county="Hennepin", class="sp", year=2000)
shapevect1<-vect(shape)
shapevect2<-vect(shape2)
shapevect<-rbind(shapevect1, shapevect2)

shape_ba<-merge(shapevect, ba, by="BLKIDFP00")

plot(buff)
plot(shape_ba, "aveage",
     cex.main = 2.5,
     plg=list( # parameters for drawing legend
       title = "Age \n (in Years)",
       title.cex = 1.5, # Legend title size
       cex = 2), add=TRUE, alpha=0.3, border = NA)
