
# NOTE: 
# Beware that environmental covariates dataframe is altered throughout. 
# Therefore, all this code should be re-run for any change made to that dataframe to prevent future errors. 
# Only using data from 2013-2017 accounting for year-to-year variation in sampling. 
# e.g. Bilsa was heavily sampled in earlier years while fragments are not.

##### Cleaning and organizing capture data #####

# 1.  Reading in the capture data
fcat.hum = read.csv(file = ".\\Captures_2004_2017_Final_(MSE).csv", na.strings=c(""," ","NA"))

# 2.  Subsetting to select useful columns
colnames(fcat.hum)
fcat.hum = dplyr::select(fcat.hum, Fragmento, Locacion, Fecha, No.Anillo, Recap, Occasion, Year, EnglName.SACC, Guild_2, Family.SACC, Notas)

# 3.  Changing column names to make them easier to work with
fcat.hum = fcat.hum %>% dplyr::rename(forest = Fragmento, site = Locacion, date = Fecha, band.number = No.Anillo, recap = Recap, samp.per = Occasion, year = Year, species = EnglName.SACC, guild = Guild_2, family = Family.SACC, notes = Notas)

# 4.  Making  "site x year" and "site x year x sampling period" unique identifiers. we'll cross-reference with sampling effort dataframe to make sure nothing's missing from either.

fcat.hum$site.year = paste(fcat.hum$site,fcat.hum$year, sep = "_")

fcat.hum$s.y.sp = paste(fcat.hum$site, fcat.hum$year, fcat.hum$samp.per, sep = "_")

# 5.  Checking species levels to make sure there aren't any repeats, typos, etc. And converting them all to caps to make them a bit easier to work with.

fcat.hum$species = toupper(fcat.hum$species)
fcat.hum$species = as.factor(fcat.hum$species)

levels(fcat.hum$species)
levels(fcat.hum$species)[levels(fcat.hum$species) == "BROWN-WINGED SCHIFFORNIS"] = "NORTHERN SCHIFFORNIS"
levels(fcat.hum$species)[levels(fcat.hum$species) == "BLUE-LORED ANTBIRD"] = "ZELEDON'S ANTBIRD"
levels(fcat.hum$species)[levels(fcat.hum$species) == "HALF-COLLARED GNATWREN"] = "TAWNY-FACED GNATWREN"
levels(fcat.hum$species)[levels(fcat.hum$species) == "RED-EYED VIREO"] = "CHIVI VIREO"
levels(fcat.hum$species)[levels(fcat.hum$species) == "CHESTNUT-CAPPED BRUSH FINCH"] = "CHESTNUT-CAPPED BRUSHFINCH"
levels(fcat.hum$species)[levels(fcat.hum$species) == "GOLDEN-FACED TYRANNULET"] = "CHOCO TYRANNULET"

# 6.  Checking the guilds to make sure they're okay.
levels(fcat.hum$guild)
levels(fcat.hum$guild)[levels(fcat.hum$guild) == "Piscívoros or ictiófagos"] = "Piscivore"
levels(fcat.hum$guild)[levels(fcat.hum$guild) == "Vegetarian"] = "Herbivore"
levels(fcat.hum$guild)[levels(fcat.hum$guild) == "Insectivore "] = "Insectivore"

# 7.  Checking families for typos. Will need to do a deeper dive at some point to make sure species are properly classified into the right families following most recent SACC updates.
levels(fcat.hum$family) # Looks good!

# 8.  Cleaning up band numbers. Looks like there's a lot of spaces in random places in these band numbers. Removing them to avoid any promblems they may create later on (particularly because we'll probably use them to create unique IDs by sampling period and spaces might be aggravating).

bands = data.frame(levels(fcat.hum$band.number)) 

# 8a. we have 6,794 band numbers, but looking at them, there are some weird entries I don't understand: E.g., "YY 160 M-C / RS-RS", "(L) B. Y (R) R", "02315 WILDLIFE BOX 8 2601 610". After looking at the data, still no idea what the first or last ones mean. They're Green and Golden-winged Manakins respectively. Things following the second pattern seem to be some sort of code used for large birds like toucans and hawks. Other weird bands: "Red on top...", "NO SE ANILLA", "muerto", things with (JK) or (GG) after the band number, "A1", "A4". These are just a handful of cases, so I'm going to change them by hand in the original .csv file and move weird extra information from the banding column to the notes column. Will then have to re-read the data back in and re-run all of this code (just one time).
# - Changed NO SE ANILLA to unbanded (blank in csv)
# - Removed (GG) and (JK) from band number and put in notes
# - Removed all code after band numbers for the YY Green Manakins and put in notes
# - Muerto was the same individual/capture event entered twice. Deleted a row, changed to unbanded, and put muerto in the notes.
# - Left the other cases mentioned above alone.

# 8b. Getting rid of the spaces. Need to make sure we still have the same number of bands. There will likely be less after. For example if there are two bands, DD 123 and DD123, in the dataset, they would get combined to DD123. However, if that's the same individual, the second instance should appear as a recapture regardless of any typos in the band number. Need to cross reference with NEW captures to find any problems.
fcat.hum$band.number = gsub('\\s+', '', fcat.hum$band.number)
fcat.hum$band.number = as.factor(fcat.hum$band.number)
bands2 = data.frame(levels(fcat.hum$band.number)) # Now we're down to 6,759 (lost 35)
rm(bands, bands2)

# 9.  Adding hummingbird Y/N and banded (Y/N) columns to easily separate our focal groups from the rest. Also turning recap into Y/N. This will come in handy later!
fcat.hum$hummer = NA
fcat.hum$hummer = ifelse(fcat.hum$family == "Trochilidae", "Y", "N")
fcat.hum$hummer = as.factor(fcat.hum$hummer)

fcat.hum$banded = NA
fcat.hum$banded = ifelse(is.na(fcat.hum$band.number), "N", "Y")
fcat.hum$banded = as.factor(fcat.hum$banded)

fcat.hum$recap = as.character(fcat.hum$recap)
fcat.hum$recap[is.na(fcat.hum$recap)] = "N"
fcat.hum$recap = as.factor(fcat.hum$recap)

# 10. We turned the Fragmento column into 'forest.' Let's change those NAs to BBS
fcat.hum$forest = factor(fcat.hum$forest, levels = c("Fragment", "BBS"))
fcat.hum$forest = fcat.hum$forest %>% replace_na("BBS")

# 11. Creating a season column to check for differences between wet and dry season.
fcat.hum$season = ifelse(grepl("Jan", fcat.hum$date) | grepl("Feb", fcat.hum$date) | grepl("Mar", fcat.hum$date) | grepl("Apr", fcat.hum$date) | grepl("May", fcat.hum$date) | grepl("Jun", fcat.hum$date) | grepl("Jul", fcat.hum$date), "wet", "dry")

# 12. Filtering out everything before 2013
fcat.hum = fcat.hum %>% 
  filter(year >= 2013) %>%
  droplevels()

dim(fcat.hum)