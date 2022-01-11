

#This code cleans and prepares habitat data by removing 
#A. sites without data
#B. undersampled fragments
#C. twice sampled sites 

#####A. Removing sites without data #####
# A.1.  Reading in the data.
env.cov = read.csv(".\\Habitat_final_2017_MSE_21Dec2019.csv", na.strings=c(""," ","NA"))

# A.2.  Creating a column for average canopy openness, removing extras, renaming, removing data before 2013
env.cov = env.cov %>% mutate(openness = (DEN.1.N + DEN.2.S + DEN.3.E + DEN.4.W)/4)

colnames(env.cov)

env.cov = env.cov %>%
  dplyr::select(Locacion, Year, ELEV, Total..10.DAP, Total..50.DAP, X..Cecropias, Dosel_metros, openness) %>%
  dplyr::rename(site = Locacion, year = Year, elevation = ELEV, dbh10 = Total..10.DAP, dbh50 = Total..50.DAP, cecropias = X..Cecropias, canopy.height = Dosel_metros) %>%
  filter(year >= 2013)

# A.3.  Averaging it all together to create one value per covariate per site.
colnames(env.cov)

env.cov = env.cov %>% 
  mutate_at(vars(elevation:openness), as.numeric) %>%
  group_by(site, year) %>% 
  dplyr::summarize_at(vars(elevation:openness), mean, na.rm = TRUE) %>%
  mutate_at(vars(elevation:openness), round, 2)

# A.4.  Checking to make sure we have covariates for all sites.
setdiff(union(fcat.hum$site, env.cov$site),
        intersect(fcat.hum$site, env.cov$site))
# Not in both datasets: "100" "101"

# A.5.  Removing the extras:
fcat.hum = fcat.hum %>%
  filter(site %!in% setdiff(union(fcat.hum$site, env.cov$site),
                            intersect(fcat.hum$site, env.cov$site))) %>%
  droplevels()

env.cov = env.cov %>%
  filter(site %!in% setdiff(union(fcat.hum$site, env.cov$site),
                            intersect(fcat.hum$site, env.cov$site))) %>%
  droplevels()

# A.6.  Adding in a site.year column
env.cov$site.year = paste(env.cov$site, env.cov$year, sep = "_")

# A.7.  Checking out differences in site.year since env.cov has 8 more sites than fcat.hum

setdiff(union(fcat.hum$site.year, env.cov$site.year),
        intersect(fcat.hum$site.year, env.cov$site.year)) # 8 site.year fields not in both datasets

setdiff(env.cov$site.year, fcat.hum$site.year) # Looks like all 8 are in env.cov. Let's remove them

env.cov = env.cov %>%
  filter(site.year %in% fcat.hum$site.year)
```

##### B. Removing undersampled fragments ######

# B.1.  Reading in the effort data
effort = read.csv("./Session_summary_2004-17_22Apr2019.csv", na.strings=c(""," ","NA"))

# B.2.  Thinning it out, renaming columns, removing data pre-2013
colnames(effort)
effort = dplyr::select(effort, Fragmento, Locacion, UTM1, UTM2, Year, Fecha, Periodo.Muestreo, Net.Hours, NOTAS) %>%
  dplyr::rename(forest = Fragmento, site = Locacion, year = Year, date = Fecha, samp.per = Periodo.Muestreo, net.hours = Net.Hours) %>%
  filter(year >= 2013)

# B.3.  Creating the same identification columns in the effort sheet
effort$site.year = paste(effort$site, effort$year, sep = "_")
effort$s.y.sp = paste(effort$site, effort$year, effort$samp.per, sep = "_")
effort$site.date = paste(effort$site, effort$date, sep = "_")

# B.4.  Creating a season column

effort$season = ifelse(grepl("Jan", effort$date) | grepl("Feb", effort$date) | grepl("Mar", effort$date) | grepl("Apr", effort$date) | grepl("May", effort$date) | grepl("Jun", effort$date) | grepl("Jul", effort$date), "wet", "dry")

# B.5.  Let's first filter out effort that doesn't match hummers's site.year column.

effort = effort %>% filter(site.year %in% fcat.hum$site.year)

# B.6.  Getting rid of capture data that doesn't match up with our effort data.

# B.6a. Firsts let's see how many banding days we don't have effort for.
anti_join(fcat.hum, effort, by = "site.date") # All clear. 

# B.6b. Next, Let's look to make sure we've got all the same sampling events in the effort and capture dataframes. This is critical because individual days AREN'T our sampling units, it's a combination of days in a sampling period.
fcat.hum$s.y.sp = as.factor(fcat.hum$s.y.sp)
effort$s.y.sp = as.factor(effort$s.y.sp)
setdiff(union(fcat.hum$s.y.sp, effort$s.y.sp),
        intersect(fcat.hum$s.y.sp, effort$s.y.sp))
# Well looky there. They don't match up because apparently sampling periods were entered differently in the two .csv files. That's unfortunate. For example, March 2015 is Occasion 1 in the effort sheet and Occasion 2 in the captures sheet. Super duper. Seems to be that "Occasion" or "Sampling period" was defined differently. In the effort .csv, 1 is probably the first sampling period of every year. In the capture .csv, 1 is the first time the site was ever banded at and numbers continue building on top of that. 

# B.6d. To fix this problem, we'll have to replace sampling period in one .csv file with sampling period from the other. We can do that by matching site.date identifiers between files, then replacing the associated sampling periods in one file with the matching sampling periods from the other.
setdiff(union(fcat.hum$site.date, effort$site.date),
        intersect(fcat.hum$site.date, effort$site.date)) # Site dates don't all match up.

setdiff(fcat.hum$site.date, effort$site.date) # All extra site dates are in effort. Let's chuck them before proceeding.

effort = effort %>% filter(site.date %in% fcat.hum$site.date)

fcat.hum$samp.per = as.factor(effort$samp.per[match(fcat.hum$site.date,                                                   effort$site.date)])

# B.6e. Now to recalculate the site_year_sampling.period column
fcat.hum$s.y.sp = as.factor(paste(fcat.hum$site, fcat.hum$year, fcat.hum$samp.per, sep = "_"))

# B.7.  Removing sites we don't have env.cov or captures for. (Should be same results regardless of using hummers or env.cov since we already made sure both contain the same data).
setdiff(union(fcat.hum$site, effort$site), intersect(fcat.hum$site, effort$site)) # I guess this was resolved in the previous step.

# B.8.  Okay, let's calculate sampling effort by s.y.sp
effort = effort %>%
  mutate(forest = factor(forest, levels = c("Fragment", "BBS"))) %>%
  mutate(forest = replace_na(forest, "BBS")) %>%
  group_by(forest, site, year, season, s.y.sp) %>%
  dplyr::summarize(nh = sum(net.hours), UTM1 = mean(UTM1, na.rm = T), UTM2 = mean(UTM2, na.rm = T))

# We're all over the place here. 148 hours is about 3 days with 10 nets and 5 hours or 3 days with 8 nets and 6 hours. We have a of days with far fewer net hours (as low as 12) or many more (as high as 185.3)

# B.9.  Reordering hummers data and dropping some columns we probably won't use
fcat.hum = fcat.hum %>% dplyr::select(forest, site, year, site.year, s.y.sp, season, species, family, guild)
```


##### C. Removing twice sampled sites ####

#In examining fcat.hum$site.year, it's clear some sites were sampled in multiple years. 
#Cannot just add those captures together because the doubling in effort for a few sites will screw with average diversity estimates and comparisons among sites, 
#so I'll have to remove extra sampling events for each site that was hit multiple times. 
#Rather than removing extra sampling events by sample coverage, we're going to do it by number of net hours.

# C.1.  Ideally, we'd like all of our sites to have 3 days worth of sampling, but if we limit ourselves to 3 days, we'll drop 15 fragments and 4 Bilsa sites. Let's create a dataframe with sites that had at least 2 days worth of sampling (i.e., >90 net hours)

effort = effort %>%
  group_by(forest, site, year) %>%
  dplyr::summarize(nh = sum(nh), UTM1 = mean(UTM1), UTM2 = mean(UTM2)) %>%
  ungroup(forest, site, year) %>%
  mutate(site.year = paste(site, year, sep = "_")) %>%
  mutate_at(vars(site.year), factor) %>%
  dplyr::filter(nh >= 90) %>%
  group_by(forest, site) %>%
  mutate(rank = rank(desc(nh))) %>% # This ranks sites by nh, with the highest nh getting a 1 and lowest a 2.
  ungroup() %>%
  filter(rank == 1) %>% # Removes sampling events with lower nh.
  dplyr::select(-rank) %>%
  droplevels()

# C.2.  Using these same data to create a dataframe with effort summed up by forest. 
effort.for = effort %>% 
  group_by(forest) %>%
  dplyr::summarize(nh = sum(nh)) %>%
  mutate_at(vars(forest), factor)

# C.3.  Removing sites from the capture database and env.cov database
fcat.hum = fcat.hum %>% filter(site.year %in% effort$site.year)

env.cov = env.cov %>% filter(site.year %in% effort$site.year)

# C.4.  Filtering out everything but hummers
fcat.hum = fcat.hum %>%
  filter(family == "Trochilidae")

dim(fcat.hum)
```