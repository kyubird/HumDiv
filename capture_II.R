
#### Dealing with recaps ####

# Birds can only be counted once per sampling period. 
# If they were all banded, this would be easy because we could just remove every band number
# that occurs more than once in a sampling period.
# Unfortunately they're not, and we can't apply unique IDs posthoc in the dataframe 
# because the same bird would receive different ID #s every time it was caught. 

# We can't remove ALL recaps either because some of them will have been banded in different sampling periods,
# and the first time they're recaptured in a sampling period still counts towards our diversity estimates
# But in that case, the recaptured birds should have a band number. 
# Since unbanded recaps should only happen in places where they ran out of bands, 
# and since JK is pretty sure they all had feathers clipped, 
# Fercho would have recognized them and either let them go or recorded them as recaps. 
# If all of that is true, then we can just remove all unbanded recaps. 
# Especially since there shouldn't be unbanded recaps from previous sampling periods 
# Sampling periods are spread across seasons, so birds should have molted.


#SO IN SUMMARY, I can  solve this problem in two steps. 
#First, all birds that meet the following criteria can be removed: 
#Recap = Y + Band.No. = NA. 
#Second, I can remove all instances of recurring band numbers 
#in a sampling period after the first occurrence 
#(which will either be the day it was first banded 
#or the first time a recap from a previous sampling period was recaught).
  
# 1. Checking to see what species are unbanded recaptures. 
unbanded.recaps = fcat.hum %>% 
  filter(recap == "Y" & is.na(fcat.hum$band.number)) %>%
  group_by(species) %>%
  dplyr::summarize(n = length(species))
# These are all hummingbirds (n=71) with one kingfisher. That means either Fercho wasn't clipping passerine tail feathers or he was and he never recaptured one (possible). This is worth considering since so many passerines went unbanded over the years, particularly in the fragments in 2015 when it seems they ran out of small bands.

# 2. Removing recaptured birds with no band number from the dataset (same ones as above).
fcat.hum = fcat.hum %>% filter(recap == "N" | banded == "Y")

# 3. If birds from previous sampling periods were never recaptured in later sampling periods, then we don't have to worry about accidentally removing them. Let's see what we're working with in terms of recaptures per specific site-year-sampling period combos.
reoccurring = fcat.hum %>% 
  filter(!is.na(band.number)) %>% 
  group_by(s.y.sp, band.number) %>% 
  dplyr::summarize(no_rows = length(band.number))

# PROBLEM. If you search the captures for some of these band numbers that keep reappearing, you'll find there seem to be duplicate entries! For example, DD 1258 was caught 4 times during 79_2014_1. THREE of those instances are marked as new captures on the same day! That shouldn't be possible; there should only be one! This is easy to confirm by searching the original for band numbers in the code above. So let's see how many of these are NOT recaptures.
reoccurring2 = fcat.hum %>% 
  filter(!is.na(band.number) & recap == "N") %>% 
  group_by(s.y.sp, band.number) %>% 
  dplyr::summarize(no_rows = length(band.number)) %>%
  filter(no_rows > 1)

# Fortunately, it looks like the vast majority of these are birds being entered into the database twice. Probably what happened was birds were caught twice on the same day and rather than being released or recaptured, they were marked as new since they were technically new to that day or sampling period. So we can just delete those! Unfortunately, if this is the case, it suggests all those unbanded birds, whose tails were theoretically clipped, were probably also double or triple counted as new birds if they were recaptured in the same sampling period. Sadly, there's not really anything we can do about that.
# Getting rid of those since we won't need them again.
rm(reoccurring, reoccurring2, unbanded.recaps)

# 4. Removing duplicate entries and the rest of the recaptures. This will take a couple steps. First, we'll split the dataset into banded and unbanded birds. Then we'll create unique IDs for the banded birds by band number and sampling period and call it a factor. We'll then create a new dataframe using just the levels of that factor since that will be equivalent to having one record of each banded individual per sampling period. This step wouldn't work for unbanded birds since there would be multiple individuals with the same unique ID, henche the splitting of the dataset Finally, we'll add the two dataframes (banded unique ID levels and unbanded birds) back together into a single dataframe.

unbanded = fcat.hum %>% filter(is.na(band.number)) # Unbanded birds. Remember, unbanded recaps have already been removed.
banded = fcat.hum %>% filter(!is.na(band.number)) # Banded birds

banded$ID = paste(banded$s.y.sp, banded$band.number, sep = "_") # Creating the unique ID.
banded$ID = as.factor(banded$ID)

# Ordering the captures by site, then year, then sampling period, then date to make sure they're all spatially and chronologically ordered (though doing it just by date may do it)
banded = banded[with(banded, order(site, year, samp.per, date)), ]
# This code came from here (it's suggested in the link for the next code below):https://stackoverflow.com/questions/1296646/how-to-sort-a-dataframe-by-multiple-columns

banded = banded[match(unique(banded$ID), banded$ID),] # This chunk of code extracts the first occurrences of each ID. So that should equate to the first time each individual was captured per sampling period. Thus it removes all recaptures from birds first banded during that sampling period and also those duplicate entries mentioned above. I found this code here: https://stackoverflow.com/questions/19944334/extract-rows-for-the-first-occurrence-of-a-variable-in-a-data-frame

# That last step dropped 481 data points! To check it's accuracy you, can rerun banded from the start and create a new dataframe from the last step (just give it a different name). Then, compare the two by ordering by ID. You'll see in the unfiltered banded dataframe, there are multiple instances of both recaps and duplicate data. It's all getting weeded out in the right order.

# Okay, we're ready to reattach the unbanded and the filtered banded dataframes back together.
unbanded$ID = NA# First, adding a column to the unbanded dataframe so it matches up with the banded one.

# Next, binding them together.
good = bind_rows(banded, unbanded)
# And reordering
good = good[with(good, order(site, year, samp.per, date, species)), ]

# Tidying up. Going to keep using caps.1415 for awhile, so we'll set that equal to good.1415 and then toss the latter. Also getting rid of banded.1415, unbanded.1415, and recaps since we won't be using them again.
fcat.hum = good
rm(banded, good, unbanded)

# More evidence that whole section went well: You can find bird XX860 is a recapture included after all the filtering. That's because it was first banded in 2010 but was recaptured during a different period in 2015. That's exactly what we were looking for.

# 5.  Now that we're sure we've got all the right captures and none of the wrong ones in our dataframe, we can get rid of all the band info since it will just be clutter from here on out. Also, before we correct for unequal sampling effort, we can tally up species captures by date or sampling period. Let's do all that now.
colnames(fcat.hum)
fcat.hum = fcat.hum %>% 
  dplyr::select(site, forest, season, date, year, site.year, s.y.sp, species, banded, recap, hummer, family, guild, notes)

# 6.  Adding a site_date column to match up with sampling effort below.
fcat.hum$site.date = paste(fcat.hum$site, fcat.hum$date, sep = "_")