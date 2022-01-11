
# This code chooses relevant site coordinates to 
# A. prepare a map
# B. calculate surrounding forest cover 
# C. calculate fragment area
# D. examine and transform environmental covariates
# E. center and scale 
# F. create species matrices for diversity analysis

##### A. Preparing the map #####

# 1.  First creating a 'sites' dataframe to manipulate
sites = dplyr::select(effort, site, year, UTM1, UTM2) %>%
  group_by(site, year) %>%
  dplyr::summarize(UTM1 = mean(UTM1), UTM2 = mean(UTM2)) %>%
  ungroup(site, year)# If we don't do this, it may force in a site column.

sites$site.year = paste(sites$site, sites$year, sep = "_")

sites = sites %>%
  dplyr::select(site.year, year, UTM1, UTM2)

# 2.  Turning it into a spatial layer.
sp::coordinates(sites) = ~UTM1+UTM2

# 3.  Setting CRS
crs(sites) = CRS("+proj=utm +zone=17N ellps=WGS84")

# 4.  Visualizing
ggplot() + layer_spatial(sites, color = "black")

# 5.  Let's check it out alongside our fragment and BBS boundaries
frag.bound = readOGR("./Hummer_Files/FCAT_2013-2017_Diversity_Fragments_26Mar2020.shp")

bbs.bound = readOGR("./Hummer_Files/Bilsa boundary_pg.shp")

ggplot() + layer_spatial(frag.bound) + layer_spatial(bbs.bound) + layer_spatial(sites)

# Not great. It looks like a number of our sites fall outside of the fragments, which means if we were to use them to calculate surrounding forest cover, our estimates would be off.

# 6.  Let's find fragment centroids and use those instead of the sites we calculated above. We'll still use the BBS points though.

centroids = gCentroid(frag.bound, byid = TRUE, id = frag.bound@data$layer)

ggplot() + layer_spatial(frag.bound) + layer_spatial(bbs.bound) + layer_spatial(centroids) # Much better

# 7.  Replacing the 'site' UTMs with the 'centroids' UTMs

# 7a. Not sure how to swap coordinates out of spatial objects, so we'll go back to the beginning.
sites = dplyr::select(effort, site, year, UTM1, UTM2) %>%
  group_by(site, year) %>%
  dplyr::summarize(UTM1 = mean(UTM1), UTM2 = mean(UTM2)) %>%
  ungroup(site, year)

sites$site.year = paste(sites$site, sites$year, sep = "_")

sites = sites %>%
  dplyr::select(site, site.year, year, UTM1, UTM2)

centroids.df = as.data.frame(centroids@coords)
centroids.df = rownames_to_column(centroids.df, "site")

sites$site = as.character(sites$site)

# NOTE: Apparently can't use %in% here because it matches by site but then replaces by row id rather than by site.

centroids.df = centroids.df %>% arrange(site) #If they're not in the same order, the next line of code won't work. It'll match values by rowid instead.

sites = sites %>%
  mutate(UTM1 = ifelse(site %in% centroids.df$site, centroids.df$x, UTM1)) %>%
  mutate(UTM2 = ifelse(site %in% centroids.df$site, centroids.df$y, UTM2))

# 8.  Making it a spatial object again
sp::coordinates(sites) = ~UTM1+UTM2
crs(sites) = CRS("+proj=utm +zone=17N ellps=WGS84")

ggplot() + layer_spatial(frag.bound) + layer_spatial(bbs.bound) + layer_spatial(sites)


##### B. Calculating forest cover #####
# 1.  Creating a clipping layer that restricts our forest cover rasters to the study area and projecting to UTM
cover.clip.layer = readOGR("./Hummer_Files/FCAT_clipping_extent.shp")

#ggplot() + layer_spatial(data = cover.clip.layer, fill = NA, colour = "black") + theme_void()

# 2.  Reading in forest cover from the year 2000. This is our base layer from which we'll have to subtract forest cover loss through our target year and then add forest cover gain, though we only have to do the latter if the cell didn't have cover in 2000, which is probably super  unlikely.
cover.2000 = raster("./Hummer_Files/Hansen_GFC-2018-v1.6_treecover2000_10N_080W.tif")
# Not gonna visualize this one cause it takes a long time.

# 3.  Clipping the 2000 forest cover raster to our study area extent.
cover.2000 = crop(cover.2000, cover.clip.layer)

#ggplot() + layer_spatial(cover.2000)

# 4.  Reading in forest loss by year, then clipping that.
# NOTE: While the Hansen et al. dataset also includes forest gain, its values are from 0 - 255, so we can't figure out what year forest recovered. I don't think there's really enough regeneration to worry about it anyway.
cover.loss = raster("./Hummer_Files/Hansen_GFC-2018-v1.6_lossyear_10N_080W.tif")

cover.loss = crop(cover.loss, cover.clip.layer)

cover.loss #V0 = no forest loss. Values 1 - 18 indicate year of loss (2001-2018)

#ggplot() + layer_spatial(cover.loss)

# 5.  Convert to UTM coordinates
cover.2000 = projectRaster(from = cover.2000, crs = CRS("+proj=utm +zone=17N ellps=WGS84"), method = "ngb")
cover.2000

cover.loss = projectRaster(from = cover.loss, crs = CRS("+proj=utm +zone=17N ellps=WGS84"), method = "ngb")
cover.loss

# 6.  Subset to just forest cover over 95% and no forest loss between 2000 and the year banding occurred. Each cell is 1 if forest cover is above 95% threshold and forest has not been lost in that year. This code creates a stack of 12 different raster layers, one for every year that sampling occurred. 

levels(as.factor(effort$year))

cover.stack = stack()

for (year in c(13,14,15,16,17)) {
  cover.temp = cover.2000 > 95 & (cover.loss == 0 | cover.loss > year)
  cover.stack = stack(cover.stack, cover.temp)
}# Code creates a ton of warnings because our rasters are converted to UTM, but it doesn't seem to be a problem. No warnings if run without the projection.

names(cover.stack) = c("2013", "2014", "2015", "2016", "2017")

# 7.  Buffering each raster layer (forest cover by year) around every site at radii of 500m and 2km and simultaneously calculating percent forest cover within that buffer.

cover.buffers = c(500,2000)
cover.list = list()
for(i in 1:length(cover.buffers)){
  cover.list[[i]] = raster::extract(cover.stack, sites, method = "simple", buffer = cover.buffers[i], fun = function(x) sum(x)/length(x))
  names(cover.list)[[i]] = paste("Buffer", cover.buffers[i], sep = "_")
}

# 8.  Combining all that forest cover data into a workable forest cover dataframe.
cover.data = as.data.frame(do.call(cbind,lapply(names(cover.list),function(x){
  res <- cover.list[[x]]
  colnames(res) <- paste(colnames(res),x,sep="_")
  res
})))

# 9.  Adding site.year and year data to the data.frame
cover.data$site.year = sites@data$site.year
cover.data$year = sites@data$year
# Note: Irrelevant warning messages might show up.

# 10. Pulling out just the relevant data for each site/year
cover.data = cover.data %>%
  gather("column", "cover", X2013_Buffer_500:X2017_Buffer_2000) %>%
  filter(unlist(Map(function(x, y) grepl(x, y), year, column))) %>%
  mutate(column = replace(column, str_detect(column, "Buffer_500"), "cover.500")) %>%
  mutate(column = replace(column, str_detect(column, "Buffer_2000"), "cover.2000")) %>%
  spread(column, cover)

# 11. Adding cover data into the env.cov data.frame. NOTE: We don't have env.cov data for 2004-2006, so you can expect some warnings here about all the NAs that will come from the full_join of dataframes with different numbers of rows.
env.cov = full_join(env.cov, cover.data, by = "site.year")

# 12. So it looks like we have a veg data for a few site.years in which we don't have any sampling effort data. Those sites were probably removed since they didn't have enough net hours. Let's remove them, drop some incomplete columns, and reorder.
colnames(env.cov)

env.cov = env.cov %>% 
  filter(!is.na(cover.500)) %>%
  ungroup() %>%
  mutate(site = sub('_.*$','', site.year)) %>% # This one creates a site column by deleting everything after the '_' in site site.year column.
  dplyr::select(site, year.y, site.year, elevation, dbh10, dbh50, cecropias, canopy.height, openness, cover.500, cover.2000) %>%
  dplyr::rename(year = year.y)

# 13. Gosh. Upon further inspection, it looks like habitat data wasn't collected until, on average, 10-12 years AFTER banding occurred for most of the sites sampled in 2004-2006. That's crazy! Maybe there's missing data somewhere? If veg data for those sites had been collected in 2008, I'd feel better about using 2008 veg data, but 10-12 years is way too long. Let's delete those rows missing data.

env.cov = env.cov %>%
  filter(!is.na(elevation)) %>%
  droplevels()


# 15. Now let's get rid of all that GIS clutter.
rm(cover.2000, cover.clip.layer, cover.data, cover.list, cover.loss, cover.stack, cover.temp, cover.buffers, i, year, centroids.df, centroids, frag.bound, bbs.bound)


##### C. Calculating fragment area #####

# 1.  Reading in the fragment covariates that Luke Browne (?) made
cov.frag = read.csv(file = "./2014-2015_frag_covariates_28Oct2019.csv")

# 2.  Adding them to our covariates dataframe
env.cov$area = as.numeric(cov.frag$area[match(env.cov$site, cov.frag$site)])

# Okay. We're missing area for 6 fragments that weren't sampled in 2014-2015. 5 of those were sampled in 2010/2012/2016, one was sampled in 2011 & 2013. They're fragments 3, 5, 6, 7, 9, 11

# Plugging in their area, found in another file on dredge
# 3 = 48
# 6 = 4.4
# 7 = 36.3
# 9 = 8.9
# 11 = 12.5

# Also, it looks like there are two fragment 8s. One is small and beside fragment 9, the other is big and covers both those smaller ones. I think the smaller one is right, so we'll replace that value with 6.75, which I found on QGIS.

env.cov$area = ifelse(env.cov$site == "3", 48, 
                      ifelse(env.cov$site == "6", 4.4,
                             ifelse(env.cov$site == "7", 36.3,
                                    ifelse(env.cov$site == "9", 8.9,
                                           ifelse(env.cov$site == "11", 12.5, 
                                                  ifelse(env.cov$site == "8", 6.75, env.cov$area))))))

rm(cov.frag)


##### Examining and transforming env. covariates ####

# 1.  Making correlation matrices 
round(cor(env.cov[,4:12]),2)

cor.chart = chart.Correlation(env.cov[,4:12], histogram=TRUE, pch=19)

# 2.  Dropping and transforming the worst offenders

env.trans = env.cov %>%
  dplyr::select(-canopy.height, -cover.2000) %>%
  mutate(cecropias = log(cecropias +1)) %>%
  mutate(dbh50 = sqrt(dbh50)) %>%
  mutate(dbh10 = sqrt(dbh10)) %>%
  mutate(area = sqrt(area))

# 3.  Checking out that correlation chart again.
cor.chart2 = chart.Correlation(env.trans[,4:10], histogram=TRUE, pch=19)
# Area and elevation are both still problems, but we'll probably leave them in.

##### Centering and scaling habitat data #####

# There are a couple ways we can go about scaling and centering: (1) we can lump all data from all areas sampled together so everything is universally scaled and centered around a single point, or (2) we can do it independently for fragments, Bilsa, and Jama-Coaque. The former approach feels right to me, so that's what we'll go with.

env.trans.cs = data.frame(env.trans$site, env.trans$year, env.trans$site.year, scale(env.trans[,c(4:10)], center = TRUE, scale = TRUE))

colnames(env.trans.cs)

env.trans.cs = env.trans.cs %>%
  dplyr::rename(site = env.trans.site, year = env.trans.year, site.year = env.trans.site.year)

##### Creating species matrices for diversity analyses #####

# Most diversity analyses will want species as columns or rows, with spatial/temporal units going in the opposite direction and each cell containing number of individuals captured. We'll aggregate this data at multiple levels to cover our bases for different analyses.

levels(fcat.hum$forest)
levels(fcat.hum$forest) <- c(levels(fcat.hum$forest), "BBS")
fcat.hum$forest[is.na(fcat.hum$forest)] <- "BBS"

# 1.  Forest species matrix
spmat.for = fcat.hum %>%
  group_by(forest, species) %>%
  dplyr::summarize(n = n()) %>% 
  pivot_wider(names_from = forest, values_from = n) %>%
  replace(is.na(.), 0) %>%
  column_to_rownames("species")

# 2.  We'll have to calculate Hill numbers adjusted by nethour with hill_taxa(), which requires our species matrices be flipped the other way around.

spmat.sy = fcat.hum %>% # Just fragment data
  group_by(forest, site, year, species) %>%
  dplyr::summarize(n = n()) %>% 
  pivot_wider(names_from = species, values_from = n) %>%
  replace(is.na(.), 0) %>%
  mutate(f.s.y = paste(forest, site, year, sep = "_")) %>%
  column_to_rownames("f.s.y") %>%
  dplyr::select(-forest,-site,-year)
```
