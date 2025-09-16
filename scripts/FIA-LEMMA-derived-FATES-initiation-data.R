## this script uses FIA and LEMMA species composition
## data to derive initialization data for FATES
## key variables to derive: PFT, DBH, stand density scaled
## to 1ha (10000 m2) plot area
library(rFIA)
library(tidyverse)

ft2acre_2_m2ha = 0.23
inch_2_cm = 2.54
acre_2_ha = 0.405
dec = 9L
home_path = path.expand("~")
data_path = file.path("~/NASA-IDS/data")
fia_path = file.path(data_path,"CA-FIA")
#lemma_ba = data.table::fread("~/NASA-IDS/data/WRF-LiveBasalArea-Sierra_2017.csv")
lemma_ba = data.table::fread("~/NASA-IDS/data/WRF-LiveBasalArea-Sierra-CanCovMasked_2017.csv")

sp_interest = c("ponderosa pine",
                "lodgepole pine",
                "Jeffrey pine",
                "sugar pine",
                "gray or California foothill pine",
                "Coulter pine",
                "limber pine",
                "incense-cedar",
                "white fir",
                "California red fir",
                "Douglas-fir",
                "blue oak",
                "canyon live oak",
                "California black oak",
                "California live oak",
                "interior live oak")

yr_include = c(2000:2017)

#### helpler function ####
assign_na = function(x){
  x = ifelse(x>1000, NA, x)
}

replace_fun = function(idx, cn){
  idx_safe = ifelse(idx %in% seq_along(cn), idx, NA_integer_)
  cn[idx_safe]
}
rename_fun = function(x){
  x = gsub(patter="X", replacement = "CN", x)
}

nth_baCN <- function(x, k) {
  hits <- str_extract_all(x, "(?<=baCN)\\d+")
  vapply(hits, function(v) if (length(v) >= k) as.integer(v[k]) else NA_integer_, 1L)
}

fia_ca = readFIA(fia_path)
id = findEVALID(fia_ca, mostRecent = FALSE, year = yr_include)
fia_2017 = clipFIA(fia_ca, mostRecent = FALSE, evalid = id)
fia_recent = clipFIA(fia_2017, mostRecent = TRUE)

tpa_byplot = tpa(fia_recent,byPlot = TRUE,
                 landType = 'forest',
                 treeType = 'live')
tpa_bypltsp = tpa(fia_recent,bySpecies = TRUE,
                  byPlot = TRUE,
                  landType = 'forest',
                  treeType = 'live')

## only keep plots where species of interest account for > 95% of total basal area
## and total stem density

tpa_byplot = tpa_byplot                 %>% 
             dplyr::select(YEAR, pltID, 
                    PLT_CN,
                    BAA, TPA)           %>% 
             rename(tpa_plt = TPA,
                    baa_plt = BAA)
tpa_bypltsp = tpa_bypltsp               %>% 
              dplyr::select(YEAR,pltID,
                     SPCD,
                     COMMON_NAME,
                     SCIENTIFIC_NAME,
                     PLT_CN,TPA,BAA)    %>% 
              rename(tpa_sp = TPA,
                     baa_sp = BAA)
tpa_df = tpa_byplot %>% left_join(tpa_bypltsp, by=c("YEAR",
                                  "PLT_CN","pltID"))
tpa_pct = tpa_df                               %>% 
         mutate(basp_pct = baa_sp/baa_plt,
                sdsp_pct = tpa_sp/tpa_plt)    %>% 
         filter(COMMON_NAME %in% sp_interest) %>% 
         group_by(pltID,PLT_CN)               %>% 
         mutate(bapct_sum = sum(basp_pct),
                sdpct_sum = sum(sdsp_pct))    %>% 
         ungroup()
         
plt_interest = tpa_pct                       %>% 
               filter(bapct_sum >= 0.95 &
                      sdpct_sum >= 0.95)   
pltcn = unique(plt_interest$PLT_CN)
spcd = unique(tpa_bypltsp$SPCD)

tpa_plot_sf = tpa(fia_recent, 
                 byPlot = TRUE,
                 landType = 'forest',
                 treeType = 'live',
                 returnSpatial = TRUE)

tpa_plot_sf = tpa_plot_sf                    %>% 
              filter(PLT_CN %in% pltcn)      %>% 
              mutate(BAA = BAA*ft2acre_2_m2ha,
                     TPA = TPA/acre_2_ha)

hist(tpa_plot_sf$BAA)
hist(tpa_plot_sf$TPA)


ba_plot = plotFIA(tpa_plot_sf,BAA,
                  plot.title = 'Total live tree basal area (m2/ha)',
)


#p1 = plotFIA(tpa_plot_sf, BAA, grp = COMMON_NAME, x = sizeClass,
#        plot.title = 'Size-class distributions of BAA by species', 
#        x.lab = 'Size Class (inches)', text.size = .75,
#        n.max = 10, color.option = 'magma')


## tree datasets for plots of interest

tree_df = filter(fia_recent$TREE, PLT_CN %in% pltcn)
tree_df = tree_df                         %>% 
          dplyr::select(CN, PLT_CN, PLOT,
                 TREE, SPCD, SPGRPCD,
                 DIA, ACTUALHT, SUBP, STATUSCD,
                 TPA_UNADJ)
## filter only live trees
tree_df = tree_df %>% filter(STATUSCD == 1)

## combine plot lon and lat
plot_df = filter(fia_recent$PLOT, CN %in% pltcn)
plot_df = plot_df                            %>% 
          dplyr::select(CN, PLOT, LAT, LON)  %>% 
          rename(PLT_CN = CN)
tree_df = tree_df %>% left_join(plot_df, by = c("PLT_CN","PLOT"))

## species name
spcd_df = filter(fia_recent$REF_SPECIES, SPCD %in% spcd)
spcd_df = spcd_df                     %>% 
          dplyr::select(SPCD, COMMON_NAME,
                SHARED_COMMON_NAME_IND,
                 GENUS, SPECIES,
                 SCIENTIFIC_NAME,
                 SPECIES_SYMBOL)
tree_df = tree_df %>% left_join(spcd_df,by=c("SPCD"))

## filter out species that are not pine, cedar, fir, and oak
# check if there are shared common names for different species
unique(tree_df$SHARED_COMMON_NAME_IND)

#since no shared common names, we just use common names for species searching
cname = unique(tree_df$COMMON_NAME)
pines = cname[grepl("pine",cname)]
cedars = cname[grepl("cedar", cname)]
firs = cname[grepl("fir", cname)]
oaks = cname[grepl("oak", cname)]
oaks = oaks[oaks[]!="tanoak"]
keep_sp = c(pines,cedars,firs,oaks)
tree_df = filter(tree_df, COMMON_NAME %in% keep_sp)

## convert dbh to cm and only keep necessary variables
tree_df = tree_df                   %>% 
          dplyr::select(CN,PLT_CN,
                       PLOT, LAT,
                       LON,TREE,
                       COMMON_NAME,
                       SCIENTIFIC_NAME,
                       TPA_UNADJ,
                       DIA)        %>% 
         mutate(DIA = DIA*inch_2_cm)

## assign functional group
tree_df = tree_df                  %>% 
          mutate(PFT = case_when(
          COMMON_NAME %in% pines ~ "pine",
          COMMON_NAME %in% cedars ~ "cedar",
          COMMON_NAME %in% firs ~ "fir",
          COMMON_NAME %in% oaks ~ "oak"
          ))
## filter out row where DIA is N/A
tree_df = filter(tree_df, !is.na(DIA))

## round expansion factor to nearest integer 
tree_df = tree_df            %>% 
          mutate(TPA_UNADJ = round(TPA_UNADJ))



### also get plot level total basal area for selected plots

tpa_byplot = filter(tpa_byplot, PLT_CN %in% pltcn)
tpa_byplot = tpa_byplot                                  %>% 
             mutate(baa_plt = baa_plt*ft2acre_2_m2ha,
                    tpa_plt = tpa_plt / acre_2_ha)

tree_df = tree_df %>% left_join(tpa_byplot, by=c("PLT_CN"))

data.table::fwrite(tree_df, file.path(data_path,"fia-tree-data-filtered-cancov-masked_2000-2017.csv"))

## search for the nearest plots for each WRF grid
## calculate mean BA of every two plots
## select the pair which has a mean BA within [0.9*lemma BA, 1.1*lemma BA]
## if multiple pairs present, randomly draw one pair

fia_xy = tree_df                       %>% 
         dplyr::select(LON,LAT,PLT_CN,
                       baa_plt) %>% 
         distinct()

nn_pt = RANN::nn2(fia_xy[, 1:2], lemma_ba[, 1:2], searchtype = 'radius',
                  radius = 0.5, k = 10)
nn_idx = nn_pt$nn.idx
nn_dis = nn_pt$nn.dists
nn_idx = data.frame(nn_idx)
nn_dis = data.frame(nn_dis)



nn_dis = nn_dis                  %>% 
         mutate_all( ~ assign_na(.))
nn_dis = nn_dis    %>% 
         mutate(mean_dis = rowMeans(.,na.rm=TRUE))
hist(nn_dis$mean_dis)

cn = as.vector(fia_xy$PLT_CN)
lemma_ba = data.frame(lemma_ba)
lemma_ba = cbind(lemma_ba, nn_idx)
lemma_ba = lemma_ba                %>% 
           rename_at(vars(X1:X10), ~rename_fun(.))

lemma_ba = lemma_ba                 %>% 
         mutate(across(c(CN1:CN10), ~replace_fun(.,cn)))

## map plot total basal area from each nearest FIA plot to WRF grids
lemma_ba = lemma_ba                                         %>%
           mutate(across(starts_with("CN"), 
                ~ fia_xy$baa_plt[match(.,fia_xy$PLT_CN)],
                .names = "ba{.col}"))

## compute mean basal area between every two 
## possible combinations of FIA plots, and create a
## new column to store the mean BA for that pair

ba_cols <- grep("^ba", names(lemma_ba), value = TRUE)
combs <- combn(ba_cols, 5, simplify = FALSE)

combs_means <- setNames(
  as.data.frame(sapply(combs, function(p) rowMeans(lemma_ba[, p], na.rm = TRUE))),
  sapply(combs, function(p) paste("mean",p[1],p[2],p[3],p[4],p[5],sep="_"))
)

lemma_ba = cbind (lemma_ba, combs_means)

relative_ba = function(x,ba){ifelse(is.na(x) | is.na(ba), NA_integer_, x/ba)}

## calculate relative difference between lemma grid BA and mean BA between two FIA plots
lemma_ba = lemma_ba                     %>% 
           mutate(across(starts_with("mean"), ~ relative_ba(.,BA), .names="rela_{.col}"))

lemma_ba = lemma_ba                        %>% 
           mutate(across(starts_with("rela"), ~ ifelse(.>= 0.9 & .<= 1.1, TRUE, FALSE)))

pat  = "^rela_mean(_baCN[0-9]+)+$"
subset_cols = grep(pat, names(lemma_ba), value = TRUE)

# 1. return pair where the first column has BA within range
lemma_ba = lemma_ba %>%
           rowwise() %>%
           mutate(first_true = {
             vals = c_across(all_of(subset_cols))
             idx  = which(vals)
             if (length(idx) > 0) subset_cols[idx[1]] else NA_character_
             })

# 2. return pair where BA is within the range and the sum
# of the numbers in pair name is lowest (meaning distance to target grid is closer), which 
# is also ranked by first and second number if sum is same, meaning we are 
# always taking pairing plots that are averagely closest to the target
# TO-DO: can directly compare mean distance to target 

m      = regexec("baCN(\\d+)_baCN(\\d+)", subset_cols)
m      = gregexpr("(?<=baCN)\\d+", subset_cols, perl = TRUE)
parts  = regmatches(subset_cols, m)
n1     = as.integer(vapply(parts, function(x) x[1], ""))
n2     = as.integer(vapply(parts, function(x) x[2], ""))
n3     = as.integer(vapply(parts, function(x) x[3], ""))
n4     = as.integer(vapply(parts, function(x) x[4], ""))
n5     = as.integer(vapply(parts, function(x) x[5], ""))
ord    = order(n1 + n2 + n3 + n4 + n5, n1, n2, n3, n4, n5)
rank   = match(seq_along(subset_cols), ord)
mat    = as.matrix(lemma_ba[subset_cols])
storage.mode(mat) = "logical"
mat[is.na(mat)] = FALSE

lemma_ba$lowest_true_col = apply(mat, 1, function(x) {
  idx <- which(x)
  if (length(idx) == 0) return(NA_character_)
  subset_cols[ idx[ which.min(rank[idx]) ] ]
})

## store selected FIA plots

#lemma_ba$pair_CN1 = as.integer(sub("^.*baCN(\\d+)_baCN\\d+$", "\\1", lemma_ba$lowest_true_col))
#lemma_ba$pair_CN2 = as.integer(sub("^.*baCN(\\d+)$", "\\1", lemma_ba$lowest_true_col))
lemma_ba = lemma_ba                   %>% 
           mutate(comb_CN1 = nth_baCN(lowest_true_col, 1),
                  comb_CN2 = nth_baCN(lowest_true_col, 2),
                  comb_CN3 = nth_baCN(lowest_true_col, 3),
                  comb_CN4 = nth_baCN(lowest_true_col, 4),
                  comb_CN5 = nth_baCN(lowest_true_col, 5))
lemma_ba = lemma_ba           %>% 
           mutate(across(starts_with("comb"), ~ paste0("CN",.)))

## replace with FIA CN number
ba_cols <- grep("^CN\\d+$", names(lemma_ba), value = TRUE)
ba_mat <- as.matrix(lemma_ba[ba_cols])
col_idx1 <- match(lemma_ba$comb_CN1, ba_cols)
col_idx2 =  match(lemma_ba$comb_CN2, ba_cols)
col_idx3 =  match(lemma_ba$comb_CN3, ba_cols)
col_idx4 =  match(lemma_ba$comb_CN4, ba_cols)
col_idx5 =  match(lemma_ba$comb_CN5, ba_cols)
row_idx <- seq_len(nrow(lemma_ba))
lemma_ba$comb_CN1 <- ba_mat[cbind(row_idx, col_idx1)]
lemma_ba$comb_CN2 <- ba_mat[cbind(row_idx, col_idx2)]
lemma_ba$comb_CN3 <- ba_mat[cbind(row_idx, col_idx3)]
lemma_ba$comb_CN4 <- ba_mat[cbind(row_idx, col_idx4)]
lemma_ba$comb_CN5 <- ba_mat[cbind(row_idx, col_idx5)]

lemma_ba = lemma_ba %>% dplyr::select(lon, lat, BA, lowest_true_col, 
                                      comb_CN1, comb_CN2, comb_CN3, 
                                      comb_CN4, comb_CN5)
lemma_ba = data.frame(lemma_ba)

final_set = lemma_ba                     %>% 
            dplyr::select(lon,lat,BA,
            lowest_true_col,comb_CN1)    %>% 
            rename(PLT_CN = comb_CN1)    %>% 
            left_join(tree_df, 
            by = "PLT_CN")   %>% 
            mutate(patch = "FIA_p1")

final_set2 = lemma_ba                    %>% 
             dplyr::select(lon,lat,BA,
             lowest_true_col,comb_CN2)   %>% 
             rename(PLT_CN = comb_CN2)   %>% 
             left_join(tree_df, 
                       by = "PLT_CN")    %>% 
             mutate(patch = "FIA_p2")

final_set3 = lemma_ba                    %>% 
             dplyr::select(lon,lat,BA,
             lowest_true_col,comb_CN3)   %>% 
             rename(PLT_CN = comb_CN3)   %>% 
             left_join(tree_df, 
             by = "PLT_CN")              %>% 
             mutate(patch = "FIA_p3")

final_set4 = lemma_ba                    %>% 
             dplyr::select(lon,lat,BA,
             lowest_true_col,comb_CN4)   %>% 
             rename(PLT_CN = comb_CN4)   %>% 
             left_join(tree_df, 
             by = "PLT_CN")              %>% 
             mutate(patch = "FIA_p4")

final_set5 = lemma_ba                    %>% 
             dplyr::select(lon,lat,BA,
             lowest_true_col,comb_CN5)   %>% 
             rename(PLT_CN = comb_CN5)   %>% 
             left_join(tree_df, 
             by = "PLT_CN")              %>% 
             mutate(patch = "FIA_p5")
final_set = rbind(final_set, final_set2, final_set3, final_set4, final_set5)
final_set$lon_f = format(final_set$lon,format="f",digits=dec)
final_set$lat_f = format(final_set$lat,format="f",digits=dec)

data.table::fwrite(final_set, file.path(data_path,"selected-FIA-plots-SierraOnWRF_2000-2017-5plots.csv"))


## create initialization files for each target grids
wrf_grids = data.table::fread(file.path(data_path,"wrf_points_creekfire.csv"))
wrf_grids = data.frame(wrf_grids)
nn_pt2 = RANN::nn2(final_set[,1:2], wrf_grids[,1:2], k =1)
nn_idx = setNames(data.frame(nn_pt2$nn.idx), c("idx1"))
wrf_grids = cbind(wrf_grids, nn_idx)

wrf_grids = wrf_grids                                   %>% 
            mutate(
                   lon_f = final_set$lon_f[idx1],
                   lat_f = final_set$lat_f[idx1],
                   lon_wrf   = final_set$lon[idx1],
                   lat_wrf   = final_set$lat[idx1],
                   plt_name = seq_along(wrf_grids$idx1)
                   )
wrf_grids = wrf_grids %>% left_join(final_set, by = c("lon_f","lat_f"))

## since some WRF grids were not mapped to FIA plots during above
## FIA-WRF paring process, so we have to drop these grids now
## later we will find the nearest WRF grids that have FIA plots assigned to them
## to those that don't have, and use the same FIA plots for them 

wrf_grids_subset = wrf_grids %>% filter(!is.na(lowest_true_col))
wrf_grids_tomap = wrf_grids                       %>% 
                  filter(is.na(lowest_true_col))  %>% 
                  dplyr::select(-idx1)


nn_pt3 = RANN::nn2(wrf_grids_subset[,1:2], wrf_grids_tomap[,1:2], k =1)
nn_idx3 = setNames(data.frame(nn_pt3$nn.idx), c("idx1"))
wrf_grids_tomap = cbind(wrf_grids_tomap, nn_idx3)
wrf_grids_tomap = wrf_grids_tomap                   %>% 
                  dplyr::select(Longitude,
                         Latitude,
                         lon_wrf,lat_wrf,
                         lon_f,lat_f,
                         plt_name,
                         idx1)                      %>% 
                  mutate(lon_f_org = lon_f,
                         lat_f_org = lat_f)

wrf_grids_tomap = wrf_grids_tomap                                   %>% 
  mutate(
    lon_f = wrf_grids_subset$lon_f[idx1],
    lat_f = wrf_grids_subset$lat_f[idx1]
  )
to_join = wrf_grids_subset %>% dplyr::select(-c(Longitude,Latitude,plt_name,idx1,lon_wrf,lat_wrf))
wrf_grids_tomap = wrf_grids_tomap %>% left_join(to_join, by = c("lon_f","lat_f"))
wrf_grids_tomap = wrf_grids_tomap           %>% 
                  dplyr::select(-c(lon_f, lat_f))  %>% 
                  rename(lon_f = lon_f_org,
                         lat_f = lat_f_org)

wrf_grids = rbind(wrf_grids_subset, wrf_grids_tomap)
creek_wrf_grids = wrf_grids %>% dplyr::select(Longitude,Latitude,lon_wrf,lat_wrf)
data.table::fwrite(creek_wrf_grids, file.path(data_path,"creek-grids-5plots.csv"))

plot_year = 2017

## first create initialization files for grids that are mapped to FIA plots
nplots = length(unique(wrf_grids$plt_name))
plot_area = 1  #1 ha
for(p in sequence(nplots)){
  plot_now = unique(wrf_grids$plt_name)[p]
  sz = wrf_grids%>% filter(plt_name==plot_now)
  plot_name = paste0(unique(sz$lon_f), "_", unique(sz$lat_f))
  sz$status  = "A"
  
  sz = sz                       %>%  
    rename(dbh = DIA,
           pft = PFT)
  
  sz = sz                          %>% 
    mutate(pft = case_when(
      pft == "pine"  ~ 1,
      pft == "cedar" ~ 2,
      pft == "fir"   ~ 3,
      pft == "oak"   ~ 4
    ))
  
  
  ### create pss file (patch structure)
  
  npatch = length(unique(sz$patch))
  time   = rep(plot_year, npatch)
  patch_df = as.data.frame(time)
  patch_df$patch = seq(npatch)
  patch_df$trk = rep(2, npatch)
  patch_df$age = rep(0, npatch)
  patch_df$area = rep((1/npatch), npatch)
  patch_df = patch_df %>% dplyr::select(time,patch, trk, age, area)
  write.table(patch_df, file.path(data_path,sprintf('%s_%i.pss',  plot_name, plot_year)), 
              row.names=FALSE, sep = " ")
  ### create css file (cohort structure)
  co_df = sz
  co_df$time = plot_year
  co_df$patch = as.numeric(sub("FIA_p(\\d+)$", "\\1", co_df$patch))
  co_df$height = -1
  patch_size=plot_area * 10000 / npatch
  co_df$nplant = 1/patch_size
  co_df = co_df %>% mutate(nplant = nplant*TPA_UNADJ*0.5) # multiply expansion factor to get stem density for 1 acre 
                                                          # and convert to 0.2 ha (patch size)
  co_df = co_df %>% dplyr::select(time, patch,dbh,height,pft,nplant)
  co_df = co_df %>% mutate(dbh = round(dbh,2))
  write.table(co_df, file.path(data_path,sprintf('%s_%i.css', plot_name, plot_year)), 
              row.names=FALSE, sep = ' ')
}


## write the control text file for FATES 
init_dir_derecho = "/glade/work/xiugao/fates-input/creek-fire/init/"
ctrl_file = wrf_grids %>% dplyr::select(Longitude,Latitude,lon_f,lat_f) %>% distinct()
ctrl_file = ctrl_file                                       %>% 
            mutate(plot_name = paste0(lon_f,"_",lat_f),
                   plot_year = 2017,
                   pss_name = paste0(plot_name,"_", plot_year,".pss"),
                   css_name = paste0(plot_name,"_", plot_year,".css"),
                   pss_path = paste0(init_dir_derecho,pss_name),
                   css_path = paste0(init_dir_derecho,css_name),
                   #lon = 360 + lon,
                   type = 1)
ctrl_file = ctrl_file                       %>% 
            dplyr::select(type,Latitude,
                   Longitude,
                   pss_path, css_path)      %>%
            rename_all(toupper)
write.table(ctrl_file, file.path(data_path,"creek-init-ctrl-5plot.txt"),row.names = FALSE,quote=FALSE)

