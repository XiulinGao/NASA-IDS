#### Prescribed ignitions over WRF domain ####
#### Description ####
#### This prescribed ignition file is generaged
#### using fire spread data with lon, lat, and date 
#### when a fire spreaded onto the grid
#### Author: Xiulin Gao ####
#### Date: 2026-02-11 ####

#### Packages ####
library(ncdf4)
library(tidyverse)
library(data.table)
library(lubridate)

#### Paths ####
main_path = file.path("~/Git-repos/NASA-IDS")
out_path = file.path(main_path,"data")
wrf_file = file.path("~/Google Drive/My Drive/9km-WRF-1980-2020/1981-01.nc")

#### Data file ####
spread_file = file.path(out_path,"detailed_wrf_grid_intersections_2020_centers.csv")


#### Processing ####
spread_df = fread(spread_file)
spread_df = spread_df                                    %>% 
            mutate(Date = gsub("PM", "", Date),
                   Date = paste(Date, "2020", sep=" "),
                   time = as.Date(Date, 
                                  format="%b %d %Y"),
                   ignition = 1)                   %>% 
            select(time, Longitude_center, Latitude_center, ignition)
spread_df = spread_df                                   %>% 
            mutate(lon = round(Longitude_center,5),
                   lat = round(Latitude_center,5))

spread_df$lon_f = factor(spread_df$lon,levels=unique(spread_df$lon))
spread_df$lat_f = factor(spread_df$lat,levels=unique(spread_df$lat))

wrf_grids = nc_open(wrf_file)
xlon = ncvar_get(wrf_grids,"LONGXY")
ylat = ncvar_get(wrf_grids,"LATIXY")
dummy = nc_close(wrf_grids)
x_vec = as.vector(xlon)
y_vec = as.vector(ylat)
n_grid = length(x_vec)

time_start = as.Date("2020-01-01")
time_end = as.Date("2020-12-31")
time_seq = seq(from = time_start, to = time_end, by="days")
time_seq_noleap = time_seq[format(time_seq, "%m-%d") != "02-29"]
n_time = length(time_seq_noleap)
tstamp = rep(time_seq_noleap, each=n_grid)


x_rep_time = rep(x_vec,time=n_time)
y_rep_time = rep(y_vec,time=n_time)
creek_ig = as.data.frame(cbind(x_rep_time,y_rep_time)) 
creek_ig$time = tstamp
creek_ig = creek_ig                            %>% 
           mutate(lon = round(x_rep_time, 5),
                  lat = round(y_rep_time, 5))


creek_ig$lon_f = factor(creek_ig$lon,levels=unique(creek_ig$lon))
creek_ig$lat_f = factor(creek_ig$lat,levels=unique(creek_ig$lat))

spread_df = spread_df %>% select(lon_f,lat_f,time,ignition)
creek_ig = creek_ig                                            %>% 
           left_join(spread_df,by=c("lon_f","lat_f","time"))   %>% 
           mutate(ignition = ifelse(is.na(ignition), 0, ignition))

#### Generate prescribed ignition file ####

ig_out = array(creek_ig$ignition, dim=c(dim(xlon)[1], dim(xlon)[2], length(time_seq_noleap)))
time_out = sequence(length(time_seq_noleap)) -1

xx  = ncdim_def( name="lon"   ,units="",vals= sequence(dim(xlon)[1])           , create_dimvar=FALSE)
yy  = ncdim_def( name="lat"   ,units="",vals= sequence(dim(xlon)[2])           , create_dimvar=FALSE)
tt  = ncdim_def( name="time"  ,units="",vals=sequence(length(time_seq_noleap)) , create_dimvar=FALSE)
nc_xy  = list   (xx,yy)
nc_xyt = list   (xx,yy,tt)
nc_t   = list   (tt)
xy     = c(dim(xlon)[1],dim(xlon)[2])
xyt    = c(dim(xlon)[1],dim(xlon)[2],length(time_seq_noleap))
curent_date = Sys.Date()
file_name = paste0("CreekFire-prescribed-ignition-",current_date,".nc")
file_name = file.path(out_path, file_name)
nc_vlist        = list()
nc_vlist$LONGXY = ncvar_def(  name      = "LONGXY"
                              , units    = "degrees_east"
                              , dim      = nc_xy
                              , missval  = undef
                              , longname = "longitude"
)#end ncvar_def
nc_vlist$LATIXY = ncvar_def( name       = "LATIXY"
                             , units    = "degrees_north"
                             , dim      = nc_xy
                             , missval  = undef
                             , longname = "latitude"
)#end ncvar_def
nc_vlist$time   = ncvar_def( name        = "time"
                             , units    = "days since 2020-01-01 00:00:00"
                             , dim      = nc_t
                             , missval  = undef
                             , longname = "observation time"
)#end ncvar_def
nc_vlist$lnfm   = ncvar_def( name        = "lnfm"
                             , units    = "number/km2/hr"
                             , dim      = nc_xyt
                             , missval  = undef
                             , longname = "number of ignition"
)#end ncvar_def

### define global attributes

att_template = list( title            = "To be replaced when looping through months"
                     , date_created   = paste0(as.character(now(tzone="UTC")), "UTC")
                     , source_code    = "Creek-fire-prescribed-ignition.R"
                     , code_notes     = "prescribed ignition for creek fire on WRF domain "
                     , code_developer = paste0( author_name
                                                ," <"
                                                , author_email
                                                ,">"
                     )#end paste0
                     , file_author    = paste0(author_name," <",author_email,">")
)#end list

nc_new <- nc_create(filename=file_name,vars=nc_vlist,verbose=FALSE)
dummy = ncvar_put(nc=nc_new,varid="LONGXY",vals=array(data=x_vec,dim=xy))
dummy = ncatt_put(nc=nc_new,varid="LONGXY",attname="mode",attval="time-invariant")
dummy = ncvar_put(nc=nc_new,varid="LATIXY", vals=array(data=y_vec, dim=xy))
dummy = ncatt_put(nc=nc_new,varid="LATIXY",attname="mode",attval="time-invariant")
dummy = ncvar_put(nc=nc_new, varid ="time",vals=time_out)
dummy = ncatt_put(nc=nc_new, varid ="time",attname="calendar",attval="noleap")
dummy = ncvar_put(nc=nc_new, varid ="lnfm",vals=ig_out)
dummy = ncatt_put(nc=nc_new, varid ="lnfm",attname="mode",attval="time-dependent")

nc_title   = "Prescribed ignition for Creek Fire on WRF domain "
att_global = modifyList( x = att_template, val = list( title = nc_title ))


# Loop through global attributes
for (l in seq_along(att_global)){
  # Current attribute information
  att_name  = names(att_global)[l]
  att_value = att_global[[l]]
  
  # Add attribute 
  dummy = ncatt_put(nc=nc_new,varid=0,attname=att_name,attval=att_value)
}#end for (l in seq_along(att_global))
nc_close(nc_new)

