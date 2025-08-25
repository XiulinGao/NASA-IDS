## this script create a land mask for specified region
## on the WRF grids 

library(tidyverse)
library(ncdf4)
library(stars)
library(lubridate)

dec = 9L
wrf_proj = "+proj=lcc +lat_1=30 +lat_0=38 
               +lon_0=-70 +lat_2=60 +R=6370000 
            +datum=WGS84 +units=m +no_defs"

author_name  = "Xiulin Gao"
author_email = "xiulingao@lbl.gov"
undef        = -9999


grid_file = file.path("~/NASA-IDS/data/creek-grids.csv")
wrf_path = file.path("~/Google Drive/My Drive/9km-WRF-1980-2020/1981-01.nc")
wrf_domain = nc_open(wrf_path)
XLONG = ncvar_get(wrf_domain,"LONGXY")
XLAT  = ncvar_get(wrf_domain, "LATIXY")
dummy = nc_close(wrf_domain)
x_vec    = as.vector(XLONG)
y_vec    = as.vector(XLAT)
wrf_df   = as.data.frame(cbind(x_vec,y_vec))

grids_active = data.table::fread(grid_file)
grids_active = data.frame(grids_active)
grids_active$mask = 1

sc = 1e9L

row_preserving_join = function(lhs, rhs, xlon, xlat, ylon, ylat, cols = NULL, scale = 1e6) {
  # build stable keys at fixed precision
  key_l <- paste0(round(lhs[[xlon]] * scale), "_", round(lhs[[xlat]] * scale))
  key_r <- paste0(round(rhs[[ylon]] * scale), "_", round(rhs[[ylat]] * scale))
  
  # row-wise index of rhs for each lhs row
  m <- match(key_l, key_r)
  
  # default: bring all non-key columns from rhs
  if (is.null(cols)) cols <- setdiff(names(rhs), c(ylon, ylat))
  
  cbind(lhs, rhs[m, cols, drop = FALSE])
}

mask = row_preserving_join(lhs = wrf_df, 
                           rhs = grids_active,
                           xlon = "x_vec",
                           xlat = "y_vec",
                           ylon = "lon_wrf",
                           ylat = "lat_wrf",
                           cols=NULL,
                           scale = sc)

fil_var   = matrix(as.factor(mask$mask), nrow=147,ncol=151)
x_arr     = matrix(mask$x_vec,nrow=147,ncol=151)
y_arr     = matrix(mask$y_vec,nrow=147,ncol=151)
wrf_star  = st_as_stars(fil_var)
wrf_star  = st_as_stars(wrf_star, curvilinear=list(X1=x_arr,X2=y_arr), crs="+init=EPSG:4326")
wrf_sf    = st_as_sf(wrf_star,as_points=FALSE,na.rm=FALSE)
ca_co     = USAboundaries::us_counties(resolution = "high", states = "CA")

###plot to see how the active domain looks like
  ggplot()                                                                      + 
  geom_sf(data=wrf_sf,colour="grey50", aes(fill=A1),lwd=0.5)                    +
  coord_sf(crs=st_crs(wrf_proj))                                                + 
  theme_bw() + theme(panel.ontop=TRUE, panel.background=element_blank())        +
  labs(x="",y="")                                                               +
  scale_fill_manual(values=c("1"="lightblue","0"= "grey50"),
                    guide="colorbar")                                           +
  geom_sf(data = ca_co, color = alpha("black", alpha=0.2),lwd=0.1,fill=NA)      


####create the mask NetCDF file
mask = mask %>% mutate(mask = ifelse(is.na(mask), 0, mask))
land_mask  = array(mask$mask,dim=c(147,151))
land_mkdif = land_mask


## create new nc file as land mask
xx  = ncdim_def( name="lon"   ,units="",vals= sequence(147)  ,create_dimvar=FALSE)
yy  = ncdim_def( name="lat"   ,units="",vals= sequence(151)  ,create_dimvar=FALSE)
nc_xy  = list   (xx,yy)
xy     = c(147,151)
file_name = file.path("~/NASA-IDS/data/creek-fire-domain-mask.nc")
nc_vlist        = list()
nc_vlist$LONGXY = ncvar_def(  name      = "lsmlon"
                              , units    = "degrees_east"
                              , dim      = nc_xy
                              , missval  = undef
                              , longname = "longitude"
)#end ncvar_def
nc_vlist$LATIXY = ncvar_def( name       = "lsmlat"
                             , units    = "degrees_north"
                             , dim      = nc_xy
                             , missval  = undef
                             , longname = "latitude"
)#end ncvar_def
nc_vlist$mask1   = ncvar_def( name      = "landmask"
                              , units    = "unitless"
                              , dim      = nc_xy
                              , missval  = undef
                              , longname = "mask for land domain, 1 being cell active"
)#end ncvar_def
nc_vlist$mask2   = ncvar_def( name      = "mod_lnd_props"
                              , units    = "unitless"
                              , dim      = nc_xy
                              , missval  = undef
                              , longname = "mask for modifying land property, 1 being active land cell"
)#end ncvar_def

### define global attributes

att_template = list( title            = "To be replaced when looping through months"
                     , date_created   = paste0(as.character(now(tzone="UTC")), "UTC")
                     , source_code    = "mask-wrf-domain.R"
                     , code_notes     = "land mask for Sierra region on WRF domain "
                     , code_developer = paste0( author_name
                                                ," <"
                                                , author_email
                                                ,">"
                     )#end paste0
                     , file_author    = paste0(author_name," <",author_email,">")
)#end list

nc_new <- nc_create(filename=file_name,vars=nc_vlist,verbose=FALSE)
dummy = ncvar_put(nc=nc_new,varid="lsmlon",vals=array(data=x_vec,dim=xy))
dummy = ncvar_put(nc=nc_new,varid="lsmlat", vals=array(data=y_vec, dim=xy))
dummy = ncvar_put(nc=nc_new, varid ="landmask",vals=land_mask)
dummy = ncvar_put(nc=nc_new, varid ="mod_lnd_props",vals=land_mkdif)


nc_title   = "Land mask for creek fire on WRF domain "
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
