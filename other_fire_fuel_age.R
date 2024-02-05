library(FireHistory)

# regular workflow --------------------------------------------------------

aoi <- user_aoi("./warren/vectors/warren_tenure_diss.shp", name = "warren")

# data <- assemble_dataNEW(fire_path = "./corp_data/DBCA_Fire_History_DBCA_060.shp",
#                          from = 1988, to = 2022, aoi = aoi, accessed_on = "24/04/23")

data <- assemble_data(fire_path = "./corp_data/DBCA_Fire_History_DBCA_060.shp",
                      FYfrom = 1950, FYto = 2023, aoi = aoi, accessed_on = "28/04/23")



# additional packages -----------------------------------------------------

library(terra)
library(fs)
library(glue)




# setup folder structure --------------------------------------------------

f1 <- "./01_by"
f2 <- "./02_yob"
f3 <- "./03_yslb"
f4 <- "./04_wfm"
f5 <- "./05_wffa"
f6 <- "./04_ofm"
f7 <- "./05_offa"

fols <- c(f1, f2, f3, f4)
fs::dir_create(fols)


# individual burn years ---------------------------------------------------

f_vecs <- terra::vect(data[["fh_alb"]]) %>%
  dplyr::arrange(fin_y)

# unique fire years
u_fyrs <- f_vecs %>%
  tidyterra::pull(var = fin_y) %>%
  unique()

template <- terra::rast(data[["aoi_alb"]], res = 30)
aoi_msk <- data[["aoi_msk"]]

for(i in seq_along(u_fyrs)){
  by <- f_vecs %>%
    tidyterra::filter(fin_y == u_fyrs[i]) %>%
    terra::rasterize(template, field = "fin_y") %>%
    terra::crop(aoi_msk, mask = TRUE)
  names(by) <- u_fyrs[i]
  nom <- glue::glue("by",  u_fyrs[i], ".tif")
  by_nom <- fs::path(f1, nom)
  terra::writeRaster(by, by_nom)
}

# find non burn years and add an infill blank year
zero_rst <- template
terra::values(zero_rst) <- 0
zero_yr <- zero_rst %>%
  terra::crop(aoi_msk, mask = TRUE)
minyr <- min(u_fyrs)
maxyr <- max(u_fyrs)
full_fyrs <- minyr:maxyr
missing <- full_fyrs[!(full_fyrs %in% u_fyrs)]

for(i in seq_along(missing)){
  nom <- glue::glue("by",  missing[i], ".tif")
  my_nom <- fs::path(f1, nom)
  terra::writeRaster(zero_yr, my_nom)
}


# year of burn ------------------------------------------------------------

yob_noms <- fs::path(f2, gsub("by", "yob", dir(f1)))

# start year same as first burn year
fs::file_copy(path = fs::dir_ls(f1)[1], new_path = yob_noms[1])

# 2nd yob
rst1 <- terra::rast(fs::dir_ls(f1)[1])
rst2 <- terra::rast(fs::dir_ls(f1)[2])
by_stk <- c(rst1, rst2)
yob <- terra::app(by_stk, fun ="max", na.rm = TRUE)
yob_nom <- yob_noms[2]
terra::writeRaster(yob, yob_nom)

# subsequent yob
yob_iter <- length(full_fyrs)-1
for(i in 2:yob_iter){
  if(i < yob_iter){
    rst1 <- terra::rast(yob_noms[i])
    rst2 <- terra::rast(fs::dir_ls(f1)[i+1])
    by_stk <- c(rst1, rst2)
    yob <- terra::app(by_stk, fun ="max", na.rm = TRUE)
    terra::writeRaster(yob, yob_noms[i+1])
  } else {
    rst1 <- terra::rast(yob_noms[i])
    rst2 <- terra::rast(fs::dir_ls(f1)[i+1])
    by_stk <- c(rst1, rst2)
    yob <- terra::app(by_stk, fun ="max", na.rm = TRUE)
    yob[yob == 0] <- NA
    terra::writeRaster(yob, yob_noms[i+1])
  }
}


# year since last burn ----------------------------------------------------

yslb_noms <- fs::path(f3, gsub("yob", "yslb", dir(f2)))

for(i in seq_along(yslb_noms)){
  yob <- terra::rast(fs::dir_ls(f2)[i])
  yslb <- full_fyrs[i] - yob
  terra::writeRaster(yslb, yslb_noms[i])
}


# fuel age area matrix ----------------------------------------------------

yslb_stk <- terra::rast(fs::dir_ls(f3))
names(yslb_stk) <- full_fyrs

# total area
tot_area <- dplyr::as_tibble(terra::freq(aoi_msk)) %>%
  dplyr::mutate(area_ha = count * 0.09) %>%
  dplyr::pull(area_ha)
# calc fuel age area stats
fuel_df <- tibble::tibble()
for(i in seq_along(full_fyrs)){
  out <- dplyr::as_tibble(terra::freq(yslb_stk[[i]])) %>%
    dplyr::mutate(layer = full_fyrs[i],
                  value = paste0("fa", value),
                  area_ha = count * 0.09) %>%
    dplyr::rename(fuel_age = value,
                  year = layer) %>%
    dplyr::select(-count)
  fuel_df <- dplyr::bind_rows(fuel_df, out)
}

# fuel age area matrix
fuel_mat <- fuel_df %>%
  tidyr::pivot_wider(names_from = fuel_age, values_from = area_ha) %>%
  replace(is.na(.), 0) %>%
  dplyr::mutate(total = rowSums(pick(where(is.numeric), -year)),
                unknown = tot_area - total,
                id = dplyr::row_number()) %>%
  dplyr::arrange(desc(id)) %>%
  dplyr::select(-id)

readr::write_csv(fuel_mat, file = "./fuel_age_area_matrix.csv")


# other fire masks --------------------------------------------------------

of_vecs <- terra::vect(data[["fh_alb"]]) %>%
  tidyterra::filter(fih_fire_t != "WF") %>%
  dplyr::arrange(fin_y)

# unique fire years
u_ofyrs <- of_vecs %>%
  tidyterra::pull(var = fin_y) %>%
  unique()

ofm_noms <- fs::path(f6, paste0("ofm", u_ofyrs, ".tif"))

for(i in seq_along(u_ofyrs)){
  of1 <- of_vecs %>%
    tidyterra::filter(fin_y == u_ofyrs[i]) %>%
    terra::rasterize(template) %>%
    terra::crop(aoi_msk, mask = TRUE)
  names(of1) <- u_ofyrs[i]
  terra::writeRaster(of1, ofm_noms[i])
}


# fuel age burnt by wild fire ---------------------------------------------

offa_noms <- fs::path(f7, gsub("ofm", "offa", dir(f6)))

# can't calculate for of if present in first year
if(u_ofyrs[1] == u_fyrs[1]){
  of_iter <- u_ofyrs[-1]
} else {
  of_iter <- u_ofyrs
}


for(yr in of_iter){
  ofm <- terra::rast(fs::path(f6, paste0("ofm", yr, ".tif")))
  of_burnt <- terra::rast(fs::path(f3, paste0("yslb", yr-1, ".tif"))) %>%
    terra::mask(ofm)
  terra::writeRaster(of_burnt, fs::path(f7, paste0("offa", yr, ".tif")))
}



offa_stk <- terra::rast(fs::dir_ls(f7))
names(offa_stk) <- of_iter


# calc fuel age area stats
offa_df <- tibble::tibble()
for(i in seq_along(of_iter)){
  out <- dplyr::as_tibble(terra::freq(offa_stk[[i]])) %>%
    dplyr::mutate(layer = paste0("other", of_iter[i]),
                  value = paste0("fa", value),
                  area_ha = count * 0.09) %>%
    dplyr::rename(fuel_age = value,
                  year = layer) %>%
    dplyr::select(-count)
  offa_df <- dplyr::bind_rows(offa_df, out)
}

# fuel age area matrix
offa_mat <- offa_df %>%
  tidyr::pivot_wider(names_from = fuel_age, values_from = area_ha) %>%
  replace(is.na(.), 0) %>%
  dplyr::mutate(of_total = rowSums(pick(where(is.numeric), -year)),
                id = dplyr::row_number()) %>%
  dplyr::arrange(desc(id)) %>%
  dplyr::select(-id)

readr::write_csv(offa_mat, file = "./other_fire_fuel_age_area_matrix.csv")
