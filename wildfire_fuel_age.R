library(FireHistory)

# regular workflow --------------------------------------------------------

district <- "albany"

district_vect <- paste0("./vectors/", district, "_ten_diss.shp")

aoi <- user_aoi(district_vect, name = district)

# data <- assemble_dataNEW(fire_path = "./corp_data/DBCA_Fire_History_DBCA_060.shp",
#                          from = 1988, to = 2022, aoi = aoi, accessed_on = "24/04/23")

data <- assemble_data(fire_path = "./corp_data/DBCA_Fire_History_DBCA_060.shp",
                      FYfrom = 1950, FYto = 2023, aoi = aoi, accessed_on = "28/04/23")



# additional packages -----------------------------------------------------

library(terra)
library(fs)




# setup folder structure --------------------------------------------------

f1 <- paste0("./", district, "/01_by")
f2 <- paste0("./", district, "/02_yob")
f3 <- paste0("./", district, "/03_tsf")
f4 <- paste0("./", district, "/04_wfm")
f5 <- paste0("./", district, "/05_wffa")

fols <- c(f1, f2, f3, f4)
fs::dir_create(fols)


# individual burn years ---------------------------------------------------

f_vecs <- terra::vect(data[["fh_alb"]]) %>%
  dplyr::arrange(fin_y)

# unique fire years
u_fyrs <- f_vecs %>%
  tidyterra::pull(var = fin_y) %>%
  unique()

#  can change resolution here
resolution <- 30
template <- terra::rast(data[["aoi_alb"]], res = resolution)
aoi_msk <- terra::rasterize(terra::vect(data[["aoi_alb"]]), template)
# aoi_msk <- data[["aoi_msk"]]

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


# time since last fire ----------------------------------------------------

tsf_noms <- fs::path(f3, gsub("yob", "tsf", dir(f2)))

for(i in seq_along(tsf_noms)){
  yob <- terra::rast(fs::dir_ls(f2)[i])
  tsf <- full_fyrs[i] - yob
  terra::writeRaster(tsf, tsf_noms[i])
}


# fuel age area matrix ----------------------------------------------------

pix_ha <- resolution^2/10000

tsf_stk <- terra::rast(fs::dir_ls(f3))
names(tsf_stk) <- full_fyrs

# total area
tot_area <- dplyr::as_tibble(terra::freq(aoi_msk)) %>%
  dplyr::mutate(area_ha = count * pix_ha) %>%
  dplyr::pull(area_ha)
# calc fuel age area stats
fuel_df <- tibble::tibble()
for(i in seq_along(full_fyrs)){
  out <- dplyr::as_tibble(terra::freq(tsf_stk[[i]])) %>%
    dplyr::mutate(layer = full_fyrs[i],
                  value = paste0("fa", value),
                  area_ha = count * pix_ha) %>%
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

readr::write_csv(fuel_mat, file = paste0("./", district, "/", district, "_fuel_age_area_matrix.csv"))


# wild fire masks ---------------------------------------------------------

wf_vecs <- terra::vect(data[["fh_alb"]]) %>%
  tidyterra::filter(fih_fire_t == "WF") %>%
  dplyr::arrange(fin_y)

# unique fire years
u_wfyrs <- wf_vecs %>%
  tidyterra::pull(var = fin_y) %>%
  unique()

wfm_noms <- fs::path(f4, paste0("wfm", u_wfyrs, ".tif"))

for(i in seq_along(u_wfyrs)){
  wf1 <- wf_vecs %>%
    tidyterra::filter(fin_y == u_wfyrs[i]) %>%
    terra::rasterize(template) %>%
    terra::crop(aoi_msk, mask = TRUE)
  names(wf1) <- u_wfyrs[i]
  terra::writeRaster(wf1, wfm_noms[i])
}


# fuel age burnt by wild fire ---------------------------------------------

wffa_noms <- fs::path(f5, gsub("wfm", "wffa", dir(f4)))
  
# can't calculate for wf if present in first year
if(u_wfyrs[1] == u_fyrs[1]){
  wf_iter <- u_wfyrs[-1]
} else {
  wf_iter <- u_wfyrs
}


for(yr in wf_iter){
  wfm <- terra::rast(fs::path(f4, paste0("wfm", yr, ".tif")))
  wf_burnt <- terra::rast(fs::path(f3, paste0("yslb", yr-1, ".tif"))) %>%
    terra::mask(wfm)
  terra::writeRaster(wf_burnt, fs::path(f5, paste0("wffa", yr, ".tif")))
}



wffa_stk <- terra::rast(fs::dir_ls(f5))
names(wffa_stk) <- wf_iter


# calc fuel age area stats
wffa_df <- tibble::tibble()
for(i in seq_along(wf_iter)){
  out <- dplyr::as_tibble(terra::freq(wffa_stk[[i]])) %>%
    dplyr::mutate(layer = paste0("wf", wf_iter[i]),
                  value = paste0("fa", value),
                  area_ha = count * 0.09) %>%
    dplyr::rename(fuel_age = value,
                  year = layer) %>%
    dplyr::select(-count)
  wffa_df <- dplyr::bind_rows(wffa_df, out)
}

# fuel age area matrix
wffa_mat <- wffa_df %>%
  tidyr::pivot_wider(names_from = fuel_age, values_from = area_ha) %>%
  replace(is.na(.), 0) %>%
  dplyr::mutate(wf_total = rowSums(pick(where(is.numeric), -year)),
                id = dplyr::row_number()) %>%
  dplyr::arrange(desc(id)) %>%
  dplyr::select(-id)

readr::write_csv(wffa_mat, file = "./wild_fire_fuel_age_area_matrix.csv")
