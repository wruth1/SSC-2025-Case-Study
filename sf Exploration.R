

library(sf)
library(data.table)
library(dplyr)
library(magrittr)
library(ggplot2)
library(lubridate)
library(purrr)      # map / possibly
library(pbapply)




data_raw = fread("Glider data/cabot_20200717_114_delayed_corrected_v4.csv")

data_renamed = data_raw %>% rename(
    temp = sea_water_temperature,
    salt = sea_water_practical_salinity,
    density = sea_water_density,
    oxygen = micromoles_of_oxygen_per_unit_mass_in_sea_water
) %>% mutate(
    elevation = -depth
)

sf_data = st_as_sf(data_renamed, coords = c("longitude", "latitude"), crs = 4326)

plot(sf_data["time"], axes=TRUE)

set.seed(1)

start_ind = 100001
num_points = 10000

sf_data_small = sf_data[start_ind:(start_ind + num_points - 1),]

plot(sf_data_small["time"], axes=TRUE)


ggplot() +
  geom_sf(data = sf_data_small,
          aes(color = time),   # map attribute to colour
          size = 0.15) +
  scale_color_viridis_c(option = "turbo") +
  coord_sf(crs = 4326) +  # WGS 84, hide graticule box
  labs(
    x = "Longitude (°E)",
    y = "Latitude (°N)"
  ) +
  theme_minimal(base_size = 11)



# ---------------------------------------------------------------------------- #
#                                    Binning                                   #
# ---------------------------------------------------------------------------- #


# ------------------------------ Construct bins ------------------------------ #

sf_data_for_binning = st_transform(sf_data, 3857)

grid_info = st_make_grid(
    sf_data_for_binning,
    cellsize = 1000,
    square = TRUE
)
grid = st_sf(id = seq_along(grid_info), grid_info)


# ------------------------- Compute summaries in bins ------------------------ #

binned_means = aggregate(sf_data_for_binning, grid, mean)




# --------------------------- Plot binned summaries -------------------------- #

ggplot(binned_means) +
  geom_sf(aes(fill = time)) +
  theme_void()



# ----------------------------- Add map underlay ----------------------------- #

library(ggspatial)

ggplot(data = filter(binned_means, !is.na(time))) +
  annotation_map_tile(zoomin = 0, cachedir = "map_tiles_ggspatial") +
  geom_sf(aes(fill = time), alpha = 0.5) +
  scale_fill_viridis_c() +
  theme_void()




# ---------------------------------------------------------------------------- #
#                               Multiple Missions                              #
# ---------------------------------------------------------------------------- #



data_large_raw = data.table()

data_names = list.files("Glider data", pattern = ".csv")

num_datasets = length(data_names)
# num_datasets = 3

for(i in 1:num_datasets){
    print(paste0("Processing dataset ", i, " of ", num_datasets))
    this_data = fread(paste0("Glider data/", data_names[i]))
    this_data$mission = as.character(i)
    data_large_raw = rbind(data_large_raw, this_data)
}


data_large_renamed = data_large_raw %>%
    rename(
        temp = sea_water_temperature,
        salt = sea_water_practical_salinity,
        density = sea_water_density,
        oxygen = micromoles_of_oxygen_per_unit_mass_in_sea_water
    ) %>% mutate(
        elevation = -depth,
        .keep = "unused"
    )



# ----------------------- Screen for impossible values ----------------------- #

# # Non-positive oxygen concentration
# data_large_dirty = data_large %>%
#     filter(salt <= 0)


# this_mission = 14
# this_oxygen = data_large %>% filter(mission == this_mission) %>% pull(oxygen)
# inds_bad = which(this_oxygen <= 0)
# inds_bad = inds_bad[c(1, 2, 9:11, 17, 24, 35, 37)]

# bad_diffs = diff(inds_bad)
# cbind(seq_along(bad_diffs), bad_diffs)

# for(j in inds_bad){
#     plot(this_oxygen[(j-1000):(j+1000)], main = paste0("Mission ", this_mission, " index ", j))
#     # plot(this_oxygen[(j-1000):(j+1000)], main = paste0("Mission ", this_mission, " index ", j), ylim = c(0, max(this_oxygen)))
#     readline(prompt="Press [enter] to continue")
# }


data_large = filter(data_large_renamed, oxygen > 0)

# Relationship between elevation, density and others


# small_data = data_large[100000:110000,]
# ggplot(small_data, aes(x = density, y = elevation, color = oxygen)) + geom_point()



# ------------------- Plot relationships between variables ------------------- #

data_small = data_large[100000:110000,]
ggplot(data_small, aes(x = oxygen, y = elevation)) + geom_point()

data_local = filter(data_large, )





# ---------------------------------------------------------------------------- #
#                                Construct Bins                                #
# ---------------------------------------------------------------------------- #


sf_data_large = st_as_sf(data_large, coords = c("longitude", "latitude"), crs = 4326) #* Med



sf_data_for_binning = st_transform(sf_data_large, 3857)                               #* Slow

grid_info = st_make_grid(                                                             #* Fast
    sf_data_for_binning,
    # cellsize = c(10000, 10000), #!!!!!!!!!!!!!!!!! Low Resulution
    cellsize = c(2500, 2500),
    # cellsize = c(1000, 1000), #!!!!!!!!!!!!!!!!! High Resulution
    square = TRUE
)
grid = st_sf(id = seq_along(grid_info), grid_info)                                    #* Fast
# plot(grid)



# ------------------------- Compute summaries in bins ------------------------ #
mean_or_name = function(x) {
  browser()
    if (is.numeric(x)) {
        mean(x, na.rm = TRUE)
    } else if (is.POSIXct(x)){
        mean(x, na.rm = TRUE)
    } else {
        all_missions = unique(x)
        if (length(all_missions) == 1) {
            all_missions
        } else {
            "Multiple"
        }
    }
}

binned_means = aggregate(sf_data_for_binning, grid, mean_or_name)                 #* Slow

test = aggregate(sf_data_for_binning, grid, function(x) list(x))                 #* Slow

binned_means$mission






# --------------------------- Plot binned summaries -------------------------- #

ggplot(binned_means) +
  geom_sf(aes(fill = time)) +
  theme_void()


ggplot(binned_means) +
  geom_sf(aes(fill = mission)) +
  theme_void()


# ----------------------------- Add an actual map ---------------------------- #

library(ggspatial)


ggplot(data = filter(binned_means, !is.na(mission))) +
  annotation_map_tile(zoomin = 0, cachedir = "map_tiles_ggspatial") +
  geom_sf(aes(fill = mission), alpha = 0.5) +
  scale_fill_viridis_d() +
  theme_void()


# Make the "Multiple" category stand-out from the other missions.

library(forcats)
library(scico)

target_mission = "Multiple"

# base_palatte = met.brewer("Isfahan1", 15)
base_palatte = scico(n = 15, palette = "batlow")
standout_colour = "#E41A1C"
NA_colour = "#808080"

mission_levels = unique(binned_means$mission)
full_palatte = setNames(
    c(base_palatte, standout_colour, NA_colour)[match(mission_levels, c(setdiff(mission_levels, c(target_mission, NA)), target_mission, NA))],
    mission_levels
)

data_bin_plot_mission = binned_means |>
                          filter(!is.na(mission)) |>
                          mutate(mission_value = ifelse(mission == target_mission, 1000, as.integer(mission)),
                                  mission = fct_reorder(mission, mission_value, .fun = unique))




pdf("Mission Map.pdf", width = 15, height = 9.5)

ggplot(data = data_bin_plot_mission) +
  annotation_map_tile(zoomin = 0, cachedir = "map_tiles_ggspatial") +
  geom_sf(aes(fill = mission), alpha = 0.5) +
  scale_fill_viridis_d() +
  scale_fill_manual(values = full_palatte, guide = guide_legend(ncol=2)) +
  theme_void() + labs(fill = "Mission") +
  theme(legend.title = element_text(face = "bold"))

dev.off()


# --------------------------- Plot other variables --------------------------- #
ggplot(data = filter(binned_means, !is.na(oxygen))) +
  annotation_map_tile(zoomin = 0, cachedir = "map_tiles_ggspatial") +
  geom_sf(aes(fill = oxygen), alpha = 0.5) +
  scale_fill_viridis_c() +
  theme_void()






# ------------------------- Identify points with bins ------------------------ #

data_bins = binned_means %>% filter(!is.na(mission)) %>% mutate(id = row_number())          #* Fast

data_large_bin_tags = st_join(sf_data_for_binning, data_bins["id"], join = st_within)       #* Slow





# ---------------------------------------------------------------------------- #
#                           Linear models within bins                          #
# ---------------------------------------------------------------------------- #

# --------------------------------- Simple lm -------------------------------- #

# Actual analysis
library(broom)
library(tidyr)

all_bin_info = data_large_bin_tags %>%
        st_drop_geometry() %>%
        nest_by(id) %>%
        mutate(
            fit = list(lm(oxygen ~ elevation, data = data)),
            stats = list(glance(fit)),
            coef = list(tidy(fit))
        ) %>%
        ungroup()


bins_R2 = all_bin_info %>%
            select(id, stats) %>%
            unnest(stats) %>%
            select(id, r.squared)

data_bins_with_R2 = left_join(data_bins, bins_R2, by = "id")



R2_plot = ggplot(data = data_bins_with_R2) +
  annotation_map_tile(zoomin = 0, cachedir = "map_tiles_ggspatial") +
  geom_sf(aes(fill = r.squared), alpha = 0.5) +
  scale_fill_viridis_c() +
  theme_minimal()


# Add a single point
extra_point_lon_lat = data.frame(lat = 48.7, lon = -62.5)   # High R2
extra_point_lon_lat = data.frame(lat = 47.4, lon = -60.33)   # Low R2
extra_point_lon_lat = data.frame(lat = 47.58, lon = -63.33)   # Another low R2
extra_point = extra_point_lon_lat %>%
                st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
                st_transform(crs = st_crs(binned_means))

R2_plot + geom_sf(data = extra_point,
          shape = 21,        # filled circle
          size  = 3,
          fill  = "red",
          colour = "black",
          stroke = 0.4)     # outline thickness


# ----------------------- Explore some interesting bins ---------------------- #

# Find bin id
# bin_data = st_transform(data_bins, 4326) %>%
#                         # mutate(across(everything(), ~ .x)) %>%
#                         st_drop_geometry()

# lon_lat = st_transform(data_bins, 4326) %>%
#             st_coordinates() %>%
#             as_tibble() %>%
#             group_by(L2) %>%
#             summarise(lon = mean(X), lat = mean(Y))
# bin_data_lon_lat = cbind(bin_data,
#                              tibble(lon = lon_lat$lon,
#                                     lat = lon_lat$lat))

library(pbapply)

all_dist_to_extra_point = pbsapply(1:nrow(data_bins), function(i) st_distance(
    data_bins[i,],
    extra_point
))

# hist(all_dist_to_extra_point)
# min(all_dist_to_extra_point)
arr_ind_min_dist = which(all_dist_to_extra_point == min(all_dist_to_extra_point))
id_min_dist = data_bins$id[arr_ind_min_dist]





# Extract and plot bin data
data_this_bin = filter(data_large_bin_tags, id == id_min_dist)
# data_high_R2 = data_this_bin
# data_low_R2 = data_this_bin
other_data_low_R2 = data_this_bin

ggplot(data = data_this_bin, aes(x = oxygen, y = elevation, colour = mission)) + geom_point()


this_fit = lm(oxygen ~ elevation, data = data_this_bin)
summary(this_fit)
plot(this_fit)



# ---------------------------------------------------------------------------- #
#                      Add fixed effect for mission to lm                      #
# ---------------------------------------------------------------------------- #


# ------------------------ Bins with multiple missions ----------------------- #

ids_mult_missions = data_bins |>
                        filter(mission == "Multiple") |>
                        pull(id)



min_n <- 50                 # need at least k+1 rows (k default = 10 → 11 pts)
default_k <- 10             # basis dimension for s()

some_ids = ids_mult_missions[11:20]


# -------------------------------- Fit models -------------------------------- #

tic()
lm_tbl <- data_large_bin_tags |>
#   filter(id %in% some_ids) |>
  filter(id %in% ids_mult_missions) |>
  st_drop_geometry() |>
  group_by(id) |>
  filter(n() >= min_n) |>
  nest() |>
  mutate(
    
    fit   = map(data,
                 possibly(~ lm(oxygen ~ elevation + mission,
                                 data = .x),
                          otherwise = NULL)),

    summ = map(fit, \(m)
                   if (is.null(m)) NULL
                   else summary(m)),
                   
    r2 = map_dbl(summ, \(s)
            if (is.null(s)) NA_real_
            else {
            this_summ = summ[[1]]
            this_summ$r.sq
            }),
    adj_r2 = map_dbl(summ, \(s)
            if (is.null(s)) NA_real_
            else {
            this_summ = summ[[1]]
            this_summ$adj.r.sq
            }),
  ) |>
  ungroup()
toc()


# ---------------------- Plot R2 for mission-adjusted lm --------------------- #

# Add mission-adjusted R2 to bin data
data_bins_mult = data_bins |>
                    left_join(bins_R2, by = "id") |>
                    filter(id %in% ids_mult_missions) |>
                    left_join(select(lm_tbl, id, r2), by = "id") |>
                    rename(simple_R2 = r.squared, mission_R2 = r2) |>
                    mutate(R2_ratio = mission_R2 / simple_R2, log_R2_ratio = log10(R2_ratio))


# Build plots

## Simple R2 (O2 ~ elevation)
ggplot(data = data_bins_mult) +
  annotation_map_tile(zoomin = 0, cachedir = "map_tiles_ggspatial") +
  geom_sf(aes(fill = simple_R2), alpha = 0.5) +
  scale_fill_viridis_c() + ggtitle("Simple R2: O2 ~ elevation") +
  theme_minimal()

## Mission-adjusted R2 (O2 ~ elevation + mission)
ggplot(data = data_bins_mult) +
    annotation_map_tile(zoomin = 0, cachedir = "map_tiles_ggspatial") +
    geom_sf(aes(fill = mission_R2), alpha = 0.5) +
    scale_fill_viridis_c() + ggtitle("Mission-Adjusted R2: O2 ~ elevation + mission") +
    theme_minimal()

## log-Ratio of R2s (mission / simple)
ggplot(data = data_bins_mult) +
    annotation_map_tile(zoomin = 0, cachedir = "map_tiles_ggspatial") +
    geom_sf(aes(fill = log_R2_ratio), alpha = 0.5) +
    scale_fill_viridis_c() + ggtitle("log-Ratio of R2s (base 10): Mission / Simple") +
    theme_minimal()


# --------------------- Explore bins with large R2 jumps --------------------- #

data_bins_mult |>
    filter(R2_ratio > 10 & R2_ratio < 500000) |>
    pull(R2_ratio) |>
    hist(breaks = 100)


id_large_jump = data_bins_mult |>
                    filter(R2_ratio > 1.9e+04 & R2_ratio < 2.1e+04) |>
                    pull(id)

this_data_jump = data_large_bin_tags |>
                    filter(id == id_large_jump)

ggplot(data = this_data_jump, mapping = aes(x = elevation, y = oxygen, colour = time, shape = mission)) + geom_point() 
ggplot(data = this_data_jump, mapping = aes(x = elevation, y = oxygen, colour = mission)) + geom_point() 

str(data_large_bin_tags)


# ---------------------------------------------------------------------------- #
#                                Explore Splines                               #
# ---------------------------------------------------------------------------- #


library(splines)
library(mgcv)
library(tictoc)

tic()
this_spline = smooth.spline(data_this_bin$elevation, data_this_bin$oxygen)
toc()

plot(data_this_bin$elevation, data_this_bin$oxygen)
lines(this_spline)


tic()
fit_GAM = gam(oxygen ~ s(elevation), data = data_this_bin, method = "GCV.Cp")
toc()

summary(fit_GAM)

plot(data_this_bin$elevation, data_this_bin$oxygen)
plot(fit_GAM, residuals = TRUE)

q = summary(fit_GAM)
q$edf



tic()
this_fit_gam = gam(oxygen ~ s(elevation), data = data_this_bin)
this_DF = summary(this_fit_gam)$edf
toc()



tic()
all_DF = sum(gam(oxygen ~ s(elevation), data = data_this_bin)$edf) -1 
toc()

summary(this_fit_gam)$s.table




# ---------------------------- Loop over all bins ---------------------------- #



min_n <- 50                 # need at least k+1 rows (k default = 10 → 11 pts)
default_k <- 10             # basis dimension for s()


tic()
gam_tbl <- data_large_bin_tags |>
  st_drop_geometry() |>
  group_by(id) |>
  filter(n() >= min_n) |>
  nest() |>
  mutate(
    
    fit   = map(data,
                 possibly(~ gam(oxygen ~ s(elevation),
                                 data = .x),
                          otherwise = NULL)),

    summ = map(fit, \(m)
                   if (is.null(m)) NULL
                   else summary(m)),
    
    edf_smooth = map_dbl(summ, \(s)
                   if (is.null(s)) NA_real_
                   else {
                     # First smooth term's EDF (oxygen ~ s(elevation))
                     this_summ = summ[[1]]
                     this_summ$s.table[1, "edf"]
                   }),
                   
    r2 = map_dbl(summ, \(s)
            if (is.null(s)) NA_real_
            else {
            this_summ = summ[[1]]
            this_summ$r.sq
            })
  ) |>
  ungroup()
toc()


spline_DFs = gam_tbl$edf_smooth
hist(spline_DFs)


# ---------------------------------- Plot DF --------------------------------- #

# Add edf_smooth to data_bins
data_bin_DF = left_join(data_bins, select(gam_tbl, id, edf_smooth), by = "id") |> rename(DF = edf_smooth)


Spline_DF_plot = ggplot(data = data_bin_DF) +
  annotation_map_tile(zoomin = 0, cachedir = "map_tiles_ggspatial") +
  geom_sf(aes(fill = DF), alpha = 0.5) +
  scale_fill_viridis_c() +
  theme_minimal()
Spline_DF_plot


# --------------------------------- Bin sizes -------------------------------- #

# Number of measurements per bin
all_bin_sizes = data_large_bin_tags |>
  st_drop_geometry() |>
  group_by(id) |>
  summarise(n = n()) 
nrow(all_bin_sizes)



# Add bin sizes to gam_tbl

gam_tbl_sizes = left_join(gam_tbl, all_bin_sizes, by = "id")
hist(gam_tbl_sizes$n)

with(gam_tbl_sizes, plot(n, edf_smooth))
with(gam_tbl_sizes, plot(n, r2))

cor(gam_tbl_sizes$n, gam_tbl_sizes$edf_smooth)






# ---------------------------------------------------------------------------- #
#                              Explore Periodicity                             #
# ---------------------------------------------------------------------------- #


# -------------------------------- One Mission ------------------------------- #

this_data = filter(data_large_bin_tags, mission == 1)
this_data_small = this_data[1:10000,]

this_data_small |> pull(time) |> month(label = TRUE)

ggplot(data = this_data_small, mapping = aes(x = time, y = oxygen)) + geom_point() 

data_large_bin_tags %>% pull(time) %>% month(label = TRUE) %>% plot()
all_months = data_large_bin_tags %>% pull(time) %>% month(label = TRUE)
table(all_months)

# Auto-correlation plot
this_O2 = this_data$oxygen
acf(this_O2, lag.max = 10000)


# Simple Fourier analysis
this_fourier = fft(this_O2)
this_n = length(this_O2)  
this_freqs = seq_len(this_n) / this_n
plot(this_freqs[-1], abs(this_fourier)[-1])

which.max(abs(this_fourier))



# Time series decomposition?
# Doesn't work. Need to specify seasonal period in frequency argument to ts().
library(forecast)

this_TS = ts(this_O2)
decomp = stl(this_TS)



# Periodogram
periodogram_1 = cpgram(this_TS)

periodogram_2 = spec.pgram(this_O2)
