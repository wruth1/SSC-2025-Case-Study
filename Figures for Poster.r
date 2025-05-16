

library(sf)
library(data.table)
library(dplyr)
library(magrittr)
library(ggplot2)
library(lubridate)
library(purrr)      # map / possibly
library(pbapply)
library(ggspatial)
library(forcats)
library(scico)
library(forecast)
library(splines)
library(mgcv)
library(tictoc)
library(pbapply)
library(broom)
library(tidyr)
library(emmeans)
library(tictoc)


# ---------------------------------------------------------------------------- #
#                                   Load Data                                  #
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




data_large = filter(data_large_renamed, oxygen > 0)

# # Optionally, subsample the large dataset
set.seed(1)
data_large = slice_sample(data_large, n = 100000)




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

# test = aggregate(sf_data_for_binning, grid, function(x) list(x))                 #* Slow

# binned_means$mission


# Confirm that binning worked correctly
ggplot(binned_means) +
  geom_sf(aes(fill = mission)) +
  theme_void()





# ------------------------- Identify points with bins ------------------------ #

data_bins = binned_means %>% filter(!is.na(mission)) %>% mutate(id = row_number())          #* Fast

data_large_bin_tags = st_join(sf_data_for_binning, data_bins["id"], join = st_within)       #* Slow



# ---------------------------------------------------------------------------- #
#                         Plot bins overlayed on a map                         #
# ---------------------------------------------------------------------------- #


# Make the "Multiple" category stand-out from the other missions.

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




pdf("Figures/Mission Map.pdf", width = 15, height = 9.5)

ggplot(data = data_bin_plot_mission) +
  annotation_map_tile(zoomin = 0, cachedir = "map_tiles_ggspatial") +
  geom_sf(aes(fill = mission), alpha = 0.5) +
  scale_fill_viridis_d() +
  scale_fill_manual(values = full_palatte, guide = guide_legend(ncol=2)) +
  theme_void() + labs(fill = "Mission") +
  theme(legend.title = element_text(size = 20), legend.text = element_text(size = 15))

dev.off()




# ---------------------------------------------------------------------------- #
#                                  Fit models                                  #
# ---------------------------------------------------------------------------- #

#? Terminology:
#?      lin vs spl: Linear vs Spline
#?      bas vs main vs int: No mission vs mission main effect vs mission interaction with elevation

check_multiple_missions = function(data) length(unique(data$mission)) > 1


min_n <- 20                 # need at least k+1 rows for a spline (k default = 10 → 11 pts)

ids_large_n = data_large_bin_tags |>
                st_drop_geometry() |>
                group_by(id) |>
                summarise(n = n()) |>
                filter(n >= min_n) |>
                pull(id)


# ------------------------------ Linear - Basic ------------------------------ #

tic()
r2_info_lin_bas <- data_large_bin_tags |>
  filter(id %in% ids_large_n) |>
  st_drop_geometry() |>
  group_by(id) |>
  # filter(n() >= min_n) |>
  nest() |>
  mutate(
    
    fit   = map(data,
                 possibly(~ lm(oxygen ~ elevation, data = .x),
                          otherwise = NULL)),

    summ = map(fit, \(m)
                   if (is.null(m)) NULL
                   else summary(m)),
                   
    r2_lin_bas = map_dbl(summ, \(s)
            if (is.null(s)) NA_real_
            else {
            this_summ = summ[[1]]
            this_summ$r.sq
            }),
    adj_r2_lin_bas = map_dbl(summ, \(s)
            if (is.null(s)) NA_real_
            else {
            this_summ = summ[[1]]
            this_summ$adj.r.sq
            }),
  ) |>
  select(-data, -fit, -summ) |>
  ungroup()
toc()


# ----------------------------- Linear - Mission ----------------------------- #


tic()
r2_info_lin_main <- data_large_bin_tags |>
#   filter(id %in% some_ids) |>
  # filter(id %in% ids_mult_missions) |>
  filter(id %in% ids_large_n) |>
  st_drop_geometry() |>
  group_by(id) |>
  # filter(n() >= min_n) |>
  nest() |>
  mutate(

    n_missions = map_dbl(data, ~ n_distinct(.x$mission)),
    
    fit   = map(data,
                 possibly(~ if(check_multiple_missions(.x)) {
                                  lm(oxygen ~ elevation + mission, data = .x)
                              } else {
                                  lm(oxygen ~ elevation, data = .x)
                              },
                          otherwise = NULL)),

    summ = map(fit, \(m)
                   if (is.null(m)) NULL
                   else summary(m)),
                   
    r2_lin_main = map_dbl(summ, \(s)
            if (is.null(s)) NA_real_
            else {
            this_summ = summ[[1]]
            this_summ$r.sq
            }),
    adj_r2_lin_main = map_dbl(summ, \(s)
            if (is.null(s)) NA_real_
            else {
            this_summ = summ[[1]]
            this_summ$adj.r.sq
            }),
  ) |>
  select(-data, -fit, -summ) |>
  ungroup()
toc()


# #* Effect on adjusted R2 of adding a mission main effect. Only plot for bins with > 1 mission
# data_lin_bas_plus_main = left_join(r2_info_lin_bas, r2_info_lin_main, by = "id") |>
#                 filter(n_missions > 1)

# with(data_lin_bas_plus_main, plot(adj_r2_lin_bas, adj_r2_lin_main))
# abline(a = 0, b = 1, col = "red")


# --------------------------- Linear - Interaction --------------------------- #

tic()
r2_info_lin_int <- data_large_bin_tags |>
#   filter(id %in% some_ids) |>
  # filter(id %in% ids_mult_missions) |>
  filter(id %in% ids_large_n) |>
  st_drop_geometry() |>
  group_by(id) |>
  # filter(n() >= min_n) |>
  nest() |>
  mutate(    
    fit   = map(data,
                 possibly(~ if(check_multiple_missions(.x)) {
                                  lm(oxygen ~ elevation * mission, data = .x)
                              } else {
                                  lm(oxygen ~ elevation, data = .x)
                              },
                          otherwise = NULL)),

    summ = map(fit, \(m)
                   if (is.null(m)) NULL
                   else summary(m)),
                   
    r2_lin_int = map_dbl(summ, \(s)
            if (is.null(s)) NA_real_
            else {
            this_summ = summ[[1]]
            this_summ$r.sq
            }),
    adj_r2_lin_int = map_dbl(summ, \(s)
            if (is.null(s)) NA_real_
            else {
            this_summ = summ[[1]]
            this_summ$adj.r.sq
            }),
  ) |>
  select(-data, -fit, -summ) |>
  ungroup()
toc()


# #* Effect on adjusted R2 of adding a mission main effect. Only plot for bins with > 1 mission
# data_lin_main_plus_int = left_join(r2_info_lin_main, r2_info_lin_int, by = "id") |>
#                 filter(n_missions > 1)

# with(data_lin_main_plus_int, plot(adj_r2_lin_main, adj_r2_lin_int))
# abline(a = 0, b = 1, col = "red")






# ------------------------------ Spline - Basic ------------------------------ #



tic()
r2_info_spl_bas <- data_large_bin_tags |>
  st_drop_geometry() |>
  filter(id %in% ids_large_n) |>
  group_by(id) |>
  nest() |>
  mutate(
    
    fit   = map(data,
                 possibly(~ gam(oxygen ~ s(elevation), data = .x),
                          otherwise = NULL)),

    summ = map(fit, \(m)
                   if (is.null(m)) NULL
                   else summary(m)),
    
    edf_spl_bas = map_dbl(summ, \(s)
                   if (is.null(s)) NA_real_
                   else {
                     # First smooth term's EDF (oxygen ~ s(elevation))
                     this_summ = summ[[1]]
                     this_summ$s.table[1, "edf"]
                   }),
                   
    adj_r2_spl_bas = map_dbl(summ, \(s)
            if (is.null(s)) NA_real_
            else {
            this_summ = summ[[1]]
            this_summ$r.sq
            })
  ) |>
  select(-data, -fit, -summ) |>
  ungroup()
toc()



# ----------------------------- Spline - Mission ----------------------------- #


tic()
r2_info_spl_main <- data_large_bin_tags |>
  st_drop_geometry() |>
  filter(id %in% ids_large_n) |>
  group_by(id) |>
  nest() |>
  mutate(
    
    fit   = map(data,
                 possibly(~ if(check_multiple_missions(.x)) {
                                  gam(oxygen ~ s(elevation) + mission, data = .x)
                              } else {
                                  gam(oxygen ~ s(elevation), data = .x)
                              },
                          otherwise = NULL)),

    summ = map(fit, \(m)
                   if (is.null(m)) NULL
                   else summary(m)),
    
    edf_spl_main = map_dbl(summ, \(s)
                   if (is.null(s)) NA_real_
                   else {
                     # First smooth term's EDF (oxygen ~ s(elevation))
                     this_summ = summ[[1]]
                     this_summ$s.table[1, "edf"]
                   }),
                   
    adj_r2_spl_main = map_dbl(summ, \(s)
            if (is.null(s)) NA_real_
            else {
            this_summ = summ[[1]]
            this_summ$r.sq
            })
  ) |>
select(-data, -fit, -summ) |>
ungroup()
toc()




# --------------------------- Spline - Interaction --------------------------- #


tic()
r2_info_spl_int <- data_large_bin_tags |>
  st_drop_geometry() |>
  filter(id %in% ids_large_n) |>
  mutate(mission = as.factor(mission)) |>
  group_by(id) |>
  nest() |>
  mutate(
    
    fit   = map(data,
                 possibly(~ if(check_multiple_missions(.x)) {
                                  gam(oxygen ~ mission + s(elevation, by = mission), data = .x)
                              } else {
                                  gam(oxygen ~ s(elevation), data = .x)
                              },
                          otherwise = NULL)),

    summ = map(fit, \(m)
                   if (is.null(m)) NULL
                   else summary(m)),

    all_edfs_spl_int = map(summ, \(s)
                   if (is.null(s)) NULL
                   else {
                     # All smoothing terms' EDFs
                     this_summ = summ[[1]]
                     this_summ$s.table[, "edf"]
                   }),
    
    total_edf_spl_int = map_dbl(all_edfs_spl_int, \(DFs)
                   if (is.null(DFs)) NA_real_
                   else {
                     # First smooth term's EDF (oxygen ~ s(elevation))
                     sum(DFs)
                   }),
    max_edf_spl_int = map_dbl(all_edfs_spl_int, \(DFs)
                   if (is.null(DFs)) NA_real_
                   else {
                     # First smooth term's EDF (oxygen ~ s(elevation))
                     max(DFs)
                   }),

    mean_edf_spl_int = map_dbl(all_edfs_spl_int, \(DFs)
                    if (is.null(DFs)) NA_real_
                    else {
                        # First smooth term's EDF (oxygen ~ s(elevation))
                        mean(DFs)
                    }),
                   
    adj_r2_spl_int = map_dbl(summ, \(s)
            if (is.null(s)) NA_real_
            else {
            this_summ = summ[[1]]
            this_summ$r.sq
            })
  ) |>
  select(-data, -fit, -summ, -all_edfs_spl_int) |>
  ungroup()
toc()



# ---------------------------------------------------------------------------- #
#                        Compile findings across models                        #
# ---------------------------------------------------------------------------- #

# ---------------------------------- Linear ---------------------------------- #

r2_info_lin = r2_info_lin_bas |>
  left_join(r2_info_lin_main, by = "id") |>
  left_join(r2_info_lin_int, by = "id")

adj_r2_lin_tidy = pivot_longer(r2_info_lin,
                               cols = c("adj_r2_lin_bas", "adj_r2_lin_main", "adj_r2_lin_int"),
                               names_to = "model", values_to = "adj_r2", names_prefix = "adj_r2_lin_") |>
                               select(-n_missions, -starts_with("r2_") )

r2_lin_tidy = pivot_longer(r2_info_lin,
                          cols = c("r2_lin_bas", "r2_lin_main", "r2_lin_int"),
                          names_to = "model", values_to = "r2", names_prefix = "r2_lin_") |>
                          select(-n_missions, -starts_with("adj_r2_"))

data_lin_tidy = left_join(adj_r2_lin_tidy, r2_lin_tidy, by = c("id", "model")) |>
                        mutate(model = case_match(model,  "bas" ~ "Basic", "main" ~ "Main_Effect", "int" ~ "Interaction"), Type = "Linear")



# ---------------------------------- Spline ---------------------------------- #

r2_info_spl = r2_info_spl_bas |>
  left_join(r2_info_spl_main, by = "id") |>
  left_join(r2_info_spl_int, by = "id")

adj_r2_spl_tidy = pivot_longer(r2_info_spl,
                               cols = c("adj_r2_spl_bas", "adj_r2_spl_main", "adj_r2_spl_int"),
                               names_to = "model", values_to = "adj_r2", names_prefix = "adj_r2_spl_") |>
                               select(-contains("edf") )

edf_spl_bas_main_tidy = pivot_longer(r2_info_spl,
                           cols = c("edf_spl_bas", "edf_spl_main"),
                           names_to = "model", values_to = "edf", names_prefix = "edf_spl_") |>
                           select(-contains("adj_r2"), -contains("int"))


edf_spl_int_total_tidy = pivot_longer(r2_info_spl,
                               cols = c("total_edf_spl_int"),
                               names_to = "model", values_to = "edf", names_prefix = "total_edf_spl_") |>
                               select(-contains("adj_r2"), -contains("spl")) 

edf_spl_int_mean_tidy = pivot_longer(r2_info_spl,
                               cols = c("mean_edf_spl_int"),
                               names_to = "model", values_to = "edf", names_prefix = "mean_edf_spl_") |>
                               select(-contains("adj_r2"), -contains("spl"))

edf_spl_int_max_tidy = pivot_longer(r2_info_spl,
                               cols = c("max_edf_spl_int"),
                               names_to = "model", values_to = "edf", names_prefix = "max_edf_spl_") |>
                               select(-contains("adj_r2"), -contains("spl"))

# edf_spl_tidy = bind_rows(edf_spl_bas_main_tidy, edf_spl_int_total_tidy)
edf_spl_tidy = bind_rows(edf_spl_bas_main_tidy, edf_spl_int_max_tidy)

data_spl_tidy = left_join(adj_r2_spl_tidy, edf_spl_tidy, by = c("id", "model")) |>
  mutate(model = case_match(model,  "bas" ~ "Basic", "main" ~ "Main_Effect", "int" ~ "Interaction"), Type = "Spline")



# ----------------------------------- Both ----------------------------------- #

data_tidy = bind_rows(data_lin_tidy, data_spl_tidy)



# ---------------------------------------------------------------------------- #
#                          Simple ANOVA of R2 and EDF                          #
# ---------------------------------------------------------------------------- #

#! Results in this section will change after re-computing everything on the full dataset

# --------------------------- (Adjusted) R-Squared --------------------------- #

#? No significant effect of interaction
fit_results_adj_r2 = lm(adj_r2 ~ Type + model, data = data_tidy)
summary(fit_results_adj_r2)

#? Both alternatives are much better than baseline (p < 1e-4). Interaction vs ME is less different (p ~ 0.0071).
ANOVA_adj_r2 = emmeans(fit_results_adj_r2, ~ model)
contrast(ANOVA_adj_r2, method = "pairwise")


# ------------------------------------ EDF ----------------------------------- #

data_tidy_edf = data_tidy |>
    filter(!is.na(edf)) |>
    select(-contains("r2"))

fit_results_edf = lm(edf ~ model, data = data_tidy_edf)
summary(fit_results_edf)

#? Adding a main effect doesn't significantly change the EDF. Adding an interaction does significantly increase the maximum within-bin EDF, relative to both basic and main effect only.
ANOVA_edf = emmeans(fit_results_edf, ~ model)
contrast(ANOVA_edf, method = "pairwise")








# ---------------------------------------------------------------------------- #
#                                  Make Plots                                  #
# ---------------------------------------------------------------------------- #


# #! For making preliminary plot
# #! Remove this

# all_R2s = r2_info_lin_bas$r2_lin_bas
# low_high_R2s = quantile(all_R2s, c(0.1, 0.9), type = 1)
# low_R2 = low_high_R2s[1]
# high_R2 = low_high_R2s[2]

# ind_low_R2 = which(all_R2s == low_R2)
# ind_high_R2 = which(all_R2s == high_R2)

# id_low_R2 = pull(r2_info_lin_bas, id)[ind_low_R2]
# id_high_R2 = pull(r2_info_lin_bas, id)[ind_high_R2]

# bin_low_R2 = filter(data_large_bin_tags, id == id_low_R2)
# bin_high_R2 = filter(data_large_bin_tags, id == id_high_R2)



# pdf(file = "Figures/low_linear_R2.pdf", width = 15, height = 9.5)
# ggplot(data = bin_low_R2, aes(x = oxygen, y = elevation)) + geom_point(aes(color = mission, shape = mission),size = 1) + 
# ggtitle(paste0("R-Squared = ", signif(low_R2, digits = 3))) + 
# theme(plot.title = element_text(hjust = 0.5, size = 40), axis.title = element_text(size = 30), legend.position = "none") +
# geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 1, formula = y ~ x) +
# geom_smooth(aes(color = mission), method = "gam", se = FALSE, linewidth = 2, formula = y ~ s(x))
# dev.off()


# pdf(file = "Figures/high_linear_R2.pdf", width = 15, height = 9.5)
# ggplot(data = bin_high_R2, aes(x = oxygen, y = elevation)) + geom_point(aes(color = mission, shape = mission),size = 1) + 
# ggtitle(paste0("R-Squared = ", signif(high_R2, digits = 3))) + 
# theme(plot.title = element_text(hjust = 0.5, size = 40), axis.title = element_text(size = 30), legend.position = "none") +
# geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 1, formula = y ~ x) +
# geom_smooth(aes(color = mission), method = "gam", se = FALSE, linewidth = 2, formula = y ~ s(x))
# dev.off()

# ---------------------- Low and High complexity regions --------------------- #

data_lin_bas = filter(data_lin_tidy, model == "Basic")
all_r2s_lin_bas = pull(data_lin_bas, r2)

low_high_r2_lin_bas = quantile(all_r2s_lin_bas, c(0.1, 0.9), type = 1)
low_r2_lin_bas = low_high_r2_lin_bas[1]
high_r2_lin_bas = low_high_r2_lin_bas[2]

ind_low_r2_lin_bas = which(all_r2s_lin_bas == low_r2_lin_bas)
ind_high_r2_lin_bas = which(all_r2s_lin_bas == high_r2_lin_bas)

id_low_r2_lin_bas = pull(data_lin_bas, id)[ind_low_r2_lin_bas]
id_high_r2_lin_bas = pull(data_lin_bas, id)[ind_high_r2_lin_bas]


bin_low_r2_lin_bas = filter(data_large_bin_tags, id == id_low_r2_lin_bas)
bin_high_r2_lin_bas = filter(data_large_bin_tags, id == id_high_r2_lin_bas)


pdf(file = "Figures/low_linear_R2.pdf", width = 15, height = 9.5)
ggplot(data = bin_low_r2_lin_bas, aes(x = oxygen, y = elevation)) + geom_point(aes(color = mission, shape = mission),size = 1) + 
ggtitle(paste0("R-Squared = ", signif(low_r2_lin_bas, digits = 3))) + 
theme(plot.title = element_text(hjust = 0.5, size = 40), axis.title = element_text(size = 30), legend.position = "none") +
geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 1, formula = y ~ x) +
geom_smooth(aes(color = mission), method = "gam", se = FALSE, linewidth = 2, formula = y ~ s(x))
dev.off()

pdf(file = "Figures/high_linear_R2.pdf", width = 15, height = 9.5)
ggplot(data = bin_high_r2_lin_bas, aes(x = oxygen, y = elevation)) + geom_point(aes(color = mission, shape = mission),size = 1) + 
ggtitle(paste0("R-Squared = ", signif(high_r2_lin_bas, digits = 3))) + 
theme(plot.title = element_text(hjust = 0.5, size = 40), axis.title = element_text(size = 30), legend.position = "none") +
geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 1, formula = y ~ x) +
geom_smooth(aes(color = mission), method = "gam", se = FALSE, linewidth = 2, formula = y ~ s(x))
dev.off()



# ---------------------------- Degrees of Freedom ---------------------------- #


data_spl_bas_edf = data_tidy |> 
                    filter(!is.na(edf), model == "Basic") |>
                    select(id, edf) 

data_plot_edf = data_bins |>
                    filter(id %in% ids_large_n) |>
                    left_join(data_spl_bas_edf, by = "id")


pdf("Figures/EDF Map.pdf", width = 14.3, height = 9.5)
ggplot(data = data_plot_edf) +
  annotation_map_tile(zoomin = 0, cachedir = "map_tiles_ggspatial") +
  geom_sf(aes(fill = edf), alpha = 0.5) +
  scale_fill_viridis_c() +
  theme_void() +  labs(fill = "DF") +
  theme(legend.title = element_text(size = 30), legend.text = element_text(size = 20))
dev.off()


# ---------------------------- Adjusted R-Squared ---------------------------- #

data_plot_r2 = full_join(filter(data_bins, id %in% ids_large_n), data_tidy, by = "id") %>%
                    select(id, model, Type, adj_r2, geometry) %>%
                    filter(model != "Main_Effect") %>%
                    mutate(model = recode_factor(model,
                                                 "Basic" = "Elevation Only",
                                                 "Interaction" = "Interaction with Mission"))


pdf("Figures/R2 Map.pdf", width = 14.4, height = 9.7)
ggplot(data = data_plot_r2) +
  annotation_map_tile(zoomin = 0, cachedir = "map_tiles_ggspatial") +
  facet_grid(rows = vars(Type), cols = vars(model), switch = "y") +
  geom_sf(aes(fill = adj_r2), alpha = 0.5) +
  scale_fill_viridis_c() +
  theme_void() +  labs(fill = "DF") +
  theme(legend.title = element_text(size = 20), legend.text = element_text(size = 15), 
        strip.text = element_text(size = 15), strip.text.y = element_text(angle = 90, vjust = 0.5))
dev.off()

























# ---------------------------------------------------------------------------- #
#                           Extract and plot bin data                          #
# ---------------------------------------------------------------------------- #

# plot_id = id_min_dist
plot_id = 1716
# plot_id = 1

data_bins |> filter(mission == "Multiple") |> pull(id)

data_this_bin = filter(data_large_bin_tags, id == plot_id)
# data_high_R2 = data_this_bin
# data_low_R2 = data_this_bin
other_data_low_R2 = data_this_bin

fit = lm(oxygen ~ elevation + mission, data = data_this_bin)

fit_spl = gam(oxygen ~ s(elevation) + mission, data = data_this_bin)
q = summary(fit_spl)

this_data = st_drop_geometry(data_this_bin) |> mutate(mission = as.factor(mission))

fit_spl_int = gam(oxygen ~ mission + s(elevation, by = mission), data = this_data)
q = summary(fit_spl_int)

sum(q$s.table[, "edf"])


ggplot(data = data_this_bin, aes(x = oxygen, y = elevation, colour = mission)) + geom_point()


this_fit = lm(oxygen ~ elevation, data = data_this_bin)
summary(this_fit)
plot(this_fit)


test = gamSim(4)
fit_test = gam(y ~ fac + s(x2, by = fac), data = test)


# ---------------------------------------------------------------------------- #
#                      Add fixed effect for mission to lm                      #
# ---------------------------------------------------------------------------- #


# ------------------------ Bins with multiple missions ----------------------- #

ids_mult_missions = data_bins |>
                        filter(mission == "Multiple") |>
                        pull(id)



min_n <- 20                 # need at least k+1 rows (k default = 10 → 11 pts)
default_k <- 10             # basis dimension for s()

ids_large_n = data_large_bin_tags |>
                group_by(id) |>
                summarise(n = n()) |>
                filter(n >= min_n) |>
                pull(id)

# some_ids = ids_mult_missions[11:20]


# -------------------------------- Fit models -------------------------------- #

check_multiple_missions = function(data) length(unique(data$mission)) > 1

tic()
lm_tbl <- data_large_bin_tags |>
#   filter(id %in% some_ids) |>
  # filter(id %in% ids_mult_missions) |>
  filter(id %in% ids_large_n) |>
  st_drop_geometry() |>
  group_by(id) |>
  # filter(n() >= min_n) |>
  nest() |>
  mutate(
    
    fit   = map(data,
                 possibly(~ if(check_multiple_missions(.x)) {
                                  lm(oxygen ~ elevation + mission, data = .x)
                              } else {
                                  lm(oxygen ~ elevation, data = .x)
                              },
                          otherwise = NULL)),

    summ = map(fit, \(m)
                   if (is.null(m)) NULL
                   else summary(m)),
                   
    r2_mission = map_dbl(summ, \(s)
            if (is.null(s)) NA_real_
            else {
            this_summ = summ[[1]]
            this_summ$r.sq
            }),
    adj_r2_mission = map_dbl(summ, \(s)
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
                    filter(id %in% ids_large_n) |>
                    left_join(bins_R2, by = "id") |>
                    # filter(id %in% ids_mult_missions) |>
                    left_join(select(lm_tbl, id, r2_mission), by = "id") |>
                    rename(simple_R2 = r2, mission_R2 = r2_mission) |>
                    mutate(R2_ratio = mission_R2 / simple_R2, log_R2_ratio = log10(R2_ratio))

simple_R2s = data_bins_mult |> select(id, simple_R2) |> st_drop_geometry()
mission_R2s = data_bins_mult |> select(id, mission_R2) |> st_drop_geometry()

#Extract lower 10th percentile of R2s
low_high_R2s =quantile(simple_R2s$simple_R2, probs = c(0.1, 0.9), type=1)
low_high_mission_R2s =quantile(mission_R2s$mission_R2, probs = c(0.1, 0.9), type=1)


id_low_R2 = simple_R2s$id[which(simple_R2s$simple_R2 == low_high_R2s[1])]
id_high_R2 = simple_R2s$id[which(simple_R2s$simple_R2 == low_high_R2s[2])]
id_low_mission_R2 = mission_R2s$id[which(mission_R2s$mission_R2 == low_high_mission_R2s[1])]
id_high_mission_R2 = mission_R2s$id[which(mission_R2s$mission_R2 == low_high_mission_R2s[2])]

data_low_R2 = data_large_bin_tags |>
                filter(id %in% id_low_R2)
data_high_R2 = data_large_bin_tags |>
                filter(id %in% id_high_R2)
data_low_mission_R2 = data_large_bin_tags |>
                filter(id %in% id_low_mission_R2)
data_high_mission_R2 = data_large_bin_tags |>
                filter(id %in% id_high_mission_R2)


#! These plots go in the poster
#! Add lines and splines?
ggplot(data = data_low_R2, aes(x = oxygen, y = elevation, colour = mission)) + geom_point() +
    ggtitle("Low R2: O2 ~ elevation + mission")
ggplot(data = data_high_R2, aes(x = oxygen, y = elevation, colour = mission)) + geom_point() +
    ggtitle("High R2: O2 ~ elevation + mission")

ggplot(data = data_low_mission_R2, aes(x = oxygen, y = elevation, colour = mission)) + geom_point() +
    ggtitle("Low R2: O2 ~ elevation + mission")
ggplot(data = data_high_mission_R2, aes(x = oxygen, y = elevation, colour = mission)) + geom_point() +
    ggtitle("High R2: O2 ~ elevation + mission")


# Build plots

## Simple R2 (O2 ~ elevation)
ggplot(data = data_bins_mult) +
  annotation_map_tile(zoomin = 0, cachedir = "map_tiles_ggspatial") +
  geom_sf(aes(fill = simple_R2), alpha = 0.5) +
  scale_fill_viridis_c() + ggtitle("Simple R2: O2 ~ elevation") +
  theme_minimal()

## Mission-adjusted R2 (O2 ~ elevation + mission)
#! This goes in the poster
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

ggplot(data = this_data_jump, mapping = aes(x = oxygen, y = elevation, colour = time, shape = mission)) + geom_point() 
ggplot(data = this_data_jump, mapping = aes(x = oxygen, y = elevation, colour = mission)) + geom_point() 

str(data_large_bin_tags)


# ---------------------------------------------------------------------------- #
#                                Explore Splines                               #
# ---------------------------------------------------------------------------- #



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



min_n <- 20                 # need at least k+1 rows (k default = 10 → 11 pts)
default_k <- 10             # basis dimension for s()


tic()
gam_tbl <- data_large_bin_tags |>
  st_drop_geometry() |>
  group_by(id) |>
  filter(n() >= min_n) |>
  nest() |>
  mutate(
    
    fit   = map(data,
                 possibly(~ if(check_multiple_missions(.x)) {
                                  gam(oxygen ~ s(elevation) + mission, data = .x)
                              } else {
                                  gam(oxygen ~ s(elevation), data = .x)
                              },
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
data_bin_DF = right_join(data_bins, select(gam_tbl, id, edf_smooth), by = "id") |> rename(DF = edf_smooth)


Spline_DF_plot = ggplot(data = data_bin_DF) +
  annotation_map_tile(zoomin = 0, cachedir = "map_tiles_ggspatial") +
  geom_sf(aes(fill = DF), alpha = 0.5) +
  scale_fill_viridis_c() +
  theme_minimal()
Spline_DF_plot


# ------------------------------ Plot spline R2 ------------------------------ #

data_bin_spline_R2 = right_join(data_bins, select(gam_tbl, id, r2), by = "id") |> rename(r2_spline = r2)


Spline_R2_plot = ggplot(data = data_bin_spline_R2) +
  annotation_map_tile(zoomin = 0, cachedir = "map_tiles_ggspatial") +
  geom_sf(aes(fill = r2_spline), alpha = 0.5) +
  scale_fill_viridis_c() +
  theme_minimal()
Spline_R2_plot


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




this_TS = ts(this_O2)
decomp = stl(this_TS)



# Periodogram
periodogram_1 = cpgram(this_TS)

periodogram_2 = spec.pgram(this_O2)
