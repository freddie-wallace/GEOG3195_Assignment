# Load packages
library(sf)
library(mgcv)
library(tidyverse)
library(gstat)
library(sp)
library(cols4all)
library(gridExtra)
library(ggplot2)
library(ggspatial)

# Load the data
msoa <- st_read("msoa.gpkg", quiet = TRUE)
df <- tibble(read.csv("msoa_vacc_data.csv"))

# Convert to numeric
df$vacc <- as.numeric(df$vacc)

# Data summary
head(df)

# Join
msoa |> left_join(df) -> msoa

# Calculate centroid coordinates
coords <- st_coordinates(st_centroid(msoa))

# Create a data frame with additional columns
df.gam <- msoa %>%
  mutate(
    Intercept = 1,  # add an Intercept column
    X = coords[, "X"] / 1000,
    Y = coords[, "Y"] / 1000
  ) %>%
  st_drop_geometry() %>%
  as_tibble()

# Drop df to save project space
rm(df)


###################### GGP GAM ######################

# Set k
k_choice <- 170

# Tune the model with pre specified k
gam.m <- gam(vacc ~ 0 + 
               Intercept + s(X, Y, bs='gp', by=Intercept, k=k_choice) + 
               o65 + s(X, Y, bs='gp', by=o65, k=k_choice) +
               l4qual + s(X, Y, bs='gp', by=l4qual, k=k_choice) +
               badhealth + s(X, Y, bs='gp', by=badhealth, k=k_choice) +
               unemp + s(X, Y, bs='gp', by=unemp, k=k_choice),
             data = df.gam)


# Plot the results
gam.check(
  gam.m,
  old.style = TRUE
)

# Display the summary
summary(gam.m)

#### Plotting #####

# Get the predicted values into a table with the spatial data
msoa |> st_drop_geometry() -> msoa

# Add the Intercept column
msoa$Intercept <- 1

# Add the coordinates
msoa$X <- coords[, "X"] / 1000
msoa$Y <- coords[, "Y"] / 1000

# Predict the values
msoa$pred <- predict(gam.m, newdata = msoa)

# Load geoms
msoa_geoms <- st_read("msoa.gpkg", quiet = TRUE)

# Join the predicted values to the spatial data
msoa_geoms |> left_join(msoa) -> msoa

# Create Function to add scale bar and north arrow to the map using ggspatial
add_scale_north <- function(map) {
  map + 
    annotation_scale(location = "br", width_hint = 0.25) +
    annotation_north_arrow(location = "tl", which_north = "true", style = north_arrow_fancy_orienteering)
}

# Map both the actual and predicted values using ggplot2 and cols4all palette
actual_map <- ggplot() +
  geom_sf(data = msoa, aes(fill = vacc), color = NA, size = 0.2) +
   scale_fill_continuous_c4a_seq(palette="carto.emrld",  limits = c(0, 1), name="Actual Vaccination Rate") +
  labs(title = "Actual Vaccination Rate") +
  theme_minimal()

predicted_map <- ggplot() +
  geom_sf(data = msoa, aes(fill = pred), color = NA, size = 0.2) +
  scale_fill_continuous_c4a_seq(palette="carto.emrld",  limits = c(0, 1), name="Predicted Vaccination Rate") +
  labs(title = "Predicted Vaccination Rate") +
  theme_minimal()

# Add scale bar and north arrow to the maps
actual_map <- add_scale_north(actual_map)
predicted_map <- add_scale_north(predicted_map)

# Display the maps side by side
grid.arrange(actual_map, predicted_map, ncol = 2)


# Plot the residuals
# Calculate residuals
msoa$residuals <- msoa$vacc - msoa$pred

# Set breaks and colors
breaks <- c(-0.2, 0, 0.2)
colors <- c("blue", "antiquewhite2", "red")

# Residuals Map
residuals_map <- ggplot() +
  geom_sf(data = msoa, aes(fill = residuals), color = NA, size = 0.2) +
  scale_fill_gradient2(
    limits = c(-0.2, 0.2),
    breaks = breaks,
    low = colors[1],
    mid = colors[2],
    high = colors[3],
    midpoint = 0,
    na.value = "gray",  # Specify a color for NA values
    name = "Residuals"
  ) +
  labs(title = "Residuals Map") +
  theme_minimal()

# Add scale bar and north arrow to the map
residuals_map <- add_scale_north(residuals_map)  

# Display the residuals map
print(residuals_map)

###################### GGP GAM SVC ######################

# Function
ggp_gam_svc <- function(model, terms, input_data) {
    n_t <- length(terms)
    input_data_copy <- input_data
    output_data <- input_data
    for (i in 1:n_t) {
        # Create a vector of zeros with length n_t
        zeros <- rep(0, n_t)
        # Set the i-th element to 1
        zeros[i] <- 1
        # Create a data frame with the i-th term as columns
        terms_df <- data.frame(matrix(rep(zeros, nrow(input_data)), ncol = n_t, byrow = TRUE))
        names(terms_df) <- terms
        # Replace columns in input_data_copy with terms_df
        input_data_copy[, terms] <- terms_df
        # Predict standard errors for the modified input_data_copy
        se.j <- predict(model, se = TRUE, newdata = input_data_copy)$se.fit
        # Predict for the modified input_data_copy
        b.j <- predict(model, newdata = input_data_copy)
        # Create expressions for assigning results to output_data
        expr1 <- paste0("b_", terms[i], " <- b.j")
        expr2 <- paste0("se_", terms[i], " <- se.j")
        # Mutate output_data with the new columns
        output_data <- output_data %>% 
            mutate(within(., !!parse(text = expr1))) %>% 
            mutate(within(., !!parse(text = expr2)))
  }
    # Add the predicted values to the output_data
    output_data$yhat <- predict(model, newdata = input_data)
    output_data
}

# Apply it
terms = c("Intercept", "o65", "l4qual", "badhealth", "unemp")
gam_svc = ggp_gam_svc(gam.m, terms, df.gam)
# Check the results
gam_svc |>
  select(starts_with("b_")) |> 
  sapply(summary) |> 
  t() |> round(1)


#### Plotting #####

# Join gam_svc to the spatial data
msoa_geoms |> left_join(gam_svc) -> gam_svc

# Map both b_Intercept and se_Intercept using ggplot2 and cols4all palette
b_Intercept_map <- ggplot() +
  geom_sf(data = gam_svc, aes(fill = b_Intercept, geometry = geom), color = NA, size = 0.2) +
  scale_fill_continuous_c4a_seq(palette="scico.lajolla", name="Smooth Term Estimate") +
  labs(title = "Intercept Smooth Term Estimate") +
  theme_minimal()

se_Intercept_map <- ggplot() +
  geom_sf(data = gam_svc, aes(fill = se_Intercept, geometry = geom), color = NA, size = 0.2) +
  scale_fill_continuous_c4a_seq(palette="tol.sunset", name="Standard Error Estimate") +
  labs(title = "Intercept Standard Error Estimate") +
  theme_minimal()

# Add scale bar and north arrow to the maps
b_Intercept_map <- add_scale_north(b_Intercept_map)
se_Intercept_map <- add_scale_north(se_Intercept_map)

# Display the maps side by side
grid.arrange(b_Intercept_map, se_Intercept_map, ncol = 2)

# Map both b_o65 and se_o65 using ggplot2 and cols4all palette
b_o65_map <- ggplot() +
  geom_sf(data = gam_svc, aes(fill = b_o65, geometry = geom), color = NA, size = 0.2) +
  scale_fill_continuous_c4a_seq(palette="scico.lajolla", name="Smooth Term Estimate") +
  labs(title = "Over 65 Smooth Term Estimate") +
  theme_minimal()

se_o65_map <- ggplot() +
  geom_sf(data = gam_svc, aes(fill = se_o65, geometry = geom), color = NA, size = 0.2) +
  scale_fill_continuous_c4a_seq(palette="tol.sunset", name="Standard Error Estimate") +
  labs(title = "Over 65 Standard Error Estimate") +
  theme_minimal()

# Add scale bar and north arrow to the maps
b_o65_map <- add_scale_north(b_o65_map)
se_o65_map <- add_scale_north(se_o65_map)

# Display the maps side by side
grid.arrange(b_o65_map, se_o65_map, ncol = 2)

# Map both b_l4qual and se_l4qual using ggplot2 and cols4all palette
b_l4qual_map <- ggplot() +
  geom_sf(data = gam_svc, aes(fill = b_l4qual, geometry = geom), color = NA, size = 0.2) +
  scale_fill_continuous_c4a_seq(palette="scico.lajolla", name="Smooth Term Estimate") +
  labs(title = "Level 4 Qualifications Smooth Term Estimate") +
  theme_minimal()

se_l4qual_map <- ggplot() +
  geom_sf(data = gam_svc, aes(fill = se_l4qual, geometry = geom), color = NA, size = 0.2) +
  scale_fill_continuous_c4a_seq(palette="tol.sunset", name="Standard Error Estimate") +
  labs(title = "Level 4 Qualifications Standard Error Estimate") +
  theme_minimal()

# Add scale bar and north arrow to the maps
b_l4qual_map <- add_scale_north(b_l4qual_map)
se_l4qual_map <- add_scale_north(se_l4qual_map)

# Display the maps side by side
grid.arrange(b_l4qual_map, se_l4qual_map, ncol = 2)

# Map both b_badhealth and se_badhealth using ggplot2 and cols4all palette
b_badhealth_map <- ggplot() +
  geom_sf(data = gam_svc, aes(fill = b_badhealth, geometry = geom), color = NA, size = 0.2) +
  scale_fill_continuous_c4a_seq(palette="scico.lajolla", name="Smooth Term Estimate") +
  labs(title = "Bad Health Smooth Term Estimate") +
  theme_minimal()

se_badhealth_map <- ggplot() +
  geom_sf(data = gam_svc, aes(fill = se_badhealth, geometry = geom), color = NA, size = 0.2) +
  scale_fill_continuous_c4a_seq(palette="tol.sunset", name="Standard Error Estimate") +
  labs(title = "Bad Health Standard Error Estimate") +
  theme_minimal()

# Add scale bar and north arrow to the maps
b_badhealth_map <- add_scale_north(b_badhealth_map)
se_badhealth_map <- add_scale_north(se_badhealth_map)

# Display the maps side by side
grid.arrange(b_badhealth_map, se_badhealth_map, ncol = 2)

# Map both b_unemp and se_unemp using ggplot2 and cols4all palette
b_unemp_map <- ggplot() +
  geom_sf(data = gam_svc, aes(fill = b_unemp, geometry = geom), color = NA, size = 0.2) +
  scale_fill_continuous_c4a_seq(palette="scico.lajolla", name="Smooth Term Estimate") +
  labs(title = "Unemployment Smooth Term Estimate") +
  theme_minimal()

se_unemp_map <- ggplot() +
  geom_sf(data = gam_svc, aes(fill = se_unemp, geometry = geom), color = NA, size = 0.2) +
  scale_fill_continuous_c4a_seq(palette="tol.sunset", name="Standard Error Estimate") +
  labs(title = "Unemployment Standard Error Estimate") +
  theme_minimal()

# Add scale bar and north arrow to the maps
b_unemp_map <- add_scale_north(b_unemp_map)
se_unemp_map <- add_scale_north(se_unemp_map)

# Display the maps side by side
grid.arrange(b_unemp_map, se_unemp_map, ncol = 2)
