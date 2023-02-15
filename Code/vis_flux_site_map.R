#************************************************************************************
# Description: The flux site distribution map
# Author: Ian R. McGregor. Edit by Xiaojie(J) Gao
# Date: 2022-03-03
#************************************************************************************
rm(list = ls())

source("Code/base.R")

library(ggplot2)
library(raster)
library(sf)



sites <- fread(file.path(gdir, 
    "sites_fluxnet2015.csv"
))


## set pch for different IGBP, then set colors
igbp_colors <- c(
    "#51B6F5", "#218A21", "#31CD31", "#9ACD31", "#97FA97",
    "#8FBB8F", "#BB8F8F", "#F5DEB3", "#DBEB9D", "#FFD600", "#EFB766",
    "#4682B2", "#FAED73", "#FF0000", "#999355", "#F5F5DC", "#BDBDBD",
    "#000000"
)
igbp_names <- c(
    "water", "enf", "ebf", "dnf", "dbf", "mixed", "closed shrubs",
    "open shrubs", "woody savannas", "savannas", "grasslands", "perm wetlands",
    "croplands", "urban", "crop/natural mosaic", "snow and ice", "barren/sparse veg",
    "unclassified"
)
igbp_abb <- c(
    "WAT", "ENF", "EBF", "DNF", "DBF", "MF", "CSH", "OSH", "WSA", "SAV",
    "GRA", "WET", "CRO", "URB", "CVM", "SNO", "BSV", "UNC"
)

igbp <- data.table(col = igbp_colors, names = igbp_names, IGBP = igbp_abb)
igbp <- igbp[IGBP %in% sites[, IGBP], .(col, IGBP)]

setkeyv(igbp, cols = c("IGBP", "col"))
setkey(sites, cols = IGBP)
sites <- sites[igbp]
setkey(sites, NULL)

## create maps
world <- ne_countries(scale = "medium", returnclass = "sf")
wgs84 <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
rob <- "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
coordinates(sites) <- ~ lon + lat
proj4string(sites) <- CRS(wgs84)
sites <- spTransform(sites, CRS(rob))
sites <- st_as_sf(sites)


map_wrl <- ggplot(world) +
    geom_sf(data = world, fill = "#CCCCCC", show.legend = FALSE) +
    geom_sf(data = sites, aes(color = IGBP), shape = 16, size = 4) +
    geom_sf(data = sites, color = "black", shape = 21, size = 4, stroke = 0.5, 
        show.legend = TRUE) +
    scale_color_manual(values = unique(sites$col)) +
    coord_sf(crs = st_crs("ESRI:54030")) +
    theme(
        legend.position = "top", legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.title = element_blank(),
        legend.text = element_text(size = 14)
    ) +
    guides(colour = guide_legend(nrow = 1))

bbeur <- bb_poly(c(-9.77, 34.004881, 34.57187, 70.16419), projection = wgs84)
bbeur <- st_transform(bbeur, crs = rob)

map_eur <- ggplot(world) +
    geom_sf(data = world, fill = "#CCCCCC", show.legend = FALSE) +
    geom_sf(data = sites, aes(color = IGBP), shape = 16, size = 4) +
    geom_sf(data = sites, color = "black", shape = 21, size = 4, 
        show.legend = TRUE) +
    scale_color_manual(values = unique(sites$col)) +
    coord_sf(
        # xlim = c(-750000, 1580000), ylim = c(3800000, 6089951),
        xlim = c(-750000, 1580000), ylim = c(3800000, 6089951),
        crs = st_crs("ESRI:54030")
    ) +
    theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none"
    )

# print map with inset
# ggdraw() +
#   draw_plot(map_wrl) +
#   draw_plot(map_eur, x=0.01, y=0.2, width=0.25, height=0.3)

# png("Output/flux_map.png", width = 1600, height = 900)

the_fig <- ggdraw() +
    draw_plot(map_wrl) +
    draw_plot(map_eur, x = 0.81, y = 0.62, width = 0.21, height = 0.25)

# fig: flux sites map
ggsave(the_fig, filename = "Output/flux_map.png", width = 16, height = 9, dpi = 150)

# dev.off()