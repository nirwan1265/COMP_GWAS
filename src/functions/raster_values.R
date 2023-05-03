# Packages required
library(sp)
library(raster)
#install.packages("rgdal")
library(rgdal)

# 

# Directory
dir_data <- "~/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/OlsenP_2023/"
dir_raster <- "~/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/OlsenP_2023/soil_info_raster/"

# Pheno file
pheno <- read.csv(paste0(dir_data,"Final_filtered_data.csv"))
pheno <- pheno[which(pheno$Continent == "Africa"), ]
colnames(pheno)



# raster files:
raster_files <- list.files(path = dir_raster, pattern = "\\.tif$", full.names = TRUE)

# Function - extracting values from a raster file for all coordinates 
extract_raster_values <- function(raster_file, df) {
  r <- raster::brick(raster_file)
  coords <- df[, c("Long", "Lat")]
  coordinates(coords) <- c("Long", "Lat")
  proj4string(coords) <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
  coords_transformed <- spTransform(coords, crs(r))
  values <- raster::extract(x = r, y = coords_transformed)
  return(as.vector(values))
}

# Looping through the files
for (raster_file in raster_files) {
  # Get the variable name from the file name
  var_name <- sub("_0-.*", "", basename(raster_file))
  
  # Extract values for all coordinates in the pheno data frame
  extracted_values <- extract_raster_values(raster_file, pheno)
  
  # Add the extracted values as a new column in the pheno data frame
  pheno[[var_name]] <- extracted_values
}

write.csv(pheno,"OlsenP_wsoilvalues.csv",row.names = F)
