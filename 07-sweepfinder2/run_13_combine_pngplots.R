# ================================================
# Combine PicMin plot with the SF2 plots
# ================================================
rm(list = ls())
library(magick)

setwd("C:/Users/u0167113/Documents/PhD research/popgen analyses/07. SweepFinder2/polarised")

IN_DIR_PIC  <- "picmin/plots"
IN_DIR_SF2  <- "plots"
OUT_DIR     <- "plots_picmin_sf2"
dir.create(OUT_DIR, showWarnings = FALSE)
pic1     <- "pops_C1C2C3C4C5_above3lines_allscafs_plot.png"
pic2     <- "pops_agriculture_above7lines_allscafs_plot.png"
pic3     <- "pops_all_theta2_above18lines_allscafs_plot.png"
outfile     <- "picmin_pops_C1C2C3C4C5_agri_all.png"

populations <- c(#"A1_AnOudin", "A2_BuSN", "A3_Mech", "A4_PL15_YEL", "A5_GenM",
    #"B1_BKN1", "B2_OM2", "B3_ZW", "B4_OHZ", "B5_DA2",
    "C1_BlfN", "C2_MO", "C3_Ter1", "C4_BW_48630", "C5_BW_36962"
    #,
    #"D1_CBOO6", "D2_LRV", "D3_BKLE5", "D4_BW_62256", "D5_BW_22050"
  )

# Read picmin plot
picmin_plot1 <- image_read(file.path(IN_DIR_PIC, pic1))
picmin_plot2 <- image_read(file.path(IN_DIR_PIC, pic2))
picmin_plot3 <- image_read(file.path(IN_DIR_PIC, pic3))

# Read SF2 plots
sf2_list <- list()

for (pop in populations) {
  
  sf2_path <- file.path(IN_DIR_SF2, pop, paste0(pop, "_allscaffs.png"))
  
  if (!file.exists(sf2_path)) {
    warning("SF2 image missing for population: ", pop, ". Skipping to next population.")
    next
  }
  
  sf2_list[[pop]] <- image_read(sf2_path)
}

# Combine all images: PicMin (first) + all SF2 plots

all_images <- do.call(c, c(#unname(sf2_list), 
                           list(picmin_plot1,
                                picmin_plot2,
                                picmin_plot3)))
combined <- image_append(all_images, stack = TRUE)

# Save 
output_file <- file.path(OUT_DIR, outfile)
image_write(combined, output_file)
  
message("Saved combined image: ", output_file)
  
