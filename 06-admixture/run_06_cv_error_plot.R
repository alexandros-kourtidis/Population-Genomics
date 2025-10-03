# Load the data
cv_data <- read.csv("cv_error.csv")

# Extract numeric values from the 'K' column
cv_data$K <- as.numeric(gsub("K=", "", cv_data$K))

# Find the index of the minimum error
min_index <- which.min(cv_data$Error)

# Create the plot
plot(cv_data$K, cv_data$Error, type = "b", pch = 19, col = "blue",
     xlab = "Number of Clusters (K)", ylab = "Cross-Validation Error",
     main = "Cross-Validation Error vs. Number of Clusters")

# Highlight the point with the lowest error in red
points(cv_data$K[min_index], cv_data$Error[min_index], col = "red", pch = 19, cex = 1.5)

# Annotate the point with the corresponding K value
text(cv_data$K[min_index], cv_data$Error[min_index], 
     labels = paste("K =", cv_data$K[min_index]), pos = 3, offset = 1.3, font = 2, cex = 1.3, col = "red")
