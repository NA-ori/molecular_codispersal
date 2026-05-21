### Stitch output files together into one CSV ###

file_location <- file.path("Data", "raw", "3-dist", "results_s1_complete", "*.txt")
filename <- file.path("Data", "3-dist", "dist_results_s1_complete.csv")


# Read raw data files
data_c <- data.frame()
files <- Sys.glob(file_location)

read_files <- function(data, files) {
    for (file in files) {
        current_data <- read.table(file, header=TRUE, sep=",")
        if (nrow(current_data) == 0) {
            print(file)
        }
        data <- rbind(data, current_data)
    }
    return(data)
}

stitched <- read_files(data_c, files)
write.csv(stitched, filename)
