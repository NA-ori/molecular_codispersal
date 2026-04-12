### Stitch output files together into one CSV ###

file_location <- file.path("untracked", "correct_results", "*.txt")
filename <- file.path("untracked", "correct_results.csv")


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
