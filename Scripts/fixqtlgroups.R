library(data.table)
library(ggplot2)

dt<-fread("peaks_10cMwindow_allinfo_nov30.txt")

#subsetting out rows where the transcript shows up more than once 
dt_subset<-dt[,.N,by = .(name)][N >1, .(name)]
dt_fix<-dt[name %in% dt_subset$name]
#sorting by transcript name
dt_fix_sorted<-dt_fix[order(name)]
#splitting data table into a list of 
list_by_transcript<- split(dt_fix_sorted, by = "name", keep.by = TRUE)
test<-as.data.table(list_by_transcript[1])





#function to fix trt group overlaps
process_data_table <- function(dt) {
  setDT(dt)
  
  # Create newtrtgroup column for grouping by markername
  dt[, newtrtgroup := .(list(trt)), by = .(markername)]
  
  # Process rows with additional criteria
  dt[, {
    # Skip processing for groups with fewer than 2 rows
    if (.N < 2) {
      list(newtrtgroup)
    } else {
      # Find pairs of rows meeting the criteria
      pairs <- combn(seq_len(.N), 2, function(idx) {
        row1 <- dt[idx[1], ]
        row2 <- dt[idx[2], ]
        row1$type == row2$type && 
          row1$chr == row2$chr && 
          row1$markername != row2$markername && 
          row1$trt != row2$trt && 
          abs(row1$pos - row2$pos) < 5
      }, simplify = FALSE)
      
      # Ensure pairs is a matrix
      if (length(pairs) == 0) {
        list(newtrtgroup)
      } else {
        pairs <- matrix(unlist(pairs), ncol = 2, byrow = TRUE)
        
        # Append trt values to newtrtgroup for pairs
        for (pair in split(pairs, col(pairs))) {
          combined_trt <- unique(c(dt[pair[1], trt], dt[pair[2], trt]))
          dt[pair, newtrtgroup := .(list(combined_trt))]
        }
        list(newtrtgroup)
      }
    }
  }, by = .(type, chr)]
  
  dt
}



processed_list_of_tables <- lapply(list_by_transcript[1:3], process_data_table)

