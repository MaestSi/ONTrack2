library("tcltk")

###Function
Run_ONTrack2 = function(demultiplexed_dir, blacklist_barcodes, ONTrack2_conf_file, ONTrack2_nf_file) {
  demultiplexed_concatenated_dir <- paste0(gsub(x = demultiplexed_dir, pattern = "/$", replacement = ""), "_concatenated")
  dir.create(demultiplexed_concatenated_dir)
  cat("#!/bin/bash\nexport NXF_DEFAULT_DSL=1\n", file = paste0(demultiplexed_concatenated_dir, "/nextflowRun.sh"))
  cat("source activate nextflow_env\n", file = paste0(demultiplexed_concatenated_dir, "/nextflowRun.sh"), append = TRUE)
  
  barcodes_dirs <- grep(x = list.dirs(path = demultiplexed_dir, recursive = TRUE, full.names = TRUE), pattern = "barcode", value = TRUE)

  if (length(blacklist_barcodes) > 0) {
    blacklist_barcodes <- gsub(x = blacklist_barcodes, pattern = "\\s*", replacement = "")
    barcodes_dirs <- barcodes_dirs[-which(grepl(pattern = paste0(blacklist_barcodes, collapse = "|"), x = barcodes_dirs))]  
  }
  
  for (i in 1:length(barcodes_dirs)) {
    sample_curr <- basename(barcodes_dirs[i])
    cat(sprintf("Concatenating reads for sample %s\n", sample_curr))
    fastq_files <- list.files(path = barcodes_dirs[i], pattern = "\\.fastq", full.names = TRUE)
    system(paste0("zless ", paste0(fastq_files, collapse = " "), " > ", paste0(demultiplexed_concatenated_dir, "/", sample_curr, ".fastq")))
  }
  
  cat(paste0("nextflow -c ", ONTrack2_conf_file, " run ", ONTrack2_nf_file, 
             " --fastq_files=\"", demultiplexed_concatenated_dir, "/*.fastq", "\"", 
             " --results_dir=", paste0(demultiplexed_concatenated_dir, "/ONTrack2_output"),
             " -w ", paste0(demultiplexed_concatenated_dir, "/ONTrack2_output/work"),
             " -profile docker -with-report report_ONTrack2.html -with-timeline timeline_ONTrack2.html"),
      file = paste0(demultiplexed_concatenated_dir, "/nextflowRun.sh"), append = TRUE)
  
  system(paste0("chmod 755 ", demultiplexed_concatenated_dir, "/nextflowRun.sh"))
  
  cat("Starting ONTrack2 pipeline...\n")
  cat(sprintf("Results will be found in folder: %s\n\n", paste0(demultiplexed_concatenated_dir, "/ONTrack2_output")))
  system(paste0("bash ", demultiplexed_concatenated_dir, "/nextflowRun.sh"))
}
###########

cat(sprintf("Starting ONTrack2 preprocessing pipeline...\n"))
cat("Is Docker Desktop running? (y/n)\n")
ans <- readLines(con = "stdin", n = 1)
if (length(grep(pattern = "y", x = ans, ignore.case = "TRUE")) > 0) {
  demultiplexed_dir <- tk_choose.dir(caption = "Select the folder containing demultiplexed fastq files\n")
  cat("Enter a comma separated list of blacklisted samples (if any, otherwise press enter):\n")
  blacklist_barcodes <- unlist(strsplit(readLines(con = "stdin", n = 1), ",")) 
  ONTrack2_conf_file <- tk_choose.files(caption = "Select ONTrack2.conf file\n", multi = FALSE)
  ONTrack2_nf_file <- tk_choose.files(caption = "Select ONTrack2.nf file\n", multi = FALSE, default = gsub(pattern = "\\.conf", replacement = ".nf", x = ONTrack2_conf_file))
  Run_ONTrack2(demultiplexed_dir = demultiplexed_dir, blacklist_barcodes = blacklist_barcodes, ONTrack2_conf_file = ONTrack2_conf_file, ONTrack2_nf_file = ONTrack2_nf_file)
} else {
  stop("Please start Docker Desktop and then rerun the script.\nExiting.")
}