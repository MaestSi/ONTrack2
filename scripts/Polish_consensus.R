#
# Copyright 2022 Simone Maestri. All rights reserved.
# Simone Maestri <simone.maestri@iit.it>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

args = commandArgs(trailingOnly=TRUE)

for(v in args)
{
  vTmp <- strsplit(v,"=")[[1]]
  assign(vTmp[[1]],vTmp[[2]])
}

#load BioStrings package
suppressMessages(library(Biostrings))

Polish_consensus <- function(draft_consensus, fastq_file, num_threads, TRP, primers_length, medaka_model) {
  num_threads <- as.numeric(num_threads)
  TRP <- as.numeric(TRP)
  target_reads_polishing <- TRP
  primers_length <- as.numeric(primers_length)
  fasta_file <- gsub(pattern = "\\.fastq", replacement = ".fasta", x = fastq_file)
  sample_dir <- dirname(fasta_file)
  sample_name <- gsub(pattern = "_decont.*", replacement = "", x = basename(fasta_file))
  sample_dir_output <- gsub(pattern = "readsClustering", replacement = paste0("consensusPolishing/", sample_name), x = dirname(dirname(fastq_file)))
  dir.create(sample_dir_output)
  logfile <- paste0(dirname(dirname(sample_dir)), "/logfile.txt")
  num_reads_mac <- as.double(system(paste0("cat ", fasta_file, " | grep \"^>\" | wc -l"), intern=TRUE))
  num_threads_medaka <- min(num_threads, 8)
  paf_file <- paste0(sample_dir_output, "/", gsub(pattern = "_decont\\.fastq$", replacement = ".paf", x = basename(fastq_file)))
  racon_consensus <- paste0(sample_dir_output, "/", sample_name, "_racon_consensus.fasta")
  seed <- 2
  consensus_untrimmed <- paste0(sample_dir_output, "/", sample_name, "_consensus_untrimmed.fasta")
  consensus <- paste0(sample_dir_output, "/", sample_name, "_consensus.fasta")
  if (num_reads_mac < target_reads_polishing) {
    target_reads_polishing <- num_reads_mac
  }
  polishing_reads_fq <- paste0(sample_dir_output, "/", sample_name, "_polishing_", target_reads_polishing, "_reads.fastq")
  system(paste0("/opt/conda/envs/ONTrack2_env/bin/seqtk sample -s ", seed , " ", fastq_file, " ",  target_reads_polishing, " > ", polishing_reads_fq))
  cat(text = paste0("Running Racon for consensus polishing of sample ", sample_name), sep = "\n")
  system(paste0("/opt/conda/envs/ONTrack2_env/bin/minimap2 -x map-ont ", draft_consensus, " ", polishing_reads_fq, " > ", paf_file))
  system(paste0("/opt/conda/envs/ONTrack2_env/bin/racon -t ", num_threads, " -m 8 -x -6 -g -8 -w 500 --no-trimming ", polishing_reads_fq, " ", paf_file, " ", draft_consensus, " > ", racon_consensus))
  if (length(which(readLines(racon_consensus) != "")) == 0) {
    cat(text = paste0("WARNING: Failed to run Racon for sample ", sample_name), sep = "\n")
    cat(text = paste0("WARNING: Failed to run Racon for sample ", sample_name),  file = logfile, sep = "\n", append = TRUE)
    system(paste0("cp ", draft_consensus, " ", racon_consensus))
  }
  cat(text = paste0("Running Medaka for consensus polishing of sample ", sample_name), sep = "\n")
  system(paste0("/opt/conda/envs/ONTrack2_env/bin/medaka_consensus -i ", polishing_reads_fq, " -d ", racon_consensus, " -m ", medaka_model, " -t ", num_threads_medaka, " -o ", sample_dir_output, "/medaka_polishing"))
  system(paste0("cp ", sample_dir_output, "/medaka_polishing/consensus.fasta ", consensus_untrimmed))
  system(paste0("/opt/conda/envs/ONTrack2_env/bin/seqtk trimfq ", consensus_untrimmed, " -b ", primers_length, " -e ", primers_length, " > ", consensus))
}

Polish_consensus(draft_consensus, fastq_file, num_threads, TRP, primers_length, medaka_model)
