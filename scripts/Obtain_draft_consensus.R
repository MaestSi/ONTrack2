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

Obtain_draft_consensus <- function(fastq_file, TRC, PLUR, num_threads, fast_alignment_flag) {
  TRC <- as.numeric(TRC)
  target_reads_consensus <- TRC
  PLUR <- as.numeric(PLUR)
  num_threads <- as.numeric(num_threads)
  fast_alignment_flag <- as.numeric(fast_alignment_flag)
  fasta_file <- gsub(pattern = "\\.fastq", replacement = ".fasta", x = fastq_file)
  sample_name <- gsub(pattern = "_decont", replacement = "", x = gsub(pattern = "\\.fasta", replacement = "", x = basename(fasta_file)))
  sample_dir <- gsub(pattern = "_decont", replacement = "", x = gsub(pattern = "readsClustering", replacement = paste0("draftConsensusCalling/", sample_name), x = dirname(dirname(fastq_file))))
  draft_consensus_tmp1 <- paste0(sample_dir, "/", sample_name, "_draft_consensus_tmp1.fasta")
  draft_consensus_tmp2 <- paste0(sample_dir, "/", sample_name, "_draft_consensus_tmp2.fasta")
  draft_consensus <- paste0(sample_dir, "/", sample_name, "_draft_consensus.fasta")
  dir.create(sample_dir)
  logfile <- paste0(dirname(dirname(sample_dir)), "/logfile.txt")
  num_reads_mac <- as.double(system(paste0("cat ", fasta_file, " | grep \"^>\" | wc -l"), intern=TRUE))
  num_reads_sample <- as.double(system(paste0("cat ", paste0(dirname(dirname(dirname(fastq_file))), "/readsFiltering/", sample_name, ".fasta"), " | grep \"^>\" | wc -l"), intern=TRUE))
  mac_sample_ratio_perc <- 100*num_reads_mac/num_reads_sample
  #create consensus sequence
  cat(text = paste0("Sample ", sample_name, ": ", num_reads_mac, " reads (", sprintf("%.2f", mac_sample_ratio_perc), "%) assigned to most abundant cluster"), sep = "\n")
  cat(text = paste0("Sample ", sample_name, ": ", num_reads_mac, " reads (", sprintf("%.2f", mac_sample_ratio_perc), "%) assigned to most abundant cluster"),  file = logfile, sep = "\n", append = TRUE)
  #if not at least 3 reads are assigned to the allele, consensus calling is skipped
  if (num_reads_mac < 3) {
    cat(text = paste0("WARNING: Only ", num_reads_mac, " reads available for sample ", sample_name), sep = "\n")
    cat(text = paste0("WARNING: Only ", num_reads_mac, " reads available for sample ", sample_name),  file = logfile, sep = "\n", append = TRUE)
    system(paste0("head -n2 ", fasta_file, " > ", draft_consensus))
  } else {
    if (num_reads_mac < target_reads_consensus) {
      target_reads_consensus <- num_reads_mac
      cat(text = paste0("WARNING: Only ", num_reads_mac, " reads available for sample ", sample_name), sep = "\n")
      cat(text = paste0("WARNING: Only ", num_reads_mac, " reads available for sample ", sample_name),  file = logfile, sep = "\n", append = TRUE)
    }
    plurality_value <- PLUR*target_reads_consensus
    sequences <- readDNAStringSet(fasta_file, "fasta")
    ws <- width(sequences)
    amplicon_length <- ceiling(mean(ws))
    draft_reads_fq <- paste0(sample_dir, "/", sample_name, "_draft_", target_reads_consensus, "_reads.fastq")
    draft_reads_fa <- paste0(sample_dir, "/", sample_name, "_draft_", target_reads_consensus, "_reads.fasta")
    seed <- 1
    system(paste0("/opt/conda/envs/ONTrack2_env/bin/seqtk sample -s ", seed , " ", fastq_file, " ",  target_reads_consensus, " > ", draft_reads_fq))
    system(paste0("/opt/conda/envs/ONTrack2_env/bin/seqtk seq -A ", draft_reads_fq, " > ", draft_reads_fa))
    mfa_file <- gsub(pattern = "\\.fasta$", replacement = ".mfa", x = draft_reads_fa)
    if (fast_alignment_flag == 1) {
      system(paste0("/opt/conda/envs/ONTrack2_env/bin/mafft --auto --thread ", num_threads, " --adjustdirectionaccurately ", draft_reads_fa, " > ", mfa_file))
    } else {
      system(paste0("/opt/conda/envs/ONTrack2_env/bin/mafft -linsi --thread ", num_threads, " --adjustdirectionaccurately ", draft_reads_fa, " > ", mfa_file))
    }
    system(paste0("/opt/conda/envs/ONTrack2_env/bin/cons -sequence ", mfa_file, " -plurality ", plurality_value, " -outseq ", draft_consensus_tmp1))
    system(paste0("sed 's/[nN]//g' ", draft_consensus_tmp1, " > ", draft_consensus_tmp2))
    DNAStringSet_obj <- readDNAStringSet(draft_consensus_tmp2, "fasta")
    DNAStringSet_obj_renamed <- DNAStringSet_obj
    original_headers <- names(DNAStringSet_obj)
    sequences <- seq(DNAStringSet_obj)
    names(DNAStringSet_obj_renamed) <- paste0("Consensus_sequence_", sample_name)
    writeXStringSet(x = DNAStringSet_obj_renamed, filepath = draft_consensus, format = "fasta", width = 20000)
  }
}

Obtain_draft_consensus(fastq_file, TRC, PLUR, num_threads, fast_alignment_flag)
