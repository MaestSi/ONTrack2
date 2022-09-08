#!/usr/bin/env nextflow
/*
========================================================================================
                         MaestSi/ONTrack2
========================================================================================
 MaestSi/ONTrack2 analysis pipeline.
 #### Homepage / Documentation
 https://github.com/MaestSi/ONTrack2
----------------------------------------------------------------------------------------
*/
def helpMessage() {
        log.info"""
    Usage:
    nextflow -c ONTrack2.conf run ONTrack2.nf --fastq_files = "/path/to/files*.fastq" --scripts_dir = "/path/to/scripts_dir" --results_dir = "/path/to/results_dir" -profile docker

    Mandatory argument:
    -profile                                                              Configuration profile to use. Available: docker, singularity
    Other mandatory arguments which may be specified in the ONTrack2.conf file
    --fastq_files                                                         Path to fastq files, use wildcards to select multiple samples
    --results_dir                                                         Path to a folder where to store results
    --scripts_dir                                                         scripts_dir is the directory containing all scripts
    --subsampling_flag                                                    subsampling_flag = true if you want to perform reads subsampling to reduce running time
    --subsampled_reads                                                    subsampled_reads is the number of subsampled reads for each sample in case subsampling_flag = true
    --minQ                                                                min Q value for reads filtering
    --minLen                                                              min read length for reads filtering
    --maxLen                                                              max read length for reads filtering
    --target_reads_consensus                                              target_reads_consensus defines the maximum number of reads used for consensus calling
    --target_reads_polishing                                              target_reads_polishing defines the maximum number of reads used for consensus polishing
    --clustering_id_threshold                                             identity threshold for clustering preliminary allele assembly
    --plurality                                                           cut-off for the number of positive matches in the multiple sequence alignment below which there is no consensus
    --fast_alignment_flag                                                 set fast_alignment_flag=1 if you want to perform fast multiple sequence alignment; otherwise set fast_alignment_flag=0
    --primers_length                                                      primers_length defines how many bases are trimmed from consensus sequences
    --medaka_model                                                        medaka model for consensus polishing
    --blast_db                                                            path to Blast-indexed database for Blasting consensus sequences
    """.stripIndent()
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

// Input of sample names, conditions, and FAST5s path.
Channel
    .fromPath(params.fastq_files)
    .map {tuple( it.name.split('\\.')[0], it )}
    .set{inputFiles_readsFiltering}

// readsFiltering
process readsFiltering {
  input:
    tuple val(sample), val(fastq) from inputFiles_readsFiltering
  
  output:
    val(sample) into readsFiltering_readsClustering
  
  script:
  if(params.readsFiltering)
  """
    mkdir -p ${params.results_dir}/readsFiltering
    export PATH=\$PATH:/opt/conda/envs/ONTrack2_env/bin/

    sample=\$(echo \$(basename \$(realpath ${fastq})) | sed \'s/\\.fa.*//\' | sed \'s/\\.fq.*//\')
    filtered_reads_fq=${params.results_dir}"/readsFiltering/"\$sample".fastq"
    filtered_reads_fa=${params.results_dir}"/readsFiltering/"\$sample".fasta"
    if ( ${params.subsampling_flag} ); then
      /opt/conda/envs/ONTrack2_env/bin/seqtk sample ${fastq} ${params.subsampled_reads} > subsampled_reads.fastq;
      cat subsampled_reads.fastq | /opt/conda/envs/ONTrack2_env/bin/NanoFilt -q ${params.minQ} --length ${params.minLen} --maxlength ${params.maxLen} > \$filtered_reads_fq;
      /opt/conda/envs/ONTrack2_env/bin/seqtk seq -A \$filtered_reads_fq > \$filtered_reads_fa;
    else
      cat ${fastq} | /opt/conda/envs/ONTrack2_env/bin/NanoFilt -q ${params.minQ} --length ${params.minLen} --maxlength ${params.maxLen} > \$filtered_reads_fq;
      /opt/conda/envs/ONTrack2_env/bin/seqtk seq -A \$filtered_reads_fq > \$filtered_reads_fa;
    fi

    ln -s \$filtered_reads_fq ./filtered_reads.fastq
    ln -s \$filtered_reads_fa ./filtered_reads.fasta
  """
  else
  """
    mkdir -p ${params.results_dir}/readsFiltering
    reads_full=\$(realpath ${fastq})
    sample=\$(echo \$(basename \$(realpath ${fastq})) | sed \'s/\\.fa.*//\' | sed \'s/\\.fq.*//\')
    filtered_reads_fq=${params.results_dir}"/readsFiltering/"\$sample".fastq"
    filtered_reads_fa=${params.results_dir}"/readsFiltering/"\$sample".fasta"
    cp ${fastq} \$filtered_reads_fq
    /opt/conda/envs/ONTrack2_env/bin/seqtk seq -A \$filtered_reads_fq > \$filtered_reads_fa
  """
}

// readsClustering
process readsClustering {
  input:
    val(sample) from readsFiltering_readsClustering
  
  output:
    val(sample) into readsClustering_draftConsensusCalling
  
  script:
  if(params.readsClustering)
  """
  mkdir -p ${params.results_dir}/readsClustering/${sample}
  clustering_dir=${params.results_dir}/readsClustering/${sample}
  
  /opt/conda/envs/ONTrack2_env/bin/vsearch --cluster_smallmem ${params.results_dir}/readsFiltering/${sample}.fasta --threads ${task.cpus} --usersort --id ${params.clustering_id_threshold} \
  --iddef 2 --clusterout_sort --fasta_width 0 --strand both --sizeout --consout \$clustering_dir/consensus_${sample}.fasta --clusters \$clustering_dir/cluster_${sample}_

  centroid_mac=\$(head -n1 \$clustering_dir/consensus_${sample}.fasta);
  n_cl=\$(ls \$clustering_dir | grep "cluster" | wc -l);
  id_centroid_tmp=\$(echo \$centroid_mac | sed \'s/centroid=//g\');
  id_centroid=\$(echo \$id_centroid_tmp | sed \'s/;seqs.*\$//g\');
  files=\$(find \$clustering_dir | grep cluster)

  for f in \$files ; do cat \$f | grep \$id_centroid && mac_file=\$f; done
  
  cat \$mac_file | grep \'^>\' | sed \'s/^>//\' | sed \'s/;.*//\' > \$clustering_dir/${sample}_ids_mac.txt

  /opt/conda/envs/ONTrack2_env/bin/seqtk subseq ${params.results_dir}/readsFiltering/${sample}.fasta \$clustering_dir/${sample}_ids_mac.txt > \$clustering_dir/${sample}_decont.fasta
  /opt/conda/envs/ONTrack2_env/bin/seqtk subseq ${params.results_dir}/readsFiltering/${sample}.fastq \$clustering_dir/${sample}_ids_mac.txt > ${params.results_dir}/readsClustering/${sample}/${sample}_decont.fastq


  """
  else
  """
  cp ${params.results_dir}/readsFiltering/${sample}.fastq \$clustering_dir/${sample}_decont.fastq
  cp ${params.results_dir}/readsFiltering/${sample}.fasta \$clustering_dir/${sample}_decont.fasta

  """
}

process draftConsensusCalling {
    input:
        val(sample) from readsClustering_draftConsensusCalling

    output:
        val(sample) into draftConsensusCalling_consensusPolishing        

    script:
    if(params.draftConsensusCalling)
    """
    mkdir -p ${params.results_dir}/draftConsensusCalling
    export PATH=\$PATH:/opt/conda/envs/ONTrack2_env/bin
    /opt/conda/envs/ONTrack2_env/bin/Rscript ${params.scripts_dir}/Obtain_draft_consensus.R fastq_file=${params.results_dir}/readsClustering/${sample}/${sample}_decont.fastq TRC=${params.target_reads_consensus} PLUR=${params.plurality} num_threads=${task.cpus} fast_alignment_flag=${params.fast_alignment_flag}
        
    """
    else
    """
        echo "Skipped."
    """
}

process consensusPolishing {
    input:
        val(sample) from draftConsensusCalling_consensusPolishing

    output:
        val(sample) into consensusPolishing_blastSearch
    
    script:
    if(params.consensusPolishing)
    """
        mkdir -p ${params.results_dir}/consensusPolishing
        export PATH=\$PATH:/opt/conda/envs/ONTrack2_env/bin/

        /opt/conda/envs/ONTrack2_env/bin/Rscript ${params.scripts_dir}/Polish_consensus.R draft_consensus=${params.results_dir}/draftConsensusCalling/${sample}/${sample}_draft_consensus.fasta fastq_file=${params.results_dir}/readsClustering/${sample}/${sample}_decont.fastq  TRP=${params.target_reads_polishing} num_threads=${task.cpus} primers_length=${params.primers_length} medaka_model=${params.medaka_model}
        
    """
    else
    """
        mkdir -p ${params.results_dir}/consensusPolishing
        /opt/conda/envs/ONTrack2_env/bin/seqtk trimfq ${params.results_dir}/draftConsensusCalling/${sample}/${sample}_draft_consensus.fasta -b ${params.primers_length} -e ${params.primers_length}  > ${params.results_dir}/consensusPolishing/${sample}/${sample}_consensus.fasta 
    """
}

process blastSearch {
    input:
          val(sample) from consensusPolishing_blastSearch

    output:
        

    script:
    if(params.blastSearch)
    """
        mkdir -p ${params.results_dir}/blastSearch/
        
        cp ${params.results_dir}/consensusPolishing/${sample}/${sample}_consensus.fasta ${params.results_dir}/blastSearch/
        /opt/conda/envs/ONTrack2_env/bin/blastn -num_threads ${task.cpus} -db ${params.blast_db} -query ${params.results_dir}/consensusPolishing/${sample}/${sample}_consensus.fasta > ${params.results_dir}/blastSearch/${sample}_blast_results.txt
        
    """
    else
    """
         mkdir -p ${params.results_dir}/blastSearch/

         cp ${params.results_dir}/consensusPolishing/${sample}/${sample}_consensus.fasta ${params.results_dir}/blastSearch/
    """
}

    
