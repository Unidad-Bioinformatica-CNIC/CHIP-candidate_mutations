workflow variant_calling{

  String sample_id
  String queue
  String accounting
  Object wf
  Object genome
  String bin_dir
  String enrichment_intervals
  String gnomad
  Boolean make_bamout
  Boolean run_ob_filter
  Array[String] reads_in
  Boolean cleanup = false 


  scatter (sample_reads_in in reads_in) {
    call sort_by_query {
      input:
        execution = wf.sort_by_query,
        queue=queue,
        accounting=accounting,
        reads_in = sample_reads_in
    }

    call uBAM_2_mappedBAM {
      input :
        execution = wf.uBAM_2_mappedBAM,
        queue=queue,
        accounting=accounting,
        genome = genome,
        reads_in = sort_by_query.reads_out
    }
  }

  call MarkDuplicates {
    input :
      sample_id = sample_id,
      execution = wf.MarkDuplicates,
      queue=queue,
      accounting=accounting,
      alignment_in =  uBAM_2_mappedBAM.alignment,
      alignment_index_in =  uBAM_2_mappedBAM.alignment_index,
      cleanup = cleanup
  }

  call GroupReadsByUmi {
    input :
      execution = wf.GroupReadsByUmi,
      queue=queue,
      accounting=accounting,
      alignment_in =  MarkDuplicates.alignment,
      bin_dir = bin_dir,
      cleanup = cleanup
  }

  call TemplateSortBam {
    input :
      execution = wf.TemplateSortBam,
      queue=queue,
      accounting=accounting,
      reads_in=GroupReadsByUmi.alignment,
      bin_dir = bin_dir
  }

  call CallMolecularConsensusReads {
    input :
      execution = wf.CallMolecularConsensusReads,
      queue=queue,
      accounting=accounting,
      reads_in=TemplateSortBam.templated_sorted,
      bin_dir = bin_dir
  }

  call mapping {
    input :
      execution = wf.mapping,
      queue=queue,
      accounting=accounting,
      genome = genome,
      reads_in = CallMolecularConsensusReads.consensus_reads
  }

  call FilterConsensusReads {
    input :
      execution = wf.FilterConsensusReads,
      queue=queue,
      accounting=accounting,
      genome = genome,
      alignment_in = mapping.alignment,
      alignment_index_in = mapping.alignment_index,
      bin_dir = bin_dir,
      cleanup = cleanup
  }

  call ClipBam {
    input :
      execution = wf.ClipBam,
      queue=queue,
      accounting=accounting,
      genome = genome,
      sample_id = sample_id,
      alignment_in = FilterConsensusReads.filtered_alignment,
      bin_dir = bin_dir
  }

  call Bam2Cram {
    input :
      sname = sample_id,
      execution = wf.Bam2Cram,
      queue = queue,
      accounting = accounting,
      genome = genome,
      alignment_in = ClipBam.clipped_alignment,
      alignment_index_in = ClipBam.clipped_alignment_index
}

  # Stats

  # uBAM mapped
  call getBarcodes {
    input :
      execution = wf.getBarcodes,
      queue=queue,
      accounting=accounting,
      alignment_in=GroupReadsByUmi.alignment
  }

  # call QualityScoreDistribution {
  #   input :
  #     execution = wf.QualityScoreDistribution,
  #     queue=queue,
  #     accounting=accounting,
  #     alignment_in=GroupReadsByUmi.alignment
  # }

  call CollectHsMetrics as raw_HsMetrics {
    input :

      execution = wf.CollectHsMetrics,
      queue=queue,
      accounting=accounting,
      genome=genome,
      alignment_in=GroupReadsByUmi.alignment,
      min_base_qual = 0,
      min_mapping_qual = 0,
      enrichment_intervals = enrichment_intervals
  }

  call CollectHsMetrics as consensus_HsMetrics {
    input :
      execution = wf.CollectHsMetrics,
      queue=queue,
      accounting=accounting,
      genome=genome,
      alignment_in=ClipBam.clipped_alignment,
      min_base_qual = 0,
      min_mapping_qual = 0,
      enrichment_intervals = enrichment_intervals
}

  call BAM_flagstat {
    input:
      execution =wf.BAM_flagstat,
      queue = queue,
      accounting = accounting,
      bamfile = ClipBam.clipped_alignment,
      bamfile_index = ClipBam.clipped_alignment_index
  }

  call SomaticVC {
    input :
      sample_id = sample_id,
      execution = wf.SomaticVC,
      queue=queue,
      accounting=accounting,
      genome = genome,
      enrichment_intervals = enrichment_intervals,
      gnomad = gnomad,
      make_bamout = make_bamout,
      run_ob_filter = run_ob_filter,
      alignment_in = Bam2Cram.cram_file,
      alignment_in_index = Bam2Cram.cram_index
}
  
  call create_gVCF {
    input :
      sample_id = sample_id,
      execution = wf.create_gVCF,
      queue = queue,
      accounting = accounting,
      genome = genome,
      enrichment_intervals = enrichment_intervals,
      alignment_in = Bam2Cram.cram_file,
      alignment_in_index = Bam2Cram.cram_index
  }

  output {
    #File alignment = uBAM_2_mappedBAM.alignment
    #File alignment_index = uBAM_2_mappedBAM.alignment_index
    Array[File] raw_aligment_log = uBAM_2_mappedBAM.aligment_log
    File markduplicates_metrics = MarkDuplicates.markduplicates_metrics 
    File umi_metrics = MarkDuplicates.umi_metrics
    File hist = GroupReadsByUmi.hist
    File consensus_alignment_log = mapping.alignment_log
    File flagstat = BAM_flagstat.flagstat
    File RX_barcode = getBarcodes.RX_barcode
    #File qual_score_dist = QualityScoreDistribution.qual_score_dist
    File raw_hs_metrics = raw_HsMetrics.hs_metrics
    File raw_target_coverage = raw_HsMetrics.target_coverage
    File raw_base_coverage = raw_HsMetrics.base_coverage
    File consensus_hs_metrics = consensus_HsMetrics.hs_metrics
    File consensus_target_coverage = consensus_HsMetrics.target_coverage
    File consensus_base_coverage = consensus_HsMetrics.base_coverage

    #File out_alignment = ClipBam.clipped_alignment
    #File out_alignment_index = ClipBam.clipped_alignment_index
    File cram_file = Bam2Cram.cram_file
    File cram_index = Bam2Cram.cram_index
    File vcf_out = SomaticVC.vcf_out
    File vcf_out_idx = SomaticVC.vcf_out_idx
    File vcf_out_stats = SomaticVC.vcf_out_stats
    File vcf_filtering_stats = SomaticVC.vcf_filtering_stats
    File gvcf_out = create_gVCF.gvcf_out
  }
}

task sort_by_query {

  Object execution
  String queue
  String accounting
  String reads_in

  String file_prefix  = basename (reads_in, ".bam")

  command {

    # Load enviroment
    source /programs/GATK/env.sh

    # job
    picard ${execution.java_args_sort_by_query} SortSam \
      TMP_DIR=$TMPDIR \
      I=${reads_in} \
      O=${file_prefix}.bam \
      SORT_ORDER="queryname"
  
    rc=$?

    # Saca el RC
    echo "ExitCode:$rc"

    exit $rc
  }
  runtime {
    backend : 'SGE_nope_long'
    cpu : execution.cpu
    memory : execution.memory
    accounting : accounting
    sge_project : queue
  }
  output {
    File reads_out = "${file_prefix}.bam"
  }
}

task uBAM_2_mappedBAM {

  Object execution
  String queue
  String accounting
  Object genome
  String reads_in

  String file_prefix  = basename (reads_in, ".bam") + ".mapped"

  String ref_fasta = "${genome.basedir}/${genome.ref_fasta}"

  String bwa_commandline="-K 100000000 -p -Y -v 3 -t ${execution.cpu} ${ref_fasta}"
  
  command {

    # Load enviroment
    source /programs/GATK/env.sh

    # job
    picard ${execution.java_args_SamToFastq} SamToFastq \
      TMP_DIR=$TMPDIR \
      I=${reads_in} \
      F=/dev/stdout \
      INTERLEAVE=true \
      INCLUDE_NON_PF_READS=false \
      CLIPPING_ATTRIBUTE=XT \
      CLIPPING_ACTION=X| \
    bwa mem ${bwa_commandline} /dev/stdin - 2>${file_prefix}.bwa.stderr.log | \
    picard ${execution.java_args_MergeBamAlignment} MergeBamAlignment \
        TMP_DIR=$TMPDIR \
        VALIDATION_STRINGENCY=SILENT \
        ALIGNED_BAM=/dev/stdin \
        UNMAPPED_BAM=${reads_in} \
        R=${ref_fasta} \
        O=${file_prefix}.bam \
        EXPECTED_ORIENTATIONS=FR \
        PAIRED_RUN=true \
        SORT_ORDER="coordinate" \
        ADD_MATE_CIGAR=true \
        MAX_INSERTIONS_OR_DELETIONS=-1 \
        ALIGNER_PROPER_PAIR_FLAGS=false \
        CREATE_INDEX=true && \
    samtools quickcheck -q ${file_prefix}.bam && \
    grep -m1 "read .* ALT contigs" ${file_prefix}.bwa.stderr.log | grep -v "read 0 ALT contigs"

    rc=$?

    # Saca el RC
    echo "ExitCode:$rc"

    exit $rc
  }
  runtime {
    backend : 'SGE_nope_long'
    cpu : execution.cpu
    memory : execution.memory
    accounting : accounting
    sge_project : queue
  }
  output {
    File alignment = "${file_prefix}.bam"
    File alignment_index = "${file_prefix}.bai"
    File aligment_log = "${file_prefix}.bwa.stderr.log"
  }
}

task MarkDuplicates {

  String sample_id
  Object execution
  String queue
  String accounting
  Array[File] alignment_in
  Array[File] alignment_index_in
  Boolean cleanup

  String file_prefix = "${sample_id}.markdup"

  command {

    # Load enviroment
    source /programs/GATK/env.sh

    picard ${execution.java_args_MarkDuplicates} UmiAwareMarkDuplicatesWithMateCigar \
      TMP_DIR=$TMPDIR \
      I=${sep=' I=' alignment_in} \
      O=${file_prefix}.bam \
      M=${file_prefix}.markduplicates_metrics \
      UMI_METRICS=${file_prefix}.umi_metrics \
      MAX_EDIT_DISTANCE_TO_JOIN=1 \
      TAGGING_POLICY=All \
      TAG_DUPLICATE_SET_MEMBERS=true \
      REMOVE_SEQUENCING_DUPLICATES=true \
      OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
      ADD_PG_TAG_TO_READS=false

    rc=$?

    if [ ${cleanup} = true ] && samtools quickcheck -q ${file_prefix}.bam; then
      rm -f ../../call-uBAM_2_mappedBAM/shard-0/execution/*.bam
    fi

    # Saca el RC
    echo "ExitCode:$rc"
    exit $rc
  }
  runtime {
    backend : 'SGE_nope_long'
    cpu : execution.cpu
    memory : execution.memory
    accounting : accounting
    sge_project : queue
  }
  output {
    File markduplicates_metrics = "${file_prefix}.markduplicates_metrics"
    File umi_metrics = "${file_prefix}.umi_metrics"
    File alignment = "${file_prefix}.bam"
  }
}

task GroupReadsByUmi {

  Object execution
  String queue
  String accounting
  File alignment_in
  String bin_dir
  Boolean cleanup

  String file_prefix  = basename (alignment_in, ".bam") + ".grouped"

  command {
    # Load enviroment
    source /programs/GATK/env.sh

    java ${execution.java_args_GroupReadsByUmi} -jar ${bin_dir}/fgbio.jar \
      --tmp-dir=$TMPDIR \
    GroupReadsByUmi \
      -i ${alignment_in} \
      -f ${file_prefix}.hist \
      -o ${file_prefix}.bam \
      --include-non-pf-reads=false \
      --strategy=${execution.grouping_strategy} \
      --edits=${execution.edit_distance}

    rc=$?

    if [ ${cleanup} = true ] && samtools quickcheck -q ${file_prefix}.bam; then
      rm -f ../../call-MarkDuplicates/execution/*.bam
    fi

    # Saca el RC
    echo "ExitCode:$rc"
      
    exit $rc
  }
  runtime {
    backend : 'SGE_nope_long'
    cpu : execution.cpu
    memory : execution.memory
    accounting : accounting
    sge_project : queue
  }
  output {
    File hist = "${file_prefix}.hist"
    File alignment = "${file_prefix}.bam"
  }
}

task TemplateSortBam {

  Object execution
  String queue
  String accounting
  File reads_in
  String bin_dir

  String file_prefix  = basename (reads_in, ".bam") + ".sorted"

  command {
    # Load enviroment
    source /programs/GATK/env.sh

    java ${execution.java_args_SortBam} -jar ${bin_dir}/fgbio.jar \
      --tmp-dir=$TMPDIR \
    SortBam \
    --input=${reads_in} \
    --sort-order=TemplateCoordinate \
    --output ${file_prefix}.bam && \
    samtools quickcheck -q ${file_prefix}.bam 
    
    # No borramos ../../call-GroupReadsByUmi porque lo usan multiples tareas
      
    rc=$?

    # Saca el RC
    echo "ExitCode:$rc"

    exit $rc
  }
  runtime {
    backend : 'SGE_nope_long'
    cpu : execution.cpu
    memory : execution.memory
    accounting : accounting
    sge_project : queue
  }
  output {
    File templated_sorted = "${file_prefix}.bam"
  }
}

task CallMolecularConsensusReads {

  Object execution
  String queue
  String accounting
  File reads_in
  String bin_dir

  String file_prefix  = basename (reads_in, ".bam") + ".consensus"

  command {
    # Load enviroment
    source /programs/GATK/env.sh

    java ${execution.java_args_CallMolecularConsensusReads} -jar ${bin_dir}/fgbio.jar \
      --tmp-dir=$TMPDIR \
    CallMolecularConsensusReads \
      -i ${reads_in} \
      -o ${file_prefix}.bam \
      --min-reads=1 
    rc=$?

    # Saca el RC
    echo "ExitCode:$rc"

    exit $rc
  }
  runtime {
    backend : 'SGE_nope_long'
    cpu : execution.cpu
    memory : execution.memory
    accounting : accounting
    sge_project : queue
  }
  output {
    File consensus_reads = "${file_prefix}.bam"
  }
}

task mapping {

  Object execution
  String queue
  String accounting
  Object genome
  File reads_in

  String ref_fasta = "${genome.basedir}/${genome.ref_fasta}"

  String bwa_commandline="-K 100000000 -p -Y -v 3 -t ${execution.cpu} ${ref_fasta}"

  String file_prefix  = basename (reads_in, ".bam") + ".mapped"


  command {

    # Load enviroment
    source /programs/GATK/env.sh

    # job
    picard ${execution.java_args_SamToFastq} SamToFastq \
      TMP_DIR=$TMPDIR \
      I=${reads_in} \
      F=/dev/stdout \
      INTERLEAVE=true \
      INCLUDE_NON_PF_READS=false| \
    bwa mem ${bwa_commandline} /dev/stdin - 2>${file_prefix}.bwa.stderr.log | \
    picard ${execution.java_args_MergeBamAlignment} MergeBamAlignment \
        TMP_DIR=$TMPDIR \
        VALIDATION_STRINGENCY=SILENT \
        ALIGNED_BAM=/dev/stdin \
        UNMAPPED_BAM=${reads_in} \
        R=${ref_fasta} \
        O=${file_prefix}.bam \
        EXPECTED_ORIENTATIONS=FR \
        PAIRED_RUN=true \
        SORT_ORDER="coordinate" \
        ADD_MATE_CIGAR=true \
        MAX_INSERTIONS_OR_DELETIONS=-1 \
        ALIGNER_PROPER_PAIR_FLAGS=false \
        CREATE_INDEX=true && \
    samtools quickcheck -q ${file_prefix}.bam && \
    grep -m1 "read .* ALT contigs" ${file_prefix}.bwa.stderr.log | grep -v "read 0 ALT contigs"

    rc=$?

    # Saca el RC
    echo "ExitCode:$rc"

    exit $rc
  }
  runtime {
    backend : 'SGE_nope_long'
    cpu : execution.cpu
    memory : execution.memory
    accounting : accounting
    sge_project : queue
  }
  output {
    File alignment = "${file_prefix}.bam"
    File alignment_index = "${file_prefix}.bai"
    File alignment_log = "${file_prefix}.bwa.stderr.log"
  }
}

task FilterConsensusReads {

  Object execution
  String queue
  String accounting
  Object genome
  File alignment_in
  File alignment_index_in
  String bin_dir
  Boolean cleanup

  String ref_fasta = "${genome.basedir}/${genome.ref_fasta}"

  String file_prefix  = basename (alignment_in, ".bam") + ".filtered"

  command {
    # Load enviroment
    source /programs/GATK/env.sh

    java ${execution.java_args_FilterConsensusReads} -jar ${bin_dir}/fgbio.jar \
      --tmp-dir=$TMPDIR \
    FilterConsensusReads \
      -i ${alignment_in} \
      -o ${file_prefix}.bam \
      --ref=${ref_fasta} \
      --min-reads=${execution.supporting_min_reads} \
      --min-base-quality=${execution.supporting_min_base_quality} \
      --reverse-per-base-tags=${execution.reverse_per_base_tags} \
      --require-single-strand-agreement=${execution.require_single_strand_agreement}

    rc=$?

    if [ ${cleanup} = true ] && samtools quickcheck -q ${file_prefix}.bam; then
          rm -f ../../call-mapping/execution/*.bam
          rm -f ../../call-TemplateSortBam/execution/*.bam
    fi
    

    # Saca el RC
    echo "ExitCode:$rc"

    exit $rc
  }
  runtime {
    backend : 'SGE_nope_long'
    cpu : execution.cpu
    memory : execution.memory
    accounting : accounting
    sge_project : queue
  }
  output {
    File filtered_alignment = "${file_prefix}.bam"
  }
}

task ClipBam {

  Object execution
  String queue
  String accounting
  Object genome
  String sample_id
  File alignment_in
  String bin_dir

  String ref_fasta = "${genome.basedir}/${genome.ref_fasta}"

  command {
    # Load enviroment
    source /programs/GATK/env.sh

    java ${execution.java_args_ClipBam} -jar ${bin_dir}/fgbio.jar \
      --tmp-dir=$TMPDIR \
    ClipBam \
      -i ${alignment_in} \
      -o ${sample_id}.bam \
      --ref=${ref_fasta} \
      --clip-overlapping-reads=${execution.clip_overlapping_reads} \
      --clipping-mode=${execution.clipping_mode} && \
    samtools quickcheck -q ${sample_id}.bam

    rc=$?

    # Saca el RC
    echo "ExitCode:$rc"

    exit $rc
  }
  runtime {
    backend : 'SGE_nope_long'
    cpu : execution.cpu
    memory : execution.memory
    accounting : accounting
    sge_project : queue
  }
  output { 
    File clipped_alignment = "${sample_id}.bam"
    File clipped_alignment_index = "${sample_id}.bai"
  }
}

task Bam2Cram {

    Object execution
    String queue
    String accounting
    Object genome
    String sname
    File alignment_in
    File alignment_index_in

    File genome_fasta = "${genome.basedir}/${genome.ref_fasta}"
    command {

    source /programs/GATK/env.sh

    #samtools
    samtools view -C -T ${genome_fasta} -o ${sname}.cram ${alignment_in}
    samtools index ${sname}.cram
    
    # No borramos call-ClipBam porque lo usan varias carreras

    exit $rc

    }
    runtime {
        backend : 'SGE_nope'
        cpu : execution.cpu
        memory : execution.memory
        accounting : accounting
        sge_project : queue
    }
    output {
    File cram_file = "${sname}.cram"
    File cram_index = "${sname}.cram.crai"
    }
}

task BAM_flagstat {

  Object execution
  String queue
  String accounting
  File bamfile
  File bamfile_index

  String file_prefix = basename( bamfile , ".bam")

  command {

    # Load enviroment
    source /programs/GATK/env.sh

    samtools flagstat ${bamfile} > ${file_prefix}.flagstat

    rc=$?

    # Saca el RC
    echo "ExitCode:$rc"

    exit $rc
  }
  runtime {
    backend : 'SGE_nope_long'
    cpu : execution.cpu
    memory : execution.memory
    accounting : accounting
    sge_project : queue
  }
  output {
  File flagstat = "${file_prefix}.flagstat"
  }
}

task getBarcodes {

  Object execution
  String queue
  String accounting
  File alignment_in

  String file_prefix = basename( alignment_in , ".bam")

  command <<<

    # Load enviroment
    source /programs/GATK/env.sh
    # 0x0040: first in pair
    # 0x0100: not primary alignment 

    samtools view -f 0x0040 -F 0x0100 ${alignment_in}|awk '{for (i=12; i<=NF; ++i) { if ($i ~ "^RX:Z:"){ s = substr($i,6,length($i)-5); }; print s }}'|sort|uniq > ${file_prefix}.RX.barcode
    rc=$?

    # Saca el RC
    echo "ExitCode:$rc"

    exit $rc
  >>>
  runtime {
    backend : 'SGE_nope_long'
    cpu : execution.cpu
    memory : execution.memory
    accounting : accounting
    sge_project : queue
  }
  output {
    File RX_barcode = "${file_prefix}.RX.barcode"
  }
}

task QualityScoreDistribution {

  Object execution
  String queue
  String accounting
  File alignment_in

  String file_prefix = basename( alignment_in , ".bam")

  command {

    # Load enviroment
    source /programs/GATK/env.sh

    picard ${execution.java_args_QualityScoreDistribution} QualityScoreDistribution \
      I=${alignment_in} \
      O=${file_prefix}.qual_score_dist \
      CHART=${file_prefix}.qual_score_dist.pdf \
      PF_READS_ONLY=true

    rc=$?

    # Saca el RC
    echo "ExitCode:$rc"
    exit $rc
  }
  runtime {
    backend : 'SGE_nope_long'
    cpu : execution.cpu
    memory : execution.memory
    accounting : accounting
    sge_project : queue
  }
  output {
    File qual_score_dist = "${file_prefix}.qual_score_dist"
  }
}

task CollectHsMetrics {

  Object execution
  String queue
  String accounting
  Object genome
  File alignment_in
  Int min_base_qual
  Int? min_mapping_qual=20
  String enrichment_intervals

  String file_prefix = basename(alignment_in, ".bam")

  String ref_fasta = "${genome.basedir}/${genome.ref_fasta}"

  command {

    # Load enviroment
    source /programs/GATK/env.sh

    picard ${execution.java_args_CollectHsMetrics} CollectHsMetrics \
      I=${alignment_in} \
      O=${file_prefix}.hs_metrics \
      R=${ref_fasta} \
      BAIT_INTERVALS=${enrichment_intervals} \
      TARGET_INTERVALS=${enrichment_intervals} \
      PER_TARGET_COVERAGE=${file_prefix}.target_coverage \
      PER_BASE_COVERAGE=${file_prefix}.base_coverage  \
      COVERAGE_CAP=${execution.covarage_cap} \
      MINIMUM_BASE_QUALITY=${min_base_qual} \
      MINIMUM_MAPPING_QUALITY=${min_mapping_qual}  

    # Saca el RC
    echo "ExitCode:$?"
    exit $rc
  }
  runtime {
    backend : 'SGE_nope_long'
    cpu : execution.cpu
    memory : execution.memory
    accounting : accounting
    sge_project : queue
  }
  output {
    File hs_metrics ="${file_prefix}.hs_metrics"
    File target_coverage = "${file_prefix}.target_coverage"
    File base_coverage = "${file_prefix}.base_coverage"
  }
}

task SomaticVC {
  String sample_id
  Object execution
  String queue
  String accounting
  Object genome
  String enrichment_intervals
  String gnomad
  Boolean make_bamout
  Boolean run_ob_filter
  File alignment_in
  File alignment_in_index
  
  String genome_fasta = "${genome.basedir}/${genome.ref_fasta}"
  
  command {
      # Load enviroment
      source /programs/GATK/env_vep.sh

      set -e

      gatk --java-options "${execution.java_args_SomaticVC}" Mutect2 \
          ${"-L " + enrichment_intervals} \
          ${"-R " + genome_fasta} \
          ${"-I " + alignment_in} \
          ${"--germline-resource " + gnomad} \
          -O "${sample_id}.raw.vcf.gz" \
          ${true='--bam-output bamout.bam' false='' make_bamout} \
          ${true='--f1r2-tar-gz f1r2.tar.gz' false='' run_ob_filter}

      #mv ${sample_id}.raw.vcf.gz.stats ${sample_id}.vcf.gz.stats

      gatk --java-options "${execution.java_args_FilterMutectCalls}" FilterMutectCalls \
        ${"-R " + genome_fasta} \
        -V ${sample_id}.raw.vcf.gz \
        --filtering-stats "${sample_id}.vcf.gz.filteringStats.tsv" \
        -O "${sample_id}.tmp.vcf.gz"

      gatk RenameSampleInVcf \
      -I ${sample_id}.tmp.vcf.gz \
      -O ${sample_id}.vcf.gz \
      --NEW_SAMPLE_NAME ${sample_id}

    bcftools index -t ${sample_id}.vcf.gz

  }
  runtime {
    backend : 'SGE_nope_long'
    cpu : execution.cpu
    memory : execution.memory
    accounting : accounting
    sge_project : queue
  }
  output {
    File vcf_out = "${sample_id}.vcf.gz"
    File vcf_out_idx = "${sample_id}.vcf.gz.tbi"
    File vcf_out_stats = "${sample_id}.raw.vcf.gz.stats"
    File vcf_filtering_stats = "${sample_id}.vcf.gz.filteringStats.tsv"
    #File bamout = "bamout.bam"
    #File bamout_idx = "bamout.bai"
    #File f1r2 = "f1r2.tar.gz"
  }
}

task GerminalVC {
  String sample_id
  Object execution
  String queue
  String accounting
  Object genome
  String enrichment_intervals
  File alignment_in
  File alignment_in_index
  
  String genome_fasta = "${genome.basedir}/${genome.ref_fasta}"
  String dbsnp = "${genome.basedir}/${genome.dbSNP_vcf}"
  String snps_hq = "${genome.basedir}/${genome.one_thousand_genomes_resource_vcf}"
  
  command <<<
      # Load enviroment
      source /programs/GATK/env_vep.sh
      set -e

      gatk --java-options "${execution.java_args_HaplotypeCaller}" HaplotypeCaller \
        ${"-I " + alignment_in} \
        ${"-R " + genome_fasta} \
        -O "${sample_id}.tmp.vcf.gz" \
        ${"-L " + enrichment_intervals} \
        --max-reads-per-alignment-start 0 \
        ${"--dbsnp " + dbsnp}

      gatk RenameSampleInVcf \
        -I ${sample_id}.tmp.vcf.gz \
        -O ${sample_id}.germline.vcf.gz \
        --NEW_SAMPLE_NAME ${sample_id}

      bcftools index -t ${sample_id}.germline.vcf.gz
      bcftools isec ${sample_id}.germline.vcf.gz ${snps_hq} -n =2 -w 1 > ${sample_id}.germline.HQ.vcf.gz

    >>>
  runtime {
    backend : 'SGE_nope_long'
    cpu : execution.cpu
    memory : execution.memory
    accounting : accounting
    sge_project : queue
  }
  output {
    File vcf_out = "${sample_id}.germline.vcf.gz"
    File vcf_out_idx = "${sample_id}.germline.vcf.gz.tbi"
    File vcf_hq_out = "${sample_id}.germline.HQ.vcf.gz"
  }
}

task create_gVCF {
  String sample_id
  Object execution
  String queue
  String accounting
  Object genome
  String enrichment_intervals
  File alignment_in
  File alignment_in_index
  
  String genome_fasta = "${genome.basedir}/${genome.ref_fasta}"
  String dbsnp = "${genome.basedir}/${genome.dbSNP_vcf}"  
  command {
      # Load enviroment
      source /programs/GATK/env.sh

      set -e

      gatk --java-options "${execution.java_args_HaplotypeCaller}" HaplotypeCaller \
        ${"-I " + alignment_in} \
        ${"-R " + genome_fasta} \
        -O "${sample_id}.tmp.g.vcf.gz" \
        ${"-L " + enrichment_intervals} \
        --max-reads-per-alignment-start 0 \
        -ERC GVCF \
        --pair-hmm-implementation ${execution.gatk_gkl_pairhmm_implementation} \
        --native-pair-hmm-threads ${execution.gatk_gkl_pairhmm_threads} \
        --smith-waterman ${execution.smith_waterman_implementation} && \
        gatk RenameSampleInVcf \
          -I ${sample_id}.tmp.g.vcf.gz \
          -O ${sample_id}.g.vcf.gz \
          --NEW_SAMPLE_NAME ${sample_id}      

  }
  runtime {
    backend : 'SGE_nope_long'
    cpu : execution.cpu
    memory : execution.memory
    accounting : accounting
    sge_project : queue
  }
  output {
    File gvcf_out = "${sample_id}.g.vcf.gz"
  }
}






