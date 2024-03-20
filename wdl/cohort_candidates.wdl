workflow cohort_candidates {

	String cohort_id
	Array[String] samples
	String queue
  	String accounting
	Object wf
	Object genome
	String references_dir

	call merge_VCF {
		input:
			cohort_id = cohort_id,
			samples = samples,
			execution = wf.merge_VCF,
      		queue = queue,
      		accounting = accounting
	}

	call normalize_VCF {
		input:
			cohort_id = cohort_id,
			vcf_in = merge_VCF.merged,
			genome = genome,
			execution = wf.normalize_VCF,
  			queue = queue,
  			accounting = accounting,
			references_dir = references_dir
	}

	call annotate_VCF {
		input:
			cohort_id = cohort_id,
			vcf_in = normalize_VCF.normalized,
			vcf_in_idx = normalize_VCF.normalized_idx,
			genome = genome,
			execution = wf.annotate_VCF,
      		queue = queue,
      		accounting = accounting,
			references_dir = references_dir
	}
	output {
		File annotated = annotate_VCF.annotated
	}
}

task merge_VCF {
	String cohort_id
	Array[String] samples
	Object execution
  	String queue
  	String accounting

	command <<<
		source /programs/GATK/env.sh

		bcftools merge -l ${write_lines(samples)} --filter-logic x --merge all -O v -o ${cohort_id}.merged.vcf.gz
	>>>
	runtime {
    	backend : 'SGE_nope_long'
    	cpu : execution.cpu
    	memory : execution.memory
    	accounting : accounting
    	sge_project : queue
	}
	output {
		File merged = "${cohort_id}.merged.vcf.gz"
	} 
}

task normalize_VCF {
	String cohort_id
	File vcf_in
	Object genome
	Object execution
  	String queue
  	String accounting
	String references_dir

	String genome_fasta = "${genome.basedir}/${genome.ref_fasta}"

	command <<<
		source /programs/GATK/env.sh
		# NOTA: este entorno tiene una version de bcftools mas nueva

		bcftools norm -m-any ${vcf_in} | bcftools norm -Oz --check-ref w -f ${genome_fasta} -o ${cohort_id}.merged.norm.vcf.gz

		bcftools index -f -t ${cohort_id}.merged.norm.vcf.gz
	>>>
	runtime {
    	backend : 'SGE_nope_long'
    	cpu : execution.cpu
    	memory : execution.memory
    	accounting : accounting
    	sge_project : queue
	}
	output {
		File normalized = "${cohort_id}.merged.norm.vcf.gz"
		File normalized_idx = "${cohort_id}.merged.norm.vcf.gz.tbi"
	} 
}

task annotate_VCF {
	String cohort_id
	File vcf_in
	File vcf_in_idx
	Object genome
	Object execution
  	String queue
  	String accounting
	String references_dir
	
	String genome_fasta = "${genome.basedir}/${genome.ref_fasta}"
	
	command <<<

	source /programs/conda_genomics/envs/vep.sh

	vep --offline --cache --refseq -species homo_sapiens --cache_version ${execution.cache_version} \
	--dir_cache ${execution.cache_dir} --dir_plugins ${execution.plugins_dir} \
	--use_transcript_ref \
	--fasta ${genome_fasta} \
	-i ${vcf_in} \
	-o ${cohort_id}.vcf.gz --vcf --compress_output gzip \
	-gene_phenotype Cancer Gene Census ClinVar COSMIC HGMD-PUBLIC NHGRI-EBI GWAS catalog \
	--check_existing \
	--everything \
	--humdiv \
	--no_stats \
	--plugin CADD,"${execution.CADD}" \
	--plugin LoFtool,"${execution.plugins_dir}/LoFtool_scores.txt"

	#SIFT4G_score,SIFT4G_converted_rankscore,SIFT4G_pred,Polyphen2_HDIV_score,Polyphen2_HDIV_rankscore,Polyphen2_HDIV_pred,GERP++_NR,GERP++_RS,GERP++_RS_rankscore,SIFT_score,SIFT_converted_rankscore,SIFT_pred,MetaLR_score,MetaLR_rankscore,MetaLR_pred,MutationAssessor_score,MutationAssessor_rankscore,MutationAssessor_pred,REVEL_score,REVEL_rankscore

	#NOTA: --use_transcript_ref; usamos el alelo de referencia no el dado por el VCF y corregido para REfSeq(ver https://www.ensembl.org/info/docs/tools/vep/script/vep_other.html#refseq_bam)
	#NOTA: Quitamos "Downstream" porque parece que da algun error
	>>>
	runtime {
    	backend : 'SGE_nope_long'
    	cpu : execution.cpu
    	memory : execution.memory
    	accounting : accounting
    	sge_project : queue
	}
	output {
		File annotated = "${cohort_id}.vcf.gz"
	}
}