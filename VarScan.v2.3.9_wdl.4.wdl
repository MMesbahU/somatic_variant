version 1.0

###################################################################################
## Version 2021-05-15
## Contact Md Mesbah Uddin <mdmesbah@gmail.com>
## Natarajan Lab MGH/Broad Institute

###################################################################################

####################################################################################
# Requires 
	# 1. VarScan jar file "VarScan.v2.3.9.jar" in the Google bucket (can be downloaded from http://varscan.sourceforge.net/ )
    # 2. GATK docker e.g. "us.gcr.io/broad-gatk/gatk:4.1.9.0"
    # 3. Reference fasta, fai and dict
    # 4. Bam file
#####    
    
### VarScan2: java -jar VarScan.v2.3.9.jar mpileup2cns --help
# OPTIONS:
#	--min-coverage	Minimum read depth at a position to make a call [8]
#	--min-reads2	Minimum supporting reads at a position to call variants [2]
#	--min-avg-qual	Minimum base quality at a position to count a read [15]
#	--min-var-freq	Minimum variant allele frequency threshold [0.01]
#	--min-freq-for-hom	Minimum frequency to call homozygote [0.75]
#	--p-value	Default p-value threshold for calling variants [99e-02]
#	--strand-filter	Ignore variants with >90% support on one strand [1]
#	--output-vcf	If set to 1, outputs in VCF format
#	--vcf-sample-list	For VCF output, a list of sample names in order, one per line
#	--variants	Report only variant (SNP/indel) positions [0]
#####################################################################################    

# VarScan
workflow VarScan2{
  input {
    File ref_fasta
    File ref_fai
    File ref_dict
    File bam
    File? bai
    
    # Samtools Parameters
    Int MapQ = 60
    Int baseQ = 60
    Int maxDepth = 1000000000
    
    # VarScan Parameters
    Int minCoverage = 10
    Int minSupportReads = 2
    Int minAvgBaseQual = 60
    Float minVAF = 0.01
    Float minVAFforHom = 0.75
    Float pvalue = 0.05
    Int strandFilter = 1
    Int reportVariant = 0
    
    #
    # File VARSCAN2="gs://fc-2800aa50-51bc-487c-b599-64f3341c6265/whi_Ref_hg19/VarScan.v2.3.9.jar"
    File VARSCAN2
    # String gatk_docker="us.gcr.io/broad-gatk/gatk:4.1.9.0"
    String gatk_docker
    Int boot_disk_size = 3
    Int preemptible_tries = 1
    Int cpu = 1
  }

call varScan2{
 		input:
        ref_fasta = ref_fasta,
        ref_fai = ref_fai,
        ref_dict = ref_dict,
        bam = bam,
        bai = bai,
        MapQ = MapQ,
        baseQ = baseQ,
        maxDepth = maxDepth,
        minCoverage =minCoverage,
        minSupportReads = minSupportReads,
        minAvgBaseQual = minAvgBaseQual,
        minVAF = minVAF,
        minVAFforHom = minVAFforHom,
        pvalue = pvalue,
        strandFilter = strandFilter,
        reportVariant = reportVariant,
        VARSCAN2 = VARSCAN2,
        docker_image = gatk_docker,
        preemptible_tries = preemptible_tries,
        cpu = cpu
        }

  output {
  	File output_vcf = varScan2.output_vcf
    File output_vcf_tbi = varScan2.output_vcf_tbi
    # File output_snp_vcf = varScan2.output_snp_vcf
    # File output_snp_tbi = varScan2.output_snp_tbi
    # File output_indel_vcf = varScan2.output_indel_vcf
    # File output_indel_tbi = varScan2.output_indel_tbi
  }
}

#Task Definitions
task varScan2{
	input {
    File ref_fasta
    File ref_fai
    File ref_dict
    File bam
    File? bai
    
    Int MapQ
    Int baseQ
    Int maxDepth
    
    Int minCoverage
    Int minSupportReads
    Int minAvgBaseQual
    Float minVAF
    Float minVAFforHom
    Float pvalue
    Int strandFilter
    Int reportVariant
    
    File VARSCAN2
    String name = basename(bam, ".bam")
    # Runtime parameters
    Int addtional_disk_size = 2
    String machine_mem_size = 5
    String docker_image
    Int preemptible_tries
    Int cpu
    }
    
	#adjust disk size
	Float input_size = size(bam, "GB")
	#Float ref_size = size(ref_fasta, "GB") + size(ref_fai, "GB")
	Float output_size = size(bam, "GB") * 0.5
    Int ref_size = ceil(size(ref_fasta, "GB") + size(ref_dict, "GB") + size(ref_fai, "GB"))
	Int disk_size = ceil(input_size + output_size ) + ref_size + addtional_disk_size
    
    #Calls samtools view to do the conversion
command <<<
        #Set -e and -o says if any command I run fails in this script, make sure to return a failure
set -e
set -o pipefail
        
# mv ~{name}.bam.bai ~{name}.bai

## SNPs and indels
samtools mpileup -d ~{maxDepth} -q ~{baseQ} -Q ~{MapQ} -f ~{ref_fasta} ~{bam} | java -jar ~{VARSCAN2} mpileup2cns --output-vcf 1 --min-var-freq ~{minVAF} --strand-filter ~{strandFilter} --variants ~{reportVariant} --p-value ~{pvalue} --min-coverage ~{minCoverage} --min-reads2 ~{minSupportReads} --min-avg-qual ~{minAvgBaseQual} --min-freq-for-hom ~{minVAFforHom} | bgzip -c > varScan2_snpindel.~{name}.vcf.gz

        ## SNPs
# samtools mpileup -d 100000000000000 -q ~{baseQ} -Q ~{MapQ} -f ~{ref_fasta} ~{bam} | java -jar ~{VARSCAN2} mpileup2snp --output-vcf 1 --min-var-freq 1e-8 --strand-filter 1 --variants 0 --p-value 1 --min-coverage 100 --min-reads2 1 --min-avg-qual 60 | bgzip -c > varScan2_mpileup2snp.~{name}.vcf.gz
         
       ## indels
# samtools mpileup -d 100000000000000 -q ~{baseQ} -Q ~{MapQ} -f ~{ref_fasta} ~{bam} | java -jar ~{VARSCAN2} mpileup2indel --output-vcf 1 --min-var-freq 1e-8 --strand-filter 1 --variants 0 --p-value 1 --min-coverage 100 --min-reads2 1 --min-avg-qual 60 | bgzip -c > varScan2_mpileup2indel.~{name}.vcf.gz
            
      ## Index vcf
tabix -p vcf varScan2_snpindel.~{name}.vcf.gz

# tabix -p vcf varScan2_mpileup2snp.~{name}.vcf.gz
      
# tabix -p vcf varScan2_mpileup2indel.~{name}.vcf.gz
            
    >>>

  runtime {
    docker: docker_image
    memory: machine_mem_size + " GB"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: preemptible_tries
    cpu: cpu
    zones: "us-central1-a us-central1-b us-central1-c us-central1-f us-east1-b us-east1-c us-east1-d us-west1-a us-west1-b us-west1-c"  
    }

    output {
        File output_vcf = "varScan2_snpindel.${name}.vcf.gz"
        File output_vcf_tbi = "varScan2_snpindel.${name}.vcf.gz.tbi"
        # File output_snp_vcf = "varScan2_mpileup2snp.${name}.vcf.gz"
        # File output_snp_tbi = "varScan2_mpileup2snp.${name}.vcf.gz.tbi"
        # File output_indel_vcf = "varScan2_mpileup2indel.${name}.vcf.gz"
        # File output_indel_tbi = "varScan2_mpileup2indel.${name}.vcf.gz.tbi"
    }
    
}


