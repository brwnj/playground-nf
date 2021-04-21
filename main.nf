params.manifest = false
params.fasta = false
params.outdir
params.gff = false
sexchroms = params.sexchroms ?: 'X,Y'
sexchroms = sexchroms.replaceAll(" ", "")

// bed intervals of interest
// params.regions

outdir = file(params.outdir)
// regions = file(params.regions)

if( !params.manifest ) { exit 1, "--manifest is not defined" }
if( !params.fasta ) { exit 1, "--fasta is not defined" }
// assumes .alt, .amb, .ann, .bwt, .pac, and .sa are present

fasta = file(params.fasta)
faidx = file("${params.fasta}.fai")
gff = params.gff ? file(params.gff) : params.gff

Channel
    .value(file("${params.fasta}.{amb,ann,bwt,pac,sa}"))
    .set { bwaidx_ch }

// ingests the sample/file manifest
// format of the CSV is: sample, group, lane, r1, r2
Channel.fromPath(params.manifest)
    .splitCsv()
    .map { row ->
        def sample_id = row[0]
        // some unique identifier in case samples are duplicated on same sequencing run
        def group_id = row[1]
        def run_id = row[2]
        def r1 = file(row[3])
        def r2 = file(row[4])
        [sample_id, group_id, run_id, r1, r2]
    }
    .set { fastq_ch }


// qc the fastqs
process fastp {
    publishDir path: "$outdir/json", pattern: "*.json"
    publishDir path: "$outdir/html", pattern: "*.html"

    input:
    set sample_id, group_id, run_id, file(r1), file(r2) from fastq_ch

    output:
    set sample_id, group_id, run_id, file("${sample_id}_R1.fastq.gz"), file("${sample_id}_R2.fastq.gz") into checked_fastq_ch
    file("${sample_id}_${run_id}.fastp.json") into fastp_report_ch
    file("${sample_id}.fastp.html")

    script:
    """
    fastp --thread ${task.cpus} --in1 $r1 --out1 ${sample_id}_R1.fastq.gz \
        --in2 $r2 --out2 ${sample_id}_R2.fastq.gz -y \
        --json ${sample_id}_${run_id}.fastp.json --html ${sample_id}.fastp.html
    """
}


process bwamem {
    tag "${sample_id}:${run_id}"

    input:
    set sample_id, group_id, run_id, file(r1), file(r2) from checked_fastq_ch
    file(fasta)
    file(bwaidx) from bwaidx_ch

    output:
    set sample_id, file("${sample_id}_${group_id}_${run_id}.bam") into bwa_ch

    script:
    rg = "@RG\\tID:${run_id}\\tPU:${run_id}\\tSM:${sample_id}\\tLB:${sample_id}\\tPL:illumina"
    """
    bwa mem -K 100000000 -R \"${rg}\" -t ${task.cpus} -M ${fasta} ${r1} ${r2} \
        | samtools sort -n --threads ${task.cpus} -m 2G --output-fmt BAM -o ${sample_id}_${group_id}_${run_id}.bam
    """
}


process markduplicates {
    tag "$sample_id"
    publishDir path: "$outdir/alignments", overwrite: true

    input:
    set sample_id, file(bam) from bwa_ch.groupTuple()

    output:
    set sample_id, file("${sample_id}.md.bam"), file("${sample_id}.md.bam.bai") into md_ch

    script:
    """
    gatk --java-options -Xmx${task.memory.toGiga()}g MarkDuplicatesSpark \
        ${bam.collect { "--input $it" }.join(" ")} \
        --output ${sample_id}.md.bam \
        --tmp-dir . \
        --spark-master \'local[*]\'
    gatk --java-options -Xmx${task.memory.toGiga()}g BuildBamIndex \
        --INPUT ${sample_id}.md.bam \
        --OUTPUT ${sample_id}.md.bam.bai \
        --TMP_DIR .
    """
}


(md_ch, indexcov_intermediate_ch) = md_ch.into(2)
// need the bam index only for indexcov/covviz
indexcov_intermediate_ch
    .map { row -> row[2] }
    .set { indexcov_ch }


process alignstats {
    publishDir path: "$outdir/json", pattern: "*.json"

    input:
    set sample_id, file(bam), file(bai) from md_ch

    output:
    file("*.json")

    script:
    """
    alignstats -C -P ${task.cpus} -i $bam -j bam -o ${sample_id}.alignstats.json
    """
}


process indexcov {
    publishDir path: "$outdir/reports/indexcov", mode: "copy"
    label 'covviz'

    input:
    file(bai) from indexcov_ch.collect()
    file faidx

    output:
    file("upward*.png")
    file("*.html")
    file("upward*.bed.gz") into bed_ch
    file("upward*.ped") into indexcov_ped_ch
    file("upward*.roc") into roc_ch

    script:
    excludepatt = params.exclude ? "--excludepatt \"${params.exclude}\"" : ""
    """
    goleft indexcov --sex $sexchroms $excludepatt --directory upward --fai $faidx $bai
    mv upward/* .
    """
}


process covviz {
    publishDir path: "$outdir/reports", mode: "copy", pattern: "*.html"
    label 'covviz'

    input:
    file ped from indexcov_ped_ch
    file bed from bed_ch
    file gff

    output:
    file("covviz_report.html")

    script:
    gff_opt = params.gff ? "--gff ${gff}" : ""
    """
    covviz --min-samples ${params.minsamples} --sex-chroms ${params.sexchroms} --exclude '${params.exclude}' \
        --z-threshold ${params.zthreshold} --distance-threshold ${params.distancethreshold} \
        --slop ${params.slop} --ped ${ped} ${gff_opt} --skip-norm ${bed}
    """
}
