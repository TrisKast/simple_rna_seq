params.reads = ""
params.transcriptome = ""
params.outdir = "results"
params.multiqc = "$baseDir/multiqc"

log.info """\
 R N A S E Q - N F   P I P E L I N E
 ===================================
 transcriptome: ${params.transcriptome}
 index        : ${params.index}
 reads        : ${params.reads}
 outdir       : ${params.outdir}
 """



multiqc_file = file(params.multiqc)

Channel
    .fromFilePairs( params.reads )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .into { read_pairs_ch; read_pairs2_ch }

if (params.index){

  Channel
      .fromPath( params.index )
      .ifEmpty { error "Cannot find any index matching: ${params.index}" }
      .into {index_ch}

} else {

    transcriptome_file = file(params.transcriptome)

    process index {
        tag "$transcriptome.simpleName"

        input:
        file transcriptome from transcriptome_file

        output:
        file 'index' into index_ch

        script:
        """
        salmon index --threads $task.cpus -t $transcriptome -i index
        """
    }

}


process quant {
    tag "$pair_id"

    input:
    file index from index_ch
    set pair_id, file(reads) from read_pairs_ch

    output:
    file(pair_id) into quant_ch

    script:
    """
    salmon quant --threads 20 --libType=A -i $index -1 ${reads[0]} -2 ${reads[1]} -o $pair_id
    """
}

process fastqc {
    tag "FASTQC on $sample_id"
    publishDir params.outdir

    input:
    set sample_id, file(reads) from read_pairs2_ch

    output:
    file("fastqc_${sample_id}_logs") into fastqc_ch


    script:
    """
    mkdir fastqc_${sample_id}_logs
    fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads}
    """
}


process multiqc {
    publishDir params.outdir, mode:'copy'

    input:
    file('*') from quant_ch.mix(fastqc_ch).collect()
    file(config) from multiqc_file

    output:
    file('multiqc_report.html')

    script:
    """
    cp $config/* .
    multiqc .
    """
}

workflow.onComplete {
	log.info ( workflow.success ? "\nDone! Open the following report in your browser --> $params.outdir/multiqc_report.html\n" : "Oops .. something went wrong" )
}
