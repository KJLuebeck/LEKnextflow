nextflow.enable.dsl = 2

params.accession = "M21012"
params.outdir = "${baseDir}/out/${params.accession}"
params.fastafiles = "${baseDir}/out/fastafiles"
params.storedir = "${baseDir}/cache"

//
process downloadfasta {
    //publishDir "${params.outdir}", mode: 'copy', overwrite: true
    storeDir "${params.storedir}"
    output:
        path "*.fasta" 
    """
    wget "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id="${params.accession}"&rettype=fasta&retmode=text" -O "${params.accession}".fasta
    """
}
//collect all fasta in one fasta
process collectfiles {
    //publishDir "${params.outdir}", mode: 'copy', overwrite: true
    input:
        path fastafiles
        path fastachannel
    output:
        path "allfasta.fasta"
    """
    cat *.fasta > allfasta.fasta
    """
}

//Mafft aligner
process mafft {
   //publishDir "${params.outdir}", mode: 'copy', overwrite: true 
   container "https://depot.galaxyproject.org/singularity/mafft%3A7.505--hec16e2b_0"
   input:
        path fastafile
   output:
        path "fertigesfasta.fasta"
    """
   mafft --auto allfasta.fasta >fertigesfasta.fasta
    """
}

//trimAl
process trimal {
   publishDir "${params.outdir}", mode: 'copy', overwrite: true 
   container "https://depot.galaxyproject.org/singularity/trimal%3A1.4.1--h9f5acd7_6"
   input:
        path fastafile
   output:
        path "trimmedfasta.fasta"
        path "trimmedfasta.html"

   """
   trimal -in fertigesfasta.fasta -out trimmedfasta.fasta -htmlout trimmedfasta.html -automated1
   """
}

workflow {
    fastachannel = downloadfasta()
    fastafiles = Channel.fromPath("${params.fastafiles}/*.fasta").collect()
    view(fastafiles)
    collectedfasta = collectfiles(fastachannel, fastafiles)
    alignment = mafft(collectedfasta)
    trimal(alignment)

}