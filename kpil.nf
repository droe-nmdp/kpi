#!/usr/bin/env nextflow

/*
 * kpil
 *
 * KIR structural interpretation for long reads. Beta version.
 *
 * usage: docker run --rm -v <pathToFastqGz>:/opt/kpi/raw/ -v <pathToOutput>:/opt/kpi/output/ droeatnmdp/kpi /opt/kpi/kpil.nf
 * 
 * Predict a (possibly-amiguous) pair of haplotypes given the
 * presence/absence (PA) genotype from one individual and a collection
 * of PA reference haplotypes.
 *
 * Individually processes all ".fastq" files in /opt/kpi/raw/. 
 * Mount your local folder here.
 * 
 * Outputs 3 interpretation files for each fastq file in /opt/kpi/output/. 
 * Mount your local folder here.
 * Three files:
 *   1. genotype: PA call for each gene
 *   2. reads: PA call for each read (middle 50-500 bases of each gene)
 *   3. interp: haplotype pair prediction
 *
 * @author Dave Roe
 * @todo By default, BBDuk has maskmiddle and rcomp set to true
 */

// things that may change per run
// here are the FASTQ files
fqNameSuffix = 'fastq.gz'          // extension on the file name
//fqNameSuffix = 'fasta'          // extension on the file name
fqDir = '/opt/kpi/raw/'
resultDir = '/opt/kpi/output'
// input: kmc probe txt files
kmcNameSuffix = 'fasta'          // extension on the file name

// things that probably won't change per run
fqPath = fqDir + '*.' + fqNameSuffix
probeFile = '/opt/kpi/input/geneHapSigMarkers_v1.txt'
haps = '/opt/kpi/tmp/input/HapSet18_v2.txt'
unmappedInd = "unmapped" // add 'unmapped' meta locus

fqs1 = Channel.fromPath(fqPath).ifEmpty { error "cannot find any fastq files matching ${fqPath}" }.map { path -> tuple(sample(path), path) }
fqs2 = Channel.fromPath(fqPath).ifEmpty { error "cannot find any fastq files matching ${fqPath}" }.map { path -> tuple(sample(path), path) }

/*
 * fq2locusBin
 *
 * Given a FASTQ file, bin the reads into separate files based on probe markers.
 * 
 * eval "bbduk.sh in=KP420443_KP420444.bwa.read1_short.fastq outm=2DL5.bin1 literal=TTGAACCCTCCATCACAGGTCCTGG k=25 maskmiddle=f"
 * 
 * Output files have an extension of 'bin1'.
 * @todo change extension to 'bin2.fastq'
 * @todo Remove the zero-length .bin1 files?
 */
process fq2locusBin {
  publishDir = resultDir
  input:
    set s, file(fq) from fqs1
  output:
    file{"*.bin1"} into bin1Fastqs

    """
    probeBinFastqs.groovy -i ${fq} -p ${probeFile} -s bin1
    """
} // fq2locusBin

/*
 * kmc2locusBin
 *
 * Given a kmc output file, bin the hit reads into separate files based on lcus.
 * 
 * e.g., ./kmc2Locus.groovy -j gonl-100a.txt -p kmers.txt -e txt -o output
 * 
 * Input files: e.g., gonl-100a.fasta
 * Output files have an extension of 'bin1'.
 * @todo change extension to 'bin2.txt', etc.
 */
process kmc2locusBin {
  publishDir resultDir, mode: 'copy', overwrite: true

  input:
    set s, file(kmc) from fqs1
  output:
    file{"${id}*.bin1"} into bin1Fastqs
    file{"${id}*.vbin1"} into vbin1Fastqs
    val id into idChannel
    val vid into vidChannel
//todo: add vilches interp output
  maxForks 1

script: 
    // e.g., gonl-100a.fasta
    // todo: document this
    (dataset, id) = kmc.name.replaceFirst(kmcNameSuffix, "").replaceFirst("\\.fasta","").split('-')
    id = id[0..-2] // remove the trailing dot
    vid = id
    """
    kmc2Locus.groovy -j ${kmc} -p ${probeFile} -e ${bin1Suffix} -j ${kmcNameSuffix}  -i ${id} -v ${vilchesFileName} -o .
    """
} // kmc2locusBin

/*
 * locusBin2ExtendedLocusBin
 * 
 * 1) Makes haplotype predictions from PA probes.
 * 2) For each gene, combine its bordering intergene reads into a single bin.
 *
 * Output files have an extension of 'bin2'.
 *
 * @todo old:   binLoci.groovy -h ${haps} -q \$fileList -p \$outFile
 */
process locusBin2ExtendedLocusBin {
  publishDir = resultDir
    // todo: add a set here?
  input:
    file(b1List) from bin1Fastqs.collect()
  output:
    file{"prediction.txt"} into predictionChannel
    file{"*.lw.fasta"} into lwChannel

    """
    outFile="prediction.txt"
    fileList=""
    ext="*bin1*"
    for bFile in $b1List; do
        if [ -s \$bFile ]; then
            if [[ \$bFile == \$ext ]]; then
                echo \$bFile
                if [ "\$fileList" == "" ]; then
                    :
                else
                    fileList+=","
                fi
                fileList+=\$bFile
            fi
        fi
    done
    echo \$fileList

    pa2Haps.groovy -h ${haps} -q \$fileList -o \$outFile
#todo: next line (change to vprediction.txt?)
#todo    assembleLocusWindow.groovy -p prediction.txt
    """
} // locusBin2ExtendedLocusBin

/*
 * kmc2locusBin
 *
 * Given a kmc output file, bin the hit reads into separate files based on lcus.
 * 
 * e.g., ./kmc2Locus.groovy -j gonl-100a.txt -p kmers.txt -e txt -o output
 * 
 * Input files: e.g., gonl-100a.fasta
 * Output files have an extension of 'bin1'.
 * @todo change extension to 'bin2.txt', etc.
 */
process kmc2locusBin {
  publishDir resultDir, mode: 'copy', overwrite: true

  input:
    set s, file(kmc) from kmcs1
  output:
    file{"${id}*.bin1"} into bin1Fastqs
    file{"${id}*.vbin1"} into vbin1Fastqs
    val id into idChannel
    val vid into vidChannel
//todo: add vilches interp output
  maxForks 1

script: 
    // e.g., gonl-100a.fasta
    // todo: document this
    (dataset, id) = kmc.name.replaceFirst(kmcNameSuffix, "").replaceFirst("\\.fasta","").split('-')
    id = id[0..-2] // remove the trailing dot
    vid = id
    """
    kmc2Locus.groovy -j ${kmc} -p ${probeFile} -e ${bin1Suffix} -j ${kmcNameSuffix}  -i ${id} -v ${vilchesFileName} -o .
    """
} // kmc2locusBin

/* 
 * assembleLocusWindow
 *
 * Assemble each unique combination of gene and its neighboring
 * genes.
 * 
 *
process assembleLocusWindow {
    publishDir resultDir, mode: 'copy', overwrite: 'true'
  input:
    file(prediction) from predictionChannel
    file(b1f2) from bin1Fastqs2
  output:
    file{"*.lw.fasta"} into finalAssembly

    """
    assembleLocusWindow.groovy -p ${prediction}
    """
} // assembleLocusWindowb
*/
/* 
 * unitigs2Final
 *
 * Second, final, level Canu assembly.
 *
 *
process unitigs2Final {
    publishDir resultDir, mode: 'copy', overwrite: 'true'  //testing(todo)
  input:
    file(b3List) from unitigFastqs.collect()
  output:
    file{"assembly/*.fasta"} into assembly
//put bak(todo)    file{"assembly/*.unitigs.fasta"} into assembly

    """
    fname='all_bin3.fasta'
    name_intervening='intervening.fasta'
    for bFile in $b3List; do
        if [ -s \$bFile ]; then
            if [ "\$bFile" == *"3DP1"* ] || [ "\$bFile" == *"2DL4"* ]; then
                cat \$bFile >> \$name_intervening
            else
                cat \$bFile >> \$fname
            fi
        fi
    done
    # add unmapped reads to both bins
#todo    cat 'unmapped.unitigs.fasta' >> \$fname
#todo    cat 'unmapped.unitigs.fasta' >> \$name_intervening

    # assemble non-intervening regions
    title='assembly'
    echo canu maxMemory=8 -p \$title -d \$title genomeSize=200k -pacbio-corrected \$fname
    canu maxMemory=8 -p \$title -d \$title genomeSize=200k -pacbio-corrected \$fname

    # assemble intervening region
    title='assembly_intervening'
    echo canu maxMemory=8 -p \$title -d \$title genomeSize=15k -pacbio-corrected \$name_intervening
    canu maxMemory=8 -p \$title -d \$title genomeSize=15k -pacbio-corrected \$name_intervening
    """
} // unitigs2Final
*/

// get the per-sample name
def sample(Path path) {
  def name = path.getFileName().toString()
  int start = Math.max(0, name.lastIndexOf('/'))
  int end = name.indexOf(fqNameSuffix)
  if ( end <= 0 ) {
    throw new Exception( "Expected file " + name + " to end in '" + fqNameSuffix + "'" );
  }
  end = end -1 // Remove the trailing '.'
  return name.substring(start, end)
} // sample

workflow.onComplete {
  println "DONE: ${ workflow.success ? 'OK' : 'FAILED' }"
}

workflow.onError {
  println "ERROR: ${workflow.errorReport.toString()}"
}
