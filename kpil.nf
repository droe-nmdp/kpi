#!/usr/bin/env nextflow

/*
 * kpil
 *
 * KIR structural interpretation for long reads. Beta version.
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
 * @todo 
 * @author Dave Roe
 */

// things that may change per run
// here are the FASTQ files
fqNameSuffix = 'fastq'          // extension on the file name
fqDir = '/opt/kpi/raw/'  //todo
resultDir = '/opt/kpi/output'

// things that probably won't change per run
fqPath = fqDir + '*.' + fqNameSuffix
probeFile = '/opt/kpi/input/locus-hap_probes_v1.txt'
probeFileFasta = '/opt/kpi/input/locus-hap_probes_v1.fasta'
haps = '/opt/kpi/input/all_haps_v2.txt'
unmappedInd = "unmapped" // add 'unmapped' meta locus

fqs1 = Channel.fromPath(fqPath).ifEmpty { error "cannot find any fastq files matching ${fqPath}" }.map { path -> tuple(sample(path), path) }
fqs2 = Channel.fromPath(fqPath).ifEmpty { error "cannot find any fastq files matching ${fqPath}" }.map { path -> tuple(sample(path), path) }

//probeChannel = Channel.create()
probeChannel = []
pf = file(probeFile)
pf.readLines().each { line ->
    (gl, probe) = line.split('\t')
//    System.err.println "adding ${gl}" //todo
//        probeChannel.add([gl, probe])
//    probeChannel[gl] = probe
//    probeChannel.add(line)
//    str = "outm=${gl}.bin1 literal=${probe}"
    str = gl
    probeChannel.add(str)
}
probeChannel.add(unmappedInd)

/*
process fq2NoHitBin1 {
  input:
    set s, file(fq) from fqs1
  output:
    file{"${s}_unmapped.bin1"} into bin1NoHitFastq

    """
    eval "bbduk.sh in=${fq} out=${s}_unmapped.bin1 ref=${probeFileFasta} k=25 maskmiddle=f"

    """
}
*/
/*
 * splitFastq
 *
 * Split the fastq file into one line per read. That's the 
 * way the sbt needs it to report per read.
 * 
 * eval "bbduk.sh in=KP420443_KP420444.bwa.read1_short.fastq outm=2DL1.bin1 literal=CTGAACCCACCAGCACAGGTCCTGG k=25 maskmiddle=f"
 * 
 * Remove the zero-length .bin1 files
 * 
 * Create the hash file here. Downstream may lead to multiple
 * processes trying to create the file.
 * @todo update documentation
 */
process fq2bin1 {
  input:
    set s, file(fq) from fqs1
    each gl from probeChannel.collect()
  output:
    file{"*.bin1"} into bin1Fastqs

    """
    if [ "${gl}" == "${unmappedInd}" ]; then
        bbduk.sh in=${fq} out=${gl}.bin1 ref=${probeFileFasta} k=25 maskmiddle=f overwrite=t
    else
        bbduk.sh in=${fq} outm=${gl}.bin1 literal=${probe} k=25 maskmiddle=f overwrite=t
    fi
    """
} // fq2bin1


/*
 * bin12bin2
 * 
 * Takes the probe-derived bin1 files and outputs the haplotype predictions
 * and *.bin2 files to be assembled.
 */
process bin12bin2 {
  publishDir = resultDir
    // todo: add a set here?
  input:
    file(b1List) from bin1Fastqs.collect()
  output:
//todo(next step?)    file{"*_bin2.fastq"} into bin2Fastqs
    file{"prediction.txt"} into prediction
    file{"*.bin2"} into bin2Fastqs

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
    binLoci.groovy -h ${haps} -q \$fileList -p \$outFile
    """
} // bin12bin2


/* 
 * assemble bin2 files into unitig files using Canu
 * 
 */
process bin22bin3 {
    publishDir resultDir, mode: 'copy', overwrite: 'true'  //testing(todo)
  input:
    file(b2List) from bin2Fastqs.collect()
  output:
    file{"*/*.unitigs.fasta"} into unitigFastqs

    """
    for bFile in $b2List; do
        if [ -s \$bFile ]; then
            name=\$(echo \$bFile | cut -f 1 -d '.')
            echo canu maxMemory=8 -p \$name -d \$name genomeSize=20k -pacbio-corrected \$bFile
            canu maxMemory=8 -p \$name -d \$name genomeSize=20k -pacbio-corrected \$bFile
        fi
    done
    """
} // bin22bin3

/* 
 * bin32bin4
 *
 * Second level Canu assembly.
 */
process bin32bin4 {
    publishDir resultDir, mode: 'copy', overwrite: 'true'  //testing(todo)
  input:
    file(b3List) from unitigFastqs.collect()
  output:
    file{"assembly/*.unitigs.fasta"} into assembly

    """
    fname='all_bin3.fasta'
    for bFile in $b3List; do
        if [ -s \$bFile ]; then
            cat \$bFile >> \$fname
        fi
    done
    title='assembly'
    echo canu maxMemory=8 -p \$title -d \$title genomeSize=200k -pacbio-corrected \$fname
    canu maxMemory=8 -p \$title -d \$title genomeSize=200k -pacbio-corrected \$fname
    """
} // bin32bin4

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
