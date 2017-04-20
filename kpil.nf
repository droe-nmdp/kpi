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
 * @author Dave Roe
 */

// things that may change per run
// here are the FASTQ files
fqNameSuffix = 'fastq'          // extension on the file name
fqDir = '/opt/kpi/raw/'  //todo?
resultDir = '/opt/kpi/output'

// things that probably won't change per run
fqPath = fqDir + '*.' + fqNameSuffix
probeFile = '/opt/kpi/input/locus-hap_probes_v1.txt'
haps = '/opt/kpi/input/all_haps_v2.txt'
unmappedInd = "unmapped" // add 'unmapped' meta locus

fqs1 = Channel.fromPath(fqPath).ifEmpty { error "cannot find any fastq files matching ${fqPath}" }.map { path -> tuple(sample(path), path) }
fqs2 = Channel.fromPath(fqPath).ifEmpty { error "cannot find any fastq files matching ${fqPath}" }.map { path -> tuple(sample(path), path) }

/*
 * fq2locusBin
 *
 * Given a FASTQ file, bin the reads into separate files based on probe markers.
 * 
 * eval "bbduk.sh in=KP420443_KP420444.bwa.read1_short.fastq outm=2DL1.bin1 literal=CTGAACCCACCAGCACAGGTCCTGG k=25 maskmiddle=f"
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
 * locusBin2ExtendedLocusBin
 * 
 * 1) Makes haplotype predictions from PA probes.
 * 2) For each gene, combine its bordering intergene reads into a single bin.
 *
 * Output files have an extension of 'bin1'.
 * @todo change extension to 'bin2.fastq'.
 */
process locusBin2ExtendedLocusBin {
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
} // locusBin2ExtendedLocusBin


/* 
 * extendedLocusBin2Unitigs
 *
 * Use Canu to assemble each extended locus bin into unitigs.
 * 
 */
process extendedLocusBin2Unitigs {
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
    cp unmapped.bin2 unmapped.unitigs.fasta
    """
} // extendedLocusBin2Unitigsb

/* 
 * unitigs2Final
 *
 * Second, final, level Canu assembly.
 *
 */
process unitigs2Final {
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
    cat 'unmapped.unitigs.fasta' >> \$fname
    title='assembly'
    echo canu maxMemory=8 -p \$title -d \$title genomeSize=200k -pacbio-corrected \$fname
    canu maxMemory=8 -p \$title -d \$title genomeSize=200k -pacbio-corrected \$fname
    """
} // unitigs2Final

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
