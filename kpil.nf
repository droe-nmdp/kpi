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
 * @todo query isn't working; build the sbt binary (currently binary)
 * @author Dave Roe
 */

// things that may change per run
// here are the FASTQ files
fqNameSuffix = 'fastq'          // extension on the file name
fqDir = '/opt/kpi/raw/'  //todo
resultDir = '/opt/kpi/tmp/output' //todo: remove the 'tmp'
hashFile = resultDir + '/bt25.hash'
probeFile = '/opt/kpi/input/locus-hap_probes_v1.txt' //todo: replace 'tmp' with 'input'
probeFileFasta = '/opt/kpi/tmp/locus-hap_probes_v1.fasta' //todo: replace 'tmp' with 'input'
haps = '/opt/kpi/input/all_haps_v2.txt'
unmappedInd = "unmapped" // add 'unmapped' meta locus
/*
probes = '/opt/kpi/input/vilches_probes.txt'
gProbes = '/opt/kpi/input/geraghty_probes.txt'
haps = '/opt/kpi/input/all_haps_v1.txt'
*/
/*testing
fqDir = '/Users/droe/data/kir/simulations/KP420443_KP420444/'
resultDir = '/Users/droe/src/docker/kpi/output'
probes = '/Users/droe/src/docker/kpi/input/vilches_probes.txt'
gProbes = '/Users/droe/src/docker/kpi/input/geraghty_probes.txt'
haps = '/Users/droe/src/docker/kpi/input/all_haps_v1.txt'
*/

// things that probably won't change per run
fqPath = fqDir + '*.' + fqNameSuffix
maxMem = "-Xms5g"
// bt count (bq2bf) info
cpuThreads = 2  // todo
cutoff = 5
bf_size = 800000000
tabArg = "\$'\t'"

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


// output is output/<file>_prediction.txt
// todo: this repeats (although cached) for each item in bin1Fastqs
// todo: input fastq file and make unmapped file

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

/* left off: assemble each bin2
process bin22bin3 {
  input:
    file(b2List) from bin2Fastqs.collect()

    """

    """
} // bin22bin3
*/
/*
 * Convert each sequence file to a bloom filter count file.
 * Works on a per-read basis due to the split in the previous step.
 
process fq2bf {
  publishDir = resultDir
  input:
    set s, file(f) from perReadFastqs
  output:
    file{"${f}.bf.bv"} into bv

    """
    bt count --cutoff ${cutoff} --threads ${cpuThreads} ${hashFile} ${bf_size} ${f} ${f}.bf.bv;
    """
} // fq2bf
*/
/*
 * Create the final bloomtree database.
 *
 * todo: change to this https://www.nextflow.io/docs/latest/faq.html#how-do-i-iterate-over-nth-files-them-within-a-process

process btBuild {
  input:
    file bvList from bv.toList()

    """
    name=bt25
    if [ ! -f ${hashFile} ];
    then
       echo "Hash file has not been created"
       exit 1
    fi
    bvFiles=${resultDir}/bvFiles.txt
    if [ -f \$bvFiles ]; then rm -f \$bvFiles;fi
    for f in $bvList;do echo \$f >> \$bvFiles;done
    
    sbuild -f ${hashFile} -p ${resultDir} -l \$bvFiles -o ${resultDir} -b \$name
    """
} // btBuild
 */
/*
 * pa2Haps
 *
 * Fit PA genotypes to (potentially ambiguous) haplotype pairs.
 *
process pa2Haps {
  input:
    set s, file(fg), file(fr) from paInterp
  output:
    set s, file{"${s}_interp.txt"} into hapInterp

    """
    pa2Haps.groovy ${haps} ${fg} ${s}_interp.txt
    cp -f ${s}_interp.txt ${resultDir}
    """
} // pa2Haps
 */

// todo(remove)
//fq2I = Channel.fromFilePairs('${resultDir}/*[.${fqNameSuffix}|_interp.txt]')

/**
 * fq2IGenotype
 *
 * Make inter-gene presence/absence calls per read and summarized.
process fq2IGenotype {
  input:
    set s, file(h) from hapInterp
//    set s, file(fq), file(h) from fq2I
//        set s, file(fq) from fqs1, file(h) from hapInterp
//    file(fi) from fqs1
  output:
    set s, file{"${s}_igenotype.txt"}, file{"${s}_ireads.txt"} into iInterp

    """
    geraghtyInterp.groovy ${gProbes} ${h} ${fqDir}/${s}.${fqNameSuffix} ${s} ${s}_igenotype.txt ${s}_ireads.txt
    cp -f ${s}_igenotype.txt ${s}_ireads.txt ${resultDir}
    """
} // fq2IGenotype
 */

/**
 * combineReadInterps
 *
 * Combine gene and inter-gene presence/absence calls per read.
 *
process combineReadInterps {
  input:
    set s, file(igen), file(ireads) from iInterp
  output:
    set s, file{"${s}_fullInterp.txt"} into fullInterp

    """
    join -t ${tabArg} -1 1 -2 1 ${resultDir}/${s}_reads.txt ${s}_ireads.txt > ${s}_fullInterp.txt
    cp -f ${s}_fullInterp.txt ${resultDir}
    rm -f ${resultDir}/${s}_reads.txt ${resultDir}/${s}_ireads.txt
    """
} // combineReadInterps
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
