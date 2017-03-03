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
fqDir = '/opt/kpi/raw/'
resultDir = '/opt/kpi/output'
probes = '/opt/kpi/input/vilches_probes.txt'
gProbes = '/opt/kpi/input/geraghty_probes.txt'
haps = '/opt/kpi/input/all_haps_v1.txt'

/* testing
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
cpuThreads = 15
cutoff = 5
bf_size = 800000000
tabArg = "\$'\t'"

fqs1 = Channel.fromPath(fqPath).ifEmpty { error "cannot find any fastq files matching ${fqPath}" }.map { path -> tuple(sample(path), path) }
fqs2 = Channel.fromPath(fqPath).ifEmpty { error "cannot find any fastq files matching ${fqPath}" }.map { path -> tuple(sample(path), path) }

/*
 * fq2Genotype
 *
 * Make gene presence/absence calls per read and summarized.
 */
process fq2Genotype {
  input:
    set s, file(fq) from fqs1
  output:
    set s, file{"${s}_genotype.txt"}, file{"${s}_reads.txt"} into paInterp

    """
    vilchesInterp.groovy ${maxMem} ${probes} ${fq} ${s} ${s}_genotype.txt ${s}_reads.txt
    cp -f ${s}_genotype.txt ${s}_reads.txt ${resultDir}
    """
} // fq2Genotype

/*
 * Convert each sequence file to a bloom filter count file.
 * Works on an individual basis: one tree per individual.
 */
process fq2bf {
  input:
    set s, file(f) from fqs2
  output:
    file{"${s}.bf.bv"} into bv

    """
    hashFile=${resultDir}/bt25.hash
    if [ ! -f "\$hashFile" ];
    then
       echo "Creating hash file"
       bt hashes --k 25 \$hashFile 1
       echo "asd" >> /tmp/a
    fi
    bt count --cutoff ${cutoff} --threads ${cpuThreads} \$hashFile ${bf_size} ${f} ${s}.bf.bv

    """

} // fq2bf

/*
 * Create the final bloomtree database.
 */
process btBuild {
  input:
    file bvList from bv.toList()
    """
    name=bt25
    hashFile=${resultDir}/\$name.hash
    if [ ! -f "\$hashFile" ];
    then
       echo "Hash file has not been created"
       exit 1
    fi
    bvFiles=${resultDir}/bvFiles.txt
    if [ -f \$bvFiles ]; then rm -f \$bvFiles;fi
    for f in $bvList;do echo \$f >> \$bvFiles;done

    bt build \$hashFile \$bvFiles \$name.bt
    bt compress \$name.bt \$name.btz
    cp \$name.bt ${resultDir}
    cp \$name.btz ${resultDir}
    """

} // btBuild

/*
 * pa2Haps
 *
 * Fit PA genotypes to (potentially ambiguous) haplotype pairs.
 *
 */
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

// todo(remove)
//fq2I = Channel.fromFilePairs('${resultDir}/*[.${fqNameSuffix}|_interp.txt]')

/**
 * fq2IGenotype
 *
 * Make inter-gene presence/absence calls per read and summarized.
 */
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

/**
 * combineReadInterps
 *
 * Combine gene and inter-gene presence/absence calls per read.
 *
 */
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
