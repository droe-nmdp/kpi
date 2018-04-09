#!/usr/bin/env nextflow

/*
 * Given 25mer hit counts, interpret KIR genotypes and
 * haplotypes.
 *
 * Input files end with kmcNameSuffix (default '_hit.txt') and are
 * located in kmcDir.
 *
 * Output files end with '_prediction.txt' in resultDir.
 *
 * @author Dave Roe
 */

// input: kmc probe txt files
kmcNameSuffix = '_hits.txt'          // extension on the file name
kmcDir = '/opt/kpi/output/'
bin1Suffix = 'bin1'
resultDir = '/opt/kpi/output/'
probeFile = '/opt/kpi/input/geneHapSigMarkers_v1.fasta'
haps = '/opt/kpi/input/HapSet18_v2.txt'

// things that probably won't change per run
kmcPath = kmcDir + '*' + kmcNameSuffix

kmcs1 = Channel.fromPath(kmcPath).ifEmpty { error "cannot find any ${kmcNameSuffix} files in ${kmcPath}" }.map { path -> tuple(sample(path), path) }
kmcs2 = Channel.fromPath(kmcPath).ifEmpty { error "cannot find any ${kmcNameSuffix} files in ${kmcPath}" }.map { path -> tuple(sample(path), path) }

/*
 * kmc2locusBin
 *
 * Given a kmc output file, bin the hit reads into separate files based on locus.
 * 
 * e.g., ./kmc2Locus2.groovy -j 100a.txt -p kmers.txt -e bin1 -o output
 * 
 * Input files: e.g., 100a.fasta
 * Output files have an extension of 'bin1'.
 */
process kmc2locusBin {
  //publishDir resultDir, mode: 'copy', overwrite: true

  input:
    set s, file(kmc) from kmcs1
  output:
	set s, file{"${s}*.bin1"} into bin1Fastqs

script: 
    // e.g., gonl-100a.fasta
    // todo: document this
	String dataset
	String id
	id = kmc.name.replaceFirst(kmcNameSuffix, "")
	String retId = id
    """
    kmc2LocusAvg2.groovy -j ${kmc} -p ${probeFile} -e ${bin1Suffix} -i ${id} -o .

if ls *.bin1 1> /dev/null 2>&1; then
    : # noop
else
    echo "snp: " > "${id}_prediction.txt"

    touch "${id}_uninterpretable.bin1"
fi
    """
} // kmc2locusBin

/*
 * locusBin2ExtendedLocusBin
 * 
 * 1) Makes haplotype predictions from PA probes.
 *
 *
 * @todo document
 */
process locusBin2ExtendedLocusBin {
  publishDir resultDir, mode: 'copy', overwrite: true
  input:
		set s, file(b1List) from bin1Fastqs
  output:
	set s, file{"*_prediction.txt"} into predictionChannel

  script:
    """
    FILES="${s}*.bin1"
    outFile="${s}"
    outFile+="_prediction.txt"
    fileList=""
    id=""
    ext="*bin1*"
    for bFile in \$FILES; do
        if [ -s \$bFile ]; then
            if [[ \$bFile == \$ext ]]; then
                id=\$(basename "\$bFile")
                # '%' Means start to remove after the next character;
                # two means non-greedy
                id="\${id%%_*}"
                if [ "\$id" == "" ]; then
                    id=\$(basename "\$bFile")
                fi
                # echo \$bFile
                if [ "\$fileList" == "" ]; then
                    :
                else
                    fileList+=","
                fi
                fileList+=\$bFile
                outFile="\${id}"
                outFile+="_prediction.txt"
            fi
        fi
    done
    pa2Haps3.groovy -h ${haps} -q "\$fileList" -o "\$outFile"
    """
} // locusBin2ExtendedLocusBin

// get the per-sample name
def sample(Path path) {
  def name = path.getFileName().toString()
  int start = Math.max(0, name.lastIndexOf('/'))
  int end = name.indexOf(kmcNameSuffix)
  if ( end <= 0 ) {
    throw new Exception( "Expected file " + name + " to end in '" + kmcNameSuffix + "'" );
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