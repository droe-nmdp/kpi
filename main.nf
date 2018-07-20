#!/usr/bin/env nextflow

/*
 * kpi
 *
 * KIR structural interpretation for raw reads.
 *
 * Predict a (possibly-amiguous) pair of haplotypes given the
 * presence/absence (PA) genotype from one individual and a collection
 * of PA reference haplotypes.
 * 
 * The input is a directory of 2-column tab-separated unix-formatted text 
 * files. Each file maps IDs to one or more files. The first column of each
 * file is the ID of a set of reads (e.g., representing
 * an individual). The second column is non-space-containing path to a 
 * optionally-gzipped fastq or fasta file.
 *
 * The output contains genotype and haplotype-pair predictions.
 *
 * One individual with 60G of gzippped fastq takes about todo minutes 
 * on an 8 CPU computer with 20G memory.
 * 
 * @author Dave Roe
 */

inputSuffix = "txt"
nfKMCForks = 1 // run this many input text files in parallel
params.input = '/opt/kpi/raw/'
params.output = '/opt/kpi/output'
geneProbes  = '/opt/kpi/input/geneHapSigMarkers_v1-wRc'
nfForks = 4 // run this many input text files in parallel
// input: kmc probe txt files
kmcNameSuffix = '_hits.txt'          // extension on the file name
bin1Suffix = 'bin1'
probeFile = '/opt/kpi/input/geneHapSigMarkers_v1.fasta'
params.haps = '/opt/kpi/input/HapSet18_v2.txt'

mapDir = params.input
resultDir = params.output
haps = params.haps
// things that probably won't change per run
mapDir = params.input
resultDir = params.output
if(!mapDir.trim().endsWith("/")) {
	mapDir += "/"
}
if(!resultDir.trim().endsWith("/")) {
	resultDir += "/"
}
mapPath = mapDir + '*' + inputSuffix

fqsIn = Channel.fromPath(mapPath).ifEmpty { error "cannot find any ${inputSuffix} files in ${mapDir}" }.map { path -> tuple(sample(path), path) }

process probeFastqs {
	container = "droeatnmdp/kpi:latest"
	//publishDir resultDir, mode: 'copy', overwrite: true
    maxForks nfKMCForks

	input: set s, file(f) from fqsIn
	output:
		set s, file('*.kmc_*') into kmcdb
	script:
		"""
        probeFastqsKMC.groovy -m ${f} -o . -w .
		"""
		
} // probeFastqs

process probeDB {
	publishDir resultDir, mode: 'copy', overwrite: true

	input: set s, file(fList) from kmcdb
	output:
		set s, file{ "*_hits.txt"} into filterdb
	
	script:
		"""
        echo ${fList}
        filterMarkersKMC2.groovy -d . -p ${geneProbes} -o . -w .
		"""
		
} // probeFastqs

/*
 * kmc2locusBin
 *
 * Given a kmc output file, bin the hit reads into separate files based on locus.
 * 
 * e.g., ./kmc2Locus2.groovy -j 100a.txt -p kmers.txt -e bin1 -o output
 * 
 * Input files: e.g., 100a.fasta
 * Output files have an extension of 'bin1'.
 * @todo improve the memory usage here
 */
process kmc2locusBin {
  publishDir resultDir, mode: 'copy', overwrite: true
  maxForks nfForks

  input:
    set s, file(hits) from filterdb // left off
  output:
	set s, file{"${s}*.bin1"} into bin1Fastqs

script: 
    // e.g., gonl-100a.fasta
    // todo: document this
	String dataset
	String id
	id = hits.name.replaceFirst(kmcNameSuffix, "")
    """
    kmc2LocusAvg2.groovy -j ${hits} -p ${probeFile} -e ${bin1Suffix} -i ${id} -o .

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
				# todo: change this to kmcNameSuffix
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
  int end = name.indexOf(inputSuffix)
  if ( end <= 0 ) {
    throw new Exception( "Expected file " + name + " to end in '" + inputSuffix + "'" );
  }
  end = end -1 // Remove the trailing '.'
  return name.substring(start, end)
} // sample
