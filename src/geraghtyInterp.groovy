#!/usr/bin/env groovy

/*
 * Make inter-gene calls on long reads using the genotyping
 * algorithm by Geraghty's team and described in 
 * Supplemental of Roe, et al. 2017.
 *
 * Input
 *   1. Geraghty probe pairs
 *   2. Haplotype pair predictions from pa2Haps.groovy
 *   3. A single fastq file
 *   4. The ID for the individual
 *
 * Output
 *   1. Genotype calls
 *   2. Read annotation
 * 
 * Load the probe pairs file. Assumes the ambiguous IUPAC codes have been resolved.
 * Loop through every read.
 * For each read, check each pair of primers.
 * Output negatives and positives
 *
 * Run via, e.g., 
 *   ./geraghtyInterp.groovy input/geraghty_probes.txt $HOME/src/groovy/kir/output/KP420443_KP420444_haps.txt $HOME/data/kir/simulations/KP420443_KP420444/KP420443_KP420444.bwa.read1.fastq KP420443_KP420444 ~/repos/dev/bioinformatics/projects/ngs/scripts/groovy/output/KP420443_KP420444_genotype.txt output/KP420443_KP420444_reads.txt
 *
 * http://biojava.org/wiki/BioJava:CookBook3:FASTQ
 * Requires
 *   biojava3-core.jar
 *   biojava3-sequencing.jar
 *   http://code.google.com/p/guava-libraries/
 * @author Dave Roe
 * @todo add command line handling
 */

import com.google.common.io.*
import java.util.*
import groovy.io.*
import java.util.zip.GZIPInputStream

// things that may change per run
debugging = false
// bp is the total distance between the two primers, including the primers
bp = 100 
distWindow = 5 // allow this much cushion on each side for distance between primers

// things that probably won't change per run
primerFileName = args[0]
hapsFileName = args[1]
fastqFile = args[2]
id = args[3]
outFileName = args[4]
outReadFileName = args[5]
err = System.err
DEFAULT_BUFFER_SIZE=100000
featureSeparator = "-" // feature separator in hapsFileName

// Map: locus -> array of primers for that locus
// each primer is an Expando that has: locus, fName, fSeq, rName, rSeq, bp
// bp is the total distance between the two primers, including the primers
primers = [:]
primers = readPrimers(primerFileName, bp)

// the file containing the haplotype pair calls
// hap name -> Set of genes/features
HashMap<String,HashSet<String>> hapCalls = readHapCalls(hapsFileName)

// open output files and print the headers
// the genotype output
outWriter = new PrintWriter(new File(outFileName).newOutputStream(), true)
outWriter.println "locus\tP/A"
// the reads output
outReadWriter = new PrintWriter(new File(outReadFileName).newOutputStream(), true)
outReadWriter.println "read\tlocus list\tlength"

if(debugging) { 
    err.println "using distance of ${bp} with a window of ${distWindow}"
}

// collection of gene names
TreeSet<String> genotype = new TreeSet() 
// read name -> collection of loci names
HashMap<String, TreeSet<String>> readHits = [:]
// read name -> length of sequence
HashMap<String, Integer> sizeMap = new HashMap()
// loop through the fastq files
def fq = new File(fastqFile)
p = fq.getParent()
if(p == null) {
    p = "./"
}
fullFileName = p + "/" + fq.getName()
//if(debugging) { 
err.println "processing ${fullFileName}...\n"
//}

BufferedReader r;
if(fq.getName().endsWith("fastq.gz")) {
    r=new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(fullFileName),DEFAULT_BUFFER_SIZE)),DEFAULT_BUFFER_SIZE)
} else {
    r=new BufferedReader(new FileReader(fullFileName),DEFAULT_BUFFER_SIZE)
}
long nLine=0L
String line;
String readName = null // read name
while((line=r.readLine())!=null) {        // each line of sequence file
    nLine++
    //err.println "${(nLine%4)} ${line}"//debugging
    if(nLine%4==1) { // sequence header on line 1/4
        readName = line.trim()
    } // read name
    if(nLine%4!=2) continue; // sequences on line 2/4
    // skip past the 3 non-sequence rows
    seq = line.trim().toUpperCase()
    if(debugging) { 
        err.println "nLine=${nLine}, read name = ${readName}, seq=${seq}"
    }
    sizeMap[readName] = seq.length()
    readSet = new TreeSet() // set of features (e.g., genes) per read
    (genotype, readSet) = interpretRead(primers, seq, genotype, readName,
                                        readSet, hapCalls)

    readHits[readName] = readSet
} // each line in fastq file

//if(debugging) { 
    err.println "done interpreting: ${readHits.keySet().size()} reads\n"
    err.println "genotype size=" + genotype.size()
//}
    
// output the genotype
primers.keySet().each { locus ->
    paStr = "N"
    if(genotype.contains(locus)) { 
        paStr = "Y"
    }
    outWriter.println "${locus}\t${paStr}"
}

if(readHits.asBoolean()) {   // write the annotation of the reads
    readHits.keySet().sort().each { id ->
        locusList = readHits[id]
        ll = locusList.join("/")
        size = sizeMap[id]
        outReadWriter.println "${id}\t${ll}\t${size}"
    }
}

outWriter.close()
outReadWriter.close()
// end main

/*
 * Interpret a read in the context of the primers.
 * 
 * If a hit is made, add it to genotype 
 *
 * @param readName the name of the read
 * @param hapCalls hap name -> Set of genes/features; if null, skip this step
 * @return modified input objects genotype and readSet
 * @todo finish documentation
 */
def ArrayList<Map> interpretRead(Map primers, String sequence,
                                 TreeSet<String> genotype, String readName,
                                 TreeSet<String> readSet,
                                 HashMap<String,HashSet<String>> hapCalls) {
    primers.each { primerLocus, primerList -> // Expando
        if(debugging) { 
            err.println "analyzing primers for ${primerLocus}"
        }
        primerList.each { primer ->
/*            if(primer.fName != "Fcon1254") {  //debugging (droe)
              return
            }
*/
            if(debugging) {
                //err.println "sequence=${sequence}"
                err.println "primer pair ${primer.locus} ${primer.fName} ${primer.rName}"
                err.println "${primer.fName}=${primer.fSeq}"
                err.println "${primer.rName}=${primer.rSeq}"
            }
            // check for existance in either direction
            fSeqIdx = sequence.indexOf(primer.fSeq)
            rSeqIdx = sequence.indexOf(reverseComplement(primer.rSeq))
            if((fSeqIdx < 0) || (rSeqIdx < 0)) {
                // try the other direction
                fSeqIdx = sequence.indexOf(reverseComplement(primer.fSeq))
                rSeqIdx = sequence.indexOf(primer.rSeq)
                if(debugging) { 
                    if((fSeqIdx > 0) && (rSeqIdx > 0)) {
                        err.println "found probes (direction 1)\t${primer.locus}"
                    }
                }
            } else {
                if(debugging) { 
                    if((fSeqIdx > 0) && (rSeqIdx > 0)) {
                        err.println "found probes (direction 2)\t${primer.locus}"
                    }
                }
            }

            if(debugging) { 
                err.println "rSeqIdx=${rSeqIdx}, fSeqIdx=${fSeqIdx}, rSeq length=${primer.rSeq.length()}, fSeq length=${primer.fSeq.length()}"
                err.println "primer.bp=${primer.bp}"
            }
            if((fSeqIdx > 0) && (rSeqIdx > 0)) { // if both primers hit
                // check distance between primers
                distance = rSeqIdx - fSeqIdx + primer.fSeq.length()
                if(fSeqIdx >= rSeqIdx) {
                    distance = fSeqIdx - rSeqIdx + primer.fSeq.length()
                }
                if(debugging) { 
                    err.println "preliminary match, distance=${distance}"
                }
                if(((primer.bp + distWindow) >= distance) && 
                   ((primer.bp - distWindow) <= distance)) { 
                    // both primers hit and within the expected distance
                    if(debugging) {
                        err.println "primers match within distance\t" +
                                    primer.locus
                    }

                    // now check expectations relative to haplotype calls
                    hapHits = false
                    if(hapCalls != null) {
                        (lfeature, rfeature) =
                            primer.locus.split(featureSeparator)
                        hapCalls.each { hap, featureSet ->
                            if(debugging) { 
                                err.println "looking for ${lfeature} and ${rfeature}"
                                err.println "in ${featureSet.join(',')}"
                            }
                            if(featureSet.contains(lfeature) &&
                               featureSet.contains(rfeature)) {
                                hapHits = true
                            }
                        }
                    }
                    if((hapCalls == null) || hapHits) {
                        if(debugging) {
                            err.println "hit for " + primer.locus +
                                        " in ${readName}"
                        }
                        genotype.add(primer.locus)
                        // record the gene hits per read
                        readSet.add(primer.locus)
                    }
                } else { // if not in distance range
                    if(debugging) {
                        err.println "not in distance"
                    }
                } 
            } else { // if both primers don't hit
                if(debugging) {
                    err.println "both primers don't match"
                }
            } 
        } // each pair of primers
    } // each locus

    return [genotype, readSet]
} //interpretRead

/*
 * readPrimers
 *
 * Read the file containing the information on the primers.
 *
 * @param primerFileName a String containing the full path to the primer file
 * @param bp the total distance between the two primers, including the primers
 * @return a Map: locus names to array of 
 *    Expandos (locus, fName, fSeq, rName, rSeq, bp)

 */ 
def Map readPrimers(String primerFileName, Integer bp) { 
    if(debugging) { 
        err.println "readPrimers(primerFileName=${primerFileName})\n"
    }
    // open primer file
    primerReader = new BufferedReader(new FileReader(primerFileName))

    locusPrimersMap = [:]
    primersCount = 0

    primerReader.eachLine { line ->
        if(line.contains("name")) { // skip header
            return
        }

        p = new Expando()
        // e.g., 3DL3-2DL3_up	TAGTACTAGAAACTCAAGCAGGAAAATTAGAATGGCTTCTTGTCACAATT	3DL3-2DL3_down	CATATGTAATGTGCAAAATGTCTAACAGGTATTATTAACATTATCAGAG

        (fName, fSeq, rName, rSeq) = line.split()
        locus = fName.split("_")[0]
        p.locus = locus
        p.fName = fName
        p.fSeq = fSeq.toUpperCase()
        p.rName = rName
        p.rSeq = rSeq.toUpperCase()
        p.bp = bp
        existingArray = locusPrimersMap[locus]
        if(!existingArray) { 
            primersArray = [p]
            locusPrimersMap[locus] = primersArray
        } else { 
            existingArray.add(p)
        }
        primersCount++
    } // each line

    primerReader.close()

    if(debugging) { 
        err.println "readPrimers: done; ${primersCount} pairs of primers\n"
    }
    return locusPrimersMap
} // readPrimers

/*
 * readHapCalls
 * Read the file containing the haplotype pair calls.
 *
 * @param hapsFileName a String containing the full path to the haplotype calls
 * @return a Map: haplotype name to Array of gene names
 *
 */ 
def HashMap<String,HashSet<String>> readHapCalls(String hapsFileName) { 
    if(debugging) { 
        err.println "readHapCalls(hapsFileName=${hapsFileName})\n"
    }
    // open haps file
    hapReader = new BufferedReader(new FileReader(hapsFileName))
    
    HashMap<String,HashSet<String>> retMap = new HashMap<String,HashSet<String>>()
    
    // the first line is the haplotype list, rest are single haplotypes
    // and their unordered genes
    //   e.g., haplotype list: cA03~tB02+cB02~tA01
    //         cA03~tB02 (unordered): 2DL3^2DP1^2DS1^3DL2^3DL3
    //         cB02~tA01 (unordered): 2DL2^2DL4^2DS2^2DS4^3DL1^3DL2^3DL3^3DP1
    hapReader.eachLine { line ->
        if(line.contains("haplotype list")) { // skip first line
            return
        }

        geneLinePattern = ~/(\S+?)\s*\S*:\s*(\S+)/
        matcher = geneLinePattern.matcher(line)
        if(debugging) { 
            err.println "matching ${geneLinePattern} against ${line}"
            err.println "matcher=${matcher}"
        }
        gl = null
        hapName = matcher[0][1]
        gl = matcher[0][2]
        if(gl == null) {
            err.println "hapName=${hapName}"
            err.println "ERROR: unexpected format in ${hapsFileName}: ${line}"
            return null
        } else {
            if(debugging) {
                err.println "hapName=${hapName}"
                err.println "gl=${gl}"
            }
        }

        gl.split("\\^").each { gene ->
            hapSet = retMap[hapName]
            if(hapSet == null) {
                hapSet = new HashSet<String>()
                retMap[hapName] = hapSet
            }
            if(debugging) {
                err.println "readHapCalls: adding ${gene} to hapSet"
            }
            hapSet.add(gene)
        }
    } // each line
    hapReader.close()

    if(debugging) { 
        err.println "readHapCalls: done"
    }

    return retMap
} // readHapCalls

// http://groovyconsole.appspot.com/script/29005
def String reverseComplement(String seq) { 
    if(debugging) { 
        //err.println "reverseComplement(seq=${seq})"
    }
    def complements = [ A:'T', T:'A', U:'A', G:'C', C:'G', Y:'R', R:'Y', S:'S', W:'W', K:'M', M:'K', B:'V', D:'H', H:'D', V:'B', N:'N' ]

    seq = seq.toUpperCase().replaceAll( /./ ) { complements."$it" ?: 'X' } 
    seq = seq.reverse()
    if(debugging) { 
        //err.println "reverseComplement: return ${seq}"
    }
    return seq
} // reverseComplement
