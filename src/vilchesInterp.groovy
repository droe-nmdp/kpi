#!/usr/bin/env groovy

/*
 * Make presence/absence genotype calls on short reads using the genotyping
 * algorithm by Vilches, et al. 2007.
 *
 * Assumes a single fastq file as input.
 *
 * Load the probe file. Assumes the ambiguous IUPAC codes have been resolved.
 * Loop through every read.
 * For each read, check each pair of primers.
 * Output negatives and positives
 *
 * Run via, e.g., 
 *   ./vilchesInterp.groovy -Xms5g input/vilches_probes.txt ccs_99.fastq ccs_99 output/ccs_99_genotype.txt output/ccs_99_reads.txt
 *
 * http://biojava.org/wiki/BioJava:CookBook3:FASTQ
 * Requires
 *   biojava3-core.jar
 *   biojava3-sequencing.jar
 *   http://code.google.com/p/guava-libraries/
 * @author Dave Roe
 */

import com.google.common.io.*
import java.util.*
import groovy.io.*
import java.util.zip.GZIPInputStream

// things that may change per run
debugging = false
distWindow = 7 // allow this much cushion on each side for distance between primers

loci = ['2DL1', '2DL2', '2DL3', '2DL4', '2DL5', '2DS1', '2DS2', '2DS3', '2DS4',
        '2DS5', '3DL1', '3DL2', '3DL3', '3DS1', '2DP1', '3DP1']
heapSpace = args[0] // e.g., -Xms5g
primerFileName = args[1]
fastqFile = args[2]
id = args[3]
outFileName = args[4]
outReadFileName = args[5]
err = System.err
DEFAULT_BUFFER_SIZE=5096

// Map: locus -> array of primers for that locus
// each primer is an Expando that has: locus, fName, fSeq, rName, rSeq, bp
// bp is the total distance between the two primers, including the primers
primers = [:]
primers = readPrimers(primerFileName)

// open output files and print the headers
outWriter = new PrintWriter(new File(outFileName).newOutputStream(), true)
outWriter.println "locus\tP/A"
outReadWriter = new PrintWriter(new File(outReadFileName).newOutputStream(), true)
outReadWriter.println "read\tlocus list\tlength"

if(debugging) { 
    err.println "using distance window of ${distWindow}"
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

//FastqReader fastqReader = new IlluminaFastqReader()
//FastqReader fastqReader = new SangerFastqReader()
BufferedReader r;
if(fq.getName().endsWith("fastq.gz")) {
    r=new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(fullFileName),DEFAULT_BUFFER_SIZE)),DEFAULT_BUFFER_SIZE)
} else {
    r=new BufferedReader(new FileReader(fullFileName),DEFAULT_BUFFER_SIZE)
}
long nLine=0L
String line;
String readName = null // read name
while((line=r.readLine())!=null) {
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
    readSet = new TreeSet()    
    (genotype, readSet) = interpretRead(primers, seq, genotype, readSet)

    readHits[readName] = readSet
} // each line in fastq file

if(debugging) { 
    err.println "done processing ${fullFileName}, ${readHits.keySet().size()} reads\n"
    err.println "genotype size=" + genotype.size()
}
    
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
        ll = locusList.join(", ")
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
 * Return two objects that are also input: genotype, and readSet
 */
def ArrayList<Map> interpretRead(Map primers, String sequence,
                                 TreeSet<String> genotype,
                                 TreeSet<String> readSet) {
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
                err.println "sequence=${sequence}"
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
                        err.println "found probes (opposite direction)\t${primer.locus}"
                    }
                }
            }

            if(debugging) { 
                err.println "rSeqIdx=${rSeqIdx}, fSeqIdx=${fSeqIdx}, rSeq length=${primer.rSeq.length()}, fSeq length=${primer.fSeq.length()}"
                err.println "primer.bp=${primer.bp}"
            }
            if((fSeqIdx > 0) && (rSeqIdx > 0)) {
                distance = rSeqIdx - fSeqIdx + primer.fSeq.length()
                if(fSeqIdx >= rSeqIdx) {
                    distance = fSeqIdx - rSeqIdx + primer.fSeq.length()
                }
                if(debugging) { 
                    err.println "distance=${distance}"
                }
                if(((primer.bp + distWindow) >= distance) && 
                   ((primer.bp - distWindow) <= distance)) { 
                    // this is a hit
                    if(debugging) {
                        err.println "hit\t" + primer.locus
                    }
                    genotype.add(primer.locus)

                    // record the gene hits per read
                    readSet.add(primer.locus)
                } // if in distance range
            } // if both primers hit
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
 * @return a Map: locus names to array of 
 *    Expandos (locus, fName, fSeq, rName, rSeq, bp)
 *      bp is the total distance between the two primers, including the primers

 */ 
def Map readPrimers(String primerFileName) { 
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
        // e.g., 2DL1	Fa517	gttggtcagatgtcatgtttgaa	Rc621	cctgccaggtcttgcg	142
        (locus, fName, fSeq, rName, rSeq, bp) = line.split()
        p.locus = locus
        p.fName = fName
        p.fSeq = fSeq.toUpperCase()
        p.rName = rName
        p.rSeq = rSeq.toUpperCase()
        //p.bp = bp.toInteger()
        // remove this?
        bp.split("/").each { length ->  
            p.bp = length.toInteger()
            existingArray = locusPrimersMap[locus]
            if(!existingArray) { 
                primersArray = [p]
                locusPrimersMap[locus] = primersArray
            } else { 
                existingArray.add(p)
            }
            primersCount++
        }
        
    } // each line

    primerReader.close()

    if(debugging) { 
        err.println "readPrimers: done; ${primersCount} pairs of primers\n"
    }
    return locusPrimersMap
} // readPrimers

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
