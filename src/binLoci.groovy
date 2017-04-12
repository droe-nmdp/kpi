#!/usr/bin/env groovy

/*
 * binLoci
 *
 * Given per-locus (gene and intergene) FASTQs, bin them into FASTQ files 
 * of new overlapping loci: each gene, its two inter-gene regions, and all 
 * unmapped regions.
 *
 * ./binLoci.groovy -h input/all_haps_v2.txt -q 2DL3.bin1,2DL3.bin1,2DL3.bin1,2DL3-2DL5B_2DL3-2DP1.bin1,2DL3-2DP1.bin1,2DP1.bin1,2DP1.bin1,2DP1.bin1,2DP1.bin1,2DP1-2DL1_2DP1-2DS1.bin1,2DS1.bin1,2DS1.bin1,2DS1.bin1,2DS1-3DL2.bin1,3DL2.bin1,3DL2.bin1,3DL3.bin1,2DL1.bin1,3DL3-2DL3.bin1,1G.bin1,2G.bin1 -p input/predictions.txt
 * 
 * Requires
 *  guava.jar: https://github.com/google/guava
 *    http://google.github.io/guava/releases/19.0/api/docs/com/google/common/collect/Table.html
 *
 * @author Dave Roe
 * @todo modularize/abstract methods with pa2Haps.groovy
 */

import groovy.io.*
import groovy.util.CliBuilder.*
import groovy.util.OptionAccessor
import com.google.common.collect.Table
import com.google.common.collect.HashBasedTable

// things that may change per run
debugging = 1 // TRACE=1, DEBUG=2, INFO=3

// thing that probably won't change per run
err = System.err

OptionAccessor options = handleArgs(args)

// open file with haplotype definitions
FileReader hapReader = new FileReader(new File(options.h))
// open file with prediction of haplotype pairs
FileReader predReader = new FileReader(new File(options.p))
// string with probe query hits
// e.g., 2DL3.bin1,2DL3.bin1,2DL3.bin1,2DL3.bin1,2DL3-2DL5B_2DL3-2DP1.bin1,2DL3-2DP1.bin1,2DP1.bin1,2DP1.bin1,2DP1.bin1,2DP1.bin1,2DP1-2DL1_2DP1-2DS1.bin1,2DS1.bin1,2DS1.bin1,2DS1.bin1,2DS1-3DL2.bin1,3DL2.bin1,3DL2.bin1,3DL3.bin1,3DL3.bin1,3DL3.bin1,3DL3.bin1,3DL3-2DL3.bin1,1G.bin1,2G.bin1
String qString = new String(options.q)

// load input files
// Table<hap names, locus names, gene count>
HashBasedTable<String, String, Boolean> hapTable =
       readReferenceHaplotypes(hapReader)
ArrayList<String> hapPredictions = readHaplotypePredictions(predReader)

TreeSet<String> genes = new TreeSet()
genes.addAll(hapTable.row(hapPredictions[0]).keySet())
genes.addAll(hapTable.row(hapPredictions[1]).keySet())
/*err.println "predicted haplotypes: ${hapPredictions.join('+')}"
err.println "predicted haplotype 1: ${hap1GeneMap.keySet().join('+')}"
err.println "predicted haplotype 2: ${hap2GeneMap.keySet().join('+')}"*/
err.println "genes: ${genes.join('+')}"

// direct the reads into '*.bin2'
processReads(qString, genes)
// end main

/*
 * processReads
 *
 * e.g., 2DL3.bin1,2DL3.bin1,2DL3.bin1,2DL3.bin1,2DL3-2DL5B_2DL3-2DP1.bin1,2DL3-2DP1.bin1,2DP1.bin1,2DP1.bin1,2DP1.bin1,2DP1.bin1,2DP1-2DL1_2DP1-2DS1.bin1,2DS1.bin1,2DS1.bin1,2DS1.bin1,2DS1-3DL2.bin1,3DL2.bin1,3DL2.bin1,3DL3.bin1,3DL3.bin1,3DL3.bin1,3DL3.bin1,3DL3-2DL3.bin1,1G.bin1,2G.bin1
 * @param hitString String containing loci that hit
 * @param List of Strings of genes present (non inter-genes)
 * @return Set of only the PA (not snp) probes that hit
 */
def void processReads(String reader, TreeSet<String> genes) {
    if(debugging <= 1) {
        err.println "processReads(reader=${reader})"
    }
    def newline = System.getProperty('line.separator')

    reader.split(',').each { multiLocus ->
        multiLocusStrip = multiLocus.replaceFirst(".bin1", "")            
        // skip unmapped and snp markers
        if((multiLocusStrip.size() > 2) && !(multiLocusStrip =~ /unmapped/)) {
            HashSet<String> processedSet = new HashSet()
            multiLocusStrip.split('_').each { locusPair ->
                locusPair.split('-').each { locus ->
                    if(!processedSet.contains(locus)) { 
                        newFileName = "${locus}.bin2"
                        File src = new File(multiLocus)
                        File dest = new File(newFileName)
                        src.eachLine { line -> // append
                            dest << line
                            dest << newline
                        }
                        processedSet.add(locus)
                    }
                }
            } // each multi locus (e.g., 2DP1-2DL1_2DP1_2DS1)
        } // each part of multi locus
    } // each comma-separated list of loci

    // append unmapped reads to *.bin2
    File src = new File("unmapped.bin1")
    err.println "droe1"//todo
    String dir = new File("").getAbsolutePath()
    err.println "dir=${dir}"//todo
    def p =  ~/.*\.bin2/
    new File(dir).eachFileMatch(p) { dest ->
        err.println "droe2: ${dest}"//todo
        src.eachLine { line -> // append
            dest << line
            dest << newline
        }
    }

    if(debugging <= 1) {
        err.println "processReads: return"
    }
    return
} // processReads

/*
 * readReferenceHaplotypes
 *
 * For now, just skip the nomenclature and freq columns.
 *
 * column headers: haplotype	nomenclature	freq	3DL3	2DS2	2DL2	2DL3	2DP1	2DL1	3DP1	2DL4	3DL1	3DS1	2DL5	2DS3	2DS5	2DS4	2DS1	3DL2
 *                 1/43/47/48	cA01~tA01	55.40%	1	0	0	1	1	1	1	1	1	0	0	0	0	1	0	1
 *
 * @param hapFile tab-separated file containing the reference haplotypes, their nomenclature name and frequencies
 * @return HashBasedTable containing just the _present_ loci (and ignore nomenclature and freq)
 *  rows headers are the haplotypes
 *  column headers are the loci
 *  cell values are the locus count
 * @todo modularize into a library (with pa2Haps.groovy)
 */
def HashBasedTable<String, String, Boolean> readReferenceHaplotypes(FileReader reader) {
    Integer startLocusIndex = 3 // haplotype	nomenclature	freq	3DL3 ...
    HashBasedTable<String, String, Boolean> hapTable = HashBasedTable.create()
    // first row are loci
    TreeSet<Integer> loci = new TreeSet()
    ArrayList header = reader.readLine().split('\t')
    header[startLocusIndex..-1].each { locus ->
        if(debugging <= 2) {
            err.println "readReferenceHaplotypes: adding ${locus} to loci list"
        }
        loci.add(locus)
    }

    // non-header rows
    reader.eachLine { line ->
        ArrayList cols = line.split('\t')
        String haplotype = cols[0].trim() // haplotype number
        if((haplotype == null) || (haplotype == "")) { // skip blank rows
            return
        }
        cols[startLocusIndex..-1].eachWithIndex { cell, i ->
            i+=+startLocusIndex
            Boolean genePresent = false
            if(cell.toInteger() > 0) {
                genePresent = true
            } else {
                return // skip non-present genes
            }
            locus = header[i]
            if(debugging < 2) {
                err.println "readReferenceHaplotypes: adding ${locus}=${genePresent} for ${haplotype}"
            }
            
            hapTable.put(haplotype, locus, genePresent)
        }
        if(debugging < 2) {
            err.println "readReferenceHaplotypes: added ${haplotype}=${hapTable.rowMap()[haplotype]}"
        }
    } // each line/haplotype

    if(debugging < 1) {
        err.println "readReferenceHaplotypes: return ${hapTable.rowKeySet().size()} rows, ${hapTable.columnKeySet().size()} columns"
    }
    return hapTable
} // readReferenceHaplotypes

/*
 * Reads the predicted pair of haplotypes and returns them in two ordered lists.
 * Returns a random pair if more than one after combining pa and snp markers.
 *
 * @param FileReader containing the prediction
 *   snps: 4+7/4+9/6+7/6+9/7+9/1/4/6/7/9/1+4/1+6/1+7/1+9
 *   pa: 13+3/1+7/3+6/35+42/13+35/14+25/12+23/13+14/9+99/17+99/7+98/14+42/7+99
 *   combined: 1+7
 * @return ArrayList containing the names of the two haplotype predictions
 */
ArrayList<String> readHaplotypePredictions(FileReader reader) {
    ArrayList<String> retlist = new ArrayList()
    (name, snpgl) = reader.readLine().split(':')
    (name, pagl) = reader.readLine().split(':')
    (name, bothgl) = reader.readLine().split(':')
    workgl = bothgl.trim()
    if(workgl == "") {
        workgl = pagl
    }
    if(debugging <= 2) { 
        err.println "readHaplotypePredictions: workgl=${workgl}"
    }
    firstgl = workgl
    // just take the first haplotype pair if multiple
    if(firstgl.contains("/")) { 
        (firstgl, restgls) = workgl.split("/", 2)
    }
    (hap1, hap2) = firstgl.split("\\+", 2)
    retlist.add(hap1)
    retlist.add(hap2)
    return retlist
} // readHaplotypePredictions

/*
 * handleArgs
 * 
 * Parse and return the input.
 *
 * @todo e.g., here
 *
 * @param args List of Strings containing command line arguments.
 * @return Option
 */
OptionAccessor handleArgs(String[] args) { 
    CliBuilder cli = new CliBuilder(usage:'binLoci.groovy [options] ', header:'Options:')
    cli.help('print this message')
    cli.h(longOpt:'haps', args:1, argName:'file', 'file with haplotype definitions',
      required: true)
    cli.q(longOpt:'qout', args:1, argName:'file', 'string with results of probe queries',
      required: true)
    cli.p(longOpt:'predictions', args:1, argName:'file', 'output file containing haplotype pair predictions',
      required: true)
    OptionAccessor options = cli.parse(args)
    return options
} // handleArgs
