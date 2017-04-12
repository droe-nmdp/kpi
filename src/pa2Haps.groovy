#!/usr/bin/env groovy

/*
 * pa2Haps
 *
 * Fit PA genotypes to (potentially ambiguous) haplotype pairs
 * 
 * e.g., ./pa2Haps.groovy -h input/all_haps_v1.txt -q .bin1,2DL3.bin1,2DL3.bin1,2DL3.bin1,2DL3-2DL5B_2DL3-2DP1.bin1,2DL3-2DP1.bin1,2DP1.bin1,2DP1.bin1,2DP1.bin1,2DP1.bin1,2DP1-2DL1_2DP1-2DS1.bin1,2DS1.bin1,2DS1.bin1,2DS1.bin1,2DS1-3DL2.bin1,3DL2.bin1,3DL2.bin1,3DL3.bin1,3DL3.bin1,3DL3.bin1,3DL3.bin1,3DL3-2DL3.bin1,1G.bin1,2G.bin1 -o prediction.txt
 * 
 * Requires
 *  guava.jar: https://github.com/google/guava
 *    http://google.github.io/guava/releases/19.0/api/docs/com/google/common/collect/Table.html
 *
 * @author Dave Roe
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
// string with probe query hits
// e.g., 2DL3.bin1,2DL3.bin1,2DL3.bin1,2DL3.bin1,2DL3-2DL5B_2DL3-2DP1.bin1,2DL3-2DP1.bin1,2DP1.bin1,2DP1.bin1,2DP1.bin1,2DP1.bin1,2DP1-2DL1_2DP1-2DS1.bin1,2DS1.bin1,2DS1.bin1,2DS1.bin1,2DS1-3DL2.bin1,3DL2.bin1,3DL2.bin1,3DL3.bin1,3DL3.bin1,3DL3.bin1,3DL3.bin1,3DL3-2DL3.bin1,1G.bin1,2G.bin1
String qString = new String(options.q)
// open output file
outWriter = new PrintWriter(new File(options.o).newOutputStream(), true)

// load input files
// Table<hap names, locus names, gene count>
HashBasedTable<String, String, Boolean> hapTable =
       readReferenceHaplotypes(hapReader)
// Table<hap pair names, locus names, gene count>
HashBasedTable<String, String, Boolean> hapPairTable =
    makeHapPairTable(hapTable)
if(debugging <= 3) {
    err.print "loaded ${hapTable.columnKeySet().size()} loci, "
    err.println "${hapTable.rowKeySet().size()} reference haplotypes, "
    err.println "${hapPairTable.rowKeySet().size()} reference haplotype pairs"
}
HashSet<String> querySet = readQueries(qString)

// _all_ P/A gene genotype from genes -> Boolean
HashMap<String,Boolean> genPAMap
// a set of genes _that hit_
HashSet<String> genHitSet
(genPAMap, genHitSet) = makeGenotypeMaps(querySet, hapTable.columnKeySet(),
                                         hapPairTable)

// haplotype pairs -> Map[locus: boolean]
HashMap<String, HashMap<String,Boolean>> interpPAMap = null
interpPAMap = interpretPA(hapPairTable, genPAMap)
err.println "${interpPAMap.size()} interpretation(s) "

// make hap pair predictions from probe hits
HashSet<String> interpHapSet = null
interpHapSet = interpretHapMarkers(genHitSet)
err.println "${interpHapSet.size()} haplotype pair interpretation(s)"

// combine pa and hap intepretations via overlap
/*todo(remove)
HashSet<String> combinedInterp = interpHapSet.intersect(interpPAMap.keySet())
err.println "${combinedInterp.size()} combined interpretation(s)"
*/

HashSet<String> paLoci = interpPAMap.keySet() 
writeOutput(outWriter, paLoci, interpHapSet)

outWriter.close()
// end main

/*
 * writeOutput
 *
 * Send the haplotype pair predictions to the output file.
 */
def void writeOutput(PrintWriter outWriter, HashSet<String> interpPASet,
                     HashSet<String> interpHapSet) {
    HashSet<String> outSet = new HashSet()
    Integer highestHitcount = 0
    Iterator e = interpPASet.iterator()
    err.println "interpPASet=${interpPASet}"//todo
    err.println "interpHapSet=${interpHapSet}"//todo
    while ( e.hasNext() ) {
        String paPair = (String)e.next().toString();
        hitcount = 0
        paPair.split("\\+").each { paHap ->
            if(interpHapSet.contains(paHap)) {
                hitcount++
            }
        }
        err.println "writeOutput: hitcount ${hitcount} for ${paPair}"//todo
        if(hitcount > highestHitcount) {
            err.println "writeOutput: setting highestHitcount"//todo
            outSet = new HashSet()
            outSet.add(paPair)
            highestHitcount = hitcount
        } else if(hitcount == highestHitcount) {
            err.println "writeOutput: adding highestHitcount"//todo
            outSet.add(paPair)
        }
    }

    diploidHapSet = interpHapSet
    // assume cA01~tA01 is homozygous
    if((interpHapSet.size() == 1) &&
       (interpHapSet.iterator()[0] == "1")) { // make homozygous
        val = interpHapSet.iterator()[0]
        val += "+${val}"
        diploidHapSet = new HashSet(1)
        diploidHapSet.add(val)
    }
    if(debugging <= 3) {
        err.println "${outSet.size()} combined predictions"
        err.println "snps: " + diploidHapSet.join('/')
        err.println "pa: " + interpPASet.join('/')
        err.println "combined: " + outSet.sort().join('/')
    }
    outWriter.println "snps: ${diploidHapSet.join('/')}"
    outWriter.println "pa: ${interpPASet.join('/')}"
    outWriter.println "combined: ${outSet.sort().join('/')}"
} // writeOutput

/*
 * makeHapPairTable
 *
 * @param hapTable single haplotype table of <hap, gene, Boolean>
 * @return Table double haplotype of <hap pair, gene, Boolean>
 */
HashBasedTable<String, String, Boolean> makeHapPairTable(Table hapTable) {
    if(debugging <= 1) { 
        err.println "makeHapPairTable()"
    }
    HashBasedTable<String, String, Boolean> hapPairTable = HashBasedTable.create()
    HashSet<String> hapList = hapTable.rowKeySet()
    if(debugging <= 2) { 
        err.println "hapList=${hapList}"
    }
    Set<String> locusList = hapTable.columnKeySet()

//    for(int i=0; hap1 = hapList[i++]; i < hapList.size()) {
    for(hap1 in hapList) {
//        for(int j=0; hap2 = hapList[j++]; j < hapList.size()) {
        for(hap2 in hapList) {
            if(hap1.compareTo(hap2) > 0) { // sort and eliminate duplicates
                continue
            }
            hapPair = "${hap1}+${hap2}"
            locusList.each { locus ->
                cellVal = false
                if(hapTable.get(hap1, locus) || hapTable.get(hap2, locus)) {
                    cellVal = true
                }
                hapPairTable.put(hapPair, locus, cellVal)
            } // each gene
        }
    }
    if(debugging <= 1) { 
        err.println "makeHapPairTable: return"
    }
    return hapPairTable
} // makeHapPairTable

/*
 * readReferenceHaplotypes
 *
 * For now, just skip the nomenclature and freq columns.
 *
 * column headers: haplotype	nomenclature	freq	3DL3	2DS2	2DL2	2DL3	2DP1	2DL1	3DP1	2DL4	3DL1	3DS1	2DL5	2DS3	2DS5	2DS4	2DS1	3DL2
 *                 1/43/47/48	cA01~tA01	55.40%	1	0	0	1	1	1	1	1	1	0	0	0	0	1	0	1
 *
 * @param hapFile tab-separated file containing the reference haplotypes, their nomenclature name and frequencies
 * @return HashBasedTable containing just the loci (ignore nomenclature and freq)
 *  rows headers are the haplotypes
 *  column headers are the loci
 *  cell values are the locus count
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
 * readProbes
 *
 * @param reader tab-separated file containing the loci and their probes; no column headers; locus in first column, probe in second
 * @return a Map of probe to gene/locus
 */
def HashMap<String,String> readProbes(FileReader reader) {
    HashMap<String,String> probeMap = new HashMap()
    reader.eachLine { line ->
        (gene, probe) = line.split('\t')
        if((gene == null) || (gene == "")) {
            return
        }
        probeMap[probe] = gene
    } // each line

    return probeMap
} // readProbes

/*
 * readQueries
 *
 * e.g., 2DL3.bin1,2DL3.bin1,2DL3.bin1,2DL3.bin1,2DL3-2DL5B_2DL3-2DP1.bin1,2DL3-2DP1.bin1,2DP1.bin1,2DP1.bin1,2DP1.bin1,2DP1.bin1,2DP1-2DL1_2DP1-2DS1.bin1,2DS1.bin1,2DS1.bin1,2DS1.bin1,2DS1-3DL2.bin1,3DL2.bin1,3DL2.bin1,3DL3.bin1,3DL3.bin1,3DL3.bin1,3DL3.bin1,3DL3-2DL3.bin1,1G.bin1,2G.bin1
 * @param hitString String containing loci that hit
 * @return Set of only the probes that hit
 */
def HashSet<String> readQueries(String reader) {
    if(debugging <= 1) {
        err.println "readQueries(reader=${reader})"
    }
    
    HashSet<String> probeSet = new HashSet()
    reader.split(',').each { multiLocus ->
        multiLocus.split('_').each { locus ->
            lIn = locus.replaceFirst(".bin1", "")
            if((lIn != "") && !(lIn =~ /unmapped/)) {
                probeSet.add(lIn)
            }
        } // each part of multi locus
    } // each comma-separated list of loci

    if(debugging <= 1) {
        err.println "readQueries: return ${probeSet.size()} loci: ${probeSet}"
    }
    return probeSet
} // readQueries

/*
 * makeGenotypeMaps
 *
 * @param qString Set containing Strings from the '.bin1' loci list
 * @param geneSet Set containing all loci 
 * @param hapPair Table<hap pair names, locus names, gene count>
 * @return List of two objects HashMap<String,Boolean> gene PA and HashSet<String> gene all
 */
List makeGenotypeMaps(HashSet<String> querySet, Set<String> geneSet,
                      HashBasedTable<String, String, Boolean> hapPairTable) {
    HashMap<String,Boolean> paProbeMap = new HashMap()
    HashSet<String> allProbeMap = new HashSet()

    // initialze paProbeMap to false for all genes
    for(String gene : geneSet) {
        paProbeMap[gene] = false
    } // initialze paProbeMap

    querySet.each { gene ->
//        if(geneSet.contains(gene)) { 
            paProbeMap[gene] = true
            allProbeMap.add(gene)
//        }
    } // each gene probe

    return [paProbeMap, allProbeMap]
} // makeGenotypeMaps

/*
 * hapPairs2Genotypes
 *
 * Generate synthetic genotypes from pairs of reference haplotypes.
 * 49 choose 2 is 1176.
 * Due to the properties of the Groovy Map '+' operation, a locus from one
 * haplotype overwrites the other, so false values must be null in the Maps.
 *
 * @param hapMap Map of haplotype names to a Map of loci and their Boolean values
 * @return a Map from haplotype pairs to a Map of loci and their Boolean values
 */
def HashMap<String, HashMap<String,Boolean>> hapPairs2Genotypes(HashMap<String,HashMap<String,Boolean>> hapMap, HashMap<String,String> nomMap) {
    estimate = hapMap.size() * 2
    HashMap<String, HashMap<String,Boolean>> genMap = new HashMap(estimate)
    hapMap.each { h1Name, h1Map ->
        hapMap.each { h2Name, h2Map ->
            if(h1Name > h2Name) { // avoid dups and sort by name
                return
            }
            (h1Name, h2Name) = [h1Name, h2Name].sort()
            genMap["${h1Name}+${h2Name}"] = h1Map + h2Map //todo quotes here?
        }
    }
    return genMap
} // hapPairs2Genotypes

/*
 * interpretPA
 *
 * Interpret a genotype in the context of haplotype pairs.
 * It assumes that the genotypes just include the genes that are present.
 *
 * @param hapPairTable Table<hap pair names, locus names, gene count>
 * @param genPAMap HashMap _all_ P/A gene genotype from genes -> Boolean
 * @return a Map from haplotype pairs to a Map of loci and their Boolean values; this is a potentially ambiiguous P/A interpretPAation
 */
HashMap<String, HashMap<String,Boolean>>
        interpretPA(HashBasedTable<String, String, Boolean> hapPairTable,
                    HashMap<String,Boolean> genPAMap) {
    if(debugging <= 1) {
        err.println "interpretPA()"
        err.println "genPAMap = " + genPAMap
    }
    // haplotype name -> Map[locus] -> Boolean
    HashMap<String, HashMap<String,Boolean>> retMap = new HashMap()
    // hap name -> Map[locus] -> Boolean
    Map<String, Map<String, Boolean>> rowMap = hapPairTable.rowMap()
    rowMap.each { haplotype, hapPairMap ->
        /*if(!(haplotype =~ /33\+4/)) { //debugging
            return
        }*/
        if(debugging <= 1) {
            err.println "interpretPA: comparing with haplotype ${haplotype} ${hapPairMap}"
        }

        Boolean allFound = true
        hapPairMap.each { locus, value ->
            if(allFound == false) { // hap pair is false if any locus is false
                return
            }
            genVal = genPAMap[locus]
            same = false
            if(((value == true) && (genVal == true)) ||
               (value == false) && ((genVal == null) || (genVal == false))) {
                    same = true
            }
            allFound = same
        }

        if(allFound) { 
            retMap[haplotype] = hapPairMap
        }
    } // each row/haplotype

    if(debugging <= 1) {
        err.println "interpretPA: return ${retMap}"
    }
    return retMap
} // interpretPA


/*
 * interpretHapMarkers
 *
 * Interpret pairs of haplotypes from the hap (not PA) markers.
 *
 * @param g HashSet of all markers _that hit_
 * @return Set of haplotype pair predictions (e.g., cA01~tA01+cB01~tA01, ...)
 */
HashSet<String> interpretHapMarkers(HashSet<String> g) {
    if(debugging <= 1) {
        err.println "interpretHapMarkers(g=${g})"
    }
    HashSet<String> retSet = new HashSet()

    if("1G" in g) {
        if("2G" in g) {
            if("3T" in g) {
                retSet.add("1")
            }
            if("3C" in g) {
                if("4T" in g) {
                    retSet.add("1")
                }
                if("4C" in g) {
                    retSet.add("1")
                }
            } // 3
        }
        if("2C" in g) {
            retSet.add("1")
        } // 2
    }
    if("1T" in g) {
        if("5A" in g) {
            if("6A" in g) {
                retSet.add("6")
            }
            if("6T" in g) {
                if("2DL1" in g) {
                    retSet.add("6")
                } else {
                    retSet.add("4")
                }
            }
        }
        if("5G" in g) {
            if("7G" in g) {
                retSet.add("4")
            }
            if("7A" in g) {
                if("2DS3" in g) {
                    retSet.add("9")
                }
                if("2DS5" in g) {
                    retSet.add("7")
                }
            } // 7
        } // 5
    } // 1

    HashSet<String> retSetPairs = new HashSet()
    retSet.each { hap1 ->
        retSet.each { hap2 ->
            if(hap1 != hap2) {
                hl = [hap1, hap2].sort().join('+')
            } else {
                hl = hap1
            }
            retSetPairs.add(hl)
        }
    }

    if(debugging <= 1) {
        err.println "interpretHapMarkers return: ${retSetPairs}"
    }
    return retSetPairs
} // interpretHapMarkers

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
    CliBuilder cli = new CliBuilder(usage:'pa2Haps.groovy [options] ', header:'Options:')
    cli.help('print this message')
    cli.h(longOpt:'haps', args:1, argName:'file', 'file with haplotype definitions',
      required: true)
    cli.q(longOpt:'qout', args:1, argName:'file', 'file with results of probe queries',
      required: true)
    cli.o(longOpt:'output', args:1, argName:'file', 'output file containing haplotype pair predictions',
      required: true)
    OptionAccessor options = cli.parse(args)
    return options
} // handleArgs
