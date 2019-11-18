#!/usr/bin/env groovy

/*
 * pa2Haps4
 *
 * Fit PA genotypes to (potentially ambiguous) haplotype pairs.
 * Explicitly doubles cA01~tA01 (1).
 * 
 * todo: update example
 * e.g., pa2Haps4.groovy -h $HOME/doc/kir/snp/all_haps_v4.txt -q 2DL2.bin1,3DL2.bin1,2DL3.bin1,3DL3.bin1,2DL4.bin1,3DP1-2DL4.bin1,2DP1.bin1,3DP1.bin1,2DS2.bin1,cA01~tA01.bin1,2DS4.bin1,cB02~tA01.bin1,3DL1.bin1,cB02~tB01.bin1 -o prediction.txt
 *
 * 
 * Requires
 *  guava.jar: https://github.com/google/guava
 *    http://google.github.io/guava/releases/19.0/api/docs/com/google/common/collect/Table.html
 *
 * @author Dave Roe
 * @todo make more cpu efficient?
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

if(debugging <= 5) {
    err.println "pa2Haps4 -h ${options.h} -q ${options.q} -o -${options.o}"
}

// open file with haplotype definitions
FileReader hapReader = new FileReader(new File(options.h))
// string with probe query hits
// e.g., 2DL3.bin1,2DL3.bin1,2DL3.bin1,2DL3.bin1,2DL3-2DL5B_2DL3-2DP1.bin1,2DL3-2DP1.bin1,2DP1.bin1,2DP1.bin1,2DP1.bin1,2DP1.bin1,2DP1-2DL1_2DP1-2DS1.bin1,2DS1.bin1,2DS1.bin1,2DS1.bin1,2DS1-3DL2.bin1,3DL2.bin1,3DL2.bin1,3DL3.bin1,3DL3.bin1,3DL3.bin1,3DL3.bin1,3DL3-2DL3.bin1,1G.bin1,2G.bin1
String qString = new String(options.q)

// load input files
// Table<hap number, locus names, gene count>
HashBasedTable<String, String, Boolean> hapTable
HashMap<String, Float> freqMap // hap number -> frequency
HashMap<String, String> nomenclatureMap // hap name -> hap number
(hapTable, freqMap, nomenclatureMap) = readReferenceHaplotypes(hapReader)
err.println "nomenclatureMap=" + nomenclatureMap//todo
// Table<hap pair numbers, locus names, gene count>
HashBasedTable<String, String, Boolean> hapPairTable =
    makeHapPairTable(hapTable)
if(debugging <= 3) {
    err.print "loaded ${hapTable.columnKeySet().size()} loci, "
    err.println "${hapTable.rowKeySet().size()} reference haplotypes, "
    err.println "${hapPairTable.rowKeySet().size()} reference haplotype pairs"
}
HashSet<String> lociSet = readQueries(qString)

// _all_ P/A gene genotype from genes -> Boolean
HashMap<String,Boolean> genPAMap
// a set of genes _that hit_
HashSet<String> genHitSet
(genPAMap, genHitSet) = makeGenotypeMaps(lociSet,
										 hapTable.columnKeySet(),
                                         hapPairTable, nomenclatureMap)

// haplotype pairs -> Map[locus: boolean]
HashMap<String, HashMap<String,Boolean>> interpPAMap = null
interpPAMap = interpretPA(hapPairTable, genPAMap)
err.println "${interpPAMap.size()} interpretation(s) "

// make hap pair predictions from probe hits
HashSet<String> interpHapSet = null
interpHapSet = interpretHapMarkers(genHitSet, nomenclatureMap.values().toSet())
err.println "${interpHapSet.size()} haplotype pair interpretation(s)"

HashSet<String> paLoci = interpPAMap.keySet() 
writeOutput(options, genPAMap, genHitSet, paLoci, interpHapSet,
			hapPairTable)
// end main

/*
 * writeOutput
 *
 * Make homozygous if the haplotype is cA01~tA01 (1).
 * 
 * Send the haplotype pair predictions to the output file.
 * @param options OptionAccessor contains option 'o' with String with output file name
 * @param genPAMap HashMap<String,Boolean> of input of _all_ loci genotype hit map
 * @param genHitSet HashSet<String> of loci that _hit_ (present)
 * @param interpPASet Set of Strings with haplotype predictions from gene interpretation
 * @param interpHapSet Set of Strings with haplotype predictions from haplotype interpretation
 * @param hapPairTable Table<hap number, locus names, gene count>
 *
 */
def void writeOutput(OptionAccessor options, Map genPAMap,
                     Set genHitSet, 
                     HashSet<String> interpPASet,
                     HashSet<String> interpHapSet,
					 HashBasedTable<String, String, Boolean> hapPairTable) {
    if(debugging <= 1) {
        err.println "writeOutput(interpPASet=${interpPASet}, interpHapSet=${interpHapSet})"
    }
    String outFileName = options.o
    // open output file
    outWriter = new PrintWriter(new File(outFileName).newOutputStream(), true)
    id = outFileName.replaceFirst("_prediction.txt", "")
	outWriter.println "id: $id"
	
    // output genotypes
    //outWriter.println genPAMap
	geneHitSet = genHitSet
	//intergeneHitSet = genHitSet.collect{ it.contains("-") ? it : null}.findAll()
    TreeSet<String> intergeneHitSet = new TreeSet()
/*    outWriter.println "genes and haplotypes: " + geneHitSet.sort().join("+")
*/
	
    // assume homozygous if only one hap call
    // (interpHapSet.iterator()[0] == "1")
    if((interpHapSet.size() == 1) { // make homozygous
        val = interpHapSet.iterator()[0]
        val += "+${val}"
        interpHapSet = new HashSet(1)
        interpHapSet.add(val)
    }

    // reduce the gene-only haplotype-pair ambiguity by
	// ranking by the number of haplotype hits that are found
	// in the haplotype-pair and keeping the best
    HashSet<String> pacombinedSet = paReduceByHap(genHitSet,
												  interpPASet,
												  interpHapSet)
    HashSet<String> hapcombinedSet = hapReduceByPA(genHitSet,
												   interpPASet,
												   interpHapSet)
	if(debugging <= 3) { 
		err.println "haplotype-reduced gene interp from ${interpPASet.size()} to ${pacombinedSet.size()}"
		err.println "gene-reduced hap interp from ${interpPASet.size()} to ${hapcombinedSet.size()}"
	}
//	outWriter.println "haplotype reduced: ${pacombinedSet.size()}"

	// get the best call for the combined prediction
	// best to worst: gene+hap, inter-gene reduced, gene only, hap only
	if((pacombinedSet.size() == 0) && (hapcombinedSet.size() == 0)) {
		// first, take the inter-gene reduced
		if(interpPASet.size() > 0) {
			pacombinedSet = interpPASet
		} else if(interpPASet.size() > 0) {
			pacombinedSet = interpPASet
		} else if(interpHapSet.size() > 0) {
			pacombinedSet = interpHapSet
		} 
	}

    // get least ambig between pa and hap-based interp
    leastAmbigSet = pacombinedSet
    if(hapcombinedSet.size() < pacombinedSet.size()) {
        leastAmbigSet = hapcombinedSet
    }
    
	if(debugging <= 3) {
        err.println "${pacombinedSet.size()} combined predictions"
		err.println "genotype: " + genHitSet.sort().join("+")
        err.println "gene: " + interpPASet.sort().join('|')
        err.println "haplotype: " + interpHapSet.sort().join('|')
        err.println "pa combined: " + pacombinedSet.sort().join('|')
        err.println "hap combined: " + hapcombinedSet.sort().join('|')
        err.println "least ambig: " + leastAmbigSet.sort().join('|')
    }
    if((options.a != null) && (options.a != "0") && (options.a != "F")) { 
        outWriter.println "${pacombinedSet.size()} pa combined predictions"
	    outWriter.println "genotype: " + genHitSet.sort().join("+")
        outWriter.println "gene: ${interpPASet.sort().join('|')}"
        outWriter.println "haplotype: ${interpHapSet.sort().join('|')}"
        outWriter.println "pa combined: ${pacombinedSet.sort().sort().join('|')}"
        outWriter.println "hap combined: ${hapcombinedSet.sort().sort().join('|')}"
        outWriter.println "all\t${id}\t${genHitSet.sort().join("+")}\t${interpPASet.sort().join('|')}\t${interpHapSet.sort().join('|')}\t${pacombinedSet.sort().sort().join('|')}\t${hapcombinedSet.sort().sort().join('|')}"
    }
    outWriter.println "${leastAmbigSet.sort().sort().join('|')}"
	outWriter.close()
} // writeOutput

/* 
 * reduceGeneByInterGene
 *
 * @param geneHitSet full set of gene hits
 * @param genHitSet full set of intergene hits
 * @param hapTable Table<hap number, locus names, gene count>
 */
HashSet<String> reduceGeneByInterGene(HashSet<String> interpPASet,
									  ArrayList<String> intergeneHitSet,
									  HashBasedTable<String, String, Boolean> hapPairTable) {
 	if(debugging <= 1) {
		err.println "reduceGeneByInterGene(${interpPASet.size()} haplotype pairs)"
	}
	HashSet<String> outSet = new HashSet()
    Integer highestHitcount = 0

    Iterator e = interpPASet.iterator()
	Map<String,String> hapPairRowMap = hapPairTable.rowMap()

	// loop through the PA hap pairs (e.g., 1+4)
	// and check if all the expected hit are found in the genotype
	// note, the opposite is not checked: all the genotype hits
	// are found in the hap pairs
    while ( e.hasNext() ) {
        String paPair = (String)e.next().toString();
		Integer paPairCount = 0
 		if(debugging <= 2) {
			err.println "reduceGeneByInterGene: paPair=$paPair"
		}
		// get the intergene regions for this haplotype
		hapPairColMap = hapPairRowMap[paPair]

		// count how many are really in the genotype
		hapPairColMap.each { geneOrInter, present ->
			if(debugging <= 2) {
				err.println "reduceGeneByInterGene: geneOrInter=$geneOrInter"
			}
			// if intergene present in hap pair
			if(geneOrInter.contains("-") && (present)) {
				if(debugging <= 2) {
					err.println "reduceGeneByInterGene: $paPair should have $geneOrInter ..."
				}
				if(intergeneHitSet.contains(geneOrInter)) {
					paPairCount++
					if(debugging <= 2) {
						err.println "reduceGeneByInterGene: and it does: total count now: $paPairCount"
					}
				}
			}
		} // each gene/intergene for this region 
		
		if(paPairCount > highestHitcount) { 
            if(debugging <= 2) {
                err.println "reduceGeneByInterGene: new highest: $paPair $paPairCount"
            }
            outSet = new HashSet()
            outSet.add(paPair)
            highestHitcount = paPairCount
        } else if(paPairCount == highestHitcount) {
            if(debugging <= 2) {
                err.println "reduceGeneByInterGene: adding $paPair to $paPairCount"
            }
            outSet.add(paPair)
		}
	} // each paPair

 	if(debugging <= 1) {
		err.println "reduceGeneByInterGene: return ${outSet.size()} haplotype pairs"
	}
    return outSet
} // reduceGeneByInterGene

/* 
 * Combine the two (gene and haplotype) interpretations
 * by going with one if the other doesn't exist, or reducing
 * the PA haplotype list by haplotype markers if possible.
 * 
 * @param genHitSet Set of loci (genes and haps) that hit
 * @param interpPASet Set of Strings with haplotype predictions from gene interpretation
 * @param interpHapSet Set of Strings with haplotype predictions from haplotype interpretation
 * @return Set of the new haplotype pair predictions; one pair per String
 */
HashSet<String> paReduceByHap(Set genHitSet,
							  HashSet<String> interpPASet,
							  HashSet<String> interpHapSet) {
    HashSet<String> outSet = new HashSet()
    Integer highestHitcount = 0
    Iterator e = interpPASet.iterator()
    // loop through each PA haplotype pair prediction
    while ( e.hasNext() ) {
        String paPair = (String)e.next().toString();
        // for this PA hap pair prediction, count how many of
        // those haplotypes typed as present
        hitcount = 0 
        paPair.split("\\+").each { paHap ->
			// don't knock out a hap that has no markers in the panel
            if(genHitSet.contains(paHap)) {
                hitcount++
			}
			// automatically rule in any pair of haplotypes
			// where at least one is a haplotype without markers
			if(!hapWithMarkers(paHap)) {
				outSet.add(paPair)
			}
        }
        if(debugging <= 2) {
            err.println "paReduceByHap: hitcount ${hitcount} for ${paPair}"
        }
        if(hitcount > highestHitcount) {
            if(debugging <= 2) {
                err.println "paReduceByHap: setting highestHitcount"
            }
            outSet = new HashSet()
            outSet.add(paPair)
            highestHitcount = hitcount
        } else if(hitcount == highestHitcount) {
            if(debugging <= 2) {
                err.println "paReduceByHap: adding highestHitcount"
            }
            outSet.add(paPair)
        }
    }
	
    return outSet
} // paReduceByHap

/* 
 * Combine the two (gene and haplotype) interpretations
 * by going with one if the other doesn't exist, or reducing
 * the haplotype marker list if possible using the PA interpreted pairs.
 * It only checks pairs right now.
 * 
 * @param genHitSet Set of loci (genes and haps) that hit
 * @param interpPASet Set of Strings with haplotype predictions from gene interpretation
 * @param interpHapSet Set of Strings with haplotype predictions from haplotype interpretation
 * @return Set of the new haplotype pair predictions; one pair per String
 * @todo split the pairs and check each haplotype
 */
HashSet<String> hapReduceByPA(Set genHitSet,
							  HashSet<String> interpPASet,
							  HashSet<String> interpHapSet) {
    HashSet<String> outSet = new HashSet()
    Integer highestHitcount = 0
    Iterator e = interpHapSet.iterator()
    // loop through each hap-marker haplotype pair prediction
    while ( e.hasNext() ) {
        String hapPair = (String)e.next().toString();
        // for this hap-marker hap pair prediction, count how many of
        // those haplotypes typed as present in pa data
        hitcount = 0 
        if(interpPASet.contains(hapPair)) {
            hitcount++
        }
		// automatically rule in any pair of haplotypes
		// where at least one is a haplotype without markers
        hapPair.split("\\+").each { hap ->
		    if(!hapWithMarkers(hap)) {
			    outSet.add(hapPair)
		    }
        }
        if(debugging <= 2) {
            err.println "hapReduceByPA: hitcount ${hitcount} for ${hapPair}"
        }
        if(hitcount > highestHitcount) {
            if(debugging <= 2) {
                err.println "hapReduceByPA: setting highestHitcount"
            }
            outSet = new HashSet()
            outSet.add(hapPair)
            highestHitcount = hitcount
        } else if(hitcount == highestHitcount) {
            if(debugging <= 2) {
                err.println "hapReduceByPA: adding highestHitcount"
            }
            outSet.add(hapPair)
        }
    } // while
	
    return outSet
} // hapReduceByPA

/*
 * Returns true if the input haplotype name/number is 
 * one for which we have markers.
 */
Boolean hapWithMarkers(String hap) {
	Boolean ret = false
	hapsWithMarkers = ['1', '3/11', '4', '6/25', '7/9', '8/17', '32', '98/99', 'cA01~tA01', 'cA01~tB01', 'cB02~tA01', 'cB01~tA01', 'cB01~tB01', 'cB02~tB01', 'cB04~tB03', 'cB01~tB05', 'cA01~tB04', 'cB05~tB01', 'cA01~tB05', 'cB05~tA01', 'cA01~tB06', 'cA01~tA02', 'cA03~tB02', 'cB03~tA01']
	if(hapsWithMarkers.contains(hap)) {
		ret = true
	}
	return ret
} // hapWithMarkers

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
            String hapPair = "${hap1}+${hap2}"
            locusList.each { locus ->
                Boolean cellVal = false
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
 * haplotype	nomenclature	freq	structure	3DL3	2DS2	2DL2	2DL3	2DP1	2DL1	3DP1	2DL4	3DL1	3DS1	2DL5	2DS3	2DS5	2DS4	2DS1	3DL2
 * 1	cA01~tA01	55.40%		1	0	0	1	1	1	1	1	1	0	0	0	0	1	0	1
 *
 * @param hapFile tab-separated file containing the reference haplotypes, their nomenclature name and frequencies
 * @return HashBasedTable containing just the loci (ignore nomenclature and freq)
 *    rows headers are the haplotypes
 *    column headers are the loci
 *    cell values are the locus count
 *  HashMap<haplotype number, Float> for the frequencies
 *  HashMap<haplotype number, String> for the nomenclature names
 */
def List readReferenceHaplotypes(FileReader reader) {
    ArrayList retList = new ArrayList()
    // hap number -> frequency
    HashMap<String, Float> freqMap = [:]
    // hap number -> hap name
    HashMap<String, String> nomenclatureMap = [:]

	// todo: don't hard-code this
    Integer startLocusIndex = 5 // haplotype	nomenclature	freq	structure	hap markers	3DL3
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
        freqMap[haplotype] = cols[2].toFloat()
        //old nomenclatureMap[haplotype] = cols[1]
        nomenclatureMap[cols[1]] = haplotype // new
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
    retList.add(hapTable)
    retList.add(freqMap)
    retList.add(nomenclatureMap)
    return retList
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
 * Returns a Set of the locus names based of the file names
 * with the '.bin1' extension.
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
            lIn = locus.replaceFirst("\\.bin1", "") //todo pass this in
			if(lIn =~ /2DL5/) {
				// convert 2DL5A and 2DL5B to 2DL5
				lIn = "2DL5"
			}
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
 * @param lociSet Set containing Strings from the '.bin1' loci list
 * @param geneSet Set containing all reference loci 
 * @param hapPair Table<hap pair names, locus names, gene count>
 * @return List of two objects HashMap<String,Boolean> gene PA and HashSet<String> gene all
 */
List makeGenotypeMaps(HashSet<String> lociSet, Set<String> geneSet,
                      HashBasedTable<String, String, Boolean> hapPairTable,
                      HashMap<String, String> nomenclatureMap) {
    if(debugging <= 1) {
        err.println "makeGenotypeMaps()"
    }
    HashMap<String,Boolean> paProbeMap = new HashMap()
    HashSet<String> allProbeMap = new HashSet()

    // initialze paProbeMap to false for all genes
    for(String gene : geneSet) {
        paProbeMap[gene] = false
    } // initialze paProbeMap

    lociSet.each { locus ->
		// if locus is an interlocus, expand for the sake
		// of getting the framework loci via the intergenes
        paProbeMap[locus] = true
        nomLocus = nomenclatureMap[locus]
        err.println "nomLocus=${nomLocus}, locus=${locus}"//todo
        if(nomLocus != null) {
            locus = nomLocus
        }
        allProbeMap.add(locus)
    } // each locus probe

    if(debugging <= 1) {
        err.println "makeGenotypeMaps: return ${paProbeMap}, ${allProbeMap}"
    }
    return [paProbeMap, allProbeMap]
} // makeGenotypeMaps

/*
 * Given a gene, return it.
 * Given an intergene,
 *   if it doesn't contain a framework gene, return it
 *   if it does contain a framework gene, return it 
 *     and the framework gene
 *
 */
ArrayList<String> checkInterGeneForFramework(String locus) {
	ArrayList<String> ret = new ArrayList()
	ret.add(locus)
	if(locus.contains("-")) {
		(gene1, gene2) = locus.split('-')
		if(gene1.contains("3DL3") || gene1.contains("3DP1") ||
		   gene1.contains("2DL4") || gene1.contains("3DL2")) {
			ret.add(gene1)
		}
		if(gene2.contains("3DL3") || gene2.contains("3DP1") ||
		   gene2.contains("2DL4") || gene2.contains("3DL2")) {
			ret.add(gene2)
		}
	} // intergene
	return ret
} // checkInterGeneForFramework

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
} // hapPairs2Genotypes (todo: remove)

/*
 * interpretPA
 *
 * Interpret a genotype in the context of haplotype pairs.
 * It assumes that the genotypes just include the genes that are present.
 *
 * @param hapPairTable Table<hap pair names, locus names, gene count>
 * @param genPAMap HashMap _all_ P/A gene genotype from genes -> Boolean
 * @return a Map from haplotype pairs to a Map of loci and their Boolean values; this is a potentially ambiiguous P/A interpretPAation
 * todo make a more elegant (e.g., vector-based) approach
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
    Map<String, Map<String, Boolean>> rowHapMap = hapPairTable.rowMap()
    rowHapMap.each { haplotype, hapColMap ->
        if(debugging <= 1) {
            err.println "interpretPA: comparing with haplotype ${haplotype}"
			//todo err.println "interpretPA: ${haplotype} ${hapColMap}"
        }

        Boolean allFound = true
        hapColMap.each { locus, value ->
            if(allFound == false) { // hap pair is false if any locus is false
                return
            }
			// don't intepret intergenes, just genes
			// also, skip half haplotypes, e.g., cA01, tB01, etc
			if(locus.contains("-") ||
			   locus.startsWith("c") || locus.startsWith("t")) {
				return
			}
            genVal = genPAMap[locus]
            same = false
            if(((value == true) && (genVal == true)) ||
               (value == false) && ((genVal == null) || (genVal == false))) {
                    same = true
            } else {
				if(debugging <= 1) {
					err.println "interpretPA: $locus not same: $value (hap) and $genVal (genotype)"
				}
			}
            allFound = same
        }

        if(allFound) { 
            if(debugging <= 1) {
                err.println "interpretPA: true"
            }
            retMap[haplotype] = hapColMap
        } else {
            if(debugging <= 1) {
                err.println "interpretPA: false"
            }
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
 * @param refHaps Set of all names of the reference haplotypes
 * @return Set of haplotype pair predictions (e.g., cA01~tA01+cB01~tA01, ...)
 */
HashSet<String> interpretHapMarkers(HashSet<String> g,
									Set<String> refHaps) {
    if(debugging <= 1) {
        err.println "interpretHapMarkers(g=${g.join(",")}, refHaps=${refHaps.join(",")}"
    }
    HashSet<String> retSet = new HashSet()

	g.each { gHit ->
		if(refHaps.contains(gHit)) {
			retSet.add(gHit)
		}
	}

    HashSet<String> retSetPairs = new HashSet()
	if(retSet.size() == 1) {
		hap = retSet.toList()[0]
		// don't homozygousify (for now)
        //retSetPairs.add([hap, hap].join('+'))
		retSetPairs.add(hap)
    } else { 
		retSet.each { hap1 ->
			retSet.each { hap2 ->
				if(hap2 <= hap1) {
					return
				}
				hl = [hap1, hap2].sort().join('+')
				retSetPairs.add(hl)
			}
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
    CliBuilder cli = new CliBuilder(usage:'pa2Haps4.groovy [options] ', header:'Options:')
    cli.help('print this message')
    cli.h(longOpt:'haps', args:1, argName:'file', 'file with haplotype definitions',
      required: true)
    cli.a(longOpt:'all', args:1, argName:'boolean', 'all output if set (debugging)',
      required: false)
    cli.q(longOpt:'qout', args:1, argName:'file', 'file with results of probe queries',
      required: true)
    cli.o(longOpt:'output', args:1, argName:'file', 'output file containing haplotype pair predictions',
      required: true)
    OptionAccessor options = cli.parse(args)
    return options
} // handleArgs

