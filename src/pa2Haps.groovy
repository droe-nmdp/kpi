#!/usr/bin/env groovy

/*
 * pa2Haps
 *
 * Fit PA genotypes to (potentially ambiguous) haplotype pairs
 * 
 * e.g., ./pa2Haps.groovy input/all_haps_v1.txt ~/repos/dev/bioinformatics/projects/ngs/scripts/groovy/output/ccs_99_vilches.txt output/ccs_99_out.txt
 * 
 * @author Dave Roe
 * @todo change the way the csv is parsed
 * @todo add command line handling
 */

import groovy.io.*
import org.supercsv.io.*
import org.supercsv.prefs.*
import org.supercsv.cellprocessor.ift.*

// things that may change per run
debugging = 5
// column names in the haplotypes file; others are the loci
hapColumnName = "haplotype"
nomColumnName = "nomenclature"
freqColumnName = "freq"
// column names in the genotype file
locusColumnName = "locus"
genotypeColumnName = "P/A"

// thing that probably won't change per run
err = System.err

// open haplotypes file
hapFile = new File(args[0])
FileReader hapReader = new FileReader(hapFile)
// open genotype file
genFile = new File(args[1])
FileReader genReader = new FileReader(genFile)
// open output file
outFile = new File(args[2])
outWriter = new PrintWriter(outFile.newOutputStream(), true)

/* 
 * Due to the properties of the Groovy Map '+' operation, a locus from one
 * haplotype overwrites the other, so false values must be null in the Maps.
 */
// map: haplotype name -> HashMap[locus:Boolean]
HashMap<String,HashMap<String,Boolean>> hapMap = null
// map: haplotype name -> nomenclature name
HashMap<String,String> nomMap = null
// map haplotype name -> frequency
HashMap<String,Float> freqMap = null

(hapMap, nomMap, freqMap) = readHaplotypes(hapReader)
genMap = readGenotypes(genReader)

err.println "${hapMap.size()} haplotypes"
err.println "${genMap.size()} loci"

// haplotype pairs -> map[locus: boolean]
HashMap<String, HashMap<String,Boolean>> hapPairMap = null
hapPairMap = hapPairs2Genotypes(hapMap, nomMap)
err.println "${hapPairMap.size()} pairs of reference haplotypes"

// haplotype pairs -> map[locus: boolean]
interpMap = interpret(genMap, hapPairMap)
err.println "${interpMap.size()} interpretation(s) "

writeOutput(outWriter, interpMap, hapMap, nomMap)

// end main

/*
 * writeOutput
 *
 * Send the haplotype pair predictions to the output file
 * and some to stderr. Can't output the order of genes, because
 * we don't know it for most of the haplotypes.
 */
def void writeOutput(outWriter, interpMap, hapMap, nomMap) {
    // output the haplotype list via haplotype names
    outWriter.print "haplotype list: "
    c = 1
    interpMap.keySet().sort().each { hapPair ->
        if(c++ == 2) {
            outWriter.println "+"
        }
        //err.println "hapPair=${hapPair}"
        (h1Name, h2Name) = hapPair.split("\\+")
        (h1Nom, h2Nom) = [nomMap[h1Name], nomMap[h2Name]]
        if(!h1Nom) {
            h1Nom = h1Name
        }
        if(!h2Nom) {
            h2Nom = h2Name
        }
        (h1Nom, h2Nom) = [h1Nom, h2Nom].sort()
        outWriter.println "${h1Nom}+${h2Nom}"
        err.println "${h1Nom}+${h2Nom}"
    } // each haplotype pair

// output the haplotype list via loci
    interpMap.keySet().sort().each { hapPair ->
        hapPair.split("\\+").each { hapName ->
            hopNom = nomMap[hapName]
            if(!hopNom) {
                hopNom = hapName
            }
            outWriter.print "${hopNom} (unordered): "
            outWriter.println hapMap[hapName].keySet().sort().join("^")
            err.println hapMap[hapName].keySet().sort().join("^")
        }
    }
} // writeOutput


/*
 * readHaplotypes
 *
 * column headers: haplotype	nomenclature	freq	3DL3	2DS2	2DL2	2DL3	2DP1	2DL1	3DP1	2DL4	3DL1	3DS1	2DL5	2DS3	2DS5	2DS4	2DS1	3DL2
 *
 * @param hapFile tab-separated file containing the reference haplotypes, their nomenclature name and frequencies
 * @return three Maps; haplotype names to Map of loci values, nomenclature, frequencies
 *   map: haplotype name -> HashMap (locus -> Boolean)
 *   map: haplotype name -> nomenclature name
 *   map: haplotype name -> frequency
 */
def List<Map> readHaplotypes(FileReader hapFile) { 
    ICsvListReader listReader = null
    listReader = new CsvListReader(hapFile, CsvPreference.TAB_PREFERENCE)
    // get all headers from the first row; arg is redundant
    String[] headers = listReader.getHeader(false)
    if(debugging <= 3) {
        err.println headers.join(',')
    }
    numLoci = headers.size()-3
    nameColIndex = headers.findIndexOf{it =~ hapColumnName}
    nomColIndex = headers.findIndexOf{it =~ nomColumnName}
    freqColIndex = headers.findIndexOf{it =~ freqColumnName}
    // map: haplotype name -> HashMap[locus:Boolean]
    HashMap<String,HashMap<String,Boolean>> hapMap = new HashMap(numLoci)
    // map: haplotype name -> nomenclature name
    HashMap<String,String> nomMap = new HashMap(numLoci)
    // map haplotype name -> frequency
    HashMap<String,Float> freqMap = new HashMap(numLoci)
    while( (line = listReader.read()) != null ) {
        hap = line[nameColIndex]
        nom = line[nomColIndex]
        if(nom != null) { 
            nomMap[hap] = nom
        }
        freqMap[hap] = line[freqColIndex]
        HashMap<String,Boolean> lociMap = new HashMap(numLoci)
        // each locus
        line[freqColIndex+1..-1].eachWithIndex { col, i ->
            i += freqColIndex+1
            locus = headers[i]
            genotype = col.toBoolean()
            if(genotype) { // only keep true ones
                lociMap[locus] = genotype
            }
        } // each column
        hapMap[hap] = lociMap
    } // each line

    return [hapMap, nomMap, freqMap]
} // readHaplotypes

/*
 * readGenotypes
 *
 * column headers: locus	P/A
 *
 * @param genFile tab-separated file containing the loci and their p/a values
 * @return a Map of loci to Boolean
 */
def HashMap<String,Boolean> readGenotypes(FileReader genFile) { 
    ICsvListReader listReader = null
    listReader = new CsvListReader(genFile, CsvPreference.TAB_PREFERENCE)
    // get all headers from the first row
    String[] headers = listReader.getHeader(false)
    if(debugging <= 3) {
        err.println headers.join(',')
    }
    locusColIndex = headers.findIndexOf{it =~ locusColumnName}
    genColIndex = headers.findIndexOf{it =~ genotypeColumnName}
    // map: locus -> Boolean
    HashMap<String,Boolean> genMap = new HashMap()

    while( (line = listReader.read()) != null ) {
        locus = line[locusColIndex]
        genotype = line[genColIndex].toBoolean()
        if(genotype) { // only keep the true ones
            genMap[locus] = genotype
        }
    } // each line

    return genMap
} // readGenotypes

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
        //err.println "h1Name=${h1Name}"//todo
        hapMap.each { h2Name, h2Map ->
            if(h1Name > h2Name) { // avoid dups and sort by name
                return
            }
            (h1Name, h2Name) = [h1Name, h2Name].sort()
            genMap["${h1Name}+${h2Name}"] = h1Map + h2Map
        }
    }
    return genMap
} // hapPairs2Genotypes

/*
 * interpret
 *
 * Interpret a genotype in the context of haplotype pairs.
 * It assumes that the genotypes just include the genes that are present.
 *
 * @param hapMap Map of haplotype names to a Map of loci and their Boolean values
 * @return a Map from haplotype pairs to a Map of loci and their Boolean values
 */
def HashMap<String, HashMap<String,Boolean>> interpret(HashMap<String,Boolean> genMap,
                                                       HashMap<String, HashMap<String,Boolean>> pairsMap) {
    HashMap<String, HashMap<String,Boolean>> retMap = new HashMap()
    gSize = genMap.size()
    //err.println "genMap = ${genMap.keySet()}"//todo
    pairsMap.each { pairsName, pairs ->
/*        if(!(pairsName =~ /cA03/)) {//debugging(todo)
           return
           }*/
        if((pairs.size() == gSize) &&
           pairs.keySet().containsAll(genMap.keySet())) {
            retMap[pairsName] = pairs
        }
    } // each pair of haplotypes
    return retMap
} // interpret
