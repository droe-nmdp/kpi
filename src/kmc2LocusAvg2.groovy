#!/usr/bin/env groovy

/*
 * kmc2LocusAvg2
 *
 * Convert a kmc output file to probe hits per gene -- via average.
 *
 * input: tsv of probe sequence to count
 * output: fasta file per gene containing sequences that hit (e.g. '2DL4.bin1')
 *   The loci with '/' in their names is converted to 'z'
 *     e.g., 7/24 -> 7z24 (e.g., 53a_7z24.bin1)
 *
 * usage: kmc2LocusAvg2.groovy [options]
 * Options:
 *  -e,--extension <extension>               extension for output files
 *  -help                                    print this message
 *  -i,--ID <id>                             ID for the individual
 *  -j,--kmc probe results <kmc>             input kmc counts
 *  -o,--directory to put the output <out>   output directory
 *  -p,--probes <probes>                     input probes
 *
 * @author Dave Roe
 *
 */

import groovy.io.*
import groovy.util.CliBuilder.*
import groovy.util.OptionAccessor
// http://apache.mirrors.pair.com//commons/math/binaries/commons-math3-3.6.1-bin.tar.gz
import org.apache.commons.math3.stat.StatUtils

// things that may change per run
debugging = 3 // TRACE=1, DEBUG=2, INFO=3
// ignore kmers with a count < this (currently not using)
//minKmers = 2 //todo (coded, but commented out)
mode2DL3 = false  // if you want to analyze/debug a single gene
kmcFasta = false   // kmc output format is fasta or text (false)
probeFasta = true   // probe file format is fasta (true) or text (false)

// things that probably won't change per run
err = System.err
fileSeparator = System.getProperty('file.separator')
OptionAccessor options = handleArgs(args)
String kmcFileName = options.j
String id = options.i
String extension = options.e
String outputDir = options.o
String probeFileName = options.p

// convert kmc output file to fasta format (.bin1 extension)
kmc2Fasta(kmcFileName, probeFileName, outputDir, id, extension,
		  kmcFasta, probeFasta)
err.println "done"

/*
 * kmc2Fasta
 *
 * Convert Kmc output to fasta format -- one file per gene containing
 * all the hit probes.
 *
 * @param probesHitsFile String containing the name of the Kmc output file: tab-separated probe DNA then count (including 0)
 * @todo modularize
 */
int kmc2Fasta(String probeHitsFile, String probeFileName, String outputDir, 
              String id, String extension, Boolean kmcFasta, Boolean probeFasta) { 
    if(debugging <= 1) { 
        err.println "kmc2Fasta: probeHitsFile=${probeHitsFile}, id=${id}, outputDir=${outputDir}, extension=${extension}"
    }

    int retval = 0
	// map: probe -> Set of loci
    HashMap<String, TreeSet<String>> locusProbeMap = loadProbeMap(probeFileName,
																  probeFasta) 
	if(debugging <= 1) {
		err.println "kmc2Fasta: ${locusProbeMap.keySet().size()} locus/probes in ${probeFileName} (locusProbeMap)"
	}
    // locus -> probe seq -> count
	// will build this as we loop through the kmc file 
   HashMap<String, HashMap<String, Integer>> locusProbeHitMap = new HashMap()
	// locus -> list of hit counts (including zeros)
	HashMap<String,ArrayList<Integer>> locusHitListMap = new HashMap()
    if(debugging <= 1) { 
        err.println "kmc2Fasta: probeHitsFile=${probeHitsFile}, outputDir=${outputDir}, extension=${extension}"
    }
	
    FileReader kmcReader = new FileReader(new File(probeHitsFile))
    Integer count = 0
	//populate the two data structures with the marker in the line (all ids)
	kmc2FastaLine(kmcReader, locusProbeMap, locusProbeHitMap,
				  locusHitListMap, kmcFasta)

    if(debugging <= 1) { 
        err.println "kmc2Fasta: ${locusProbeHitMap.size()} loci in locusProbeHitMap"
    }
	// add the zeros back to locusHitListMap
	locusProbeMap.keySet().each { probe ->
		locusSet = locusProbeMap[probe] // set of loci defined for this probe
		zeroCount = 0
		locusSet.each { loc ->
			Map pMap = locusProbeHitMap[loc]
			Integer pcount = 0
			Integer prcCount = 0
			if(pMap != null) {
				pcount = pMap[probe]
				prcCount = pMap[reverseComplement(probe)]
			}
			ArrayList plist = locusHitListMap[loc]
			if(plist == null) {
				plist = new ArrayList()
			}
			if(pcount == null) {
				pcount = 0
			}
			if(prcCount == null) {
				prcCount = 0
			}
			if((pcount == 0) && (prcCount == 0)) {
				plist.add(0)
				zeroCount++
			}
		} // each locus defined for this probe
		if((zeroCount > 0) && (debugging <= 1)){ 
			err.println "kmc2Fasta: added $zeroCount zeros for $probe in " +
				locusSet.join(",")
		}
	} 

    // output hits per locus with 'extension' extension
	err.println "outputting ${extension} files ..."
	loci = locusProbeHitMap.keySet()
	err.println "${loci.size()} loci in locusProbeHitMap"
    loci.each { loc ->
		/*if((loc != "2DL3") && (loc != "2DL1")) {
			return //todo(remove)
		}*debugging*/
		//err.println loc //todo
		// replace slash in locus (e.g., 7/24) with 'z'
		escapedLoc = loc.replaceAll("/", "z")
        fullFileName = "${outputDir}${fileSeparator}${id}_${escapedLoc}.${extension}"
        ArrayList plist = locusHitListMap[loc]
        if(plist == null) {
			if(debugging <= 2) {
				err.println "no hit list for ${loc}"
			}
			return
		}

		if(debugging <= 2) {
			err.println "kmc2Fasta: $loc plist=" + plist
		}
		double[] pListD = plist.toArray()
		double[] avgList = StatUtils.mode(pListD)
		if(debugging <= 2) {
			err.println "avgList=" + avgList
		}
        if((avgList == null)  || (avgList.size() == 0)){
			return
		}
		Float avg = 0
		if(avgList.size() > 1) {
			// mean the modes just to get one number
			avg = StatUtils.mean(pListD) 
			err.println "reducing the mode to one number: mean=${avg}"
		} else {
			avg = avgList[0]
			if(debugging <= 2) {
				err.println "kmc2Fasta: setting to first avgList: " + avgList[0]
			}
		}

		/* this was for some initial debugging using 2DL3;
		 * it is not needed anymore */
		// if the mode is zero, but over 200 hits, take the mode of the hits
		if((loc == "2DL3") && (mode2DL3)) {
			// err.println "in 2DL3 mode"
			// 
			if((avg == 0) && (preZero2DL3List.size() > 200)) {
				avgList = StatUtils.mode(preZero2DL3List)
			}
			if((avgList == null)  || (avgList.size() == 0)){
				return
			}
			avg = 0
			if(avgList.size() > 1) {
				avg = StatUtils.mean(preZero2DL3List)
				err.println "mean=${avg}"
			} else {
				avg = avgList[0]
				err.println "setting to first avgList(mode): " + avgList[0]
			}
            //err.println "droe1 " + plist  //todo
		}
		
		if(avg == 0) {
			return
		}

        outWriter = new PrintWriter(new File(fullFileName).newOutputStream(),
									true)
		// non-null plist from here down
        plist.each { hitCount ->
            if(debugging <= 1) { 
                err.println "kmc2Fasta: writing ${loc} ${hitCount}"
            }
			// this is in the file name
            //outWriter.println ">${loc}"
            outWriter.println "${hitCount}"
        }
        outWriter.close()
    } // each locus

    if(debugging <= 1) { 
        err.println "kmc2Fasta: return ${retval}"
    }
    return retval
} // kmc2Fasta

/*
 * kmc2FastaLine
 *
 * Process each line of the kmc input file. Populate the two data structures.
 * 
 * @param line a line from the kmc input file
 *   e.g. marker  label   pvalue 100a    100b    100c ...
 * @param locusProbeMap Map: probe -> Set of loci
 * @param locusProbeHitMap locus -> probe seq -> count
 * @param locusHitListMap Map: locus -> list of hit counts (including zeros)
 */
def kmc2FastaLine(FileReader kmcReader, HashMap<String, TreeSet<String>> locusProbeMap,
				  HashMap<String, HashMap<String, Integer>> locusProbeHitMap,
				  HashMap<String,ArrayList<Integer>> locusHitListMap,
				  Boolean kmcFasta) {
	if(debugging <= 1) {
		err.println "kmc2FastaLine()"
	}

	// todo: modularize this
	// fasta format
	// e.g., >3DP1-2DL4
	//        AACATAAGCCAGTAGAATAGCATCT
	String probe = null
	Integer count = 0
    while(line = kmcReader.readLine()) {
		if(debugging <= 1) {
			err.println "kmc2FastaLine: line=$line"
		}
		if(!line || (line == "")) {
			continue
		}
		if(kmcFasta == true) {
			locus = line[1..-1].trim()
			if(debugging <= 1) {
				err.println "kmc2FastaLine: setting locus to $locus"
			}

			line = kmcReader.readLine()  // read the sequence row
			probe = line.trim()
			if(debugging <= 2) { 
				err.println "kmc2FastaLine: $locus $probe"
			}
			count = 10
		} else {
			// kmc text format
			(probe, countStr) = line.split('\t')
			count = countStr.toInteger()
			if(debugging <= 2) { 
				err.println "kmc2FastaLine: non-fasta: $probe $count"
			}
		}
		// here is where the kmc output is standardized to the
		// complementarity in the probe input file
		probeRc = reverseComplement(probe)
		locusSet = locusProbeMap[probe]
		if(locusSet == null) {		
			// check reverse complement
			locusSet = locusProbeMap[probeRc]
			if(locusSet == null) {
				if(debugging <= 2) { 
					err.println "kmc2FastaLine: no locus set for $probe"
				}
				continue
			} else {
				probe = probeRc
			}
		}
		if(debugging <= 2) { 
			err.println "kmc2FastaLine: $probe: " + locusSet.join(", ")
		}
		locusSet.each { locus ->
			if(debugging <= 1) {
				err.println "kmc2FastaLine: hit ${probe}"
			}

			// list of probes and their hits per locus
			ArrayList hitList = locusHitListMap[locus]
			Map probeHitMap = locusProbeHitMap[locus]
			/* this can get huge
			if(debugging <= 2) {
				err.println "hitList=$hitList"
				(huge) err.println "probeHitMap=" + probeHitMap
			}
  		    */
			if(probeHitMap == null) { 
				probeHitMap = new HashMap()
				hitList = new ArrayList()
			}
			currentCount = probeHitMap[probe]
			if(currentCount != null) {
				count = count + currentCount // forward and reverse complement
			}
			hitList.add(count)
			locusHitListMap[locus] = hitList
			
			probeHitMap[probe] = count
			locusProbeHitMap[locus] = probeHitMap
		} // each locus for a probe
	} // each line
	if(debugging <= 1) {
		err.println "kmc2FastaLine: return"
	}
} // kmc2FastaLine

/*
 * readPrimers
 *
 * Read the file containing the information on the primers.
 *
 * @param primerFileName a String containing the full path to the primer file
 * @return a Map: locus names to array of 
 *    Expandos (locus, fName, fSeq, rName, rSeq, bp)
 *      bp is the total distance between the two primers, including the primers
 * @todo remove
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
        // skip header and blank lines
        if((line == null) || line.contains("name") || (line.length() == 0)) {
            return
        }
		//err.println line //todo
        p = new Expando()
        // e.g., 2DL1	Fa517	gttggtcagatgtcatgtttgaa	Rc621	cctgccaggtcttgcg	142
        (fName, fSeq, rName, rSeq) = line.split()
        (locus, rest) = fName.split("_")
        p.locus = locus
        p.fName = fName
        p.fSeq = fSeq.toUpperCase()
        p.rName = rName
        p.rSeq = rSeq.toUpperCase()
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

/*
 * loadProbeMap
 *
 * @param probeListFileName a reduced and sorted version of the exact output matrix: marker  label   pvalue 100a    100b    100c ...
 * @return Map: probe -> Set of loci
 */
HashMap<String,TreeSet<String>> loadProbeMap(String probeListFileName,
											 Boolean probeFasta) {
	if(debugging <= 1) {
		err.println "loadProbeMap(probeListFileName=${probeListFileName}, probeFasta=$probeFasta"
	}
    // open file with probes
    FileReader probeReader = new FileReader(new File(probeListFileName))
	// map: probe -> Set of loci
    HashMap<String,TreeSet<String>> probeMap = new HashMap()
 	String locus = null
	String probe = null
    probeReader.eachLine { line ->
		if(probeFasta == true) {
			start = line.startsWith('>')
			if(start == true) {
				locus = line[1..-1].trim()
			} else {
				probe = line.trim()
				if(debugging <= 2) { 
					err.println "fasta: $locus $probe"
				}
				addToProbeMap(probeMap, locus, probe)
				locus = null
				probe = null
			}
		} else {
			// kmc text format
			(locus, probe) = line.split('\t')
			if(debugging <= 2) { 
				err.println "non-fasta: $locus $probe"
			}
			addToProbeMap(probeMap, locus, probe)
		}
    } // each line of file

    probeReader.close()
	if(debugging <= 3) {
		err.println "loadProbeMap: return ${probeMap.keySet().size()} probe/regions"
	}

	//err.println probeMap //todo
    return probeMap
} // loadProbeMap

void addToProbeMap(HashMap<String,TreeSet<String>> probeMap,
				   String locus, String probe) {
	s = probeMap[probe]
	if(s == null) {
		s = new TreeSet()
		probeMap[probe] = s
	}
	s.add(locus)

} // addToProbeMap
/*
 * handleArgs
 * 
 * Parse and return the input.
 *
 * @param args List of Strings containing command line arguments.
 * @return Option
 */
OptionAccessor handleArgs(String[] args) { 
    CliBuilder cli = new CliBuilder(usage:'kmc2LocusAvg2.groovy [options] ',
      header:'Options:')
    cli.help('print this message')

    cli.j(longOpt:'kmc probe results', args:1, argName:'kmc', 'input kmc counts',
      required: true)
    cli.o(longOpt:'directory to put the output', args:1, argName:'out', 
      'output directory', required: true)
    cli.i(longOpt:'ID', args:1, argName:'id', 'ID for the individual',
      required: true)
    cli.e(longOpt:'extension', args:1, argName:'extension', 'extension for output files',
      required: true)
    cli.p(longOpt:'probes', args:1, argName:'probes', 'input probes',
      required: true)

    OptionAccessor options = cli.parse(args)
    return options
} // handleArgs
