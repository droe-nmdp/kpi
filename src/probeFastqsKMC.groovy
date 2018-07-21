#!/usr/bin/env groovy

/*
 * probeFastqsKMC
 *
 * Given a set of FASTQ files, create a KMC 3 database for each ID.
 *
 * e.g., probeFastqsKMC.groovy -m samples_map.txt -p 25mers.fasta -o . -w work
 *
 * -m can be a directory with fq.gz files or file with two columns: name to file map
 * 
 * Requires KMC 3.
 *   http://sun.aei.polsl.pl/REFRESH/index.php?page=projects&project=kmc
 *   'kmc' to build database
 *   e.g., kmc -k25 -ci2 -fq @./gonl-52b-cmd.txt ./gonl-52b work3
 *
 * @author Dave Roe
 */

import groovy.io.*
import groovy.io.FileType
import groovy.util.CliBuilder.*
import groovy.util.OptionAccessor

// things that may change per run
debugging = 3 // TRACE=1, DEBUG=2, INFO=3
kmerSize = "25"
minKmers = "3" // reads hit less than this will be ignored

// things that probably won't change per run
err = System.err
fileSeparator = System.getProperty('file.separator')

OptionAccessor options = handleArgs(args)
if(debugging <= 4) {
    err.println "kmerSize=${kmerSize} minKmers=${minKmers}"
}    

// make list of fastq files for every individual
HashMap<String,ArrayList<String>> fqMap = loadFqMap(options.m,
													options.p)
if(debugging <= 2) {
    err.println "${fqMap.keySet().size()} IDs in the fastq map"
    firstKey = fqMap.keySet().iterator().next()
    err.println "${fqMap[firstKey].size()} fastq files for ${firstKey}"
}

probeHits = probeReads(options.o, options.w, fqMap, kmerSize, minKmers)
err.println "done"
// end main

/*
 * loadFqMap
 *
 * @param fqMapFileName file containing a tab-delimited mapping between ids and fastq files
 *
 * if input is dir, process every file as a single id; take id
 * from first file name up to the first '_' or '.' if no underscore
 */
HashMap<String,ArrayList<String>> loadFqMap(String fqMapFileName,
											String fpath) { 
	// return value
    HashMap<String,ArrayList<String>> fqMap = new HashMap()

    // open file with probes
	f = new File(fqMapFileName)

	String id = null
	if(f.isDirectory()) {
		f.eachInFileRecurse InFileType.INFILES,  { inFile ->
			if(inFile.name.endsWith(".fq") || inFile.name.endsWith(".fq.gz") ||
			   inFile.name.endsWith(".fastq") || inFile.name.endsWith(".fastq.gz") ||
			   inFile.name.endsWith(".fa") || inFile.name.endsWith(".fa.gz") ||
			   inFile.name.endsWith(".fasta") || inFile.name.endsWith(".fasta.gz")) {
				if(id == null) {
					id = substring(0, names.indexOf('_'))
					if(id == null) {
						id = substring(0, names.indexOf('.'))
					}
				}
			}
		}
	} else {
		FileReader probeReader = new FileReader(f)
		probeReader.eachLine { line ->
			if(debugging <= 1) {
				err.println line
			}
			(id, shortFileName) = line.split('\t')
			fileName = fpath + fileSeparator + shortFileName
			idList = fqMap[id]
			if(idList != null) { 
				idList.add(fileName)
			} else { 
				ArrayList<String> l = new ArrayList()
				l.add(fileName)
				fqMap[id] = l
			}    
		} // each line of file
		
		probeReader.close()
	}
    return fqMap
} // loadFqMap

/*
 * probeReads
 *
 * e.g., kmc -k25 -ci2 -fq @output/gonl-157a-cmd.txt output/gonl-157a work
 * 
 * For each individual, build a database from their fastq files.
 * 
 */
def void probeReads(String outputDir, String workDir, 
                    HashMap<String,ArrayList<String>> fqMap,
                    String kmerSize, String minKmers) {
    if(debugging <= 1) {
        err.println "probeReads(workDir=${workDir}, outputDir=${outputDir})"
    }

    fqMap.keySet().sort().each { id ->
        outFile = outputDir + fileSeparator + id
        // make the command file
        defFileName = outputDir + fileSeparator + id + "-cmd.txt"
        outWriter = new PrintWriter(new File(defFileName).newOutputStream(),
                                    true)
        fqMap[id].each { fq ->
            outWriter.println fq
        }
        outWriter.close()

		//todo: change back to fq
		//https://github.com/refresh-bio/KMC/issues/40
        cmd = ["kmc", "-k${kmerSize}", "-ci${minKmers}", "-fm",
               "@${defFileName}", outFile, workDir]
        if(debugging <= 3) {
            err.println cmd
        }
        ret = cmd.execute()
        ret.waitFor()
        retVal = ret.exitValue()
        if(retVal) {
            err.println "makeKMCdb: *error* running ${cmd}"
        }

        if(debugging <= 2) {
            err.println "probeReads: returned ${retVal}"
        }
    } // each id
    if(debugging <= 1) {
        err.println "probeReads: return"
    }
} // probeReads

        
/*
 * handleArgs
 * 
 * Parse and return the input.
 *
 * @param args List of Strings containing command line arguments.
 * @return Option
 */
OptionAccessor handleArgs(String[] args) { 
    CliBuilder cli = new CliBuilder(usage:'probeFastqsKMC.groovy [options] ',
      header:'Options:')
    cli.help('print this message')

    cli.m(longOpt:'fastq name map or directory name (one ID only)', args:1,
		  argName:'map', 'fastq map', required: true)
    cli.p(longOpt:'path to the sequences files', args:1,
		  argName:'path', 'fastq path', required: true)
    cli.o(longOpt:'directory to put the output', args:1, argName:'out', 
		  'output directory', required: true)
    cli.w(longOpt:'work directory', args:1, argName:'work', 
          'work directory', required: true)

    OptionAccessor options = cli.parse(args)
    return options
} // handleArgs
