#!/usr/bin/env groovy

/*
 * probeFastqsKMC
 *
 * Given a file for a single individual, which points to a set of FASTQ files,
 * create a KMC 3 database for that ID.
 *
 * e.g., probeFastqsKMC.groovy -m samples_map.txt -p 25mers.fasta -o . -w work
 *
 * -m can be a directory with fq.gz files or file with two columns: name to file map
 * 
 * Requires KMC 3.
 *   http://sun.aei.polsl.pl/REFRESH/index.php?page=projects&project=kmc
 *   'kmc' to build database
 * Options:
 *  -help                                                       print this
 *                                                              message
 *  -m,--fastq name map or directory name (one ID only) <map>   fastq map
 *  -o,--directory to put the output <out>                      output
 *                                                              directory
 *  -w,--work directory <work>                                  work
 *                                                              directory
 * kmerSize=25 minKmers=3
 *
 *    e.g., kmc -k25 -ci2 -fq @./gonl-52b-cmd.txt ./gonl-52b work3
 *
 * @author Dave Roe
 * @todo add g option for multiple IDs per file or folder
 */

import groovy.io.*
import groovy.io.FileType
import groovy.util.CliBuilder.*
import groovy.util.OptionAccessor

// things that may change per run
debugging = 3 // TRACE=1, DEBUG=2, INFO=3
kmerSize = "25"
minKmers = "3" // < this will be ignored

// things that probably won't change per run
err = System.err
fileSeparator = System.getProperty('file.separator')

OptionAccessor options = handleArgs(args)
if(debugging <= 4) {
    err.println "kmerSize=${kmerSize} minKmers=${minKmers}"
}    

// make list of fastq files for every individual
HashMap<String,ArrayList<String>> fqList
if(options.m != null) { 
	fqList = loadFqMap(options.m)
} else if(options.g != null) {
	// todo
	fqList = loadFqMap(options.m))
}
HashMap<String,ArrayList<String>> fqList = loadFqMap(options.m)
if(debugging <= 2) {
    err.println "${fqList.keySet().size()} IDs in the fastq map"
    firstKey = fqList.keySet()[0]
    err.println "${fqList[firstKey].size()} fastq files for ${firstKey}"
}

probeHits = probeReads(options.o, options.w, fqList, kmerSize, minKmers)
err.println "done"
// end main

/*
 * loadFqMap
 *
 * @param fqListFileName file containing a tab-delimited mapping between ids and fastq files
 *
 * if input is dir, process every file as a single id; take id
 * from first file name up to the first '_' or '.' if no underscore
 */
HashMap<String,ArrayList<String>> loadFqMap(String fqListFileName) { 
	// return value
    HashMap<String,ArrayList<String>> fqMap = new HashMap()

    // open file with probes
	f = new File(fqListFileName)

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
			(id, fileName) = line.split('\t')
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
                    HashMap<String,ArrayList<String>> fqList,
                    String kmerSize, String minKmers) {
    if(debugging <= 1) {
        err.println "probeReads(workDir=${workDir}, outputDir=${outputDir})"
    }

    fqList.keySet().sort().each { id ->
        outFile = outputDir + fileSeparator + id
        // make the command file
        defFileName = outputDir + fileSeparator + id + "-cmd.txt"
        outWriter = new PrintWriter(new File(defFileName).newOutputStream(),
                                    true)
        fqList[id].each { fq ->
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
		  argName:'map', 'fastq map', required: false)
//    cli.g(longOpt:'fastq(.gz) name (assumes one file per id)', args:1,
//		  argName:'map', 'fastq map', required: false)
    cli.o(longOpt:'directory to put the output', args:1, argName:'out', 
		  'output directory', required: true)
    cli.w(longOpt:'work directory', args:1, argName:'work', 
          'work directory', required: true)

    OptionAccessor options = cli.parse(args)
    return options
} // handleArgs
