#!/usr/bin/env groovy

/*
 * pa2Haps
 *
 * Given a FASTQ file, bin the reads into separate files based on probe markers.
 *
 * Requires BBTools' bbduk.
 * 
 * @author Dave Roe
 */

import groovy.io.*
import groovy.util.CliBuilder.*
import groovy.util.OptionAccessor

// things that may change per run
debugging = 1 // TRACE=1, DEBUG=2, INFO=3

// things that probably won't change per run
err = System.err
newline = System.getProperty('line.separator')


OptionAccessor options = handleArgs(args)

// direct the reads into '*.bin2'
binReads(options.i, options.p, options.s)
// end main

/*
 * binReads
 *
 * For each probe, use bbduk to bin (extract) fastq reads per probe.
 * 
 * e.g., bbduk.sh	in=/Users/droe/data/kir/simulations/KP420443_KP420444/KP420443_KP420444.fastq	outm=6T.bin1	literal=TAGAACCAGGGATGGAGAGAGATTA	k=25 maskmiddle=f overwrite=t
 * bbduk.sh in=/Users/droe/data/kir/simulations/KP420443_KP420444/KP420443_KP420444.fastq out=unmapped.bin1 ref=/Users/droe/Dropbox/doc/kir/snp/locus-hap_probes_v1.fasta k=25 maskmiddle=f overwrite=t
 */
def void binReads(String fqFileName, String probeFileName, String suffix) {
    if(debugging <= 1) {
        err.println "binReads(fqFileName=${fqFileName}, suffix=${suffix})"
    }

    // open file with probes
    FileReader probeReader = new FileReader(new File(probeFileName))

    probeReader.eachLine { line ->
        (pname, pseq) = line.split('\t')
        cmd = ["bbduk.sh", "in=${fqFileName}", "outm=tmp.bin1",
               "literal=${pseq}", "k=25 maskmiddle=f overwrite=t fastawrap=300000"]
        if(debugging <= 3) {
            err.println cmd.join(' ')
        }
        ret = cmd.execute()
        ret.waitFor()
        bbdRetVal = ret.exitValue()

        // append the output to the correct file
        File src = new File("tmp.bin1")
        String dir = new File("").getAbsolutePath()
        File dest = new File("${pname}.${suffix}")
        src.eachLine { sline -> // append
            dest << sline
            dest << newline
        }

    } // each probe
    probeNameFastaFile = probeFileName.replace(".txt", ".fasta")
    cmd = ["bbduk.sh", "in=${fqFileName}", "out=unmapped.${suffix}",
           "ref=${probeNameFastaFile}", "k=25 maskmiddle=f overwrite=t fastawrap=300000"]
    if(debugging <= 3) {
        err.println cmd.join(' ')
    }
    ret = cmd.execute()
    ret.waitFor()
    bbdRetVal = ret.exitValue()

    // clean up
    cmd = ["rm", "-f", "tmp.bin1"]
    ret = cmd.execute()
    ret.waitFor()

    if(debugging <= 1) {
        err.println "binReads: return"
    }
    return
} // binReads

/*
 * handleArgs
 * 
 * Parse and return the input.
 *
 * @param args List of Strings containing command line arguments.
 * @return Option
 */
OptionAccessor handleArgs(String[] args) { 
    CliBuilder cli = new CliBuilder(usage:'probeBinFastqs.groovy [options] ',
      header:'Options:')
    cli.help('print this message')
    cli.i(longOpt:'in', args:1, argName:'fastq', 'input FASTQ file',
      required: true)
    cli.s(longOpt:'suffix', args:1, argName:'suffix', 'suffix for output files',
      required: true)
    cli.p(longOpt:'probes', args:1, argName:'probes', 'input probes',
      required: true)
    OptionAccessor options = cli.parse(args)
    return options
} // handleArgs
