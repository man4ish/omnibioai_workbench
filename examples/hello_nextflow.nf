#!/usr/bin/env nextflow

/*
 * Minimal Hello-World Nextflow workflow for OmnibioAI demo
 *
 * This workflow prints a message to stdout.
 */

process hello {
    /*
     * Standard output is captured and returned
     */
    output:
    stdout

    script:
    """
    echo "Hello from OmnibioAI workflow service"
    """
}

workflow {
    hello()
}

