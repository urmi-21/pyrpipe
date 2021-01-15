#!/usr/bin/env nextflow

params.outdir = 'myout'
params.infile='runids.txt'
projectdir = projectDir

srrids= new File("runids.txt").collect {it}
//create a list of outfile paths
def qflist = []
srrids.each { val -> qflist.add("$projectdir/$params.outdir/$val/salmon_out/quant.sf")}
qflist.each { println it }

process rnaseq {

    publishDir "$projectdir/$params.outdir"

    input:
    val srr from srrids
    output:
    file "${srr}*" into result
    
    script:
    """
    #!/usr/bin/env python3
    from pyrpipe import sra,mapping,quant
    salmon=quant.Salmon(index="$projectdir/human_data/salmon_index")
    thisid='${srr}'
    print('Processing:',thisid)
    sra.SRA(thisid).quant(salmon).delete_fastq()
    """
}

process merge {
    input:
    file("*") from result.collect()
    output:
    stdout res

    """
    echo ${qflist}
    python ${projectdir}/scripts/merge.py ${qflist.join(",")} ${projectdir}/${params.outdir}
    """
}

res.subscribe { println it }