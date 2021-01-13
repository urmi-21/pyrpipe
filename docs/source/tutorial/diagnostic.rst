pyrpipe_diagnostic
-------------------

The pyrpipe_diagnostic command allows users to easily examine pyrpipe logs.
It provides several options.

1. report: This command can generate a summary or detailed report of the analysis.
2. shell: This command creates a bash file containing all the commands executed via pyrpipe
3. benchmark: This command can generate benchmarks to compare walltimes of the different commands in the pipeline.
4. multiqc: This command uses the `MultiQC <https://multiqc.info/>`_  tool to generate a report using pyrpipe logs  and other logs geenrated by the pipeline.