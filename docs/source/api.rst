************
pyrpipe API
************

pyrpipe is highly modular system with multiple modules: `sra`, `assembly`, `mapping`, `quant`, `qc`, and `tools`.
Each module represents a type of RNA-Seq process and contains several classes representing tools. These modules are described below:

sra
=======================

The :py:mod:`sra` module contains the SRA class. This class contains method to access RNA-Seq data from NCBI-SRA or disk.
This class depends on the `NCBI SRA-TOOLS` program.




assembly
=========
The :py:mod:`assembly` module contains classes for transcript assembly tools. Each class represent a transcript assembly tool and provides API to it.
Currently, available tools are StringTie, Cufflinks, and Trinity


mapping
=========
The :py:mod:`mapping` module contains classes for read alignment.
Currently available tools are Hisat2, STAR, and BowTie2


quant
=========
The :py:mod:`quant` module contains classes for transcript quantification tools.
Currently available tools are Salmon and Kallisto

qc
=========
The :py:mod:`qc` module contains classes for fastq quality control tools.
The `qc` objects can be directly use with `SRA` objects to perform quality filtering.
Currently available tools are Trim Galore and BBDuk


tools
=========
The :py:mod:`tools` module contains classes providing APIs for various RNA-Seq tools.
Currently available tools are samtools, mikado, portcullis, and Ribocode
