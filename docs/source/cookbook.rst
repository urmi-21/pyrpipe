======================
Cookbook
======================

.. contents::

Using SRA objects
----------------------

.. code-block::

    from pyrpipe.sra import SRA #imports the SRA class
    
    #create an SRA object using a valid run accession
    """
    this checks if fastq files already exist in the directory,
    otherwise downloads the fastq files and stores the path in the object
    """
    myob=SRA('SRR1168424',directory='./results')
    
    #create an SRA object using fastq paths
    myob=SRA('SRR1168424',fastq='./fastq/SRR1168424_1.fastq',fastq2='./fastq/SRR1168424_2.fastq')
    
    #create an SRA object using sra path
    myob=SRA('SRR1168424',sra='./sra/SRR1168424.sra')
    
    #accessing fastq files
    print(myob.fastq,myob.fastq2)
    
    #check if fastq files are present
    print (myob.fastq_exists())
    
    #check sra file
    print (myob.sra_exists())
    
    #delete fastq
    myob.delete_fastq()
    
    #delete sra
    myob.delete_sra()
    
    #download fastq 
    myob.download_fastq()
    
    #trim fastq files
    myob.trim(qcobject)
    
    
Using Mapping objects
----------------------

.. code-block::

    from pyrpipe.mapping import Star
    
    #create a star object
    star=Star(index='path_to_star_index')
    
    #perform alignment using SRA object
    bam=star.perform_alignment(sraobject)
    #or
    sraobject.align(star)
    bam=sraobject.bam_path
    
    #execute STAR with any arguments and parameters
    kwargs={'--outFilterType' : 'BySJout',
            '--runThreadN': '6',
            '--outSAMtype': 'BAM SortedByCoordinate',
            '--readFilesIn': 'SRR3098744_1.fastq SRR3098744_2.fastq'
            }
    star.run(**kwargs)
    
    
    
    
    
    
    
    
    
    
    
