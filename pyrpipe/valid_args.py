#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  7 14:15:50 2020

@author: usingh

This module contains a list of valid arguments for the tools
"""

_args_FASTERQDUMP=['-e','-f','-t','-s','-N','-X','-a','-p','-c','-o','-O','-h','-V',
                   '-L','-v','-q','-b','-m','-x','-S','-3','-P','-M',
                   '-B','--option-file','--strict','--table','--include-technical',
                   '--skip-technical','--concatenate-reads']

_args_TRIM_GALORE=['--cores','-v','-q','--phred33','--phred64','--fastqc','--fastqc_args','-a','-a2',
                   '--illumina','--nextera','--small_rna','--consider_already_trimmed',
                   '--max_length','--stringency','-e','--gzip','--dont_gzip','--length',
                   '--max_n','--trim-n','-o','--no_report_file','--suppress_warn',
                   '--clip_R1','--clip_R2','--three_prime_clip_R1','--three_prime_clip_R2',
                   '--2colour','--path_to_cutadapt','--basename','-j','--hardtrim5','--hardtrim3',
                   '--clock','--polyA','--rrbs','--non_directional','--keep','--paired','-t',
                   '--retain_unpaired','-r1','-r2']


_args_BBDUK=['in','in2','ref','literal','touppercase','interleaved','qin','reads','copyundefined',
             'samplerate','samref','out','out2','outm','outm2','outs','stats','refstats','rpkm',
             'dump','duk','nzo','overwrite','showspeed','ziplevel','fastawrap','qout','statscolumns',
             'rename','refnames','trd','ordered','maxbasesout','maxbasesoutm','','json','bhist','qhist',
             'qchist','aqhist','bqhist','lhist','phist','gchist','ihist','gcbins','maxhistlen','histbefore',
             'ehist','qahist','indelhist','mhist','idhist','idbins','varfile','vcf','ignorevcfindels',
             'k','rcomp','maskmiddle','minkmerhits','minkmerfraction','mincovfraction','hammingdistance',
             'qhdist','editdistance','hammingdistance2','qhdist2','editdistance2','forbidn','removeifeitherbad',
             'trimfailures','findbestmatch','skipr1','skipr2','ecco','recalibrate','sam','le.','amino',
             'threads','prealloc','monitor','minrskip','maxrskip','rskip','qskip','speed','ktrim','kmask',
             'maskfullycovered','ksplit','mink','qtrim','trimq','trimclip','minlength','mlf','maxlength',
             'minavgquality','maqb','minbasequality','maxns','mcb','ottm','tp','tbo','strictoverlap',
             'minoverlap','mininsert','tpe','forcetrimleft','forcetrimright','forcetrimright2',
             'forcetrimmod','restrictleft','restrictright','mingc','maxgc','gcpairs','tossjunk',
             'swift','chastityfilter','barcodefilter','barcodes','xmin','ymin','xmax','ymax','trimpolya',
             'trimpolygleft','trimpolygright','trimpolyg','filterpolyg','pratio','plen','entropy','entropywindow',
             'entropyk','minbasefrequency','entropytrim','entropymask','entropymark','cardinality',
             'cardinalityout','loglogk','loglogbuckets','-Xmx','-eoom','-da']

_args_STAR=['--parametersFiles','--sysShell','--runMode','--runThreadN','--runDirPerm','--runRNGseed','--quantMode','--quantTranscriptomeBAMcompression','--quantTranscriptomeBan','--twopassMode','--twopass1readsN',
            '--genomeDir','--genomeLoad','--genomeFastaFiles','--genomeChrBinNbits','--genomeSAindexNbases','--genomeSAsparseD','--genomeSuffixLengthMax','--genomeChainFiles','--genomeFileSizes',
            '--sjdbFileChrStartEnd','--sjdbGTFfile','--sjdbGTFchrPrefix','--sjdbGTFfeatureExon','--sjdbGTFtagExonParentTranscript','--sjdbGTFtagExonParentGene','--sjdbOverhang','--sjdbScore','--sjdbInsertSave',
            '--inputBAMfile','--readFilesIn','--readFilesCommand','--readMapNumber','--readMatesLengthsIn','--readNameSeparator','--clip3pNbases','--clip5pNbases','--clip3pAdapterSeq','--clip3pAdapterMMp','--clip3pAfterAdapterNbases',
            '--limitGenomeGenerateRAM','--limitIObufferSize','--limitOutSAMoneReadBytes','--limitOutSJoneRead','--limitOutSJcollapsed','--limitBAMsortRAM ','--limitSjdbInsertNsj','--outFileNamePrefix','--outTmpDir','--outTmpKeep',
            '--outStd','--outReadsUnmapped','--outQSconversionAdd','--outMultimapperOrder','--outSAMtype','--outSAMmode','--outSAMstrandField','--outSAMattributes','--outSAMattrIHstart','--outSAMunmapped','--outSAMorder',
            '--outSAMprimaryFlag','--outSAMreadID','--outSAMmapqUnique','--outSAMflagOR','--outSAMflagAND','--outSAMattrRGline','--outSAMheaderHD','--outSAMheaderPG','--outSAMheaderCommentFile','--outSAMfilter','--outSAMmultNmax',
            '--outBAMcompression','--outBAMsortingThreadN','--bamRemoveDuplicatesType','--bamRemoveDuplicatesMate2basesN','--outWigType','--outWigStrand','--outWigReferencesPrefix','--outWigNorm','--outFilterType',
            '--outFilterMultimapScoreRange','--outFilterMultimapNmax','--outFilterMismatchNmax','--outFilterMismatchNoverLmax','--outFilterMismatchNoverReadLmax','--outFilterScoreMin','--outFilterScoreMinOverLread',
            '--outFilterMatchNmin','--outFilterMatchNminOverLread','--outFilterIntronMotifs','--outSJfilterReads','--outSJfilterOverhangMin','--outSJfilterCountUniqueMin','--outSJfilterCountTotalMin','--outSJfilterDistToOtherSJmin',
            '--outSJfilterIntronMaxVsReadN','--scoreGap','--scoreGapNoncan','--scoreGapGCAG ','--scoreGapATAC','--scoreGenomicLengthLog2scale','--scoreDelOpen','--scoreDelBase','--scoreInsOpen','--scoreInsBase','--scoreStitchSJshift',
            '--seedSearchStartLmax','--seedSearchStartLmaxOverLread','--seedSearchLmax','--seedMultimapNmax','--seedPerReadNmax','--seedPerWindowNmax','--seedNoneLociPerWindow','--alignIntronMin','--alignIntronMax','--alignMatesGapMax',
            '--alignSJoverhangMin','--alignSJstitchMismatchNmax','--alignSJDBoverhangMin','--alignSplicedMateMapLmin','--alignSplicedMateMapLminOverLmate','--alignWindowsPerReadNmax','--alignTranscriptsPerWindowNmax','--alignTranscriptsPerReadNmax',
            '--alignEndsType','--alignEndsProtrude','--alignSoftClipAtReferenceEnds','--winAnchorMultimapNmax','--winBinNbits','--winAnchorDistNbins','--winFlankNbins','--winReadCoverageRelativeMin','--winReadCoverageBasesMin',
            '--chimOutType','--chimSegmentMin','--chimScoreMin','--chimScoreDropMax','--chimScoreSeparation','--chimScoreJunctionNonGTAG','--chimJunctionOverhangMin','--chimSegmentReadGapMax','--chimFilter','--chimMainSegmentMultNmax']

_args_HISAT2=['-x','-1','-2','-U','--sra-acc','-S','-q','--qseq','-f','-r','-c','-s',
              '-u','-5','-3','--phred33','--phred64','--int-quals',
              '--sra-acc','--n-ceil','--ignore-quals','--nofw','--norc','--pen-cansplice',
              '--pen-noncansplice','--pen-canintronlen','--pen-noncanintronlen','--min-intronlen'
              ,'--max-intronlen','--known-splicesite-infile','--novel-splicesite-outfile',
              '--novel-splicesite-infile','--no-temp-splicesite','--no-spliced-alignment',
              '--rna-strandness','--tmo','--dta','--dta-cufflinks','--avoid-pseudogene',
              '--no-templatelen-adjustment','--mp','--sp','--no-softclip','--np','--rdg',
              '--rfg','--score-min','-k','-I','-X','--fr','--rf','--ff','--no-mixed',
              '--no-discordant','-t','--un','--al','--un-conc','--al-conc','--un-gz',
              '--summary-file','--new-summary','--quiet','--met-file','--met-stderr',
              '--met','--no-head','--no-sq','--rg-id','--rgit-sec-seq','-o','-p',
              '--reorder','--mm','--qc-filter','--seed','--non-deterministic',
              '--remove-chrname','--add-chrname','--version']

_args_BOWTIE2=['-x','-1','-2','-U','--interleaved','-S','-b','-q','--tab5','--tab6','--qseq','-f','-r','-F','-c','-s','-u','-5','-3',
                '--trim-to','--phred33','--phred64','--int-quals','--very-fast','--fast',
                '--sensitive','--very-sensitive','--very-fast-local','--fast-local',
                '--sensitive-local','--very-sensitive-local','-N','-L','-i','--n-ceil',
                '--dpad','--gbar','--ignore-quals','--nofw','--norc','--no-1mm-upfront',
                '--end-to-end','--local','--ma','--mp','--np','--rdg','--rfg','--score-min',
                '-k','-a','-D','-R','-I','-X','--fr','--rf','--ff','--no-mixed','--no-discordant',
                '--dovetail','--no-contain','--no-overlap','--align-paired-reads','--preserve-tags',
                '-t','--un','--al','--un-conc','--al-conc','--un-gz','--quiet','--met-file',
                '--met-stderr','--met','--no-unal','--no-head','--no-sq','--rg-id','--rg',
                '--omit-sec-seq','--sam-no-qname-trunc','--xeq','--soft-clipped-unmapped-tlen',
                '-p','--threads','--reorder','--mm','--qc-filter','--seed','--non-deterministic',
                ]

_args_STRINGTIE=['-G','--version','--conservative','--rf','--fr','-o','-l',
                '-f','-L','-m','-a','-j','-t','-c','-s','-v','-g','-M',
                '-p','-A','-B','-b','-e','-x','-u','-h','--merge','-F','-T','-i']

_args_CUFFLINKS=['-o','--output-dir','-p','--num-threads','--seed','-G','--GTF','-g','--GTF-guide','-M','--mask-file','-b','--frag-bias-correct','-u','--multi-read-correct','--library-type','--library-norm-method',
                '-m','--frag-len-mean','-s','--frag-len-std-dev','--max-mle-iterations','--compatible-hits-norm','--total-hits-norm','--num-frag-count-draws','--num-frag-assign-draws','--max-frag-multihits','--no-effective-length-correction',
                '--no-length-correction','-N','--upper-quartile-norm','--raw-mapped-norm','-L','--label','-F','--min-isoform-fraction','-j','--pre-mrna-fraction','-I','--max-intron-length','-a','--junc-alpha','-A','--small-anchor-fraction',
                '--min-frags-per-transfrag','--overhang-tolerance','--max-bundle-length','--max-bundle-frags','--min-intron-length','--trim-3-avgcov-thresh','--trim-3-dropoff-frac','--max-multiread-fraction','--overlap-radius',
                '--no-faux-reads','--3-overhang-tolerance','--intron-overhang-tolerance','-v','--verbose','-q','--quiet','--no-update-check']


_args_TRINITY=['--seqType','--max_memory','--left','--right','--single','--SS_lib_type','--CPU','--min_contig_length',
              '--long_reads','--genome_guided_bam','--jaccard_clip','--trimmomatic','--normalize_reads','--no_distributed_trinity_exec',
              '--output','--full_cleanup','--cite','--verbose','--version','--show_full_usage_info','--KMER_SIZE','--prep','--no_cleanup',
              '--no_version_check','--min_kmer_cov','--inchworm_cpu','--no_run_inchworm','--max_reads_per_graph','--min_glue','--no_bowtie',
              '--no_run_chrysalis','--bfly_opts','--PasaFly','--CuffFly','--group_pairs_distance','--path_reinforcement_distance','--no_path_merging',
              '--min_per_id_same_path','--max_diffs_same_path','--max_internal_gap_same_path','--bflyHeapSpaceMax','--bflyHeapSpaceInit',
              '--bflyGCThreads','--bflyCPU','--bflyCalculateCPU','--bfly_jar','--quality_trimming_params','--normalize_max_read_cov',
              '--normalize_by_read_set','--genome_guided_max_intron','--genome_guided_min_coverage','--genome_guided_min_reads_per_partition',
              '--grid_conf','--grid_node_CPU','--grid_node_max_memory']

_args_KALLISTO=['-i','--index','-o','--output-dir','--bias','-b','--bootstrap-samples',
                '--seed','--plaintext','--fusion','--single','--fr-stranded','--rf-stranded',
                '-l','--fragment-length','-s','--sd','-t','--threads','--pseudobam']

_args_SALMON=['--help-reads','-i','--index','-l','--libType','-r','--unmatedReads',
              '-1','--mates1','-2','--mates2','-o','--output','--discardOrphansQuasi',
              '--allowOrphansFMD','--seqBias','--gcBias','-p','--threads','--incompatPrior',
              '-g','--geneMap','-z','--writeMappings','--meta','--alternativeInitMode',
              '--auxDir','-c','--consistentHits','--dumpEq','-d','--dumpEqWeights',
              '--fasterMapping','--minAssignedFrags','--reduceGCMemory','--biasSpeedSamp',
              '--strictIntersect','--fldMax','--fldMean','--fldSD','-f','--forgettingFactor',
              '-m','--maxOcc','--initUniform','-w','--maxReadOcc','--noLengthCorrection',
              '--noEffectiveLengthCorrection','--noFragLengthDist','--noBiasLengthThreshold',
              '--numBiasSamples','--numAuxModelSamples','--numPreAuxModelSamples','--useVBOpt',
              '--rangeFactorizationBins','--numGibbsSamples','--numBootstraps','--thinningFactor',
              '-q','--perTranscriptPrior','--vbPrior','--writeOrphanLinks','--writeUnmappedNames',
              '-x','--quasiCoverage']


_args_SAMTOOLS=None