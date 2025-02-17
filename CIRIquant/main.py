#! /usr/bin/env python
# -*- encoding:utf-8 -*=
import os
import sys
import re
import argparse


def main():
    from version import __version__
    import circ
    import pipeline
    import pandas as pd
    from logger import get_logger
    from utils import check_file, check_dir, check_config, get_thread_num
    from utils import CIRCparser, TOOLS
    from fastq import fastq_cat

    # Init argparser
    parser = argparse.ArgumentParser(prog='CIRIquant')

    # required arguments
    parser.add_argument('-m', '--metadata', dest='metadata', metavar='FILE', default=None,
                        help='Input metadata', )    
    parser.add_argument('--config', dest='config_file', metavar='FILE',
                        help='Config file in YAML format', )
    # parser.add_argument('-1', '--read1', dest='mate1', metavar='MATE1',
    #                     help='Input mate1 reads (for paired-end data)', )
    # parser.add_argument('-2', '--read2', dest='mate2', metavar='MATE2',
    #                     help='Input mate2 reads (for paired-end data)', )

    # optional arguments
    parser.add_argument('-o', '--out', dest='output', metavar='DIR', default=None,
                        help='Output directory, default: ./', )
    parser.add_argument('-t', '--threads', dest='cpu_threads', default=4, metavar='INT',
                        help='Number of CPU threads, default: 4', )
    parser.add_argument('-a', '--anchor', dest='anchor', default=5, metavar='INT',
                        help='Minimum anchor length for junction alignment, default: 5', )
    parser.add_argument('-l', '--library-type', dest='library_type', metavar='INT', default=0,
                        help='Library type, 0: unstranded, 1: read1 match the sense strand,'
                             '2: read1 match the antisense strand, default: 0', )
    parser.add_argument('--ciri3', dest='ciri3', default=False, action='store_true',
                        help='Use CIRI3 for quantification', )

    parser.add_argument('-v', '--verbose', dest='verbosity', default=False, action='store_true',
                        help='Run in debugging mode', )
    parser.add_argument('--version', action='version',
                        version='%(prog)s {version}'.format(version=__version__))
    parser.add_argument('-e', '--log', dest='log_file', default=None, metavar='LOG',
                        help='Log file, default: out_dir/prefix.log', )

    # provide pre-defined list of circRNAs
    parser.add_argument('--bed', dest='bed', metavar='FILE', default=None,
                        help='bed file for putative circRNAs (optional)', )
    parser.add_argument('--circ', dest='circ', metavar='FILE', default=None,
                        help='circRNA prediction results from other softwares', )
    parser.add_argument('--tool', dest='tool', metavar='TOOL', default=None,
                        help='circRNA prediction tool, required if --circ is provided', )

    # # when provide RNase R result, do RNase R correction
    # parser.add_argument('--RNaseR', dest='rnaser', metavar='FILE', default=None,
    #                     help='CIRIquant result of RNase R sample', )

    # skip hisat2 alignment for RNA-seq data
    parser.add_argument('--bam', dest='bam', metavar='BAM', default=None,
                        help='hisat2 alignment to reference genome', )

    # skip stringtie prediction
    parser.add_argument('--no-gene', dest='gene_exp', default=False, action='store_true',
                        help='Skip stringtie estimation for gene abundance', )

    # skip FSJ calculation
    parser.add_argument('--no-fsj', dest='no_fsj', default=False, action='store_true',
                        help='Skip FSJ extraction to reduce run time', )
    parser.add_argument('--bsj-file', dest='bsj_read_file', metavar='FILE', default=None,
                        help='output BSJ read IDs to file (optional)')
    

    args = parser.parse_args()

    """Check required parameters"""
    # check input reads
    if args.metadata:
        metadata_file = check_file(args.metadata)
    else:
        sys.exit('No input metadata specified, please see manual for detailed information')

    try:
        lib_type = int(args.library_type)
    except ValueError:
        sys.exit('Wrong library type, please check your command.\nSupported types:\n0 - unstranded;\n'
                 '1 - read1 match the sense strand;\n2 - read1 match the antisense strand;')

    if lib_type not in [0, 1, 2]:
        sys.exit('Wrong library type, please check your command.\nSupported types:\n0 - unstranded;\n'
                 '1 - read1 match the sense strand;\n2 - read1 match the antisense strand;')

    # check configuration
    if args.config_file:
        config = check_config(check_file(args.config_file))
    else:
        sys.exit('A config file is needed, please see manual for detailed information.')

    """Check optional parameters"""
    # use circRNA bed file if provided
    bed_file = check_file(args.bed) if args.bed else None
    # metadata_file = check_file(args.metadata) if args.metadata else None
    circ_file = check_file(args.circ) if args.circ else None
    circ_tool = args.tool

    # # user provided RNase R CIRIquant results
    # rnaser_file = check_file(args.rnaser) if args.rnaser else None

    # pre aligned hisat2 bam
    hisat_bam = check_file(args.bam) if args.bam else None

    # # Output prefix
    # if args.prefix is None:
    #     try:
    #         prefix = re.search(r'(\S+)[_/-][12]', os.path.basename(reads[0])).group(1)
    #     except AttributeError:
    #         sys.exit('Ambiguous sample name, please manually select output prefix')
    # else:
    #     prefix = args.prefix

    # check output dir
    outdir = './Superdepth' if args.output is None else args.output
    outdir = check_dir(outdir)
    
    ################################ Change this later
    superdepth_outdir = outdir + '/circ'
    superdepth_outdir = check_dir(superdepth_outdir)


    # Parse arguments
    log_file = os.path.abspath(args.log_file) if args.log_file else '{}/CIRIquant.log'.format(outdir) ##############
    verbosity = args.verbosity
    logger = get_logger('CIRIquant', log_file, verbosity)

    # Add lib to PATH
    lib_path = os.path.dirname(os.path.split(os.path.realpath(__file__))[0]) + '/libs'
    os.environ['PATH'] = lib_path + ':' + os.environ['PATH']
    ciri2_exec = lib_path + '/CIRI2.pl'
    ciri3_exec = lib_path + '/CIRI3.jar'

    """Start Running"""
    os.chdir(outdir)
    logger.info("Input metadata: {}".format(os.path.basename(args.metadata)))

    if lib_type == 0:
        lib_name = 'unstranded'
    elif lib_type == 1:
        lib_name = 'ScriptSeq'
    elif lib_type == 2:
        lib_name = 'TAKARA SMARTer'
    else:
        sys.exit('Unsupported library type, please check the manual for instructions.')

    logger.info('Library type: {}'.format(lib_name))
    logger.info('Output directory: {}'.format(outdir))
    logger.info('Config: {} Loaded'.format(config))

    thread = get_thread_num(int(args.cpu_threads))
    anchor = int(args.anchor)


    # Step1: CircRNA prep

    ## Preliminary check of metadata
    Samples = pd.read_csv(metadata_file)

    # Define required columns
    required_columns = {'sample', 'fastq_1', 'fastq_2'}

    # Check if required columns are present
    if not required_columns.issubset(Samples.columns):
        missing_columns = required_columns - set(Samples.columns)
        logger.error("Error: The following required columns are missing: {}".format(missing_columns))
        sys.exit(1)

    # Check if the 'bam' column exists
    if 'bam' in Samples.columns:
        # Store the 'bam' column in a specific variable
        bam_list = Samples['bam'].dropna().tolist()  # Drop NaN values and convert to list
        logger.info("BAM column found. Storing BAM samples.")
    else:
        bam_list = [None] * len(Samples)  # Initialize an empty list if 'bam' column is not present
        logger.info("No BAM column found.")

    # Store the other required columns in respective variables
    prefix_list = Samples['sample'].tolist()
    fastq_1_list = Samples['fastq_1'].tolist()
    fastq_2_list = Samples['fastq_2'].tolist()  

    # set SD prefix
    prefix_SD = 'Superdepth'

    ## 
    if bed_file:
        logger.info('Using user-provided circRNA bed file: {}'.format(bed_file))
        logger.info('Reminder: Ensure circRNAs are derived from a Superdepth-based run from any circRNA detection tool')
    else:
        if circ_file or circ_tool:
            if circ_file and circ_tool:
                logger.info('Using predicted circRNA results from {}: {}'.format(circ_tool, circ_file))
                logger.info('Reminder: Ensure circRNAs are derived from a Superdepth-based run run from any circRNA detection tool')
                circ_parser = CIRCparser(circ_file, circ_tool)
            else:
                sys.exit('--circ and --tool must be provided in the same time!')
        elif metadata_file:
            ## Prepare SD reads alignment and circRNA detection
            logger.info('No circRNA information provided, run Superdepth-based CIRI2 for junction site prediction ..')

            ## Cocatenate fastq files
            logger.info('Concatenating reads 1 ..')
            fastq_cat(fastq_1_list, thread, superdepth_outdir + '/Superdepth_1.fastq.gz')

            logger.info('Concatenating reads 2 ..')
            fastq_cat(fastq_2_list, thread, superdepth_outdir + '/Superdepth_2.fastq.gz')

            SD_reads = [superdepth_outdir + '/Superdepth_1.fastq.gz', superdepth_outdir + '/Superdepth_2.fastq.gz']
            
            # run bwa alignment
            bwa_sam = pipeline.run_bwa(log_file, thread, SD_reads, outdir, prefix_SD)
            if args.ciri3:
                ciri_file = pipeline.run_ciri3(log_file, thread, bwa_sam, outdir, prefix_SD, ciri3_exec)
            else:
                ciri_file = pipeline.run_ciri2(log_file, thread, bwa_sam, outdir, prefix_SD, ciri2_exec)
            circ_parser = CIRCparser(ciri_file, 'CIRI2')

        bed_file = '{}/Superdepth.bed'.format(outdir)
        circ_parser.convert(bed_file)


    # Step 1.2 Generate fasta and hisat2 index using the predicted circRNAs
    circ_info, denovo_index = circ.circ_index(log_file, thread, bed_file, outdir)


    # Step 2 & 3: HISAT2 mapping and Gene Abundance Estimation
    # if hisat_bam is None:

    for prefix, fastq_1, fastq_2, hisat_bam in zip(prefix_list, fastq_1_list, fastq_2_list, bam_list):
        reads = [fastq_1, fastq_2]

        outdir_prefix = outdir + '/' + prefix
        outdir_prefix = check_dir(outdir_prefix)

        if hisat_bam:
            logger.info('Skipping alignment for {}, BAM file exists: {}'.format(prefix, hisat_bam))
        else:
            logger.info('Aligning RNA-seq reads to reference genome of {}'.format(prefix))
            hisat_bam = pipeline.align_genome(log_file, thread, reads, outdir_prefix, prefix)
        
        if args.gene_exp:
            logger.info('Skipping gene abundance estimation for {}'.format(prefix))
        else:
            logger.info('Estimate gene abundance for {}'.format(prefix))    
            pipeline.gene_abundance(log_file, thread, outdir_prefix, prefix, hisat_bam)


        # Step4: estimate circRNA expression level
        if args.no_fsj:
            logger.info('Skipping FSJ reads extraction')
        out_file = circ.proc(log_file, thread, circ_info, denovo_index, hisat_bam, reads, outdir_prefix, prefix, anchor, lib_type,
                            args.no_fsj, args.bsj_read_file)
        
        # Remove temporary file
        pipeline.clean_tmp(outdir_prefix, prefix)

        logger.info('circRNA Expression profile: {}'.format(os.path.basename(out_file)))
        



    # else:
    #     logger.info('HISAT2 alignment bam provided, skipping alignment step ..')
    # logger.debug('HISAT2 bam: {}'.format(os.path.basename(hisat_bam)))

    # # Step3: Estimate Gene Abundance
    # if args.gene_exp:
    #     logger.info('Skipping gene abundance estimation')
    # else:
    #     for sample in range(len(Samples)):
    #         sample = Samples.iloc[sample]
    #         fastq_1 = sample['fastq_1']
    #         fastq_2 = sample['fastq_2']
    #         prefix = sample['sample']
    #         reads = [fastq_1, fastq_2]
    #         pipeline.gene_abundance(log_file, thread, outdir, prefix, hisat_bam)


    # # Step4: estimate circRNA expression level
    # if args.no_fsj:
    #     logger.info('Skipping FSJ reads extraction')
    
    # for hisat_bam, reads, outdir_sample, prefix in zip(hisat_bam_list, read, outdir_prefix_list, prefix_list):
    #     out_file = circ.proc(log_file, thread, circ_info, denovo_index, hisat_bam, reads, outdir_sample, prefix, anchor, lib_type,
    #                         args.no_fsj, args.bsj_read_file)

    # # Remove temporary files
    # pipeline.clean_tmp(outdir, prefix)

    # logger.info('circRNA Expression profile: {}'.format(os.path.basename(out_file)))

    logger.info('Finished!')


if __name__ == '__main__':
    main()
