import os

import pypeliner
import pypeliner.managed as mgd

from wgs.utils import helpers

def calc_bam_metrics(config, bamfile, 
                    picard_insert_metrics, picard_insert_pdf, flagstat_metrics,
                    picard_gc, picard_gc_summary, picard_gc_pdf, 
                    picard_wgs, picard_wgs_params,
                    tempdir
):
    workflow = pypeliner.workflow.Workflow()
    ref = config["reference"]

    workflow.transform(
        name = "calc_picard_insert_metrics",
        func = bam_collect_insert_metrics,
        args = (
            bamfile,
            flagstat_metrics,
            picard_insert_metrics,
            picard_insert_pdf,
            tempdir,
            docker_image = config["docker"]["picard_insert"]
        )
    )
     workflow.transform(
        name = "calc_picard_gc_metrics",
        func = bam_collect_gc_metrics,
        args = (
            bamfile,
            ref,
            picard_gc,
            picard_gc_summary,
            picard_gc_pdf,
            tempdir,
            docker_image = config["docker"]["picard_gc"]
        )
    )
     workflow.transform(
        name = "calc_picard_wgs_metrics",
        func = bam_collect_wgs_metrics,
        args = (
            bamfile,
            ref,
            picard_wgs,
            picard_wgs_params,
            tempdir,
            docker_image = config["docker"]["picard_wgs"]
        )
    )
 



    

def index_and_flagstat(bamfile, indexfile, flagstatfile, docker_image=None):
    cmd = ['samtools', 'index', bamfile, indexfile]
    pypeliner.commandline.execute(*cmd, docker_image=docker_image)

    cmd = ['samtools', 'flagstat', bamfile, '>', flagstatfile]
    pypeliner.commandline.execute(*cmd, docker_image=docker_image)


def flagstat(bamfile, flagstatfile, docker_image=None):
    cmd = ['samtools', 'flagstat', bamfile, '>', flagstatfile]
    pypeliner.commandline.execute(*cmd, docker_image=docker_image)


def index(bamfile, indexfile, docker_image=None):
    cmd = ['samtools', 'index', bamfile, indexfile]
    pypeliner.commandline.execute(*cmd, docker_image=docker_image)


def markdups(input, output, metrics, tempdir, mem="2G", docker_image=None):
    cmd = ['picard', '-Xmx' + mem, '-Xms' + mem,
           '-XX:ParallelGCThreads=1',
           'MarkDuplicates',
           'INPUT=' + input,
           'OUTPUT=' + output,
           'METRICS_FILE=' + metrics,
           'REMOVE_DUPLICATES=False',
           'ASSUME_SORTED=True',
           'VALIDATION_STRINGENCY=LENIENT',
           'TMP_DIR=' + tempdir,
           'MAX_RECORDS_IN_RAM=150000'
           ]

    pypeliner.commandline.execute(*cmd, docker_image=docker_image)


def picard_merge_bams(inputs, output, mem="2G", **kwargs):
    if isinstance(inputs, dict):
        inputs = inputs.values()

    cmd = ['picard', '-Xmx' + mem, '-Xms' + mem,
           '-XX:ParallelGCThreads=1',
           'MergeSamFiles',
           'OUTPUT=' + output,
           'SORT_ORDER=coordinate',
           'ASSUME_SORTED=true',
           'VALIDATION_STRINGENCY=LENIENT',
           'MAX_RECORDS_IN_RAM=150000'
           ]

    for bamfile in inputs:
        cmd.append('I=' + os.path.abspath(bamfile))

    pypeliner.commandline.execute(*cmd, **kwargs)


def bam_index(infile, outfile, **kwargs):
    pypeliner.commandline.execute(
        'samtools', 'index',
        infile,
        outfile,
        **kwargs)


def merge_bams(inputs, output, picard_docker_image=None, samtools_docker_image=None):
    output_index = output + '.bai'
    picard_merge_bams(inputs, output, docker_image=picard_docker_image)
    bam_index(output, output_index, docker_image=samtools_docker_image)


def bwa_mem_paired_end(fastq1, fastq2, output,
                         reference, readgroup,
                       numthreads,
                         **kwargs):
    """
    run bwa aln on both fastq files,
    bwa sampe to align, and convert to bam with samtools view
    """

    if not numthreads:
        numthreads=1

    if not readgroup:
        pypeliner.commandline.execute(
            'bwa', 'mem', '-M',
            '-t', numthreads,
            reference, fastq1, fastq2,
            '>', output,
            **kwargs)
    else:
        try:
            readgroup_literal = '"' + readgroup + '"'
            pypeliner.commandline.execute(
                'bwa', 'mem', '-M', '-R', readgroup_literal,
                '-t', numthreads,
                reference, fastq1, fastq2,
                '>', output,
                **kwargs)
        except pypeliner.commandline.CommandLineException:
            pypeliner.commandline.execute(
                'bwa', 'mem', '-M', '-R', readgroup,
                '-t', numthreads,
                reference, fastq1, fastq2,
                '>', output,
                **kwargs)


def get_readgroup(read_group_info, sample_id, lane_id):
    if read_group_info:
        rg_id = read_group_info['ID'].format(sample_id=sample_id, lane_id=lane_id)
        read_group = ['@RG', 'ID:{0}'.format(rg_id)]
        for key, value in sorted(read_group_info.items()):
            if key == 'ID':
                continue
            value = value.format(sample_id=sample_id, lane_id=lane_id)
            read_group.append(':'.join((key, value)))
        read_group = '\\t'.join(read_group)
    elif sample_id or lane_id:
        sample_id = sample_id if sample_id else ''
        lane_id = lane_id if lane_id else ''
        ids = '-'.join([sample_id, lane_id])
        read_group = ['@RG', 'ID:{0}'.format(ids)]
        read_group = '\\t'.join(read_group)
    else:
        read_group = None

    return read_group


def samtools_sam_to_bam(samfile, bamfile,
                         **kwargs):
    pypeliner.commandline.execute(
        'samtools', 'view', '-bSh', samfile,
        '>', bamfile,
        **kwargs)


def align_bwa_mem(
        read_1, read_2, ref_genome, aligned_bam, threads, tempdir,
        sample_id=None, lane_id=None, read_group_info=None,
        docker_config=None
):

    readgroup = get_readgroup(read_group_info, sample_id, lane_id)

    helpers.makedirs(tempdir)

    bwa_mem_output = os.path.join(tempdir, "bwa_mem.sam")
    bwa_mem_paired_end(
        read_1, read_2, bwa_mem_output, ref_genome,
        readgroup, threads, docker_image=docker_config['bwa']
    )

    samtools_sam_to_bam(
        bwa_mem_output, aligned_bam, docker_image=docker_config['samtools']
    )


def bam_sort(bam_filename, sorted_bam_filename, threads=1, mem="2G", docker_image=None):

    pypeliner.commandline.execute(
        'samtools', 'sort', '-@', threads, '-m', mem,
        bam_filename,
        '-o',
        sorted_bam_filename,
        docker_image=docker_image)

#taken from single cell
def bam_collect_wgs_metrics(bam_filename, ref_genome, metrics_filename,
                            config, tempdir, mem="2G", docker_image=None):
    if not os.path.exists(tempdir):
        makedirs(tempdir)

    pypeliner.commandline.execute(
        'picard', '-Xmx' + mem, '-Xms' + mem,
        '-XX:ParallelGCThreads=1',
        'CollectWgsMetrics',
                  'INPUT=' + bam_filename,
                  'OUTPUT=' + metrics_filename,
                  'REFERENCE_SEQUENCE=' + ref_genome,
                  'MINIMUM_BASE_QUALITY=' +
                  str(config['min_bqual']),
                  'MINIMUM_MAPPING_QUALITY=' +
                  str(config['min_mqual']),
        'COVERAGE_CAP=500',
        'VALIDATION_STRINGENCY=LENIENT',
                  'COUNT_UNPAIRED=' +
                  ('True' if config['count_unpaired'] else 'False'),
                  'TMP_DIR=' + tempdir,
        'MAX_RECORDS_IN_RAM=150000',
        docker_image=docker_image
    )


def bam_collect_gc_metrics(bam_filename, ref_genome, metrics_filename,
                           summary_filename, chart_filename, tempdir,
                           mem="2G", docker_image=None):
    if not os.path.exists(tempdir):
        makedirs(tempdir)

    pypeliner.commandline.execute(
        'picard', '-Xmx' + mem, '-Xms' + mem,
        '-XX:ParallelGCThreads=1',
        'CollectGcBiasMetrics',
                  'INPUT=' + bam_filename,
                  'OUTPUT=' + metrics_filename,
                  'REFERENCE_SEQUENCE=' + ref_genome,
                  'S=' + summary_filename,
                  'CHART_OUTPUT=' + chart_filename,
        'VALIDATION_STRINGENCY=LENIENT',
                  'TMP_DIR=' + tempdir,
        'MAX_RECORDS_IN_RAM=150000',
        docker_image=docker_image
    )


def bam_collect_insert_metrics(bam_filename, flagstat_metrics_filename,
                               metrics_filename, histogram_filename, tempdir,
                               mem="2G", docker_image=None):
    # Check if any paired reads exist
    has_paired = None
    with open(flagstat_metrics_filename) as f:
        for line in f:
            if 'properly paired' in line:
                if line.startswith('0 '):
                    has_paired = False
                else:
                    has_paired = True

    if has_paired is None:
        raise Exception('Unable to determine number of properly paired reads from {}'.format(
            flagstat_metrics_filename))

    if not has_paired:
        with open(metrics_filename, 'w') as f:
            f.write('## FAILED: No properly paired reads\n')
        with open(histogram_filename, 'w'):
            pass
        return

    if not os.path.exists(tempdir):
        makedirs(tempdir)

    pypeliner.commandline.execute(
        'picard', '-Xmx' + mem, '-Xms' + mem,
        '-XX:ParallelGCThreads=1',
        'CollectInsertSizeMetrics',
                  'INPUT=' + bam_filename,
                  'OUTPUT=' + metrics_filename,
                  'HISTOGRAM_FILE=' + histogram_filename,
        'ASSUME_SORTED=True',
        'VALIDATION_STRINGENCY=LENIENT',
                  'TMP_DIR=' + tempdir,
        'MAX_RECORDS_IN_RAM=150000',
        docker_image=docker_image
    )