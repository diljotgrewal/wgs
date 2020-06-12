'''
Created on Feb 21, 2018

@author: pwalters
'''
import os

import pypeliner
import pypeliner.managed as mgd
from wgs.utils import helpers
from wgs.config import config


def create_consensus_workflow(
        museq_germline,
        museq_snv,
        strelka_snv,
        strelka_indel,
        somatic_calls,
        somatic_snpeff,
        somatic_ma,
        somatic_ids,
        indel_calls,
        indel_snpeff,
        indel_ma,
        indel_ids,
        germline_calls,
        germline_snpeff,
        germline_ma,
        germline_ids,
        refdir
):
    params = config.default_params('variant_calling')
    chromosomes = config.refdir_data(refdir)['params']['chromosomes']

    # germline_snpeff_annotations = os.path.join(
    #     outdir, '{}_germline_snpeff_annotations.csv.gz'.format(sample_id)
    # )
    # indel_snpeff_annotations = os.path.join(
    #     outdir, '{}_indel_snpeff_annotations.csv.gz'.format(sample_id)
    # )
    # somatic_snpeff_annotations = os.path.join(
    #     outdir, '{}_somatic_snpeff_annotations.csv.gz'.format(sample_id)
    # )
    #
    # germline_ma_annotations = os.path.join(
    #     outdir, '{}_germline_ma_annotations.csv.gz'.format(sample_id)
    # )
    # indel_ma_annotations = os.path.join(
    #     outdir, '{}_indel_ma_annotations.csv.gz'.format(sample_id)
    # )
    # somatic_ma_annotations = os.path.join(
    #     outdir, '{}_somatic_ma_annotations.csv.gz'.format(sample_id)
    # )
    #
    # germline_ids_annotations = os.path.join(
    #     outdir, '{}_germline_ids_annotations.csv.gz'.format(sample_id)
    # )
    # indel_ids_annotations = os.path.join(
    #     outdir, '{}_indel_ids_annotations.csv.gz'.format(sample_id)
    # )
    # somatic_ids_annotations = os.path.join(
    #     outdir, '{}_somatic_ids_annotations.csv.gz'.format(sample_id)
    # )

    workflow = pypeliner.workflow.Workflow()

    workflow.transform(
        name='parse_museq_germlines',
        ctx=helpers.get_default_ctx(
            memory=15,
            walltime='8:00', ),
        func='wgs.workflows.variant_calling_consensus.tasks.parse_vcf',
        args=(
            mgd.InputFile(museq_germline, extensions=['.csi', '.tbi']),
            mgd.OutputFile(germline_calls, extensions=['.yaml']),
            mgd.OutputFile(germline_snpeff, extensions=['.yaml']),
            mgd.OutputFile(germline_ma, extensions=['.yaml']),
            mgd.OutputFile(germline_ids, extensions=['.yaml']),
            params["parse_museq"],
            chromosomes,
            mgd.TempSpace("tempdir_parse_germlines")
        ),
    )

    workflow.transform(
        name='parse_strelka_indel',
        ctx=helpers.get_default_ctx(
            memory=15,
            walltime='8:00', ),
        func='wgs.workflows.variant_calling_consensus.tasks.parse_vcf',
        args=(
            mgd.InputFile(strelka_indel, extensions=['.csi', '.tbi']),
            mgd.OutputFile(indel_calls, extensions=['.yaml']),
            mgd.OutputFile(indel_snpeff, extensions=['.yaml']),
            mgd.OutputFile(indel_ma, extensions=['.yaml']),
            mgd.OutputFile(indel_ids, extensions=['.yaml']),
            params["parse_strelka"],
            chromosomes,
            mgd.TempSpace("tempdir_strelka_indel")
        ),
    )

    workflow.transform(
        name='parse_museq_snv',
        ctx=helpers.get_default_ctx(
            memory=15,
            walltime='8:00', ),
        func='wgs.workflows.variant_calling_consensus.tasks.parse_vcf',
        args=(
            mgd.InputFile(museq_snv, extensions=['.csi', '.tbi']),
            mgd.TempOutputFile('museq_snv.csv', extensions=['.yaml']),
            mgd.TempOutputFile('museq_snpeff.csv', extensions=['.yaml']),
            mgd.TempOutputFile('museq_ma.csv', extensions=['.yaml']),
            mgd.TempOutputFile('museq_ids.csv', extensions=['.yaml']),
            params["parse_museq"],
            chromosomes,
            mgd.TempSpace("tempdir_parse_museq_snv")
        ),
    )

    workflow.transform(
        name='parse_strelka_snv',
        ctx=helpers.get_default_ctx(
            memory=15,
            walltime='8:00', ),
        func='wgs.workflows.variant_calling_consensus.tasks.parse_vcf',
        args=(
            mgd.InputFile(strelka_snv, extensions=['.csi', '.tbi']),
            mgd.TempOutputFile('strelka_snv.csv', extensions=['.yaml']),
            mgd.TempOutputFile('strelka_snv_snpeff.csv', extensions=['.yaml']),
            mgd.TempOutputFile('strelka_snv_ma.csv', extensions=['.yaml']),
            mgd.TempOutputFile('strelka_snv_ids.csv', extensions=['.yaml']),
            params["parse_strelka"],
            chromosomes,
            mgd.TempSpace("tempdir_parse_strelka_snv")
        ),
    )

    workflow.transform(
        name='merge_snvs',
        ctx=helpers.get_default_ctx(
            memory=15,
            walltime='8:00', ),
        func='wgs.workflows.variant_calling_consensus.tasks.merge_overlap',
        args=(
            [mgd.TempInputFile('strelka_snv.csv', extensions=['.yaml']),
             mgd.TempInputFile('museq_snv.csv', extensions=['.yaml'])],
            mgd.OutputFile(somatic_calls, extensions=['.yaml']),
        ),
    )

    workflow.transform(
        name='merge_snpeff',
        ctx=helpers.get_default_ctx(
            memory=15,
            walltime='8:00', ),
        func='wgs.workflows.variant_calling_consensus.tasks.merge_overlap',
        args=(
            [mgd.TempInputFile('strelka_snv_snpeff.csv', extensions=['.yaml']),
             mgd.TempInputFile('museq_snpeff.csv', extensions=['.yaml'])],
            mgd.OutputFile(somatic_snpeff, extensions=['.yaml']),
        ),
        kwargs={'on': ['chrom', 'pos']}
    )

    workflow.transform(
        name='merge_ma',
        ctx=helpers.get_default_ctx(
            memory=15,
            walltime='8:00', ),
        func='wgs.workflows.variant_calling_consensus.tasks.merge_overlap',
        args=(
            [mgd.TempInputFile('strelka_snv_ma.csv', extensions=['.yaml']),
             mgd.TempInputFile('museq_ma.csv', extensions=['.yaml'])],
            mgd.OutputFile(somatic_ma, extensions=['.yaml']),
        ),
        kwargs={'on': ['chrom', 'pos']}
    )

    workflow.transform(
        name='merge_ids',
        ctx=helpers.get_default_ctx(
            memory=15,
            walltime='8:00', ),
        func='wgs.workflows.variant_calling_consensus.tasks.merge_overlap',
        args=(
            [mgd.TempInputFile('strelka_snv_ids.csv', extensions=['.yaml']),
             mgd.TempInputFile('museq_ids.csv', extensions=['.yaml'])],
            mgd.OutputFile(somatic_ids, extensions=['.yaml']),
        ),
        kwargs={'on': ['chrom', 'pos', 'value', 'type']}
    )

    return workflow
