import matplotlib.pyplot as plt
import pandas as pd
import pypeliner
from wgs.utils import helpers



def merge_cna_tables(amps, dels, labels, output, cohort):
    cna = {d1: [a,d] for (d1, a), (d2, d) in zip(amps.items(), dels.items()) }
    number=0
    for s, files in cna.items():
        
        label = labels[(cohort, s)]
        for f in files:
        
            data = pd.read_csv(f, usecols=["cn_type", "gene_name", "pass_filter"])

            data=data[data.pass_filter == True]

            data=data.rename(columns={"cn_type": "CN", "gene_name":"Gene"})

            data["Sample_name"] = [label] * len(data)

            data = data[["Gene", "Sample_name", "CN"]]

            if number==0:
                header=True
            else:
                header=False

            if not data.empty:
                number+=1
                data.to_csv(output, index=False, mode='a', header=header, sep="\t")


def classify_remixt(sample_label, remixt, gtf, output_dir, amps, dels, docker_image=None):

    cmd = [
        "classifycopynumber", gtf, output_dir, sample_label, amps, dels, "--remixt_h5_filename", remixt, "--plot", False
    ]

    pypeliner.commandline.execute(*cmd, docker_image=docker_image)



def merge_mafs(mafs, merged_maf, write_header=True):
    for m in mafs:

        maf = pd.read_csv(m, sep="\t", chunksize=10e6)
        for chunk in maf:

            chunk.to_csv(merged_maf, sep="\t", index=False, header=write_header, mode='a')
        #only write the first header
        write_header=False   


def merge_relable_mafs(mafs, merged_maf, labels, class_label=None, write_header=True):

    mafs = list(mafs.values())
    labels = list(labels.values())

    #only write the first header
    write_header = write_header
    #TODO better v names
    for label, maf in zip(labels, mafs):
        maf = pd.read_csv(maf, sep="\t", chunksize=10e6, skiprows=1)
        for chunk in maf:
            if labels:
                chunk["Tumor_Sample_Barcode"] = [label] * len(chunk.index)
            if class_label:
                chunk["Variant_Classification"] = chunk.Variant_Classification.apply(lambda c: c+class_label)
            chunk.to_csv(merged_maf, sep="\t", index=False, header=write_header, mode='a')

        #only write the first header
        write_header=False
    

def annotate_maf_with_oncokb(
        maf, api_key, tmpspace, annotated_maf, docker_image=None
):
    helpers.makedirs(tmpspace)

    cmd = [
        "MafAnnotator.py", "-i", maf, "-o", annotated_maf, "-b", api_key
    ]
    cmd = [
        "/juno/work/shah/abramsd/oncokb-annotator/MafAnnotator.py", "-i", maf, "-o", annotated_maf, "-b", api_key
    ]

    pypeliner.commandline.execute(*cmd, docker_image=docker_image)


def filter_annotated_maf(annotated_maf, filtered_maf):
    annotated_maf = pd.read_csv(annotated_maf, sep="\t")

    # filt_maf = annotated_maf[
    #     (annotated_maf.ONCOGENIC == "Oncogenic")
    #     | (annotated_maf.ONCOGENIC == "Likely Oncogenic")
    #     ]
    filt_maf = annotated_maf[
        (annotated_maf.oncogenic == "Oncogenic")
        | (annotated_maf.oncogenic == "Likely Oncogenic")
        ]
    filt_maf.to_csv(filtered_maf, sep="\t")


def plot_mutation_burden(maf, burden_plot_path):
    maf = pd.read_csv(maf, sep="\t").drop_duplicates()
    data = maf.groupby("Tumor_Sample_Barcode").size().sort_values(ascending=False)

    fig, axis = plt.subplots(figsize=(15, 5))
    nums = range(len(data.index))
    axis.bar(nums, data, align='center', color="black")
    plt.xticks(nums, data.index, rotation='vertical')
    axis.set_ylabel("Number of mutations")
    fig.savefig(burden_plot_path, format="png")
    plt.close()


def make_R_cohort_plots(
        cohort_maf, cntable, oncoplot_path, somatic_interactions,
        mafsummary, filtered_maf=None
):
    if not filtered_maf:
        filtered_maf = cohort_maf

    plots_cmd = [
        "maftools_plots.R", cohort_maf, cntable, filtered_maf,
        oncoplot_path, somatic_interactions, mafsummary
    ]

    pypeliner.commandline.execute(*plots_cmd)


def make_report(cohort_label, oncoplot, somatic_interactions, mafsummary, burden_plot, report_path):
    cmd = [
        "run_report.sh", report_path, cohort_label, oncoplot,
        somatic_interactions, mafsummary, burden_plot
    ]
    pypeliner.commandline.execute(*cmd)


c("Frame_Shift_Del_somatic", "Frame_Shift_Ins_somatic", "Splice_Site_somatic", "Translation_Start_Site_somatic",
"Nonsense_Mutation_somatic", "Nonstop_Mutation_somatic", "In_Frame_Del_somatic","In_Frame_Ins_somatic", "Missense_Mutation_somatic",
"Frame_Shift_Del_germline", "Frame_Shift_Ins_germline", "Splice_Site_germline", "Translation_Start_Site_germline",
"Nonsense_Mutation_germline", "Nonstop_Mutation_germline", "In_Frame_Del_germline","In_Frame_Ins_germline", "Missense_Mutation_germline")