import random
import vcf
import gzip
import os
import consensus
import pypeliner
from shutil import copyfile
import csv
import pandas as pd

def get_test_line(v1, v2, n=10000):
    #choose a random chrom
    not_matched = False
    reader1 = vcf.Reader(filename=v1)
    reader2 = vcf.Reader(filename=v2)
    chrom = random.choice(list(map(str, list(range(1, 23)) + ["X"])))
    lines1 = reader1.fetch(chrom)
    line = random.choice([next(lines1) for _ in range(n)])
    while not_matched == False:
        try:
            next(reader2.fetch(line.CHROM, start=line.POS-1, end=line.POS))
            line = random.choice([next(lines1) for _ in range(n)])

        except StopIteration:
            not_matched=True
            return line


def record_to_str(record, f, tmpdir):

    dummyfile = os.path.join(tmpdir, "dummy")
    pypeliner.commandline.execute("scp", f, dummyfile)

    r = vcf.Reader(open(dummyfile, "rt"))
    w = vcf.Writer(open(dummyfile, "at"), r)  

    template_record = next(r)

    ffs = w._map(str, [record.CHROM, record.POS, template_record.ID, record.REF] \
            + [w._format_alt(record.ALT), template_record.QUAL or '.', w._format_filter(template_record.FILTER),
                w._format_info(template_record.INFO)])

    if template_record.FORMAT:
        ffs.append(template_record.FORMAT)

    samples = [w._format_sample(template_record.FORMAT, sample)
        for sample in template_record.samples]   

    w.close()

    #delete dummy file
    pypeliner.commandline.execute("rm", dummyfile)

    output = "\t".join(ffs+samples)
    return output


def add_entry_to_vcf(record, v, tmpdir):
    '''
    you CANNOT append to a file wiht pyvcf, as it automatically writes out a vcf header with the first call to
    write_record. the work around is the make a writer on a dummy tmp vcf, use that writer to  get the formatted
    str version of thevariant, than normal io append that string to an unzipped v, then delete the dummy and recipe and index v
    '''
    #unzip
    pypeliner.commandline.execute("gunzip", v)
    gunzipped = v.replace(".gz", "")

    recordstr = record_to_str(record, gunzipped, tmpdir)

    #write
    with open(gunzipped, "at") as f:
        f.write(recordstr + "\n")
    f.close()

    #rezip
    s = gunzipped+"_s"
    pypeliner.commandline.execute("vcf-sort", gunzipped ,">", s)
    pypeliner.commandline.execute("mv", s, gunzipped)
    pypeliner.commandline.execute("bgzip", gunzipped)

    #now reindex
    pypeliner.commandline.execute("tabix", "-f", "-p", "vcf", v)


def run_consensus_test(vcfs, test1, test2, tmpdir):

    testfile1 = os.path.join(tmpdir, os.path.basename(vcfs[test1]))
    testfile2 = os.path.join(tmpdir, os.path.basename(vcfs[test2]))

    copyfile(vcfs[test1], testfile1)
    copyfile(vcfs[test1] + ".tbi", testfile1 + ".tbi")
    copyfile(vcfs[test2], testfile2)

    copyfile(vcfs[test2] + ".tbi", testfile2 + ".tbi")

    vcfs[test1] = testfile1
    vcfs[test2] = testfile2

    #make test
    record = get_test_line(testfile1, testfile2)
    
    add_entry_to_vcf(record, testfile2, tmpdir) 

    reader1 = vcf.Reader(filename=testfile1)
    assert next(reader1.fetch(record.CHROM, start=record.POS-1, end=record.POS))
    reader2 = vcf.Reader(filename=testfile2)
    assert next(reader2.fetch(record.CHROM, start=record.POS-1, end=record.POS))

    #run consensus
    chroms = list(map(str, list(range(1, 23)) + ["X"]))
    consensus_vcf = os.path.join(tmpdir, "consensus.vcf")
    counts = os.path.join(tmpdir, "counts.tsv")

    consensus.main(vcfs["museq"], 
        vcfs["freebayes"], 
        vcfs["rtg"], 
        vcfs["samtools"], 
        consensus_vcf, 
        counts, 
        chroms    
    )

    #read in consensus vcf and make sure the variant is included
    consensus_reader = pd.read_csv(consensus_vcf, sep="\t")
    assert not consensus_reader[(consensus_reader.POS==record.POS) & consensus_reader["#CHROM"]==record.CHROM].empty


def run_normalize_test(vcfs, test1, test2, a1, r1, a2, r2, tmpdir):

    testfile1 = os.path.join(tmpdir, os.path.basename(vcfs[test1]))
    testfile2 = os.path.join(tmpdir, os.path.basename(vcfs[test2]))

    copyfile(vcfs[test1], testfile1)
    copyfile(vcfs[test1] + ".tbi", testfile1 + ".tbi")
    copyfile(vcfs[test2], testfile2)

    copyfile(vcfs[test2] + ".tbi", testfile2 + ".tbi")

    vcfs[test1] = testfile1
    vcfs[test2] = testfile2

    #make test
    record = get_test_line(testfile1, testfile2)

    record.REF = r1
    record.ALT = a1
    add_entry_to_vcf(record, testfile1, tmpdir) 

    record.REF = r2
    record.ALT = a2
    add_entry_to_vcf(record, testfile2, tmpdir) 

    reader1 = vcf.Reader(filename=testfile1)
    assert next(reader1.fetch(record.CHROM, start=record.POS-1, end=record.POS))
    reader2 = vcf.Reader(filename=testfile2)
    assert next(reader2.fetch(record.CHROM, start=record.POS-1, end=record.POS))

    #run consensus
    chroms = list(map(str, list(range(1, 23)) + ["X"]))
    consensus_vcf = os.path.join(tmpdir, "consensus.vcf")
    counts = os.path.join(tmpdir, "counts.tsv")

    consensus.main(vcfs["museq"], 
        vcfs["freebayes"], 
        vcfs["rtg"], 
        vcfs["samtools"], 
        consensus_vcf, 
        counts, 
        chroms    
    )

    #read in consensus vcf and make sure the variant is included
    consensus_reader = pd.read_csv(consensus_vcf, sep="\t")
    assert not consensus_reader[(consensus_reader.POS==record.POS) & consensus_reader["#CHROM"]==record.CHROM].empty
    


def test_museq_freebayes(freebayes, museq, rtg, samtools, testdir):
    vcfs = {"freebayes": freebayes, "museq": museq, "rtg": rtg, "samtools": samtools}
    run_consensus_test(vcfs, "museq", "freebayes", testdir)

def test_museq_rtg(freebayes, museq, rtg, samtools, testdir):
    vcfs = {"freebayes": freebayes, "museq": museq, "rtg": rtg, "samtools": samtools}
    run_consensus_test(vcfs, "rtg", "museq", testdir)

def test_museq_samtools(freebayes, museq, rtg, samtools, testdir):
    vcfs = {"freebayes": freebayes, "museq": museq, "rtg": rtg, "samtools": samtools}
    run_consensus_test(vcfs, "samtools", "museq", testdir)

def test_samtools_rtg(freebayes, museq, rtg, samtools, testdir):
    vcfs = {"freebayes": freebayes, "museq": museq, "rtg": rtg, "samtools": samtools}
    run_consensus_test(vcfs, "samtools", "rtg", testdir)

def test_samtools_freebayes(freebayes, museq, rtg, samtools, testdir):
    vcfs = {"freebayes": freebayes, "museq": museq, "rtg": rtg, "samtools": samtools}
    run_consensus_test(vcfs, "samtools", "freebayes", testdir)

def test_rtg_freebayes(freebayes, museq, rtg, samtools, testdir):
    vcfs = {"freebayes": freebayes, "museq": museq, "rtg": rtg, "samtools": samtools}
    run_consensus_test(vcfs, "rtg", "freebayes", testdir)

def test_normalization_1(freebayes, museq, rtg, samtools, testdir):
    # A-> AT and AT->ATT
    vcfs = {"freebayes": freebayes, "museq": museq, "rtg": rtg, "samtools": samtools}   
    r1 =  "A"
    a1 = "AT"
    r2 = "AT"
    a2 = "ATT"
    consensus_variant =  run_normalize_test(vcfs, "rtg", "samtools",  a1, r1, a2, r2, testdir)

    assert consensus_variant.REF == "A"
    assert consensus_variant.ALT == "AT"

def test_normalization_2(freebayes, museq, rtg, samtools, testdir):
    #  GTA -> G and GTATA -> GTA
    vcfs = {"freebayes": freebayes, "museq": museq, "rtg": rtg, "samtools": samtools}   
    r1 =  "GTA"
    a1 = "G"
    r2 = "GTATA"
    a2 = "GTA"
    consensus_variant =  run_normalize_test(vcfs, "rtg", "samtools",  a1, r1, a2, r2, testdir)

    assert consensus_variant.REF == "GTA"
    assert consensus_variant.ALT == "G"

def test_normalization_3(freebayes, museq, rtg, samtools, testdir):
    #   CA -> CTA  and C -> CT 
    vcfs = {"freebayes": freebayes, "museq": museq, "rtg": rtg, "samtools": samtools}   
    r1 =  "CA"
    a1 = "CTA"
    r2 = "C"
    a2 = "CT"
    consensus_variant =  run_normalize_test(vcfs, "rtg", "samtools",  a1, r1, a2, r2, testdir)

    assert consensus_variant.REF == "C"
    assert consensus_variant.ALT == "CT"

def test_normalization_4(freebayes, museq, rtg, samtools, testdir):
    #  T -> TAC and T -> TACACAC
    vcfs = {"freebayes": freebayes, "museq": museq, "rtg": rtg, "samtools": samtools}   
    r1 =  "T"
    a1 = "TAC"
    r2 = "T"
    a2 = "TACACAC"
    consensus_variant =  run_normalize_test(vcfs, "rtg", "samtools",  a1, r1, a2, r2, testdir)

    assert consensus_variant.REF == "T"
    assert consensus_variant.ALT == "TAC"

freebayes = "/juno/work/shah/isabl_data_lake/analyses/63/79/6379/results/germline/SA1181_N/SA1181_N_freebayes_germline.vcf.gz"
museq = "/juno/work/shah/isabl_data_lake/analyses/63/79/6379/results/germline/SA1181_N/SA1181_N_museq_single_annotated.vcf.gz"
rtg = "/juno/work/shah/isabl_data_lake/analyses/63/79/6379/results/germline/SA1181_N/SA1181_N_rtg_germline.vcf.gz"
samtools = "/juno/work/shah/isabl_data_lake/analyses/63/79/6379/results/germline/SA1181_N/SA1181_N_samtools_germline.vcf.gz"

test_museq_freebayes(freebayes, museq, rtg, samtools, "/juno/work/shah/abramsd/CODE/TEST")

# make_test("/juno/work/shah/abramsd/CODE/TEST/SA1181_N_freebayes_germline.vcf.gz",
# "/juno/work/shah/abramsd/CODE/TEST/SA1181_N_museq_single_annotated.vcf.gz")