from collections import defaultdict
import os
import pandas
from snakemake.utils import validate


def parse_samplelist(f, config, PREPROCESS):
    df = pandas.read_csv(f, sep="\t", dtype={"sampleID": str, "runID": int, "assemblyGroup": str,
                                         "fileName": str, "pair": str, "interleaved": str})
    validate(df, "config/samples.schema.yaml")
    dict_shell = lambda: defaultdict(dict_shell)  # Dictionary with arbitrary number of levels
    assemblyGroups = dict_shell()
    samples = dict_shell()
    df.fillna("", inplace=True)

    for i in list(df.index):
        sample = df.iloc[i]["sampleID"]
        runID = str(df.iloc[i]["runID"])
        R1 = df.iloc[i]["fileName"]
        groups = []
        r2 = False

        # Initiate keys for all assembly group values
        if "assemblyGroup" in df.columns:
            groups = df.iloc[i]["assemblyGroup"].split(",")
            # Remove empty assembly groups
            groups = [g for g in groups if g!= ""]
            for g in groups:
                if g not in assemblyGroups.keys() and g != "":
                    assemblyGroups[g] = dict_shell()

        if "interleaved" in df.columns and df.iloc[i]["interleaved"]:
            # If interleaved fastq is provided, add filepaths to split fastq files and later produce these using the
            # deinterleave_fastq rule in preprocessing.rules.
            inter = R1
            R1 = os.path.join(config["intermediate_path"], "deinterleaved", "{}_{}_R1.fastq.gz".format(sample, runID))
            R2 = os.path.join(config["intermediate_path"], "deinterleaved", "{}_{}_R2.fastq.gz".format(sample, runID))
            samples[sample][runID]["interleaved"] = inter
            samples[sample][runID]["R1"] = R1
            samples[sample][runID]["R2"] = R2
            for g in groups:
                assemblyGroups[g][sample][runID]["R1"] = [os.path.join(config["intermediate_path"], "preprocess",
                                         "{}_{}_R1{}.fastq.gz".format(sample, runID, PREPROCESS))]
                assemblyGroups[g][sample][runID]["R2"] = [os.path.join(config["intermediate_path"], "preprocess",
                                         "{}_{}_R2{}.fastq.gz".format(sample, runID, PREPROCESS))]
            continue

        # Handling of paired and/or single end sequence files
        # If the sample annotation file has a 'pair' column, add the read files as 'R1' and 'R2'
        if "pair" in df.columns and df.iloc[i]["pair"]:
            R2 = df.iloc[i]["pair"]
            r2 = True
            samples[sample][runID]["R1"] = R1
            samples[sample][runID]["R2"] = R2
            # Add filepaths to preprocessed output files for each of the read files in each of the assembly groups
            # This will be the initial input to the assembly rule
            for g in groups:
                if r2:
                    assemblyGroups[g][sample][runID]["R1"] = [os.path.join(config["intermediate_path"], "preprocess",
                                         "{}_{}_R1{}.fastq.gz".format(sample, runID, PREPROCESS))]
                    assemblyGroups[g][sample][runID]["R2"] = [os.path.join(config["intermediate_path"], "preprocess",
                                         "{}_{}_R2{}.fastq.gz".format(sample, runID, PREPROCESS))]
                else:
                    assemblyGroups[g][sample][runID]["se"] = [os.path.join(config["intermediate_path"], "preprocess",
                                        "{}_{}_se{}.fastq.gz".format(sample, runID, PREPROCESS))]

        # If there is no 'pair' column, add the single file path as 'se'
        else:
            samples[sample][runID]["se"] = R1
            for g in groups:
                assemblyGroups[g][sample][runID]["se"] = [os.path.join(config["intermediate_path"], "preprocess",
                                         "{}_{}_se{}.fastq.gz".format(sample, runID, PREPROCESS))]
    return samples, assemblyGroups

def check_sequencing_type(samples):
  types = []
  for sample in samples.keys():
    for run in samples[sample].keys():
      if "R2" in samples[sample][run].keys():
        types.append("paired")
      else:
        types.append("single")
  if len(set(types))>1: return 'mixed'
  else: return types[0]
