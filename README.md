# AFLA Demux

Sorts sequencing reads into one set of FASTQ files per sample, for **Illumina** or **Nanopore** data, using a sample sheet.
Runs in [EPI2ME Desktop](https://epi2me.nanoporetech.com/). **Education and research use only. Not for clinical use.**

## Before you run
Install EPI2ME Desktop with Docker. In EPI2ME choose **Import workflow** and give this repository's address.

## Illumina
Starts from FASTQ files such as `Undetermined_S0_L001_R1_001.fastq.gz`. Indexes are read from the read names, or from `_I1` and
`_I2` files next to the reads. Converting raw BCL files needs Illumina's own software.

Sample sheet: an Illumina `SampleSheet.csv`, or a CSV like this.

```text
sample,index,index2
Patient01,ATCACGTT,GGCTACAG
Patient02,CGATGTAA,TTAGGCAT
```

- One mismatch per index is allowed by default. The sheet is refused if two samples could be confused.
- The direction of index 2 differs between instruments and is detected from the reads.
- Reads whose two indexes belong to different samples are reported as possible index hopping and are never assigned.

## Nanopore
For custom or third-party barcodes. For standard ONT barcoding kits, MinKNOW already demultiplexes. Starts from FASTQ.

```text
sample,start,end
Patient01,ACGT...full sequence before the insert...,TGCA...full sequence after the insert...
```

`end` is optional; if left out, it is taken as the reverse complement of `start`. `start` is the barcode with any fixed sequence
around it. The barcode is searched for at both ends of each read, and reads too close to two barcodes stay unclassified.

## Results
`<sample>/` holds the FASTQ files, `demux_report.html` the per-sample counts and the reasons reads were not assigned. These sample
folders can be used directly as input for the AFLA CFTR workflows.

## Testing status
Illumina: tested on a simulated dual-index run including index errors, unknown indexes and index hopping. Not yet tested on a real
run. Nanopore: tested on real amplicon reads with simulated custom barcodes and barcode errors added. Research use: check the
per-sample counts in the report for every run.
