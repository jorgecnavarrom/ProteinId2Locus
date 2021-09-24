# ProteinId2Locus

A simple script to retrieve the genetic context around a given protein-encoding gene. Specifically, it queries the identical protein groups (IPG) database.

## Requirements

* Python 3
* Biopython

A conda environment is recommended

## Usage

```
usage: protein_id_2_genomic_context.py [-h] [-l LIST] [-p PROTEINID] [-e EXTRA] [-o OUTPUTFOLDER]

optional arguments:
  -h, --help            show this help message and exit
  -l LIST, --list LIST  File with protein accession list
  -p PROTEINID, --proteinid PROTEINID
                        A single protein id to query
  -e EXTRA, --extra EXTRA
                        Extra kbps at either side of the target protein-coding gene to download. Default=20
  -o OUTPUTFOLDER, --outputfolder OUTPUTFOLDER
                        Where to put retrieved files. Default=output
```

## Input

The input is either a single protein id, or a list of protein ids. If using the file, each line must only contain a single protein id.

## Output

* GenBank files. Naming is `[original protein name]_[best result in IPG]_locus.gbk`
* `ipg_results.txt`. IPG query results, per entry
* `accession_map.tsv`. A tab-separated file linking the original protein id with the resulting GenBank filename
