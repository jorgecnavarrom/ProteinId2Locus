#!/usr/bin/env python

"""
ProtoTuna: Protein to neighborhood retriever.

Given a protein id, or a list of them, retrieve a number of kbs upstream and
downstream from NCBI.

To do this, ProtoTuna queries NCBI's identical protein database and chooses one
of the entries there to retrieve the context.

Input:
- A single protein id or
- A list with protein ids (one for each line)

Output:
- GenBank files with the context
- A report of the query, 'ipg_results.txt'. This file has the following info:
  * The title of the original protein id query, and length
  * A list of the database hits: accession, source, genomic accession, start and
    stop coordinates. 
    Example:
XP_014084778.1 (670)

XP_014084778.1	RefSeq	NW_014024938.1	1511240	1513715
XP_014084778.1	RefSeq	XM_014229303.1	1	2013
EMD96011.1	INSDC	KB445570.1	1532759	1535234
ENI10869.1	INSDC	KB733444.1	1511240	1513715

---

A note about NCBI accession prefixes:
https://www.ncbi.nlm.nih.gov/books/NBK21091/table/ch18.T.refseq_accession_numbers_and_mole/
"""

from pathlib import Path
import sys
import os
import argparse
from time import sleep
from Bio import Entrez


__author__ = "Jorge Navarro"
__version__ = "1.0"
__maintainer__ = "Jorge Navarro"
__email__ = "jorge.c.navarro.munoz@gmail.com"


class acc_data:
    def __init__(self, na="", src="", nuca="", begin=0, end=0):
        self.new_acc = na
        self.source = src
        self.nuc_acc = nuca
        self.start = begin
        self.stop = end

        return


def parameter_parser():
    def_extra = 20
    def_o = Path("./output")

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-l", "--list", help="File with protein accession list", type=Path
    )
    parser.add_argument(
        "-p", "--proteinid", help="A single protein id to query", type=str
    )
    parser.add_argument(
        "-e",
        "--extra",
        help=f"Extra kbps at either side of \
        the target protein-coding gene to download. Default={def_extra}",
        default=def_extra,
        type=int,
    )
    parser.add_argument(
        "-o",
        "--outputfolder",
        help=f"Where to put retrieved \
        files. Default={def_o}",
        type=Path,
        default=def_o,
    )

    return parser.parse_args()


def get_record_size(acc):
    e1 = Entrez.read(Entrez.epost("nuccore", id=acc))
    we = e1["WebEnv"]
    qk = e1["QueryKey"]
    x = Entrez.read(Entrez.esummary(db="nuccore", WebEnv=we, query_key=qk))
    return int(x[0]["Length"])


def get_best(current_best, candidate):
    # for the initial case (empty current_best)
    if current_best.new_acc == "":
        return candidate

    # discard RefSeq's standalone protein entries (starting with e.g. "XM_")
    if candidate.source == "RefSeq" and candidate.nuc_acc[0] == "X":
        return current_best
    # alternative
    if current_best.source == "RefSeq" and current_best.nuc_acc[0] == "X":
        return candidate

    # e.g. Swiss-Prot vs INSDC/RefSeq
    if current_best.nuc_acc == "" and candidate.nuc_acc != "":
        return candidate
    # alternative
    if candidate.nuc_acc == "" and current_best != "":
        return current_best

    # INSDC vs RefSeq/INSDC
    # easiest is a RefSeq version against its INSDC counterpart: they have
    # the same coordinates
    if current_best.start == candidate.start and current_best.stop == candidate.stop:
        if candidate.source == "RefSeq":
            return candidate
        else:
            return current_best

    # INSDC vs RefSeq: prefer RefSeq always (but probably not _always_ the best?)
    if current_best.source == "INSDC" and candidate.source == "RefSeq":
        return candidate
    # alt.
    if current_best.source == "RefSeq" and candidate.source == "INSDC":
        return current_best

    # Both are INSDC. Tricky one. Best strategy is to compare contig size...
    candidate_contig_l = get_record_size(candidate.nuc_acc)
    current_best_contig_l = get_record_size(current_best.nuc_acc)
    if candidate_contig_l >= current_best_contig_l:
        return candidate

    return current_best


def download_genbank(output, original_acc, best_candidate, extra):
    if best_candidate.nuc_acc == "":
        # got standalone protein, use entry accession instead
        filename = f"{original_acc}_{best_candidate.new_acc}_protein.fasta"
        handle = Entrez.efetch(
            db="protein", id=best_candidate.new_acc, rettype="fasta", retmode="text"
        )

        with open(output / filename, "w") as gb:
            gb.write(handle.read())

        # get, and return, header
        fasta = handle.read()
        header = fasta.split("\n")[0][1:]
        # TODO this is empty apparently...
        return f"\t{header}"
    else:
        filename = f"{original_acc}_{best_candidate.new_acc}_locus.gbk"
        start = max(0, best_candidate.start - (extra * 1000))
        # we don't have the length of the contig, but I think this should work (same above?)
        end = best_candidate.stop + (extra * 1000)
        handle = Entrez.efetch(
            db="nucleotide",
            id=best_candidate.nuc_acc,
            seq_start=str(start),
            seq_stop=str(end),
            rettype="gb",
            retmode="text",
        )
        with open(output / filename, "w") as gb:
            gb.write(handle.read())
        return f"{filename}"


def parse_accessions(args):
    """
    Reads arguments and parses the list of accession ids
    """

    # get protein ids
    accession_list = []
    accession_set = set()
    if args.list:
        acc_list_file = args.list
        if not acc_list_file.is_file():
            sys.exit("Provided list is not a valid file")
        # suppose the list is nicely formatted
        with open(acc_list_file) as alf:
            accession_list = [acc.strip() for acc in alf.readlines()]
        accession_set = set(accession_list)
    if args.proteinid:
        if args.proteinid not in accession_set:
            accession_list.append(args.proteinid)
        else:
            print(f"Warning: got {args.proteinid} as individual argument, but it's already included in the list")

    print(f"Got {len(accession_list)} accessions")
    accession_list = list(sorted(set(accession_list)))
    print(f"Got {len(accession_list)} unique accessions")

    return accession_list


def main():
    args = parameter_parser()

    Entrez.email = ""
    Entrez.api_key = ""
    Entrez.tool = "ProtoTuna"

    # Minimal checking of input
    if not (args.list or args.proteinid):
        sys.exit("No input. Stop")

    o = args.outputfolder
    if not o.is_dir():
        os.makedirs(o, exist_ok=True)

    # Get accessions
    accession_list = parse_accessions(args)

    # Validation: doesn't work too well as each 'post' generates a different 
    # query_key and then Entrez.efetch will only get results for a single query
    # Use first accession to get session stuff
    acc_to_qk = dict()
    e1 = Entrez.read(Entrez.epost("protein", id=accession_list[0]))
    we = e1["WebEnv"]
    qk = e1["QueryKey"]
    acc_to_qk[accession_list[0]] = qk

    # Query everything else in the same 'post' + try to detect anomalies with 
    # accessions
    for acc in accession_list[1:]:
        try:
            e = Entrez.read(Entrez.epost("protein", id=acc, Webenv=we))
        except RuntimeError:
            print(f"Problem with {acc}")
        else:
            acc_to_qk[acc] = e["QueryKey"]

    with open(o / "ipg_results.txt", "w") as ipg_res, \
        open(o / "accession_map.tsv", "w") as am:
        for acc, qk in acc_to_qk.items():
            iden_prots = Entrez.efetch(
                db="protein", rettype="ipg", retmode="xml", WebEnv=we, query_key=qk
            )
            sleep(0.34)
            try:
                ipgs = Entrez.read(iden_prots)
            except Entrez.Parser.ValidationError:
                # some of these don't work through Entrez, but do work on the website (?)
                print(
                    f"ERROR! Can't read efetch results for {acc}. Try manually: https://www.ncbi.nlm.nih.gov/ipg/?term={acc}"
                )
                continue
            except RuntimeError:
                print(
                    f"RuntimeError for {acc}. Try manually: https://www.ncbi.nlm.nih.gov/ipg/?term={acc}"
                )
                continue

            for report in ipgs.values():
                original_acc = report.attributes["product_acc"]
                slen = report["Product"].attributes["slen"]
                ipg_res.write(f"{original_acc} ({slen})\n\n")

                current_best = acc_data()

                # options: we only have a single, Swiss-Prot result
                for prot in report["ProteinList"]:
                    source = prot.attributes["source"]
                    acc = prot.attributes["accver"]
                    nuc_acc = ""
                    seq_start = ""
                    seq_stop = ""

                    if "CDSList" in prot:
                        for cds in prot["CDSList"]:
                            nuc_acc = cds.attributes["accver"]
                            seq_start = int(cds.attributes["start"])
                            seq_stop = int(cds.attributes["stop"])
                            ipg_res.write(
                                f"{acc}\t{source}\t{nuc_acc}\t{seq_start}\t{seq_stop}\n"
                            )

                            candidate = acc_data(
                                na=acc,
                                src=source,
                                nuca=nuc_acc,
                                begin=seq_start,
                                end=seq_stop,
                            )
                            current_best = get_best(current_best, candidate)
                    else:
                        # e.g. Swiss-Prot
                        ipg_res.write(f"{acc}\t{source}\t\t\t\n")

                        candidate = acc_data(na=acc, src=source, nuca=nuc_acc)
                        current_best = get_best(current_best, candidate)

                # gets the actual genbank/fasta
                filename_currentacc = download_genbank(
                    o, original_acc, current_best, args.extra
                )
                am.write(f"{original_acc}\t{filename_currentacc}\n")

                ipg_res.write("\n\n")

    print("\nDone")


if __name__ == "__main__":
    main()
    