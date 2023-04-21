import logging
import os
import pyfastx


def infer_bin_name_from_bin_inputs(input_bins):
    # remove common prefix of bin dir to get nicer label
    commonprefix_len = len(os.path.commonprefix(input_bins))
    reversed_strings = ["".join(reversed(s)) for s in input_bins]
    commonsufix_len = len(os.path.commonprefix(reversed_strings))

    bin_name_to_bin_dir = {d[commonprefix_len : len(d) - commonsufix_len]: d for d in input_bins}

    logging.debug(f"input bins : {input_bins}")
    logging.debug(f"commonprefix {os.path.commonprefix(input_bins)}")
    logging.debug(f"commonsuffix {os.path.commonprefix(reversed_strings)}")
    logging.debug(f"bin_name_to_bin_dir {bin_name_to_bin_dir}")

    return bin_name_to_bin_dir


def write_bin_info(bins, output):

    header = ["bin_id", "origin", "name", "completeness", "contamination", "score", "size", "N50", "contig_count"]
    with open(output, "w") as fl:
        fl.write("\t".join(header) + "\n")
        for bin_obj in bins:

            line = [
                bin_obj.id,
                bin_obj.origin,
                bin_obj.name,
                bin_obj.completeness,
                bin_obj.contamination,
                bin_obj.score,
                bin_obj.length,
                bin_obj.N50,
                len(bin_obj.contigs),
            ]

            fl.write("\t".join((str(e) for e in line)) + "\n")


def write_bin_info_debug(bins, output):

    header = [
        "origin",
        "name",
        "id",
        "completeness",
        "contamination",
        "score",
        "size",
        "N50",
        "contig_count",
        "contigs",
    ]

    with open(output, "w") as fl:
        fl.write("\t".join(header) + "\n")
        for bin_obj in bins:

            line = [
                bin_obj.origin,
                bin_obj.name,
                bin_obj.id,
                bin_obj.completeness,
                bin_obj.contamination,
                bin_obj.score,
                bin_obj.length,
                bin_obj.N50,
                len(bin_obj.contigs),
                ";".join((str(c) for c in bin_obj.contigs)),
            ]

            fl.write("\t".join((str(e) for e in line)) + "\n")


def write_bins_fasta(selected_bins, contigs_fasta, outdir):
    fa = pyfastx.Fasta(contigs_fasta, build_index=True)

    for sbin in selected_bins:
        outfile = os.path.join(outdir, f"bin_{sbin.id}.fa")

        with open(outfile, "w") as outfl:
            sequences = (f">{c}\n{fa[c]}" for c in sbin.contigs)
            outfl.write("\n".join(sequences) + "\n")


def check_contig_consistency(contigs_from_assembly, contigs_from_elsewhere, assembly_file, elsewhere_file):
    logging.debug("check_contig_consistency.")
    are_contigs_consistent = len(set(contigs_from_elsewhere) | set(contigs_from_assembly)) <= len(
        set(contigs_from_assembly)
    )

    issue_countigs = len(set(contigs_from_elsewhere) - set(contigs_from_assembly))
    message = f"{issue_countigs} contigs found in file {elsewhere_file} \
                were not found in assembly_file ({assembly_file})."
    assert are_contigs_consistent, message
