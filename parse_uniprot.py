#! /usr/bin/env python3

import logging
import sys
import re
import urllib
from Bio import ExPASy
from Bio import SwissProt


def parse_entry(uniprot_accession_number):
    """
    Parse the uniprot entry.

    :param uniprot_accession_number: the uniprot accession number for the protein
    :type uniprot_accession_number: str
    :return: the dictionary describing the protein
    :rtype: dict
    """
    logging.info(f"Parsing uniprot entry: {uniprot_accession_number}")

    # Load the Uniprot entry
    # uniprot_file_path = os.path.join(out_dir, "{}.dat".format(uniprot_accession_number))
    try:
        handle = ExPASy.get_sprot_raw(uniprot_accession_number)
        uniprot = SwissProt.read(handle)
        pdb_acc_num = None
        cross_ref = None
        for cross_ref in uniprot.cross_references:
            if cross_ref[0] == "PDB":
                pdb_acc_num = cross_ref[1]
                break
    except urllib.error.HTTPError as http_err:
        if http_err.code == 404:
            msg = (f"{http_err} {uniprot_accession_number} accession number does not exists in UniProt database, check "
                   f"https://www.uniprot.org/ to get the correct accession number.")
        elif http_err.code == 504:
            msg = f"{http_err} UniProt database web site unavailable."
        else:
            msg = http_err
        logging.error(msg)
        sys.exit(1)

    # Parse the Uniprot entry
    # get the organism and entry name in UniProt
    organism = uniprot.organism
    organism = organism.split(" (")[0]
    organism = re.sub("[^A-Za-z0-9]+", "_", organism)

    entry_name = uniprot.entry_name
    uniprot_sequence = uniprot.sequence

    # create a dictionary for the protein
    protein = {"accession_number": uniprot_accession_number, "seq": uniprot_sequence, "organism": organism,
               "entry_name": entry_name}

    # if the PDB ID was found, add it
    if pdb_acc_num:
        protein["PDB"] = pdb_acc_num

    # create a dictionary to get the structure positions => type
    structures = {}
    # parse the UniProt Features 
    # (see https://www.uniprot.org/help/?query=*&fil=category%3A%22PTM+processing%22+AND+section%3Amanual&columns=title)
    for feature in uniprot.features:
        # structure information
        if feature.type in ["HELIX", "STRAND", "TURN"]:
            structures[feature.location.start] = {"type": feature.type, "end": feature.location.end}
        # signal peptide
        if feature.type == "SIGNAL":
            if "signal_peptide" in protein.keys():
                protein["signal_peptide"].append((feature.location.start, feature.location.end))
            else:
                protein["signal_peptide"] = [(feature.location.start, feature.location.end)]
        # Modified residue: Modified residues excluding lipids, glycans and protein cross-links
        if feature.type == "MOD_RES":
            modif = feature.qualifiers["note"].rstrip(".;,")
            if "modified_residue" not in protein.keys():
                protein["modified_residue"] = {modif: [(feature.location.start, feature.location.end)]}
            else:
                if modif in protein["modified_residue"].keys():
                    protein["modified_residue"][modif].append((feature.location.start, feature.location.end))
                else:
                    protein["modified_residue"][modif] = [(feature.location.start, feature.location.end)]
        # Glycosylation: Covalently attached glycan group(s)
        if feature.type == "CARBOHYD":
            glyco = feature.qualifiers["note"].rstrip(".;,")
            if "glycosylation" not in protein.keys():
                protein["glycosylation"] = {glyco: [(feature.location.start, feature.location.end)]}
            else:
                if glyco in protein["glycosylation"].keys():
                    protein["glycosylation"][glyco].append((feature.location.start, feature.location.end))
                else:
                    protein["glycosylation"][glyco] = [(feature.location.start, feature.location.end)]
        # Site: Any interesting single amino acid site on the sequence
        if feature.type == "SITE":
            site = feature.qualifiers["note"].rstrip(".;,")
            if "site" not in protein:
                protein["site"] = {site: [(feature.location.start, feature.location.end)]}
            else:
                if site in protein["site"]:
                    protein["site"][site].append((feature.location.start, feature.location.end))
                else:
                    protein["site"][site] = [(feature.location.start, feature.location.end)]
        # Pro-peptide: Part of a protein that is cleaved during maturation or activation
        if feature.type == "PROPEP":
            propep = feature.qualifiers["note"].rstrip(".;,")
            if "propeptide" not in protein:
                protein["propeptide"] = {propep: [(feature.location.start, feature.location.end)]}
            else:
                if propep in protein["propeptide"].keys():
                    protein["propeptide"][propep].append((feature.location.start, feature.location.end))
                else:
                    protein["propeptide"][propep] = [(feature.location.start, feature.location.end)]
        # Disulfure bond, get the two linked positions on the same chain (no interchain)
        if feature.type == "DISULFID":
            if "note" in feature.qualifiers and not feature.qualifiers["note"].startswith("Interchain"):
                if "disulfid" in protein.keys():
                    protein["disulfid"].append((feature.location.start, feature.location.end))
                else:
                    protein["disulfid"] = [(feature.location.start, feature.location.end)]

    # print("****************************************")
    # for start in structures:
    #     print("{}: start {}, end {}".format(structures[start]["type"], start, structures[start]["end"]))
    # print("****************************************")

    # update protein with the structures
    sequence_length = len(protein["seq"])
    structure_last_position = 0
    for position in sorted(structures.keys()):
        if structure_last_position == 0 and position == 1:
            if "structure" in protein.keys():
                protein["structure"][position] = structures[position]["type"]
            else:
                protein["structure"] = {position: structures[position]["type"]}
            structure_last_position = structures[position]["end"]
        elif position >= structure_last_position + 1:
            if "structure" in protein.keys():
                protein["structure"][structure_last_position + 1] = "FREE"
            else:
                protein["structure"] = {structure_last_position + 1: "FREE"}
            protein["structure"][position] = structures[position]["type"]
            structure_last_position = structures[position]["end"]
    if structure_last_position < sequence_length:
        if "structure" in protein.keys():
            protein["structure"][structure_last_position + 1] = "FREE"
        else:
            protein["structure"] = {structure_last_position + 1: "FREE"}

    logging.info(f"\tProtein: {protein['entry_name']}")
    logging.info(f"\tOrganism: {protein['organism']}")
    if "PDB" in protein:
        logging.info(f"\tPDB accession number: {protein['PDB']}\t{cross_ref}")
    else:
        logging.info("\tPDB: No accession number in Uniprot entry")

    ### TOREMOVE
    # print("#####################################################")
    # print("Keys in protein dictionary: {}".format(protein.keys()))
    # for k, v in protein.items():
    #     print("{}:".format(k))
    #     if k == "organism":
    #         print("\t{}".format(v))
    #     if k == "entry_name":
    #         print("\t{}".format(v))
    #     if k == "seq":
    #         print("\tlength: {}\n\t{}".format(len(v), v))
    #     if k == "structure":
    #         for start_pos, structure in v.items():
    #             print("\t{}: {}".format(structure, start_pos))
    #     if k == "modified_residue":
    #         for k2, v2 in v.items():
    #             print("\t{}".format(k2))
    #             for pos in v2:
    #                 print("\t\t{}".format(pos))
    #     elif k == "glycosylation":
    #         for k2, v2 in v.items():
    #             print("\t{}".format(k2))
    #             for pos in v2:
    #                 print("\t\t{}".format(pos))
    #     elif k == "site":
    #         for k2, v2 in v.items():
    #             print("\t{}".format(k2))
    #             for pos in v2:
    #                 print("\t\t{}".format(pos))
    #     elif k == "propeptide":
    #         for k2, v2 in v.items():
    #             print("\t{}".format(k2))
    #             for pos in v2:
    #                 print("\t\t{}".format(pos))
    #     elif k == "disulfid":
    #         for bond in v:
    #             print("\t{} disulfid bond to {}".format(bond[0], bond[1]))
    #     elif k == "signal_peptide":
    #         for signal_peptide in v:
    #             print("\tsignal peptide from {} to {}".format(signal_peptide[0], signal_peptide[1]))
    # print("#####################################################\n\n")

    return protein
