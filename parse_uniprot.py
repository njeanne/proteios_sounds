#! /usr/bin/env python3

import sys
import re
import urllib
from Bio import ExPASy
from Bio import SwissProt
from Bio import SeqIO

def parse_entry(uniprot_AN, out_dir, logger):
    '''
    Parse the uniprot entry.
    :param int uniprot_AN: the uniprot accession number for the protein
    :param str out_dir: the path to the output directory
    :param logger logger: the logger
    :return: the dictionary describing the protein
    :rtype: dict
    '''
    ### Load the Uniprot entry
    # uniprot_file_path = os.path.join(out_dir, '{}.dat'.format(uniprot_AN))
    try:
        handle = ExPASy.get_sprot_raw(uniprot_AN)
        uniprot = SwissProt.read(handle)
    except urllib.error.HTTPError as err:
        msg = '{}: {} accession number does not exists in UniProt database, check https://www.uniprot.org/ for the correct accession number.'.format(err, uniprot_AN)
        print('ERROR: {}'.format(msg))
        logger.error(msg)
        sys.exit(1)
    #TODO:
    #coder erreur uniprot indisponible
    # except urllib.error.ATROUVER:
    #     msg = 'Uniprot web site unavailable (https://www.uniprot.org/)'
    #     if not os.path.exists(uniprot_file_path):
    #         msg = '{}, try again later!'
    #         print('ERROR: {}'.format(msg))
    #         logger.error(msg)
    #         sys.exit(1)
    #     else:
    #         msg = '{}, using previously downloaded UniProt data: {}'.format(msg, uniprot_file_path)
    #         uniprot = SwissProt.read(uniprot_file_path)

    ### Parse the Uniprot entry
    # get the organism and entry name in UniProt
    organism = uniprot.organism
    if '(' in organism:
        organism = organism.split('(')[1].split(')')[0]
    organism = organism.rstrip('.,;:').replace(' ', '_')
    entry_name = uniprot.entry_name
    uniprot_sequence = uniprot.sequence

    # create a dictionary for the protein
    protein = {'seq': uniprot_sequence, 'organism': organism, 'entry_name': entry_name}

    # create a dictionary to get the structures positions => type
    structures = {}
    # parse the UniProt Features (see https://www.uniprot.org/help/?query=*&fil=category%3A%22PTM+processing%22+AND+section%3Amanual&columns=title)
    # for feature in uniprot.features:
    #     print(feature)
    # print('###########################################\n\n')
    for feature in uniprot.features:
        # strucutre information
        if feature[0] == 'HELIX' or feature[0] == 'STRAND' or feature[0] == 'TURN':
            structures[feature[1]] = {'type': feature[0], 'end': feature[2]}
        # signal peptide
        if feature[0] == 'SIGNAL':
            if 'signal_peptide' in protein.keys():
                protein['signal_peptide'].append((feature[1], feature[2]))
            else:
                protein['signal_peptide'] = [(feature[1], feature[2])]
        # Modified residue: Modified residues excluding lipids, glycans and protein crosslinks
        if feature[0] == 'MOD_RES':
            modif = re.split(' {', feature[3])[0].rstrip('.;,')
            if 'modified_residue' not in protein.keys():
                protein['modified_residue'] = {modif: [(feature[1], feature[2])]}
            else:
                if modif in protein['modified_residue'].keys():
                    protein['modified_residue'][modif].append((feature[1], feature[2]))
                else:
                    protein['modified_residue'][modif] = [(feature[1], feature[2])]
        # Glycosylation: Covalently attached glycan group(s)
        if feature[0] == 'CARBOHYD':
            glyco = re.split(' ', feature[3])[0].rstrip('.;,')
            if 'glycosylation' not in protein.keys():
                protein['glycosylation'] = {glyco: [(feature[1], feature[2])]}
            else:
                if glyco in protein['glycosylation'].keys():
                    protein['glycosylation'][glyco].append((feature[1], feature[2]))
                else:
                    protein['glycosylation'][glyco] = [(feature[1], feature[2])]
        # Site: Any interesting single amino acid site on the sequence
        if feature[0] == 'SITE':
            site = re.split(' {', feature[3])[0].rstrip('.;,')
            if 'site' not in protein:
                protein['site'] = {site: [(feature[1], feature[2])]}
            else:
                if site in protein['site']:
                    protein['site'][site].append((feature[1], feature[2]))
                else:
                    protein['site'][site] = [(feature[1], feature[2])]
        # Propeptide: Part of a protein that is cleaved during maturation or activation
        if feature[0] == 'PROPEP':
            propep = re.split(' {', feature[3])[0].rstrip('.;,')
            if 'propeptide' not in protein:
                protein['propeptide'] = {propep: [(feature[1], feature[2])]}
            else:
                if propep in protein['propeptide'].keys():
                    protein['propeptide'][propep].append((feature[1], feature[2]))
                else:
                    protein['propeptide'][propep] = [(feature[1], feature[2])]
        # Disulfid bond, get the two linked positions on the same chain (no interchain)
        if feature[0] == 'DISULFID':
            if not feature[3].startswith('Interchain'):
                if 'disulfid' in protein.keys():
                    protein['disulfid'].append((feature[1], feature[2]))
                else:
                    protein['disulfid'] = [(feature[1], feature[2])]


    # update protein with the structures
    sequence_length = len(protein['seq'])
    structure_last_position = 0
    for position in sorted(structures.keys()):
        if structure_last_position == 0 and position == 1:
            if 'structure' in protein.keys():
                protein['structure'][position] = structures[position]['type']
            else:
                protein['structure'] = {position: structures[position]['type']}
            structure_last_position = structures[position]['end']
        elif position >= structure_last_position + 1:
            if 'structure' in protein.keys():
                protein['structure'][structure_last_position + 1] = 'FREE'
            else:
                protein['structure'] = {structure_last_position + 1: 'FREE'}
            protein['structure'][position] = structures[position]['type']
            structure_last_position = structures[position]['end']
    if structure_last_position < sequence_length:
        if 'structure' in protein.keys():
            protein['structure'][structure_last_position + 1] = 'FREE'
        else:
            protein['structure'] = {structure_last_position + 1: 'FREE'}

    # ### TOREMOVE
    # print('#####################################################\n\n')
    # print('structures :')
    # for structure in structures:
    #     print(structure)
    # print('#####################################################\n\n')
    # for k, v in protein.items():
    #     if k == 'modified_residue':
    #         print(k)
    #         for k2, v2 in v.items():
    #             print('\t{}'.format(k2))
    #             for pos in v2:
    #                 print('\t\t{}'.format(pos))
    #     elif k == 'glycosylation':
    #         print(k)
    #         for k2, v2 in v.items():
    #             print('\t{}'.format(k2))
    #             for pos in v2:
    #                 print('\t\t{}'.format(pos))
    #     elif k == 'site':
    #         print(k)
    #         for k2, v2 in v.items():
    #             print('\t{}'.format(k2))
    #             for pos in v2:
    #                 print('\t\t{}'.format(pos))
    #     elif k == 'propeptide':
    #         print(k)
    #         for k2, v2 in v.items():
    #             print('\t{}'.format(k2))
    #             for pos in v2:
    #                 print('\t\t{}'.format(pos))
    #     elif k == 'disulfid':
    #         print(k)
    #         for bond in v:
    #             print('\t{} disulfid bond to {}'.format(bond[0], bond[1]))
    #     elif k == 'signal_peptide':
    #         print(k)
    #         for signal_peptide in v:
    #             print('\tsignal peptide from {} to {}'.format(signal_peptide[0], signal_peptide[1]))
    # #################

    return protein
