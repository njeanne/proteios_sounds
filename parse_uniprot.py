#! /usr/bin/env python3

import sys
import logging
import re
import urllib
from Bio import ExPASy
from Bio import SwissProt


def midi_keys_from_seq(sequence, protein_dict, list_keys_midi):
    '''Add the ordered MIDI keys, depending on the frequency of the amino acids
    and the protein sequence, to the protein dictionary.

    :param str sequence: the protein amino acid sequence
    :param dict protein_dict: the representation of the protein as a dictionary
    :param list list_keys_midi: the list MIDI values with the initial order
    :return: the protein dictionary
    :rtype: dict
    '''
    protein_dict['seq'] = {}
    logging.debug('AA sequence ({} AA): {}'.format(len(sequence), sequence))
    for i in range(len(sequence)):
        protein_dict['seq'][i] = sequence[i]
    # frequence of AA in the sequence
    set_AA = set(''.join(sequence))
    proportion_AA = {}
    for aa in set_AA:
        proportion_AA[aa] = sequence.count(aa) / len(sequence)
    # sort by decreasing frequency
    proportion_AA = sorted(proportion_AA.items(),
                           key=lambda kv: kv[1],
                           reverse=True)

    protein_dict['midi_keys'] = {}
    for idx, aa_proportion in enumerate(proportion_AA):
        protein_dict['midi_keys'][aa_proportion[0]] = list_keys_midi[idx]

    return protein_dict


def get_uniprot_features(entry, protein_dict):
    '''Update the protein dictionary with the protein structures.

    :param dict entry: the uniprot dictionary.
    :param dict protein_dict: the dictionary describing the protein.
    :return: the protein dictionary
    :rtype: dict
    '''
    for feature in entry.features:
        # signal peptide
        if feature[0] == 'SIGNAL':
            if 'signal_peptide' in protein_dict.keys():
                protein_dict['signal_peptide'].append((feature[1], feature[2]))
            else:
                protein_dict['signal_peptide'] = [(feature[1], feature[2])]
        # Modified residue: Modified residues excluding lipids, glycans and protein_dict crosslinks
        if feature[0] == 'MOD_RES':
            modif = re.split(' {', feature[3])[0].rstrip('.;,')
            if 'modified_residue' not in protein_dict.keys():
                protein_dict['modified_residue'] = {modif: [(feature[1], feature[2])]}
            else:
                if modif in protein_dict['modified_residue'].keys():
                    protein_dict['modified_residue'][modif].append((feature[1], feature[2]))
                else:
                    protein_dict['modified_residue'][modif] = [(feature[1], feature[2])]
        # Glycosylation: Covalently attached glycan group(s)
        if feature[0] == 'CARBOHYD':
            glyco = re.split(' ', feature[3])[0].rstrip('.;,')
            if 'glycosylation' not in protein_dict.keys():
                protein_dict['glycosylation'] = {glyco: [(feature[1], feature[2])]}
            else:
                if glyco in protein_dict['glycosylation'].keys():
                    protein_dict['glycosylation'][glyco].append((feature[1], feature[2]))
                else:
                    protein_dict['glycosylation'][glyco] = [(feature[1], feature[2])]
        # Site: Any interesting single amino acid site on the sequence
        if feature[0] == 'SITE':
            site = re.split(' {', feature[3])[0].rstrip('.;,')
            if 'site' not in protein_dict:
                protein_dict['site'] = {site: [(feature[1], feature[2])]}
            else:
                if site in protein_dict['site']:
                    protein_dict['site'][site].append((feature[1], feature[2]))
                else:
                    protein_dict['site'][site] = [(feature[1], feature[2])]
        # Propeptide: Part of a protein that is cleaved during maturation or activation
        if feature[0] == 'PROPEP':
            propep = re.split(' {', feature[3])[0].rstrip('.;,')
            if 'propeptide' not in protein_dict:
                protein_dict['propeptide'] = {propep: [(feature[1], feature[2])]}
            else:
                if propep in protein_dict['propeptide'].keys():
                    protein_dict['propeptide'][propep].append((feature[1], feature[2]))
                else:
                    protein_dict['propeptide'][propep] = [(feature[1], feature[2])]
        # Disulfid bond, get the two linked positions on the same chain (no interchain)
        if feature[0] == 'DISULFID':
            if not feature[3].startswith('Interchain'):
                if 'disulfid' in protein_dict.keys():
                    protein_dict['disulfid'].append((feature[1], feature[2]))
                else:
                    protein_dict['disulfid'] = [(feature[1], feature[2])]
    return protein_dict


def get_uniprot_structures(entry, protein_dict):
    '''Update the protein dictionary with the protein structures.

    :param dict entry: the uniprot dictionary.
    :param dict protein_dict: the dictionary describing the protein.
    :return: the protein dictionary
    :rtype: dict
    '''
    # create a dictionary to get the structures positions => type
    structures = {}
    # parse the UniProt Features (see https://www.uniprot.org/help/?query=*&fil=category%3A%22PTM+processing%22+AND+section%3Amanual&columns=title)
    for feature in entry.features:
        # strucutre information
        if feature[0] == 'HELIX' or feature[0] == 'STRAND' or feature[0] == 'TURN':
            # get START as key and a dictionay TYPE: END
            structures[feature[1]] = {'type': feature[0], 'end': feature[2]}

    sequence_length = len(protein_dict['seq'])
    structure_last_position = 0
    for position in sorted(structures.keys()):
        if structure_last_position == 0 and position == 1:
            if 'structure' in protein_dict.keys():
                protein_dict['structure'][position] = structures[position]['type']
            else:
                protein_dict['structure'] = {position: structures[position]['type']}
            structure_last_position = structures[position]['end']
        elif position >= structure_last_position + 1:
            if 'structure' in protein_dict.keys():
                protein_dict['structure'][structure_last_position + 1] = 'FREE'
            else:
                protein_dict['structure'] = {structure_last_position + 1: 'FREE'}
            protein_dict['structure'][position] = structures[position]['type']
            structure_last_position = structures[position]['end']
    if structure_last_position < sequence_length:
        if 'structure' in protein_dict.keys():
            protein_dict['structure'][structure_last_position + 1] = 'FREE'
        else:
            protein_dict['structure'] = {structure_last_position + 1: 'FREE'}
    return protein_dict


def parse_entry(uniprot_AN, midi_keys):
    '''
    Parse the uniprot entry.

    :param int uniprot_AN: the uniprot accession number for the protein
    :param list midi_keys: the MIDI keys list
    :return: the dictionary describing the protein
    :rtype: dict
    '''

    # Load the Uniprot entry
    logging.info('UniProt accession number: {}'.format(uniprot_AN))
    try:
        handle = ExPASy.get_sprot_raw(uniprot_AN)
        uniprot = SwissProt.read(handle)
        pdb_acc_num = None
        for cross_ref in uniprot.cross_references:
            if cross_ref[0] == 'PDB':
                pdb_acc_num = cross_ref[1]
                # print('{}: {}\t{}'.format(cross_ref[0], cross_ref[1], cross_ref))
                break
    except urllib.error.HTTPError as http_err:
        if http_err.code == 404:
            err_msg = '{} {} accession number does not exists in UniProt database, check https://www.uniprot.org/ to get the correct accession number.'.format(http_err, uniprot_AN)
        elif http_err.code == 504:
            err_msg = '{} UniProt database web site unavailable.'.format(http_err)
        else:
            err_msg = http_err
        print('ERROR: {}'.format(err_msg))
        logging.error(err_msg)
        sys.exit(1)

    # Parse the Uniprot entry
    # get the organism and entry name in UniProt
    entry_name = uniprot.entry_name
    logging.info('Protein: {}'.format(entry_name))
    organism = uniprot.organism
    organism = organism.split(' (')[0]
    organism = re.sub('[^A-Za-z0-9]+', '_', organism)
    logging.info('Organism: {}'.format(organism))

    # create a dictionary for the protein
    protein = {'accession_number': uniprot_AN,
               'seq': uniprot.sequence,
               'organism': organism,
               'entry_name': entry_name}

    # if the PDB ID was found add it
    if pdb_acc_num:
        protein['PDB'] = pdb_acc_num
        logging.info('PDB: {} (Protein DataBase accession number, see http://www.rcsb.org/)'.format(pdb_acc_num))
    else:
        logging.info('PDB: No accession number in Uniprot entry')

    # update protein with the structures
    protein = get_uniprot_structures(uniprot, protein)

    # add features of the protein
    protein = get_uniprot_features(uniprot, protein)

    # add the MIDI keys ordered depending on the sequence and AA frequency
    protein = midi_keys_from_seq(uniprot.sequence, protein, midi_keys)

    ### TOREMOVE
    # print('#####################################################')
    # print('Keys in protein dictionary: {}'.format(protein.keys()))
    # for k, v in protein.items():
    #     print('{}:'.format(k))
    #     if k == 'organism':
    #         print('\t{}'.format(v))
    #     if k == 'entry_name':
    #         print('\t{}'.format(v))
    #     if k == 'seq':
    #         print('\tlength: {}\n\t{}'.format(len(v), v))
    #     if k == 'structure':
    #         for start_pos, structure in v.items():
    #             print('\t{}: {}'.format(structure, start_pos))
    #     if k == 'modified_residue':
    #         for k2, v2 in v.items():
    #             print('\t{}'.format(k2))
    #             for pos in v2:
    #                 print('\t\t{}'.format(pos))
    #     elif k == 'glycosylation':
    #         for k2, v2 in v.items():
    #             print('\t{}'.format(k2))
    #             for pos in v2:
    #                 print('\t\t{}'.format(pos))
    #     elif k == 'site':
    #         for k2, v2 in v.items():
    #             print('\t{}'.format(k2))
    #             for pos in v2:
    #                 print('\t\t{}'.format(pos))
    #     elif k == 'propeptide':
    #         for k2, v2 in v.items():
    #             print('\t{}'.format(k2))
    #             for pos in v2:
    #                 print('\t\t{}'.format(pos))
    #     elif k == 'disulfid':
    #         for bond in v:
    #             print('\t{} disulfid bond to {}'.format(bond[0], bond[1]))
    #     elif k == 'signal_peptide':
    #         for signal_peptide in v:
    #             print('\tsignal peptide from {} to {}'.format(signal_peptide[0], signal_peptide[1]))
    # print('#####################################################\n\n')

    return protein
