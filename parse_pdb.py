#! /usr/bin/env python3

import os
import sys
import logging
# add pymol to the python path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'ressources', 'lib',
                                'python3.7', 'site-packages'))
import pymol


def hamming(s1, s2):
    '''Calculate the Hamming distance between two strings.

    :param str s1: the first string.
    :param str s2: the second string.
    :return: the number of differences.
    :rtype: int
    '''

    assert len(s1) == len(s2)
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))


def get_common_coordinates(pdb_seq, uniprot_seq, pdb_coord_matching):
    '''Look for the common part between the PDB and uniprot sequences.

    :param dict pdb_coord_matching: the dictionary
    :param str pdb_seq: the PDB sequence.
    :param str uniprot_seq: the UniProt sequence.
    :return: the start, stop, the PDB chain, number of differences and if the PDB sequence is longer than the uniprot.
    :rtype: dictionary.
    '''

    if pdb_coord_matching['uniprot_longer_than_pdb']:
        short_seq = pdb_seq
        long_seq = uniprot_seq
    else:
        short_seq = uniprot_seq
        long_seq = pdb_seq
    # look for the matching part between the two sequences
    for i in range(0, len(long_seq) - len(short_seq) + 1):
        ham_dist = hamming(short_seq, long_seq[i:i+len(short_seq)])
        if ham_dist < pdb_coord_matching['nb_discrepancies']:
            pdb_coord_matching['nb_discrepancies'] = ham_dist
            pdb_coord_matching['PDB_start_on_uniprot'] = i
            if pdb_coord_matching['nb_discrepancies'] == 0:
                break
    pdb_coord_matching['PDB_stop_on_uniprot'] = pdb_coord_matching['PDB_start_on_uniprot'] + len(short_seq)
    return pdb_coord_matching


def get_pdb_info(prot_dict, pdb_dir):
    '''Get informations about PDB sequences.

    :param dictionary prot_dict: the protein dictionary.
    :param str pdb_dir: the path of the pdb data directory.
    :return: the PDB data
    :rtype: dictionary
    '''

    # get the UniProt sequence
    uniprot_seq = ''.join([aa for aa in prot_dict['seq'].values()])
    # open pymol and retrieve the protein with PDB accession number
    pymol.finish_launching(['pymol', '-qc'])  # Pymol: quiet and no GUI
    # set the path to download the PDB data
    pymol.cmd.set('fetch_path', pymol.cmd.exp_path(pdb_dir), quiet=1)
    logging.info('Fetching PDB file with accession number: {}'.format(prot_dict['PDB']))
    pymol.cmd.fetch(prot_dict['PDB'])
    pymol.cmd.disable('all')
    pymol.cmd.enable(prot_dict['PDB'])

    # create the dictionary to retrieve the best match between PDB and uniprot seq
    pdb_data = {'PDB_start_on_uniprot': 0,
                'PDB_stop_on_uniprot': 0,
                'uniprot_longer_than_pdb': True,
                'nb_discrepancies': float('inf'),
                'chain': None}

    # iterate over the PDB chains
    for chain in pymol.cmd.get_chains(prot_dict['PDB']):
        if pdb_data['nb_discrepancies'] != 0:
            pdb_data['chain'] = chain
            # get the PDB sequence
            pymol.stored_list = []
            pymol.cmd.iterate('(name ca) and (chain {})'.format(chain),
                              'pymol.stored_list.append((resi, oneletter))')
            pdb_seq = ''.join([tuple_aa[1] for tuple_aa in pymol.stored_list])
            # check if pdb seq is longer than uniprot
            if len(pdb_seq) > len(uniprot_seq):
                pdb_data['uniprot_longer_than_pdb'] = False
            pdb_data = get_common_coordinates(pdb_seq, uniprot_seq, pdb_data)
            if pdb_data['nb_discrepancies'] == 0:
                break

    pymol.cmd.quit()

    # set the frames index of the positions to color
    frames_idx = [i + pdb_data['PDB_start_on_uniprot'] for i in range(pdb_data['PDB_stop_on_uniprot'] - pdb_data['PDB_start_on_uniprot'] + 1)]
    pdb_data['frames_idx'] = frames_idx

    return pdb_data
