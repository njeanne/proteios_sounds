#! /usr/bin/env python3

import os
import time
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
        print('UNIPROT longer than PDB: uniprot seq {} AA for pdb seq {} AA\n'.format(len(uniprot_seq),
                                                                                      len(pdb_seq)))
    else:
        short_seq = uniprot_seq
        long_seq = pdb_seq
        print('UNIPROT shorter than PDB: uniprot seq {} AA for pdb seq {} AA\n'.format(len(uniprot_seq),
                                                                                       len(pdb_seq)))
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

def get_pdb_info(prot_dict, pdb_dir, logger):
    '''Get informations about PDB sequences.

    :param dictionary prot_dict: the protein dictionary.
    :param str pdb_dir: the path of the pdb data directory.
    :param logger logger: the logger.
    :return: the PDB data
    :rtype: dictionary
    '''

    # get the UniProt sequence
    uniprot_seq = ''.join([aa for aa in prot_dict['seq'].values()])
    print('UNI seq complete:\n{}\nlength: {}'.format(uniprot_seq, len(uniprot_seq)))

    # open pymol and retrieve the protein with PDB accession number
    pymol.finish_launching(['pymol', '-qc']) # Pymol: quiet and no GUI
    logger.info('creating the {} movie'.format(prot_dict['PDB']))
    # set the path to download the PDB data
    pymol.cmd.set('fetch_path', pymol.cmd.exp_path(pdb_dir), quiet=1)
    print('Fetching PDB accession number: {}'.format(prot_dict['PDB']))
    pymol.cmd.fetch(prot_dict['PDB'])
    pymol.cmd.disable('all')
    pymol.cmd.enable(prot_dict['PDB'])

    print('Chains: {}\n'.format(pymol.cmd.get_chains(prot_dict['PDB'])))
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
            print('PDB seq chain {}: {}'.format(chain, pdb_seq))
            pdb_data = get_common_coordinates(pdb_seq, uniprot_seq, pdb_data)
            if pdb_data['nb_discrepancies'] == 0:
                break

    print('\n\n\nPDB data: {}'.format(pdb_data))

    pymol.cmd.hide('all')
    pymol.cmd.show('cartoon')
    pymol.cmd.set('ray_opaque_background', 1)

    # no AA colored frames for signal peptide
    if pdb_data['PDB_start_on_uniprot'] > 0:
        img_path = os.path.join(pdb_dir,
                                'frames',
                                '{}_no-idx.png'.format(prot_dict['PDB']))
        print('Frame without index [missing in PDB]: {}'.format(img_path))
        pymol.cmd.png(img_path, quiet=0)

    # wait for pymol frame creation
    while not os.path.exists(img_path):
        time.sleep(0.5)

    pymol.cmd.quit()


    # color the amino acids and save the frames
    frames_idx = [i + pdb_data['PDB_start_on_uniprot'] for i in range(pdb_data['PDB_stop_on_uniprot'] - pdb_data['PDB_start_on_uniprot'] + 1)]

    pdb_data['frames_idx'] = frames_idx
    # with open(yaml_path, 'w') as outfile:
    #     yaml.dump(pdb_data, outfile, default_flow_style=True)

    return pdb_data
