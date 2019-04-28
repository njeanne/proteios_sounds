#! /usr/bin/env python3

import os
import pymol
import __main__

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

    short_seq = uniprot_seq
    long_seq = pdb_seq
    pdb_longer = True
    if len(uniprot_seq) > len(pdb_seq):
        pdb_longer = False
        print('PDB shorter: pdb seq {} AA for uniprot seq {} AA'.format(len(pdb_seq),
                                                                        len(uniprot_seq)))
        short_seq = pdb_seq
        long_seq = uniprot_seq
    else:
        print('UNIPROT shorter: uniprot seq {} AA for pdb seq {} AA'.format(len(uniprot_seq),
                                                                            len(pdb_seq)))
    # look for the matching part between the two sequences
    for i in range(0, len(long_seq) - len(short_seq) + 1):
        ham_dist = hamming(short_seq, long_seq[i:i+len(short_seq)])
        if ham_dist < pdb_coord_matching['nb_diff']:
            pdb_coord_matching['nb_diff'] = ham_dist
            pdb_coord_matching['start'] = i
            pdb_coord_matching['pdb_longer_than_uniprot'] = pdb_longer
            if pdb_coord_matching['nb_diff'] == 0:
                break
    pdb_coord_matching['stop'] = pdb_coord_matching['start'] + len(short_seq)
    return pdb_coord_matching


def create_molecule_movie(prot_dict, durations, output_dir, logger):
    '''Create a movie of the common part between the UniProt and PDB sequences.

    :param dictionary prot_dict: the UniProt data dictionary.
    :param list durations: the list of the durations of the music notes.
    :param str output_dir: the path of the output directory.
    :param logger logger: the logger.
    :return: the movie.
    :rtype:
    '''

    print(durations)
    print('length duration: {}'.format(len(durations)))
    __main__.pymol_argv = ['pymol', '-qc'] # Pymol: quiet and no GUI
    pymol.finish_launching()
    logger.info('creating the {} movie'.format(prot_dict['PDB']))
    # set the path to download the PDB data
    pymol.cmd.set('fetch_path', pymol.cmd.exp_path(output_dir), quiet=1)
    print('Fetching PDB accession number: {}'.format(prot_dict['PDB']))
    pymol.cmd.fetch(prot_dict['PDB'])
    img_path = os.path.join(output_dir, '{}_{}.png'.format(prot_dict['accession_number'], prot_dict['PDB']))
    pymol.cmd.disable('all')
    pymol.cmd.enable(prot_dict['PDB'])

    # get the UniProt sequence
    uniprot_seq = ''.join([aa for aa in prot_dict['seq'].values()])
    print('UNI seq complete:\n{}\nlength: {}'.format(uniprot_seq, len(uniprot_seq)))
    # remove signal peptide
    if 'signal_peptide' in prot_dict:
        uniprot_seq = uniprot_seq[prot_dict['signal_peptide'][0][1]:]
    print('UNI seq no peptide signal:\n{}\nlength: {}'.format(uniprot_seq, len(uniprot_seq)))

    # get the PDB sequence
    pymol.stored_list = []
    pymol.cmd.iterate('(name ca)', 'pymol.stored_list.append((resi, resn))')

    print('Chains: {}'.format(pymol.cmd.get_chains(prot_dict['PDB'])))
    # create the dictionary to retrieve the best match between PDB and uniprot seq
    pdb_coord_data = {'start': 0,
                      'stop': 0,
                      'pdb_longer_than_uniprot': True,
                      'nb_diff': float('inf'),
                      'pdb_chain': None}
    # iterate over the PDB chains
    for chain in pymol.cmd.get_chains(prot_dict['PDB']):
        if pdb_coord_data['nb_diff'] != 0:
            pdb_coord_data['pdb_chain'] = chain
            pymol.cmd.iterate('(name ca) and (chain {})'.format(chain),
                              'pymol.stored_list.append((resi, oneletter))')
            pdb_seq = ''.join([tuple_aa[1] for tuple_aa in pymol.stored_list])
            print('PDB seq chain {}: {}'.format(chain, pdb_seq))
            pdb_coord_data = get_common_coordinates(pdb_seq, uniprot_seq, pdb_coord_data)
            if pdb_coord_data['nb_diff'] == 0:
                break
    print(pdb_coord_data)


    pymol.cmd.hide('all')
    pymol.cmd.show('cartoon')
    pymol.cmd.set('ray_opaque_background', 1)
    pymol.cmd.color('red', 'ss h')
    pymol.cmd.color('yellow', 'ss s')
    pymol.cmd.png(img_path, quiet=1)
    pymol.cmd.quit()

# if __name__ == '__main__':
#     seqA = 'ABCDEF'
#     seqB = 'CEE'

    # less_seq_dif(seqA, seqB)
