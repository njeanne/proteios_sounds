#! /usr/bin/env python3

import os
import shutil
import time
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
        print('PDB shorter: pdb seq {} AA for uniprot seq {} AA\n'.format(len(pdb_seq),
                                                                          len(uniprot_seq)))
        short_seq = pdb_seq
        long_seq = uniprot_seq
    else:
        print('UNIPROT shorter: uniprot seq {} AA for pdb seq {} AA\n'.format(len(uniprot_seq),
                                                                              len(pdb_seq)))
    # look for the matching part between the two sequences
    for i in range(0, len(long_seq) - len(short_seq) + 1):
        ham_dist = hamming(short_seq, long_seq[i:i+len(short_seq)])
        if ham_dist < pdb_coord_matching['nb_diff']:
            pdb_coord_matching['nb_diff'] = ham_dist
            pdb_coord_matching['start_on_uniprot_without_signal_peptide'] = i
            pdb_coord_matching['pdb_longer_than_uniprot'] = pdb_longer
            if pdb_coord_matching['nb_diff'] == 0:
                break
    pdb_coord_matching['stop_on_uniprot_without_signal_peptide'] = pdb_coord_matching['start_on_uniprot_without_signal_peptide'] + len(short_seq)
    return pdb_coord_matching


def create_molecule_movie(prot_dict, durations, output_dir, logger):
    '''Create a movie of the common part between the UniProt and PDB sequences.

    :param dictionary prot_dict: the UniProt data dictionary.
    :param list durations: the list of the durations of the music notes.
    :param str output_dir: the path of the output directory.
    :param logger logger: the logger.
    :return: the movie path.
    :rtype:
    '''

    print(durations)
    print('length duration: {}\n\n'.format(len(durations)))

    # get the UniProt sequence
    uniprot_seq = ''.join([aa for aa in prot_dict['seq'].values()])
    print('UNI seq complete:\n{}\nlength: {}'.format(uniprot_seq, len(uniprot_seq)))
    # remove signal peptide
    if 'signal_peptide' in prot_dict:
        uniprot_seq = uniprot_seq[prot_dict['signal_peptide'][0][1]:]
    print('UNI seq no peptide signal:\n{}\nlength: {}\n\n'.format(uniprot_seq, len(uniprot_seq)))

    __main__.pymol_argv = ['pymol', '-qc'] # Pymol: quiet and no GUI
    pymol.finish_launching()
    logger.info('creating the {} movie'.format(prot_dict['PDB']))
    # set the path to download the PDB data
    pymol.cmd.set('fetch_path', pymol.cmd.exp_path(output_dir), quiet=1)
    print('Fetching PDB accession number: {}'.format(prot_dict['PDB']))
    pymol.cmd.fetch(prot_dict['PDB'])
    pymol.cmd.disable('all')
    pymol.cmd.enable(prot_dict['PDB'])

    # pymol.cmd.iterate('(name ca)', 'pymol.stored_list.append((resi, resn))')

    print('Chains: {}\n'.format(pymol.cmd.get_chains(prot_dict['PDB'])))
    # create the dictionary to retrieve the best match between PDB and uniprot seq
    pdb_data = {'start_on_uniprot_without_signal_peptide': 0,
                'stop_on_uniprot_without_signal_peptide': 0,
                'pdb_longer_than_uniprot': True,
                'uniprot_length': len(uniprot_seq),
                'nb_diff': float('inf'),
                'chain': None}
    # iterate over the PDB chains
    for chain in pymol.cmd.get_chains(prot_dict['PDB']):
        if pdb_data['nb_diff'] != 0:
            pdb_data['chain'] = chain
            # get the PDB sequence
            pymol.stored_list = []
            pymol.cmd.iterate('(name ca) and (chain {})'.format(chain),
                              'pymol.stored_list.append((resi, oneletter))')
            pdb_seq = ''.join([tuple_aa[1] for tuple_aa in pymol.stored_list])
            print('PDB seq chain {}: {}'.format(chain, pdb_seq))
            pdb_data = get_common_coordinates(pdb_seq, uniprot_seq, pdb_data)
            if pdb_data['nb_diff'] == 0:
                break
    pdb_data['pymol_stored_list'] = pymol.stored_list
    # pymol.cmd.quit()

    print('\n\n\nPDB data: {}'.format(pdb_data))




    ## TODO: Faire un multithreading ou multi processing (voir https://realpython.com/python-concurrency/)
    ## pour diminuer le temps de création des images. Gérer le fait que chaque processus ne produit que 2
    ## images afin de ne pas avoir le bandeau licence expirée.

    pymol.cmd.hide('all')
    pymol.cmd.show('cartoon')
    pymol.cmd.set('ray_opaque_background', 1)

    # create the frames
    frames_dir = os.path.join(output_dir, 'frames')
    if os.path.exists(frames_dir):
        shutil.rmtree(frames_dir)
    os.mkdir(frames_dir)

    # no AA colored frames for signal peptide
    idx_frame = 0
    img_path_0 = os.path.join(frames_dir,
                              '{}_{}_{}.png'.format(prot_dict['accession_number'],
                                                    prot_dict['PDB'],
                                                    idx_frame))
    print('Frame {} [signal peptide first frame]: {}'.format(idx_frame,
                                                             img_path_0))
    pymol.cmd.png(img_path_0, quiet=1)
    # wait for pymol frame creation
    while not os.path.exists(img_path_0):
        time.sleep(1)
    # copy the frame to represent the signal peptide
    for i in range(1, prot_dict['signal_peptide'][0][1]):
        idx_frame += 1
        img_path = os.path.join(frames_dir,
                                '{}_{}_{}.png'.format(prot_dict['accession_number'],
                                                      prot_dict['PDB'],
                                                      idx_frame))
        print('Frame {} [signal peptide]: {}'.format(idx_frame,
                                                     img_path))
        shutil.copyfile(img_path_0, img_path)

    # copy the frame to represent the AA not present on pdb
    if pdb_data['start_on_uniprot_without_signal_peptide'] > -1:
        for i in range(pdb_data['start_on_uniprot_without_signal_peptide']):
            idx_frame += 1
            img_path = os.path.join(frames_dir,
                                    '{}_{}_{}.png'.format(prot_dict['accession_number'],
                                                          prot_dict['PDB'],
                                                          idx_frame))
            print('Frame {} [missing in PDB]: {}'.format(idx_frame,
                                                         img_path))
            shutil.copyfile(img_path_0, img_path)

    # color the amino acids
    for i in range(pdb_data['stop_on_uniprot_without_signal_peptide'] - pdb_data['start_on_uniprot_without_signal_peptide']):
        idx_frame += 1
        # add 1 because pdb index is starts on 1
        i = i + 1
        pymol.cmd.color('red', 'resi {}'.format(i))
        if not i == 1:
            pymol.cmd.color('green', 'resi {}'.format(i-1))
        img_path = os.path.join(frames_dir,
                                '{}_{}_{}.png'.format(prot_dict['accession_number'],
                                                      prot_dict['PDB'],
                                                      idx_frame))
        print('Frame {} [in PDB]: {}'.format(idx_frame,
                                             img_path))
        pymol.cmd.png(img_path, quiet=1)

    # no AA colored frames for uniprot AA outside pdb and after pdb AA
    #TODO: tester si pdb va moins loin que uniprot, si oui copier les frames
    # en incrémentant les numéros

    pymol.cmd.quit()
