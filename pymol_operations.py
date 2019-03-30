#! /usr/bin/env python3

import os
import pymol
import __main__

def hamming(s1, s2):
    '''Calculate the Hamming distance between two strings'''
    assert len(s1) == len(s2)
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))

def less_difference_seq(pdb_seq, uniprot_seq):
    short_seq = pdb_seq
    long_seq = uniprot_seq
    if len(pdb_seq) > len(uniprot_seq):
        print('short uniprot: {}'.format(len(pdb_seq)))
        short_seq = uniprot_seq
        long_seq = pdb_seq
    else:
        print('short pdb')
    min_diff = len(short_seq)
    idx_start = 0
    for i in range(0, len(long_seq) - len(short_seq) + 1):
        ham_dist = hamming(short_seq, long_seq[i:i+len(short_seq)])
        if ham_dist < min_diff:
            min_diff = ham_dist
            idx_start = i
            if min_diff == 0:
                break
    return idx_start


def create_molecule_movie(prot_dict, durations, output_dir, logger):
    print(durations)
    print('length duration: {}'.format(len(durations)))
    __main__.pymol_argv = ['pymol', '-qc'] # Pymol: quiet and no GUI
    pymol.finish_launching()
    logger.info('creating the {} movie'.format(prot_dict['PDB']))
    # set the path to download the PDB data
    pymol.cmd.set('fetch_path', pymol.cmd.exp_path(output_dir), quiet=1)
    pymol.cmd.fetch(prot_dict['PDB'])
    img_path = os.path.join(output_dir, '{}_{}.png'.format(prot_dict['accession_number'], prot_dict['PDB']))
    pymol.cmd.disable('all')
    pymol.cmd.enable(prot_dict['PDB'])

    # get the PDB sequence
    pymol.stored_list = []
    if len(pymol.cmd.get_chains(prot_dict['PDB'])) > 1:
        # more than one chain, select the Chain A
        pymol.cmd.iterate('(name ca) and (chain A)', 'pymol.stored.list.append((resi, oneletter))')
    else:
        pymol.cmd.iterate('(name ca)', 'pymol.stored.list.append((resi, oneletter))')
    pdb_seq = ''.join([tuple_aa[1] for tuple_aa in pymol.stored.list])
    print('PDB seq: {}'.format(pdb_seq))

    uniprot_seq = ''.join([aa for aa in prot_dict['seq'].values()])
    print('UNI seq: {}'.format(uniprot_seq))

    # # remove signal peptide
    # if 'signal_peptide' in prot_dict:
    #     uniprot_seq = uniprot_seq[prot_dict['signal_peptide'][0][1]:]
    # print('UNI seq: {}'.format(uniprot_seq))
    idx_start = less_difference_seq(pdb_seq, uniprot_seq)
    idx_end = idx_start + len(pdb_seq)


    pymol.cmd.hide('all')
    pymol.cmd.show('cartoon')
    pymol.cmd.set('ray_opaque_background', 1)
    # pymol.cmd.color('red', 'ss h')
    # pymol.cmd.color('yellow', 'ss s')
    pymol.cmd.png(img_path, quiet=1)
    pymol.cmd.quit()

if __name__ == '__main__':
    seqA = 'ABCDEF'
    seqB = 'CEE'

    less_seq_dif(seqA, seqB)
