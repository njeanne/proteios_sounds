#! /usr/bin/env python3

__author__ = 'Nicolas JEANNE'
__copyright__ = 'GNU General Public License'
__version__ = '1.0.0'
__email__ = 'nicolas.jeanne@ntymail.com'

import argparse
import os
import pymol


if __name__ == '__main__':
    descr = '''
    {} v.{}

    Created by {}.
    Contact: {}
    {}

    Creates a frames from PDB data.
    '''.format(os.path.basename(__file__), __version__, __author__, __email__, __copyright__)

    # Parse arguments
    parser = argparse.ArgumentParser(description=descr, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-p', '--pdb_dir', required=True, help='path to the pdb data folder.')
    parser.add_argument('-c', '--chain', required=True, help='the chain id.')
    parser.add_argument('-n', '--frame_nb', required=True, help='the frame number')
    parser.add_argument('-i', '--idx_aa', required=True, help='the AA idx in PDB')
    parser.add_argument('pdb_AN', help='the PDB accession number')
    args = parser.parse_args()

    # open pymol and retrieve the protein with PDB accession number
    pymol.finish_launching(['pymol', '-qc']) # Pymol: quiet and no GUI
    print('THREADING creating the {} movie'.format(args.pdb_AN))
    # set the path to download the PDB data
    # pymol.cmd.set('fetch_path', pymol.cmd.exp_path(fn_arg['pdb_dir']), quiet=1)
    # print('Fetching PDB accession number: {}'.format(fn_arg['pdb_AN']))
    pymol.cmd.load(os.path.join(args.pdb_dir, '{}.cif'.format(args.pdb_AN.lower())))
    pymol.cmd.disable('all')
    pymol.cmd.enable(args.pdb_AN)
    pymol.stored_list = []
    pymol.cmd.iterate('(name ca) and (chain {})'.format(args.chain),
                      'pymol.stored_list.append((resi, oneletter))')
    pymol.cmd.hide('all')
    pymol.cmd.show('cartoon')
    pymol.cmd.set('ray_opaque_background', 1)

    print(args.idx_aa)
    pymol.cmd.color('red', 'resi {}'.format(args.idx_aa))
    img_path = os.path.join(args.pdb_dir,
                            'frames',
                            '{}_{}.png'.format(args.pdb_AN, args.frame_nb))
    print('Frame {} [in PDB]: {}'.format(args.frame_nb,
                                         img_path))
    pymol.cmd.png(img_path, quiet=1)
    pymol.cmd.quit()
