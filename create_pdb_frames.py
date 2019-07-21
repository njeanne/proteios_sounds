#! /usr/bin/env python3

__author__ = 'Nicolas JEANNE'
__copyright__ = 'GNU General Public License'
__version__ = '1.0.0'
__email__ = 'nicolas.jeanne@ntymail.com'

import argparse
import os
import sys
# add pymol to the python path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'lib/python3.7/site-packages'))
import pymol


if __name__ == '__main__':
    descr = '''
    {} v.{}

    Created by {}.
    Contact: {}
    {}

    Creates a frames from PDB data.
    '''.format(os.path.basename(__file__), __version__, __author__, __email__,
               __copyright__)

    # Parse arguments
    parser = argparse.ArgumentParser(description=descr, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-p', '--pdb_dir', required=True, help='path to the pdb data folder.')
    parser.add_argument('-c', '--chain', required=True, help='the chain id.')
    parser.add_argument('-n', '--frame_nb', required=True, help='the frame number')
    parser.add_argument('-i', '--idx_aa', required=True, help='the AA idx in PDB')
    parser.add_argument('--color_aa', required=False, action='store_true', help='color in red the AA specified by "--idx_aa".')
    parser.add_argument('pdb_AN', help='the PDB accession number')
    args = parser.parse_args()

    # open pymol and retrieve the protein with PDB accession number
    pymol.finish_launching(['pymol', '-qc'])  # Pymol: quiet and no GUI
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

    if args.color_aa:
        pymol.cmd.color('red', 'resi {}'.format(args.idx_aa))
    img_path = os.path.join(args.pdb_dir,
                            'frames',
                            '{}_{}.png'.format(args.pdb_AN, args.frame_nb))
    print('[Pymol] Frame {} (in PDB file): {}'.format(args.frame_nb,
                                                      img_path))
    pymol.cmd.png(img_path, width=800, height=600, quiet=1)
    pymol.cmd.quit()
