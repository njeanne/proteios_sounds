#! /usr/bin/env python3

__author__ = 'Nina Verstraete, Jacques TOEN & Nicolas JEANNE'
__copyright__ = 'GNU General Public License'
__version__ = '1.0.0'
__email__ = 'nicolas.jeanne@ntymail.com'

import argparse
import sys
import os
import logging
import shutil
import subprocess
import concurrent.futures
import threading
import parse_uniprot
import midi_operations
import parse_pdb


def restricted_tempo(tempo_value):
    '''
    Check range for tempo argument, must be between 60 and 150.
    param str tempo_value: value of the tempo argument in BPM.
    return: the tempo.
    rtype: int
    '''
    tempo_value = int(tempo_value)
    if tempo_value < 60 or tempo_value > 250:
        raise argparse.ArgumentTypeError('{} not in range 60 to 250.'.format(tempo_value))
    return tempo_value

if __name__ == '__main__':
    prg_id = os.path.splitext(os.path.basename(__file__))[0]
    descr = '''
    {} v.{}

    Created by {}.
    Contact: {}
    {}

    Create a MIDI file from a protein entry of the UniProt database (https://www.uniprot.org/).
    '''.format(prg_id, __version__, __author__, __email__, __copyright__)

    # Parse arguments
    parser = argparse.ArgumentParser(description=descr,
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-o', '--out', required=True, help='path to the results directory.')
    parser.add_argument('-s', '--score', required=False, action='store_true',
                        help='''use musescore software to create the score
                        corresponding to the MIDI file.''')
    parser.add_argument('-p', '--play', required=False, action='store_true',
                        help='play the music with Timidity, just for tests.')
    parser.add_argument('-t', '--tempo', required=False, type=restricted_tempo,
                        help='set the tempo in BPM. Value between 60 and 250.')
    parser.add_argument('-i', '--instruments', required=False, nargs=3,
                        help='''set channel 0, 1 and 2 instruments,
                        restricted to 3 values between 0 and 127
                        separated by spaces. Default is 0:  Acoustic Grand,
                        42: Cello and 65: Alto Sax.
                        See: http://www.pjb.com.au/muscript/gm.html#patch for details.''')
    parser.add_argument('-d', '--debug', required=False, action='store_true',
                        help='''debug mode, create a log file which details each
                         entry of the MIDI file.''')
    parser.add_argument('uniprot_accession_number',
                        help='''the protein accession number in the UniProt
                        database. Example: Human Interleukin-8 > P10145''')
    args = parser.parse_args()


    # check if instruments are between 0 and 127
    if args.instruments:
        for i in range(len(args.instruments)):
            instru = int(args.instruments[i])
            if instru < 0 or instru > 127:
                raise argparse.ArgumentTypeError('{} should be 3 integers between 0 and 127.'.format(args.instruments))
            else:
                args.instruments[i] = instru
        instrus = args.instruments
    else:
        instrus = [0, 42, 65]

    # tempo
    if args.tempo:
        tempo = int(args.tempo)
    else:
        tempo = 100  # In BPM

    ### midi notes on major mode correspondance with AA sorted by decreasing molecular weight
    # keys are set as DO (48, 60, 72) degrees I, SOL (55, 67) degrees V, FA (53, 65) degrees IV, RE (50, 62) degrees II,
    # MI (52, 64) degrees III, LA (57, 69) degrees VI and SI (59, 71) degrees VII. Finally, we add 7 alterations #
    # following the ascending quint (54, 66, 49, 61, 56, 68, 51)
    initial_midi_keys = [48, 60, 72, 55, 67, 53, 65, 50, 62, 52, 64, 57, 69, 59, 71, 54, 66, 49, 61, 56, 68, 51]
    midi_keys = {}

    ### Physico-chemical properties of AA
    AA_PHY_CHI = {'A': {'hybrophobic', 'small'},
                  'R': {'polar', 'pos_charged'},
                  'N': {'polar', 'small'},
                  'D': {'polar', 'small', 'neg_charged'},
                  'C': {'hydrophobic', 'polar', 'small'},
                  'E': {'polar', 'neg_charged'},
                  'Q': {'polar'},
                  'G': {'hydrophobic', 'small'},
                  'H': {'hydrophobic', 'polar', 'pos_charged', 'aromatic'},
                  'I': {'hydrophobic', 'aliphatic'},
                  'L': {'hydrophobic', 'aliphatic'},
                  'K': {'hydrophobic', 'polar', 'pos_charged'},
                  'M': {'hydrophobic'},
                  'F': {'hydrophobic', 'aromatic'},
                  'P': {'small'},
                  'S': {'polar', 'small'},
                  'T': {'hydrophobic', 'polar', 'small'},
                  'W': {'hydrophobic', 'polar', 'aromatic'},
                  'Y': {'hydrophobic', 'polar', 'aromatic'},
                  'V': {'hydrophobic', 'small', 'aliphatic'}}

    # create the output directory
    out_dir = os.path.abspath(args.out)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # create the log file
    log_path = os.path.join(out_dir, '{}.log'.format(prg_id))
    if os.path.exists(log_path):
        os.remove(log_path)
    logging.basicConfig(filename=log_path,
                        level=logging.DEBUG,
                        format='%(asctime)s\t%(levelname)s:\t%(message)s',
                        datefmt='%Y/%m/%d %H:%M:%S')
    logger = logging.getLogger(__name__)
    logger.info(' '.join(sys.argv))

    logger.info('Output directory: {}'.format(out_dir))
    logger.info('Tempo: {} BPM'.format(tempo))
    logger.info('Instruments: {} (general MIDI patch numbers, see: http://www.pjb.com.au/muscript/gm.html#patch)'.format(', '.join(map(str, instrus))))
    logger.info('Create score: {}'.format(args.score))

    # parsing of uniprot entry
    protein = parse_uniprot.parse_entry(args.uniprot_accession_number, logger)
    logger.info('UniProt accession number: {}'.format(args.uniprot_accession_number))
    logger.info('Protein: {}'.format(protein['entry_name']))
    logger.info('Organism: {}'.format(protein['organism']))
    if 'PDB' in protein:
        logger.info('PDB: {} (Protein DataBase accession number)'.format(protein['PDB']))
    else:
        logger.info('PDB: No accession number in Uniprot entry')

    sequence = protein['seq']
    sequence_length = len(sequence)
    protein['seq'] = {}
    if args.debug:
        logger.info('AA sequence ({} AA): {}'.format(sequence_length, sequence))
    for i in range(sequence_length):
        protein['seq'][i] = sequence[i]
    # frequence of AA in the sequence
    set_AA = set(''.join(sequence))
    proportion_AA = {}
    for aa in set_AA:
        proportion_AA[aa] = sequence.count(aa) / sequence_length
    # sort by decreasing frequency
    proportion_AA = sorted(proportion_AA.items(), key=lambda kv: kv[1], reverse=True)

    for idx, aa_proportion in enumerate(proportion_AA):
        midi_keys[aa_proportion[0]] = initial_midi_keys[idx]

    # create the MIDI file
    midi_file_path, keys_duration = midi_operations.create_midi(args.uniprot_accession_number,
                                                                protein,
                                                                midi_keys,
                                                                tempo,
                                                                instrus,
                                                                out_dir,
                                                                AA_PHY_CHI,
                                                                logger,
                                                                args.debug)
    print('MIDI file for {} {} ({}) created: {}'.format(protein['entry_name'],
                                                        protein['organism'],
                                                        args.uniprot_accession_number,
                                                        midi_file_path))

    if 'PDB' in protein:
        # create the directories for PDB data and frames
        pdb_dir = os.path.join(args.out, 'pdb', '{}_{}'.format(protein['accession_number'], protein['PDB']))
        if not os.path.exists(pdb_dir):
            os.makedirs(pdb_dir)
        frames_dir = os.path.join(pdb_dir, 'frames')
        if os.path.exists(frames_dir):
            shutil.rmtree(frames_dir)
        os.mkdir(frames_dir)

        # get data from the PDB file
        pdb_data = parse_pdb.get_pdb_info(protein,
                                          pdb_dir,
                                          logger)


        # pdb_data = None
        # with open(yaml_path, 'r') as stream:
        #     try:
        #         pdb_data = yaml.safe_load(stream)
        #         print(pdb_data)
        #         # change dictionary to list of tuples to retrieve pymol format
        #         tuples_list_pdb = []
        #         for item in sorted(pdb_data['pymol_stored_list']):
        #             tuples_list_pdb.append((str(item), pdb_data['pymol_stored_list'][item]))
        #         pdb_data['pymol_stored_list'] = tuples_list_pdb
        #         print(pdb_data)
        #     except yaml.YAMLError as yaml_exc:
        #         print(yaml_exc)


        args_list = []
        for i in range(10):
            cmd = './test.py -x {} -y 2 '.format(i)
            args_list.append(cmd)
        subprocess.run(cmd, shell=True)
        with concurrent.futures.ThreadPoolExecutor(max_workers=10) as executor:
            print('executor')
            executor.map(subprocess.run, args_list)

        # args_list = []
        # for aa_idx, frame_idx in pdb_data['frames_idx'].enumerate():
        #     cmd = './create_pdb_frames.py -p {} -c {} -n {} -i {} {}'.format(os.path.abspath(pdb_dir),
        #                                                                   pdb_data['chain'],
        #                                                                   frame_idx,
        #                                                                   aa_idx + 1,
        #                                                                   protein['PDB'])
        #     # subprocess.call(cmd, shell=True)
        #     args_list.append((cmd, eval('shell=True')))
        #     # args_list.append(cmd)
        #
        # for i in args_list:
        #     print(i)
        #
        # with concurrent.futures.ThreadPoolExecutor(max_workers=10) as executor:
        #     print('executor')
        #     executor.map(subprocess.call, args_list)
        #
        # args_list1 = [(1, 2, 'deux'), (2, 2, 'deux'), (3, 2, 'deux'), (4, 2, 'deux'), (5, 2, 'deux'), (6, 2, 'deux')]
        # def toto(argument):
        #     print('{} pour {}'.format(argument[0]*argument[1], argument[2]))
        # print('\n\n\ntoto')
        # with concurrent.futures.ThreadPoolExecutor(max_workers=5) as executor:
        #     executor.map(toto, args_list1)
        #
        # args_list2 = [1,2,3,4,5,6,7,8,9]
        # def toto2(idx):
        #     print('{} pour {}'.format(idx*2, idx))
        # print('\n\n\ntoto2')
        # with concurrent.futures.ThreadPoolExecutor(max_workers=5) as executor:
        #     executor.map(toto2, args_list2)
        #
        # args_list3 = [{'idx': 1, 'nb': 2, 'str': 'deux'},
        #               {'idx': 2, 'nb': 2, 'str': 'deux'},
        #               {'idx': 3, 'nb': 2, 'str': 'deux'},
        #               {'idx': 4, 'nb': 2, 'str': 'deux'},
        #               {'idx': 5, 'nb': 2, 'str': 'deux'},
        #               {'idx': 6, 'nb': 2, 'str': 'deux'}]
        # def toto3(argument):
        #     print('{} pour {}'.format(argument['idx']*argument['nb'], argument['str']))
        # print('\n\n\ntoto3')
        # with concurrent.futures.ThreadPoolExecutor(max_workers=5) as executor:
        #     executor.map(toto3, args_list3)



    # get the score
    if args.score:
        print('Creating the score:')
        score_basename = '{}_{}_{}_{}bpm_score.pdf'.format(args.uniprot_accession_number,
                                                           protein['entry_name'],
                                                           protein['organism'],
                                                           tempo)
        score_output = os.path.join(args.out, score_basename)
        cmd = 'mscore -o {} {}'.format(score_output, midi_file_path)
        subprocess.run(cmd, shell=True)
        print('Score created at {}'.format(score_output))

    # play the file with timidity if asked
    if args.play:
        cmd = 'timidity {}'.format(midi_file_path)
        subprocess.run(cmd, shell=True)
