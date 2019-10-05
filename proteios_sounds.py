#! /usr/bin/env python3

__author__ = 'Nina VERSTRAETE, Jacques TOEN & Nicolas JEANNE'
__copyright__ = 'GNU General Public License'
__version__ = '1.0.0'
__email__ = 'nicolas.jeanne@ntymail.com'

import argparse
import sys
import os
import logging
import subprocess
import concurrent.futures
import time
import parse_uniprot
import midi_operations
import parse_pdb
import audio_video


def restricted_tempo(tempo_value):
    '''
    Check range for tempo argument, must be between 60 and 150.
    param str tempo_value: value of the tempo argument in BPM.
    return: the tempo.
    rtype: int
    '''
    tempo_value = int(tempo_value)
    if tempo_value < 60 or tempo_value > 250:
        raise argparse.ArgumentTypeError(('{} not in range '
                                          '60 to 250.').format(tempo_value))
    return tempo_value


if __name__ == '__main__':
    prg_id = os.path.splitext(os.path.basename(__file__))[0]
    descr = '''
    {} v.{}

    Created by {}.
    Contact: {}
    {}

    Create a MIDI file and from a protein entry of the UniProt database
    (https://www.uniprot.org/).
    If the data are available in the UniProt entry, a movie file of the 3D
    representation of the protein will also be created.
    '''.format(prg_id, __version__, __author__, __email__, __copyright__)

    # Parse arguments
    parser = argparse.ArgumentParser(description=descr,
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-o', '--out',
                        required=True,
                        help='path to the results directory.')
    parser.add_argument('-s', '--score', required=False, action='store_true',
                        help=('use musescore software to create the score '
                              'corresponding to the MIDI file.'))
    parser.add_argument('-t', '--tempo', required=False, type=restricted_tempo,
                        help='set the tempo in BPM. Value between 60 and 250.')
    parser.add_argument('-i', '--instruments', required=False, nargs=3,
                        help=('set channel 0, 1 and 2 instruments, restricted '
                              'to 3 values between 0 and 127 separated by '
                              'spaces. Default is 0:  Acoustic Grand, '
                              '42: Cello and 65: Alto Sax. See: '
                              'http://www.pjb.com.au/muscript/gm.html#patch '
                              'for details.'))
    parser.add_argument('--log_level', required=False,
                        choices=['DEBUG', 'INFO', 'WARNING',
                                 'ERROR', 'CRITICAL'], default='INFO',
                        help=('set the logging level, if not set, '
                              'default is INFO.'))
    parser.add_argument('uniprot_AN',
                        help=('the protein Accession Number in the UniProt '
                              'database. Example: Human Interleukin-8 > P10145'))
    args = parser.parse_args()

    # check if instruments are between 0 and 127
    if args.instruments:
        for i in range(len(args.instruments)):
            instru = int(args.instruments[i])
            if instru < 0 or instru > 127:
                raise argparse.ArgumentTypeError(('{} should be 3 integers between '
                                                  '0 and 127.').format(args.instruments))
            args.instruments[i] = instru
        instrus = args.instruments
    else:
        instrus = [0, 42, 65]

    # tempo
    if args.tempo:
        tempo = int(args.tempo)
    else:
        tempo = 100  # In BPM

    # midi notes on major mode correspondance with AA sorted by decreasing
    # molecular weight keys are set as DO (48, 60, 72) degrees I,
    # SOL (55, 67) degrees V, FA (53, 65) degrees IV, RE (50, 62) degrees II,
    # MI (52, 64) degrees III, LA (57, 69) degrees VI and
    # SI (59, 71) degrees VII. Finally, we add 7 alterations '#' following the
    # ascending quint (54, 66, 49, 61, 56, 68, 51)
    initial_midi_keys = [48, 60, 72, 55, 67, 53, 65, 50, 62, 52, 64, 57, 69,
                         59, 71, 54, 66, 49, 61, 56, 68, 51]
    midi_keys = {}

    # Physico-chemical properties of AA
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
                        level=args.log_level,
                        format='%(asctime)s\t%(levelname)s:\t%(message)s',
                        datefmt='%Y/%m/%d %H:%M:%S')
    logger = logging.getLogger(__name__)

    logger.info(' '.join(sys.argv))
    logger.info('Output directory: {}'.format(out_dir))
    logger.info('Tempo: {} BPM'.format(tempo))
    logger.info(('Instruments: {} (general MIDI patch numbers, see: '
                 'http://www.pjb.com.au/muscript/gm.html#patch)').format(', '.join(map(str,
                                                                                       instrus))))
    logger.info('Create score: {}'.format(args.score))

    # parsing of uniprot entry
    protein = parse_uniprot.parse_entry(args.uniprot_AN, initial_midi_keys)

    # set the result files base name
    file_base_name = '{}_{}_{}_{}bpm_instrus'.format(args.uniprot_AN,
                                                     protein['entry_name'],
                                                     protein['organism'],
                                                     tempo)
    for instru in instrus:
        file_base_name = '{}-{}'.format(file_base_name, instru)

    # create the MIDI file
    midi_file_path = os.path.join(out_dir, '{}.midi'.format(file_base_name))
    keys_duration = midi_operations.create_midi(midi_file_path, protein, tempo,
                                                instrus, AA_PHY_CHI)
    msg_template = 'file for {} {} ({}) created:'.format(protein['entry_name'],
                                                         protein['organism'],
                                                         args.uniprot_AN)
    msg = 'MIDI {} {}'.format(msg_template, midi_file_path)
    logger.info(msg)
    print(msg)

    # convert MIDI to flac
    flac_file_path = audio_video.convert_midi_to_flac(midi_file_path)
    msg = 'FLAC {} {}'.format(msg_template, flac_file_path)
    logger.info(msg)
    print(msg)

    if 'PDB' in protein:
        # create the directories for PDB data and frames
        pdb_dir = os.path.join(os.path.abspath(args.out),
                               'pdb',
                               '{}_{}'.format(protein['accession_number'],
                                              protein['PDB']))
        frames_dir = os.path.join(pdb_dir, 'frames')
        if not os.path.exists(frames_dir):
            os.makedirs(frames_dir)

        # get data from the PDB file
        pdb_data = parse_pdb.get_pdb_info(protein, pdb_dir)

        # create a frame without colored AA for all AA outside the PDB data
        existing_frames = sorted([png for png in os.listdir(frames_dir)])
        if '{}_no-idx.png'.format(protein['PDB']) not in existing_frames:
            print(('\nCreating {} ({}) protein frame, '
                   'please wait..').format(protein['entry_name'],
                                           protein['PDB']))
            logger.info('Creating {} ({}) protein frame.'.format(protein['entry_name'],
                                                                 protein['PDB']))
            cmd_no_color = ('./create_pdb_frames.py -p {} -c {} '
                            '-n {} -i {} {}').format(pdb_dir,
                                                     pdb_data['chain'],
                                                     'no-idx',
                                                     1,
                                                     protein['PDB'])
            logger.debug(cmd_no_color)
            sub = subprocess.run(cmd_no_color, shell=True,
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE)
            # capturing the output
            if sub.stdout:
                logger.debug('{}: {}'.format(cmd_no_color,
                                             sub.stdout.decode('utf-8')))
            if sub.stderr:
                logger.error('{}: {}'.format(cmd_no_color,
                                             sub.stderr.decode('utf-8')))

        # create the commands for the python script which generates the pymol
        # pictures with colored AA
        cmd_list = []
        for aa_idx, frame_idx in enumerate(pdb_data['frames_idx']):
            if '{}_{}.png'.format(protein['PDB'],
                                  frame_idx) not in existing_frames:
                cmd = ('./create_pdb_frames.py -p {} -c {} '
                       '-n {} -i {} --color_aa {}').format(pdb_dir,
                                                           pdb_data['chain'],
                                                           frame_idx,
                                                           aa_idx + 1,
                                                           protein['PDB'])
                cmd_list.append(cmd)

        # threading to run the commands
        if cmd_list:
            nb_threads_to_do = len(cmd_list)
            nb_threads_done = 0
            errors = 0
            msg = ('Creating {} ({}) protein colored AA frames, '
                   'please wait..').format(protein['entry_name'],
                                           protein['PDB'])
            print('\n{}'.format(msg))
            logger.info(msg)
            with concurrent.futures.ThreadPoolExecutor(max_workers=5) as executor:
                for cmd in cmd_list:
                    logger.info(cmd)
                    thread = executor.submit(subprocess.run,
                                             cmd,
                                             shell=True,
                                             stdout=subprocess.PIPE,
                                             stderr=subprocess.PIPE)
                    # capturing the output
                    if thread.result().stdout:
                        logger.debug('{}: {}'.format(cmd,
                                                     thread.result().stdout))
                        nb_threads_done += 1
                    if thread.result().stderr:
                        logger.error('{}: {}'.format(cmd,
                                                     thread.result().stderr))
                        nb_threads_done += 1
                        errors += 1
                    print('{}/{} threads ({} errors)'.format(nb_threads_done,
                                                             nb_threads_to_do,
                                                             errors))
        # check if all frames are created else wait
        while len(os.listdir(frames_dir)) != (len(pdb_data['frames_idx']) + 1):
            time.sleep(1)
        # create the movie
        movie_path = os.path.join(out_dir, '{}.avi'.format(file_base_name))
        if not os.path.exists(movie_path):
            audio_video.create_movie(movie_path,
                                     frames_dir,
                                     keys_duration,
                                     flac_file_path)
        else:
            msg = 'Movie file already exists: {}'.format(movie_path)
            print(msg)
            logger.info(msg)
    else:
        msg = ('No movie created, no PDB entry (3D image reference) '
               'in the Uniprot entry for {}').format(args.uniprot_AN)
        logger.info(msg)
        print(msg)

    # create the score
    if args.score:
        msg = 'Creating the score:'
        logger.info(msg)
        print(msg)
        score_basename = '{}_{}_{}_{}bpm_score.pdf'.format(args.uniprot_AN,
                                                           protein['entry_name'],
                                                           protein['organism'],
                                                           tempo)
        score_output = os.path.join(args.out, score_basename)
        cmd = 'mscore -o {} {}'.format(score_output, midi_file_path)
        subprocess.run(cmd, shell=True)
        msg = 'Score created at {}'.format(score_output)
        logger.info(msg)
        print(msg)
