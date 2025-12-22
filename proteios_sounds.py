#! /usr/bin/env python3

__author__ = "Nina VERSTRAETE, Jacques TOEN & Nicolas JEANNE"
__copyright__ = "GNU General Public License"
__version__ = "1.1.0"
__email__ = "mesmeraf@gmail.com"

import argparse
import sys
import os
import logging
import multiprocessing
import subprocess
import time
import parse_uniprot
import midi_operations
import parse_pdb
import protein_movie

# add pymol to the python path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "lib/python3.7/site-packages"))
import pymol


def restricted_tempo(tempo_value):
    """
    Check range for tempo argument, must be between 60 and 150.
    :param tempo_value: value of the tempo argument in BPM.
    :type tempo_value: str
    :return: the tempo.
    :rtype: int
    """
    tempo_value = int(tempo_value)
    if tempo_value < 60 or tempo_value > 250:
        raise argparse.ArgumentTypeError("f{tempo_value} not in range 60 to 250.")
    return tempo_value


def create_log(path, level):
    """Create the log as a text file and as a stream.

    :param path: the path of the log.
    :type path: str
    :param level: the level og the log.
    :type level: str
    :return: the logging:
    :rtype: logging
    """

    log_level_dict = {"DEBUG": logging.DEBUG,
                      "INFO": logging.INFO,
                      "WARNING": logging.WARNING,
                      "ERROR": logging.ERROR,
                      "CRITICAL": logging.CRITICAL}

    if level is None:
        log_level = log_level_dict["INFO"]
    else:
        log_level = log_level_dict[level]

    if os.path.exists(path):
        os.remove(path)

    logging.basicConfig(format="%(asctime)s %(levelname)s:\t%(message)s",
                        datefmt="%Y/%m/%d %H:%M:%S",
                        level=log_level,
                        handlers=[logging.FileHandler(path), logging.StreamHandler()])
    return logging


def create_pdb_frames(pdb_accession_number, chain, idx_aa, pdb_directory, frame_nb, color_aa):
    """
    Creates the frames from PDB data.

    :param pdb_accession_number: the PDB accession number.
    :type pdb_accession_number: str
    :param chain: the chain id.
    :type chain: str
    :param idx_aa: the AA idx in the PDB.
    :type idx_aa: int
    :param pdb_directory: the path to the pdb data folder.
    :type pdb_directory: str
    :param frame_nb: the frame number.
    :type frame_nb: str
    :param color_aa: color in red the current AA.
    :type color_aa: bool
    """
    # open pymol and retrieve the protein with PDB accession number
    pymol.finish_launching(["pymol", "-qc"])  # Pymol: quiet and no GUI
    # set the path to download the PDB data
    pymol.cmd.load(os.path.join(pdb_directory, f"{pdb_accession_number.lower()}.cif"))
    pymol.cmd.disable("all")
    pymol.cmd.enable(pdb_accession_number)
    pymol.stored_list = []
    pymol.cmd.iterate(f"(name ca) and (chain {chain})", "pymol.stored_list.append((resi, oneletter))")
    pymol.cmd.hide("all")
    pymol.cmd.show("cartoon")
    pymol.cmd.set("ray_opaque_background", 1)

    if color_aa:
        pymol.cmd.color("red", f"resi {idx_aa}")
    img_path = os.path.join(pdb_directory, "frames", f"{pdb_accession_number}_{frame_nb}.png")
    logging.info(f"[Pymol] Frame {frame_nb} (in PDB file): {img_path}")
    pymol.cmd.png(img_path, width=800, height=600, quiet=1)
    pymol.cmd.quit()


if __name__ == "__main__":
    descr = f"""
    {os.path.splitext(os.path.basename(__file__))[0]} v.{__version__}

    Created by {__author__}.
    Contact: {__email__}
    {__copyright__}

    Create a MIDI file and from a protein entry of the UniProt database (https://www.uniprot.org/).
    If the data are available in the UniProt entry, a movie file of the 3D representation of the protein will also be 
    created.
    """

    # Parse arguments
    parser = argparse.ArgumentParser(description=descr,
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-o", "--out", required=True, help="path to the results directory.")
    parser.add_argument("-s", "--score", required=False, action="store_true",
                        help="use musescore software to create the score corresponding to the MIDI file.")
    parser.add_argument("-t", "--tempo", required=False, type=restricted_tempo,
                        help="set the tempo in BPM. Value between 60 and 250.")
    parser.add_argument("-i", "--instruments", required=False, nargs=3,
                        help="set channel 0, 1 and 2 instruments, restricted to 3 values between 1 and 128 separated "
                             "by spaces. Default is 1:  Acoustic Grand, 43: Cello and 66: Alto Sax. "
                             "See: https://en.wikipedia.org/wiki/General_MIDI for details.")
    parser.add_argument("-f", "--force", required=False, action="store_true",
                        help="if the video file exists, force to recreate it.")
    parser.add_argument("-l", "--log", required=False, type=str,
                        help="the path for the log file. If this option is skipped, the log file is created in the "
                             "output directory.")
    parser.add_argument("--log-level", required=False, type=str,
                        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
                        help="set the log level. If the option is skipped, log level is INFO.")
    parser.add_argument("uniprot",
                        help="the protein Accession Number in the UniProt database. Example: Human Interleukin-8 > "
                             "P10145")
    args = parser.parse_args()

    # check if instruments are between 0 and 127
    if args.instruments:
        for i in range(len(args.instruments)):
            instrument = int(args.instruments[i]) - 1
            if instrument < 0 or instrument > 127:
                raise argparse.ArgumentTypeError(f"{args.instruments} should be 3 integers between 1 and 128.")
            args.instruments[i] = instrument
        instruments = args.instruments
    else:
        instruments = [0, 42, 65]

    # tempo
    if args.tempo:
        tempo = int(args.tempo)
    else:
        tempo = 100  # In BPM

    # MIDI keys on major mode correspondance with AA sorted by decreasing molecular weight are set as "C" or "DO" in
    # French (48, 60, 72) degrees I, "G" or "SOL" in French (55, 67) degrees V, "F" or "FA" in French (53, 65) degrees
    # IV, "D" or "RE" in French (50, 62) degrees II, "E" or "MI" in French (52, 64) degrees III, "A" or "LA" in French
    # (57, 69) degrees VI and "B" or "SI" in French (59, 71) degrees VII. Finally, we add 7 alterations "#" following
    # the ascending quint (54, 66, 49, 61, 56, 68, 51)
    initial_midi_keys = [48, 60, 72, 55, 67, 53, 65, 50, 62, 52, 64, 57, 69,
                         59, 71, 54, 66, 49, 61, 56, 68, 51]
    midi_keys = {}

    # Physico-chemical properties of AA
    AA_PHY_CHI = {"A": {"hydrophobic", "small"},
                  "R": {"polar", "pos_charged"},
                  "N": {"polar", "small"},
                  "D": {"polar", "small", "neg_charged"},
                  "C": {"hydrophobic", "polar", "small"},
                  "E": {"polar", "neg_charged"},
                  "Q": {"polar"},
                  "G": {"hydrophobic", "small"},
                  "H": {"hydrophobic", "polar", "pos_charged", "aromatic"},
                  "I": {"hydrophobic", "aliphatic"},
                  "L": {"hydrophobic", "aliphatic"},
                  "K": {"hydrophobic", "polar", "pos_charged"},
                  "M": {"hydrophobic"},
                  "F": {"hydrophobic", "aromatic"},
                  "P": {"small"},
                  "S": {"polar", "small"},
                  "T": {"hydrophobic", "polar", "small"},
                  "W": {"hydrophobic", "polar", "aromatic"},
                  "Y": {"hydrophobic", "polar", "aromatic"},
                  "V": {"hydrophobic", "small", "aliphatic"}}

    # create output directory if necessary
    os.makedirs(args.out, exist_ok=True)
    # create the logger
    if args.log:
        log_path = args.log
    else:
        log_path = os.path.join(args.out, f"{os.path.splitext(os.path.basename(__file__))[0]}.log")
    create_log(log_path, args.log_level)

    logging.info(f"version: {__version__}")
    logging.info(f"CMD: {' '.join(sys.argv)}")

    logging.info(f"\tOutput directory: {args.out}")
    logging.info(f"\tTempo: {tempo} BPM")
    logging.info(f"\tInstruments: {', '.join(map(str, [x + 1 for x in instruments]))} (general MIDI patch numbers, "
                 "see: https://en.wikipedia.org/wiki/General_MIDI)")
    logging.info(f"\tCreate score: {args.score}")

    # parsing of uniprot entry
    protein = parse_uniprot.parse_entry(args.uniprot)

    sequence = protein["seq"]
    sequence_length = len(sequence)
    protein["seq"] = {}
    logging.debug(f"AA sequence ({sequence_length} AA): {sequence}")
    for i in range(sequence_length):
        protein["seq"][i] = sequence[i]
    # frequency of AA in the sequence
    set_AA = set(''.join(sequence))
    proportion_AA = {}
    for aa in set_AA:
        proportion_AA[aa] = sequence.count(aa) / sequence_length
    # sort by decreasing frequency
    proportion_AA = sorted(proportion_AA.items(), key=lambda kv: kv[1], reverse=True)

    for idx, aa_proportion in enumerate(proportion_AA):
        midi_keys[aa_proportion[0]] = initial_midi_keys[idx]

    # set the result files base name
    file_base_name = f"{args.uniprot}_{protein['entry_name']}_{protein['organism']}_{tempo}bpm_intrus"
    for instrument in instruments:
        file_base_name = f"{file_base_name}-{instrument}"

    # create the MIDI file
    midi_file_path = os.path.join(args.out, f"{file_base_name}.midi")
    keys_duration = midi_operations.create_midi(midi_file_path, protein, midi_keys, tempo, instruments, AA_PHY_CHI)

    if "PDB" in protein:
        multiprocessing.set_start_method("spawn")
        # create the directories for PDB data and frames
        pdb_dir = os.path.join(os.path.abspath(args.out), "pdb", f"{protein['accession_number']}_{protein['PDB']}")
        frames_dir = os.path.join(pdb_dir, "frames")
        os.makedirs(frames_dir, exist_ok=True)

        # get data from the PDB file
        pdb_data = parse_pdb.get_pdb_info(protein, pdb_dir)

        # create a frame without colored AA for all AA outside the PDB data
        existing_frames = sorted([png for png in os.listdir(frames_dir)])
        if f"{protein['PDB']}_no-idx.png" not in existing_frames:
            logging.info(f"Creating {protein['entry_name']} ({protein['PDB']}) protein frame, please wait..")
            processes = []
            process = multiprocessing.Process(target=create_pdb_frames,
                                              args=(protein["PDB"], pdb_data["chain"], 1, pdb_dir, "no-idx", False))
            processes.append(process)
            process.start()
            for process in processes:
                process.join()

        # create the commands for the python script which generates the pymol pictures with colored AA
        amino_acids_indexes = []
        frames_indexes = []
        for aa_idx, frame_idx in enumerate(pdb_data["frames_idx"]):
            if f"{protein['PDB']}_{frame_idx}.png" not in existing_frames:
                amino_acids_indexes.append(aa_idx)
                frames_indexes.append(frame_idx)

        if amino_acids_indexes:
            processes = []
            for idx in range(len(amino_acids_indexes)):
                process = multiprocessing.Process(target=create_pdb_frames,
                                                  args=(protein["PDB"], pdb_data["chain"], frames_indexes[idx], pdb_dir,
                                                        amino_acids_indexes[idx], True))
                processes.append(process)
                process.start()
            for process in processes:
                process.join()

        # check if all frames are created else wait
        while len(os.listdir(frames_dir)) != (len(pdb_data["frames_idx"]) + 1):
            time.sleep(1)
        # create the movie
        movie_path = os.path.join(args.out, f"{file_base_name}.mp4")
        if args.force or not os.path.exists(movie_path):
            protein_movie.create_movie(movie_path, frames_dir, keys_duration, midi_file_path)
        else:
            logging.info(f"Movie file already exists: {movie_path}")

    # create the score
    if args.score:
        logging.info("Creating the score:")
        score_basename = f"{args.uniprot}_{protein['entry_name']}_{protein['organism']}_{tempo}bpm_score.pdf"
        score_output = os.path.join(args.out, score_basename)
        cmd = f"mscore -o {score_output} {midi_file_path}"
        subprocess.run(cmd, shell=True)
        logging.info(f"\tScore created at {score_output}")

    # play the file with timidity if asked
    if args.play:
        cmd = f"timidity {midi_file_path}"
        subprocess.run(cmd, shell=True)
