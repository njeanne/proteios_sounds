#! /usr/bin/python3

__author__ = 'Nina Verstraete, Jacques TOEN & Nicolas JEANNE'
__copyright__ = 'GNU General Public License'
__version__ = '1.0'
__email__ = 'nicolas.jeanne@ntymail.com'

import argparse
import sys
import os
import re
import logging
import subprocess
from datetime import datetime
from Bio import ExPASy
from Bio import SwissProt
from midiutil import MIDIFile


# Check range for tempo argument, must be between 60 and 150.
# @param x:	[str] value of tempo argument in BPM.
# @return:		[int] the tempo.
def restricted_tempo(x):
	x = int(x)
	if x < 60 or x > 150:
		raise argparse.ArgumentTypeError('{} not in range 60 to 150.'.format(x))
	return x

descr = '''
proteios_sounds.py v.{}

Created by {}.
Contact: {}
{}

Create a MIDI file from a protein entry of the UniProt database (https://www.uniprot.org/).
'''.format(__version__, __author__, __email__, __copyright__)

# Parse arguments
parser = argparse.ArgumentParser(description=descr, formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('-o', '--out', required=True, help='path to the results directory.')
parser.add_argument('-m', '--mode', required=True, choices=['major', 'mixolydian', 'dorian', 'blues'], help='The mode to apply: major, mixolydian, dorian or blues.')
parser.add_argument('-p', '--play', required=False, action='store_true', help='play the music with Timidity, just for tests.')
parser.add_argument('-t', '--tempo', required=False, type=restricted_tempo, help='set the tempo in BPM. Value between 60 and 150.')
parser.add_argument('-d', '--debug', required=False, action='store_true', help='debug mode, create a log file which details each entry of the MIDI file.')
parser.add_argument('uniprot_accession_number', help='the protein accession number in the UniProt database. Example: Human Interleukin-8 > P10145')
args = parser.parse_args()


### midi keys correspondance with AA
# Major
MIDI_KEYS ={'major':{'A':48, 'C':50, 'D':52, 'E':53, 'F':55, 'G':57, 'H':59, 'I':60, 'K':62, 'L':64, 'M':65, 'N':67, 'P':69, 'Q':71,
'R':[48,52,55,59], 'S':[48,50,53,57], 'T':[52,55,59,62], 'U':[52,53,57,60], 'V':[50,53,55,59], 'W':[48,52,55,57], 'Y':[50,53,57,59]}}
# Mixolydian
MIDI_KEYS['mixolydian'] = MIDI_KEYS['major']
MIDI_KEYS['mixolydian']['H']=58
MIDI_KEYS['mixolydian']['R'][3]=58
MIDI_KEYS['mixolydian']['T'][2]=58
MIDI_KEYS['mixolydian']['V'][3]=58
MIDI_KEYS['mixolydian']['Y'][3]=58
MIDI_KEYS['mixolydian']['Q']=70
# Dorian
MIDI_KEYS['dorian'] = MIDI_KEYS['mixolydian']
MIDI_KEYS['dorian']['D'] = 51
MIDI_KEYS['dorian']['R'][1] = 51
MIDI_KEYS['dorian']['T'][0] = 51
MIDI_KEYS['dorian']['U'][0] = 51
MIDI_KEYS['dorian']['W'][1] = 51
MIDI_KEYS['dorian']['L'] = 63
# Blues
MIDI_KEYS['blues'] = MIDI_KEYS['dorian']
MIDI_KEYS['blues']['E'] = 52
MIDI_KEYS['blues']['S'][2] = 52
MIDI_KEYS['blues']['U'][1] = 52 # donne 'U':[52,52,57,60] ??????
MIDI_KEYS['blues']['V'][1] = 52
MIDI_KEYS['blues']['Y'][1] = 52
MIDI_KEYS['blues']['M'] = 64

# create the output directory
out_dir = os.path.abspath(args.out)
if not os.path.exists(out_dir):
	os.makedirs(out_dir)

if args.debug:
	# create the log file
	logPath = os.path.join(out_dir, '{}_{}_{}.log'.format(args.uniprot_accession_number, args.mode, datetime.now().strftime('%Y-%m-%d_%H-%M-%S')))
	logging.basicConfig(filename=logPath, level=logging.DEBUG, format='%(asctime)s\t%(levelname)s:\t%(message)s', datefmt='%Y/%m/%d %H:%M:%S')
	logger = logging.getLogger(__name__)
	logger.info(' '.join(sys.argv))

# Parse the Uniprot entry
handle = ExPASy.get_sprot_raw(args.uniprot_accession_number)
uniprot = SwissProt.read(handle)
organism = uniprot.organism.split('(')[1].split(')')[0]
entry_name = uniprot.entry_name

# create a dictionary for the protein
protein = {'seq': {}}
sequence_length = len(uniprot.sequence)
for i in range(sequence_length):
	protein['seq'][i+1] = uniprot.sequence[i]
#print(dir(uniprot))

# create a dictionary to get the structures positions => type
structures = {}
for feature in uniprot.features:
	print(feature)
	if feature[0] == 'MOD_RES':
		if 'modification' in protein.keys():
			protein['modification'][feature[1]] = re.split('\W+', feature[3])[0]
		else:
			protein['modification'] = {feature[1]: re.split('\W+', feature[3])[0]}
	if feature[0] == 'HELIX' or feature[0] == 'STRAND' or feature[0] == 'TURN':
		structures[feature[1]] = {'type': feature[0], 'end': feature[2]}

# update protein with the structures
structure_last_position = 0
for position in sorted(structures.keys()):
	if structure_last_position == 0 and position == 1:
		if 'structure' in protein.keys():
			protein['structure'][position] = structures[position]['type']
		else:
			protein['structure'] = {position: structures[position]['type']}
		structure_last_position = structures[position]['end']
	elif position >= structure_last_position + 1:
		if 'structure' in protein.keys():
			protein['structure'][structure_last_position + 1] = 'FREE'
		else:
			protein['structure'] = {structure_last_position + 1: 'FREE'}
		protein['structure'][position] = structures[position]['type']
		structure_last_position = structures[position]['end']
if structure_last_position < sequence_length:
	if 'structure' in protein.keys():
		protein['structure'][structure_last_position + 1] = 'FREE'
	else:
		protein['structure'] = {structure_last_position + 1: 'FREE'}

# create the MIDI file
mode_name = args.mode
file_base_name = '{}_{}_{}_{}'.format(args.uniprot_accession_number, entry_name, organism, mode_name)
midi_file_path = os.path.join(out_dir, '{}.midi'.format(file_base_name))

with open(midi_file_path, 'wb') as  midiFile:

	track    = 0
	time     = 0   # In beats
	if args.tempo:
		tempo = int(args.tempo)
	else:
		tempo    = 100  # In BPM
	# a channel is defined by an instrument nbr and a volume (0-127, as per the MIDI standard)
	channels = {0: {'instrument': 1, 'vol': 100}, 1: {'instrument': 42, 'vol': 40}, 2: {'instrument': 65, 'vol': 60}}

	MyMIDI = MIDIFile(numTracks=1, adjust_origin=False) # One track, defaults to format 1 (tempo track automatically created)
	MyMIDI.addTempo(track, time, tempo)
	# add the channels (1 per instrument)
	for channelNbr in channels.keys():
		MyMIDI.addProgramChange(track, channel=channelNbr, time=time, program=channels[channelNbr]['instrument'])

	mode = MIDI_KEYS[mode_name]
	for i in range(1, len(protein['seq'])):
		AA = protein['seq'][i]
		duration = 1

		if 'structure' in protein.keys():
			if i in protein['structure'].keys():
				logging.debug('{}:'.format(protein['structure'][i]))
				if protein['structure'][i] == 'HELIX':
					channels[0]['vol'] = 40
					channels[1]['vol'] = 100
					channels[2]['vol'] = 60
				elif protein['structure'][i] == 'STRAND':
					channels[0]['vol'] = 60
					channels[1]['vol'] = 40
					channels[2]['vol'] = 100
				elif protein['structure'][i] == 'TURN':
					channels[0]['vol'] = 100
					channels[1]['vol'] = 40
					channels[2]['vol'] = 60

		# test if key is a chord
		if isinstance(mode[AA], list):
			for channelNbr in channels.keys():
				for key in mode[AA]:
					MyMIDI.addNote(track, channel=channelNbr, pitch=key, time=time, duration=duration, volume=channels[channelNbr]['vol'])
		# just a key
		else:
			key = mode[AA]
			for channelNbr in channels.keys():
				MyMIDI.addNote(track, channel=channelNbr, pitch=key, time=time, duration=duration, volume=channels[channelNbr]['vol'])
		logging.debug('\tAA: {}\tkey: {}\ttime: {}\tchannel: {}\tvolume: {}\tduration: {}'.format(AA, mode[AA], time, channelNbr, channels[channelNbr]['vol'], duration))
		time = time + duration

	MyMIDI.writeFile(midiFile)

print('Done!\nMIDI file for {} {} ({}) created in: {}'.format(entry_name, organism, args.uniprot_accession_number, midi_file_path))


# play the file with timidity if asked
if args.play:
	cmd = 'timidity {}'.format(midi_file_path)
	subprocess.run(cmd, shell=True)
