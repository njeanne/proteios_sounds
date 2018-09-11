#! /usr/bin/env python3

__author__ = 'Nina Verstraete, Jacques TOEN & Nicolas JEANNE'
__copyright__ = 'GNU General Public License'
__version__ = '1.0'
__email__ = 'nicolas.jeanne@ntymail.com'

import argparse
import sys
import os
import re
import logging
import urllib
import subprocess
from datetime import datetime
from Bio import ExPASy
from Bio import SwissProt
from midiutil import MIDIFile


# Check range for tempo argument, must be between 60 and 150.
# @param x:	[str] value of tempo argument in BPM.
# @return: [int] the tempo.
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
parser.add_argument('-i', '--instruments', required=False, nargs=3, help='set channel 0, 1 and 2 instruments, restricted to 3 values between 0 and 127 separated by spaces. Default is 0:  Acoustic Grand, 42: Cello and 65: Alto Sax. See: http://www.pjb.com.au/muscript/gm.html#patch for details.')
parser.add_argument('-d', '--debug', required=False, action='store_true', help='debug mode, create a log file which details each entry of the MIDI file.')
parser.add_argument('uniprot_accession_number', help='the protein accession number in the UniProt database. Example: Human Interleukin-8 > P10145')
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

### midi notes correspondance with AA
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
if args.debug:
	logPath = os.path.join(out_dir, '{}_{}_{}bpm_{}.log'.format(args.uniprot_accession_number, args.mode, tempo, datetime.now().strftime('%Y-%m-%d_%H-%M-%S')))
	logging.basicConfig(filename=logPath, level=logging.DEBUG, format='%(asctime)s\t%(levelname)s:\t%(message)s', datefmt='%Y/%m/%d %H:%M:%S')
	logger = logging.getLogger(__name__)
	logger.info(' '.join(sys.argv))


if args.debug:
	logger.info('tempo: {} BPM'.format(tempo))

# Parse the Uniprot entry
try:
	handle = ExPASy.get_sprot_raw(args.uniprot_accession_number)
except urllib.error.HTTPError as err:
	msg = '{}{} accession number does not exists in UniProt database, check https://www.uniprot.org/ for the correct accession number.'.format(err, args.uniprot_accession_number)
	print('ERROR: {}'.format(msg))
	if args.debug:
		logger.error(msg)
	sys.exit(0)

uniprot = SwissProt.read(handle)

### TOREMOVE
print(uniprot.sequence)
for item in uniprot.features:
	print(item)
##############

# get the organism and entry name in UniProt
organism = uniprot.organism
if '(' in organism:
	 organism = organism.split('(')[1].split(')')[0]
organism = organism.rstrip('.,;:').replace(' ', '_')
entry_name = uniprot.entry_name

# create a dictionary for the protein
protein = {'seq': {}, 'modified_residue': {}, 'structure': {}, 'glycosylation': {}, 'site': {}, 'propeptide': {}}
sequence_length = len(uniprot.sequence)
for i in range(sequence_length):
	protein['seq'][i+1] = uniprot.sequence[i]


# create a dictionary to get the structures positions => type
structures = {}
# parse the UniProt Features (see https://www.uniprot.org/help/?query=*&fil=category%3A%22PTM+processing%22+AND+section%3Amanual&columns=title)
for feature in uniprot.features:
	if feature[0] == 'HELIX' or feature[0] == 'STRAND' or feature[0] == 'TURN':
		structures[feature[1]] = {'type': feature[0], 'end': feature[2]}
	# Modified residue: Modified residues excluding lipids, glycans and protein crosslinks
	if feature[0] == 'MOD_RES':
		modif = re.split(' {', feature[3])[0].rstrip('.;,')
		if modif in protein['modified_residue'].keys():
			protein['modified_residue'][modif].append((feature[1], feature[2]))
		else:
			protein['modified_residue'][modif] = [(feature[1], feature[2])]
	# Glycosylation: Covalently attached glycan group(s)
	if feature[0] == 'CARBOHYD':
		glyco = re.split(' ', feature[3])[0].rstrip('.;,')
		if glyco in protein['glycosylation'].keys():
			protein['glycosylation'][glyco].append((feature[1], feature[2]))
		else:
			protein['glycosylation'][glyco] = [(feature[1], feature[2])]
	# Site: Any interesting single amino acid site on the sequence
	if feature[0] == 'SITE':
		site = re.split(' {', feature[3])[0].rstrip('.;,')
		if site in protein['site'].keys():
			protein['site'][site].append((feature[1], feature[2]))
		else:
			protein['site'][site] = [(feature[1], feature[2])]
	# Propeptide: Part of a protein that is cleaved during maturation or activation
	if feature[0] == 'PROPEP':
		propep = re.split(' {', feature[3])[0].rstrip('.;,')
		if propep in protein['propeptide'].keys():
			protein['propeptide'][propep].append((feature[1], feature[2]))
		else:
			protein['propeptide'][propep] = [(feature[1], feature[2])]


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

### TOREMOVE
for k,v in protein.items():
	if k == 'modified_residue':
		print(k)
		for k2,v2 in v.items():
			print('\t{}'.format(k2))
			for pos in v2:
				print('\t\t{}'.format(pos))
	elif k == 'glycosylation':
		print(k)
		for k2,v2 in v.items():
			print('\t{}'.format(k2))
			for pos in v2:
				print('\t\t{}'.format(pos))
	elif k == 'site':
		print(k)
		for k2,v2 in v.items():
			print('\t{}'.format(k2))
			for pos in v2:
				print('\t\t{}'.format(pos))
	elif k == 'propeptide':
		print(k)
		for k2,v2 in v.items():
			print('\t{}'.format(k2))
			for pos in v2:
				print('\t\t{}'.format(pos))
########################

# create the MIDI file
mode_name = args.mode
file_base_name = '{}_{}_{}_{}_{}bpm_instrus'.format(args.uniprot_accession_number, entry_name, organism, mode_name,tempo)
for instru in instrus:
	file_base_name = '{}-{}'.format(file_base_name, instru)
midi_file_path = os.path.join(out_dir, '{}.midi'.format(file_base_name))

with open(midi_file_path, 'wb') as  midiFile:

	track = 0
	time = 0   # In beats

	# a channel is defined by an instrument nbr and a volume (0-127, as per the MIDI standard, see: http://www.pjb.com.au/muscript/gm.html)
	# channel 9 is for percussions (see https://pjb.com.au/muscript/gm.html#perc)
	channels = {0: {'instrument': instrus[0], 'vol': 100}, 
				1: {'instrument': instrus[1], 'vol': 40}, 
				2: {'instrument': instrus[2], 'vol': 60}}
	
	if args.debug:
		logger.info('Instrument number by channel, see: http://www.pjb.com.au/muscript/gm.html for instruments number correspondance:')
		for channel_nb in channels.keys():
			logger.info('\tchannel {}: instrument {}'.format(channel_nb, channels[channel_nb]['instrument']))

	MyMIDI = MIDIFile(numTracks=1, adjust_origin=False) # One track, defaults to format 1 (tempo track automatically created)
	MyMIDI.addTempo(track, time, tempo)
	# add the channels (1 per instrument)
	for channel_nbr in channels.keys():
		MyMIDI.addProgramChange(track, channel=channel_nbr, time=time, program=channels[channel_nbr]['instrument'])

	mode = MIDI_KEYS[mode_name]

	for i in range(1, len(protein['seq']) + 1):
		AA = protein['seq'][i]
		if i == 1:
			next_AA = protein['seq'][sequence_length]
		elif i == len(protein['seq']):
			next_AA = protein['seq'][1]
		else:
			next_AA = protein['seq'][i+1]

		# set the duration depending on the number of shared properties
		shared_properties = len(set.intersection(AA_PHY_CHI[AA], AA_PHY_CHI[next_AA]))
		if shared_properties == 0:
			duration = 1
		elif shared_properties == 1:
			duration = 1.5
		elif shared_properties == 2:
			duration = 2
		else:
			duration = 4

		# change the volume of each instrument depending on the structure
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
				else:
					channels[0]['vol'] = 100
					channels[1]['vol'] = 60
					channels[2]['vol'] = 40

		# test if note is a chord
		if isinstance(mode[AA], list):
			for channel_nbr in channels.keys():
				for note in mode[AA]:
					MyMIDI.addNote(track, channel=channel_nbr, pitch=note, time=time, duration=duration, volume=channels[channel_nbr]['vol'])
		# just a note
		else:
			note = mode[AA]
			for channel_nbr in channels.keys():
				MyMIDI.addNote(track, channel=channel_nbr, pitch=note, time=time, duration=duration, volume=channels[channel_nbr]['vol'])

		if args.debug:
			if isinstance(mode[AA], list):
				played_note = '[{}]'.format(' , '.join(str(v) for v in mode[AA]))
			else:
				played_note = mode[AA]
			logging.debug('\tposition: {:<5}\tAA: {:<5}\tnote: {:<20}\ttime: {:<5}\tduration: {:<5}'.format(i, AA, played_note, time, duration))
			for channel_nbr in channels.keys():
				logging.debug('\t\tchannel: {:<5}\tvolume: {:<5}'.format(channel_nbr, channels[channel_nbr]['vol']))

		time = time + duration

	MyMIDI.writeFile(midiFile)

print('Done!\nMIDI file for {} {} ({}) created in: {}'.format(entry_name, organism, args.uniprot_accession_number, midi_file_path))


# play the file with timidity if asked
if args.play:
	cmd = 'timidity {}'.format(midi_file_path)
	subprocess.run(cmd, shell=True)
