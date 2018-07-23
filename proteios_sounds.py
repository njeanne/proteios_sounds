#! /usr/bin/python3

import argparse
import os
import re
from Bio import ExPASy
from Bio import SwissProt
from Bio import SeqIO
from midiutil import MIDIFile
import subprocess

# Parse arguments
parser = argparse.ArgumentParser(description='Create a MIDI file from a protein entry of the UniProt database (https://www.uniprot.org/)')
parser.add_argument('-o', '--out', required=True, help='path to the results directory.')
parser.add_argument('-m', '--mode', required=True, choices=['major', 'mixolydian', 'dorian', 'blues', 'persian'], help='The mode to apply: major, mixolydian, dorian or blues.')
parser.add_argument('-p', '--play', required=False, help='play the music.')
parser.add_argument('-t', '--tempo', required=False, type=int, help='tempo selection in BPM. Value between 60 and 150.')
parser.add_argument('uniprot_accession_number', help='the protein accession number in UniProt database. Example: Hemoglobin sub-unit A1 > P69905')
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



outDir = os.path.abspath(args.out)
if not os.path.exists(outDir):
	os.makedirs(outDir)
modeName = args.mode

# Parse the Uniprot entry
handle = ExPASy.get_sprot_raw(args.uniprot_accession_number)
uniprot = SwissProt.read(handle)
organism = uniprot.organism.split('(')[1].split(')')[0]
entry_name = uniprot.entry_name

protein = {'seq': {}}
for i in range(0,len(uniprot.sequence)):
	protein['seq'][i+1] = uniprot.sequence[i]
#print(dir(uniprot))

for feature in uniprot.features:
	if feature[0] == 'MOD_RES':
		if 'modification' in protein.keys():
			protein['modification'][feature[1]] = re.split('\W?', feature[3])[0]
		else:
			protein['modification'] = {feature[1]: re.split('\W?', feature[3])[0]}
	if feature[0] == 'HELIX' or feature[0] == 'STRAND' or feature[0] == 'TURN':
		if 'structure' in protein.keys():
			protein['structure'][feature[1]] = feature[0]
		else:
			protein['structure'] = {feature[1]: feature[0]}


# create the MIDI file
file_base_name = '{}_{}_{}_{}'.format(args.uniprot_accession_number, entry_name, organism, modeName)
midiFilePath = os.path.join(outDir, '{}.midi'.format(file_base_name))
debugFilePath = os.path.join(outDir, '{}.csv'.format(file_base_name))	
with open(midiFilePath, 'wb') as  midiFile:
	
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

	debugFile = open(debugFilePath, 'w')
	debugFile.write('AA\tkey\ttime\tchannel\tvolume\tduration')
	mode = MIDI_KEYS[modeName]
	for i in range(1, len(protein['seq'])):
		AA = protein['seq'][i]
		duration = 1
		
		if 'structure' in protein.keys():
			if i in protein['structure'].keys():
				debugFile.write('\n{}'.format(protein['structure'][i]))
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
				debugFile.write('\n{};{};{};{};{};{}'.format(AA, mode[AA], time, channelNbr, channels[channelNbr]['vol'], duration))
				for key in mode[AA]:
					MyMIDI.addNote(track, channel=channelNbr, pitch=key, time=time, duration=duration, volume=channels[channelNbr]['vol'])		
		# just a key
		else:
			key = mode[AA]
			for channelNbr in channels.keys():
				MyMIDI.addNote(track, channel=channelNbr, pitch=key, time=time, duration=duration, volume=channels[channelNbr]['vol'])
				debugFile.write('\n{};{};{};{};{};{}'.format(AA, key, time, channelNbr, channels[channelNbr]['vol'], duration))
		time = time + duration
		
	MyMIDI.writeFile(midiFile)
	debugFile.close()


# play the file with timidity if asked
if args.play:
	cmd = 'timidity {}'.format(midiFilePath)
	subprocess.run(cmd, shell=True)

