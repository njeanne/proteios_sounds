#! /usr/bin/env python3

import os
import logging
from midiutil import MIDIFile

def create_midi(uniprot_accession_number,
                protein,
                midi_keys,
                tempo,
                instrus,
                output_directory,
                aa_phy_chi,
                logger,
                debug=False):

    file_base_name = '{}_{}_{}_{}bpm_instrus'.format(uniprot_accession_number, protein['entry_name'], protein['organism'], tempo)

    for instru in instrus:
        file_base_name = '{}-{}'.format(file_base_name, instru)
    midi_file_path = os.path.join(output_directory, '{}.midi'.format(file_base_name))

    with open(midi_file_path, 'wb') as  midiFile:

        track = 0
        time = 0   # In beats

        # a channel is defined by an instrument nbr and a volume (0-127, as per the MIDI standard, see: http://www.pjb.com.au/muscript/gm.html)
        # channel 9 is for percussions (see https://pjb.com.au/muscript/gm.html#perc)
        channels = {0: {'instrument': instrus[0], 'vol': 100},
                    1: {'instrument': instrus[1], 'vol': 40},
                    2: {'instrument': instrus[2], 'vol': 60}}

        if debug:
            logger.info('Instrument number by channel, see: http://www.pjb.com.au/muscript/gm.html for instruments number correspondance:')
            for channel_nb in channels.keys():
                logger.info('\tchannel {}: instrument {}'.format(channel_nb, channels[channel_nb]['instrument']))

        MyMIDI = MIDIFile(numTracks=1, adjust_origin=False) # One track, defaults to format 1 (tempo track automatically created)
        MyMIDI.addTempo(track, time, tempo)
        # add the channels (1 per instrument)
        for channel_nbr in channels.keys():
            MyMIDI.addProgramChange(track, channel=channel_nbr, time=time, program=channels[channel_nbr]['instrument'])

        sequence_length = len(protein['seq'])

        for i in range(0, sequence_length):
            AA = protein['seq'][i]
            if i == 0:
                prev_AA = protein['seq'][sequence_length - 1]
                next_AA = protein['seq'][i + 1]
            elif i == sequence_length - 1:
                prev_AA = protein['seq'][i - 1]
                next_AA = protein['seq'][0]
            else:
                prev_AA = protein['seq'][i - 1]
                next_AA = protein['seq'][i + 1]

            # set the duration of the key (current AA) depending on the number of shared properties with the next AA

            shared_properties_current_next = len(set.intersection(aa_phy_chi[AA], aa_phy_chi[next_AA]))
            if shared_properties_current_next == 0:
                duration = 1
            elif shared_properties_current_next == 1:
                duration = 1.5
            elif shared_properties_current_next == 2:
                duration = 2
            else:
                duration = 4


            # set the chords depending on number of shared properties between current AA and the previous AA
            ## TODO: codage des accords
            # utiliser un accord avec 3 notes en sautant un note à chaque fois (ex sol-si-re). utiliser le tableau
            # de la TODO_list sals les altérations pour faire la liste des notes
            # si la note initiale est une altération, prendre le numéro de l'altération -1 et rajouter les notes de l'accord
            # et si les autres notes de l'accord sortent de l'index du tableau prendre celles de l'octave du dessous pour
            # celles qui dépassent (voir avec modulo de 7)
            shared_properties_current_previous = len(set.intersection(aa_phy_chi[AA], aa_phy_chi[prev_AA]))
            if shared_properties_current_previous == 2:
                # accord 2 notes
                pass
            elif shared_properties_current_previous == 3:
                # accord 3 notes
                pass
            elif shared_properties_current_previous == 4:
                # accord 4 notes
                pass
            else:
                pass

            # change the volume of each instrument depending on the structure
            if 'structure' in protein.keys():
                if i in protein['structure'].keys():
                    logger.debug('{}:'.format(protein['structure'][i]))
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

            if debug:
                logger.debug('\tposition: {:<5}\tAA: {:<5}\tnote: {:<20}\ttime: {:<5}\tduration: {:<5}'.format(i, AA, midi_keys[AA], time, duration))
            for channel_nbr in channels.keys():
                MyMIDI.addNote(track, channel=channel_nbr, pitch=midi_keys[AA], time=time, duration=duration, volume=channels[channel_nbr]['vol'])
                if debug:
                    logger.debug('\t\tchannel: {:<5}\tvolume: {:<5}'.format(channel_nbr, channels[channel_nbr]['vol']))

            time = time + duration
        MyMIDI.writeFile(midiFile)

    return midi_file_path
