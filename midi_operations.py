#! /usr/bin/env python3

import os
from midiutil import MIDIFile

def create_chord(pitch_list, keys_in_chord, idx_key, KEYS_OCTAVE_ONLY):
    '''
    :param list pitch_list: the list of MIDI keys integers
    :param int keys_in_chord: the number of keys in the chord
    :param int idx_key: the index in the list of KEYS_OCTAVE_ALTERATIONS of the first key of the chord
    :param list KEYS_OCTAVE_ONLY: the list of keys MIDI without the alterations
    :return: the list of MIDI keys for the chord
    :rtype: list of integers
    '''
    added_keys = 1
    while added_keys < keys_in_chord:
        idx_key = idx_key + 2
        # if index out of bound, substract the length of the list
        # to change the index
        if idx_key >= len(KEYS_OCTAVE_ONLY):
            idx_key = idx_key - len(KEYS_OCTAVE_ONLY)
        pitch_list.append(KEYS_OCTAVE_ONLY[idx_key])
        added_keys += 1
    return pitch_list

def create_midi(path_midi, protein, midi_keys, tempo, instrus, aa_phy_chi,
                logger, debug=False):
    '''
    Creates the MIDI file from the protein data.
    :param str path_midi: the path to the MIDI file.
    :param dict protein: the dictionary descri the protein.
    :param dict midi_keys: the dictionary of the keys
    :param list initial_keys: the list of the initial keys
    :param int tempo: the tempo in BPM
    :param list instrus: the list of the MIDI instrument numbers
    :param dict aa_phy_chi: dictionary of the physico-chemical attributes of the amino acids
    :param logger logger: the logger
    :param boolean debug: the enable debug mode
    :return: the list of keys durations
    :rtype: list of floats
    '''

    # octaves DO, RE, MI, FA, SOL, LA, SI. 2 octaves and 1 more SOL,
    # the remaining 7 keys are altérations (#)
    KEYS_OCTAVE_ALTERATIONS = [48, 50, 52, 53, 55, 57, 59, 60, 62, 64, 65, 67,
                               69, 71, 72, 54, 66, 49, 61, 56, 68, 51]
    KEYS_OCTAVE_ONLY = KEYS_OCTAVE_ALTERATIONS[:15]

    with open(path_midi, 'wb') as  midiFile:

        track = 0
        time = 0   # In beats

        # a channel is defined by an instrument nbr and
        # a volume (0-127, as per the MIDI standard,
        # see: http://www.pjb.com.au/muscript/gm.html)
        # channel 9 is for percussions
        # (see https://pjb.com.au/muscript/gm.html#perc)
        channels = {0: {'instrument': instrus[0], 'vol': 100},
                    1: {'instrument': instrus[1], 'vol': 40},
                    2: {'instrument': instrus[2], 'vol': 60}}

        if debug:
            logger.info('Instrument number by channel, see: http://www.pjb.com.au/muscript/gm.html for instruments number correspondance:')
            for channel_nb in channels:
                logger.info('\tchannel {}: instrument {}'.format(channel_nb,
                                                                 channels[channel_nb]['instrument']))

        MyMIDI = MIDIFile(numTracks=1, adjust_origin=False) # One track, defaults to format 1 (tempo track automatically created)
        MyMIDI.addTempo(track, time, tempo)
        # add the channels (1 per instrument)
        for channel_nbr in channels:
            MyMIDI.addProgramChange(track, channel=channel_nbr,
                                    time=time,
                                    program=channels[channel_nbr]['instrument'])

        sequence_length = len(protein['seq'])

        durations_list = []

        for i in range(0, sequence_length):
            AA = protein['seq'][i]
            pitch_list = [midi_keys[AA]]

            if i == 0:
                prev_AA = protein['seq'][sequence_length - 1]
                next_AA = protein['seq'][i + 1]
            elif i == sequence_length - 1:
                prev_AA = protein['seq'][i - 1]
                next_AA = protein['seq'][0]
            else:
                prev_AA = protein['seq'][i - 1]
                next_AA = protein['seq'][i + 1]

            # set the duration of the key (current AA) depending on the number
            # of shared properties with the next AA
            if AA == 'X' or next_AA == 'X': # non determined AA
                shared_properties_current_next = 0
            else:
                shared_properties_current_next = len(set.intersection(aa_phy_chi[AA],
                                                                      aa_phy_chi[next_AA]))

            if shared_properties_current_next == 0:
                duration = 1
            elif shared_properties_current_next == 1:
                duration = 1.5
            elif shared_properties_current_next == 2:
                duration = 2
            else:
                duration = 4
            durations_list.append(float(duration))

            # set the chords depending on number of shared properties between
            # current AA and the previous AA
            if AA == 'X' or prev_AA == 'X': # non determined AA
                shared_properties_current_previous = 0
            else:
                shared_properties_current_previous = len(set.intersection(aa_phy_chi[AA],
                                                                          aa_phy_chi[prev_AA]))

            if shared_properties_current_previous == 2:
                # 2 keys chord
                keys_in_chord = 2
                idx_key = KEYS_OCTAVE_ALTERATIONS.index(midi_keys[AA])
                pitch_list = create_chord(pitch_list, keys_in_chord, idx_key,
                                          KEYS_OCTAVE_ONLY)
            elif shared_properties_current_previous == 3:
                # 3 keys chord
                keys_in_chord = 3
                idx_key = KEYS_OCTAVE_ALTERATIONS.index(midi_keys[AA])
                pitch_list = create_chord(pitch_list, keys_in_chord, idx_key,
                                          KEYS_OCTAVE_ONLY)
            elif shared_properties_current_previous >= 4:
                # 4 keys chord
                keys_in_chord = 4
                idx_key = KEYS_OCTAVE_ALTERATIONS.index(midi_keys[AA])
                pitch_list = create_chord(pitch_list, keys_in_chord, idx_key,
                                          KEYS_OCTAVE_ONLY)

            # change the volume of each instrument depending on the structure
            if 'structure' in protein.keys():
                if i in protein['structure'].keys():
                    logger.debug('{}: {}'.format(protein['structure'][i], i))
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
                logger.debug('position: {}'.format(i))
                logger.debug('AA: {}'.format(AA))
                logger.debug('pitch: {}'.format(pitch_list))
                logger.debug('time: {}'.format(time))
                logger.debug('duration: {}'.format(duration))
            for channel_nbr in channels:
                for pitch in pitch_list:
                    MyMIDI.addNote(track,
                                   channel=channel_nbr,
                                   pitch=pitch,
                                   time=time,
                                   duration=duration,
                                   volume=channels[channel_nbr]['vol'])

            time = time + duration
        MyMIDI.writeFile(midiFile)

    return durations_list
