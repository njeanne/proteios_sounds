#! /usr/bin/env python3

import logging
from midiutil import MIDIFile


def create_chord(pitch_list, keys_in_chord, idx_key, keys_octave_only):
    """
    :param pitch_list: the list of MIDI keys integers.
    :type pitch_list: list
    :param keys_in_chord: the number of keys in the chord.
    :type keys_in_chord: int
    :param idx_key: the chord's first key index in the list of keys_octave_alterations.
    :type idx_key: int
    :param keys_octave_only: the list of keys' MIDI without the alterations.
    :type keys_octave_only: list
    :return: the list of MIDI keys for the chord.
    :rtype: list of integers
    """
    added_keys = 1
    while added_keys < keys_in_chord:
        idx_key = idx_key + 2
        # if index out of bound, subtract the length of the list to change the index
        if idx_key >= len(keys_octave_only):
            idx_key = idx_key - len(keys_octave_only)
        pitch_list.append(keys_octave_only[idx_key])
        added_keys += 1
    return pitch_list


def create_midi(path_midi, protein, midi_keys, tempo, instruments, aa_phy_chi, debug=False):
    """
    Creates the MIDI file from the protein data.
    :param path_midi: the path to the MIDI file.
    :type path_midi: str
    :param protein: the dictionary describing the protein.
    :type protein: dict
    :param midi_keys: the dictionary of the keys.
    :type midi_keys: dict
    :param tempo: the tempo in BPM.
    :type tempo: int
    :param instruments: the MIDI instrument numbers list.
    :type instruments: list
    :param aa_phy_chi: amino acids physico-chemical attributes dictionary.
    :type aa_phy_chi: dict
    :param debug: the enable debug mode.
    :type debug: boolean
    :return: the list of keys durations.
    :rtype: list of floats.
    """
    logging.info("MIDI file creation:")

    # octaves DO, RE, MI, FA, SOL, LA, SI. 2 octaves and 1 more SOL,
    # the remaining 7 keys are altérations (#)
    keys_octave_alterations = [48, 50, 52, 53, 55, 57, 59, 60, 62, 64, 65, 67, 69, 71, 72, 54, 66, 49, 61, 56, 68, 51]
    keys_octave_only = keys_octave_alterations[:15]

    with open(path_midi, "wb") as midiFile:
        track = 0
        time = 0   # In beats

        # a channel is defined by an instrument nbr and
        # a volume (0-127, as per the MIDI standard,
        # see: http://www.pjb.com.au/muscript/gm.html)
        # channel 9 is for percussions
        # (see https://pjb.com.au/muscript/gm.html#perc)
        channels = {0: {"instrument": instruments[0], "vol": 100},
                    1: {"instrument": instruments[1], "vol": 40},
                    2: {"instrument": instruments[2], "vol": 60}}

        if debug:
            logging.debug("Instrument number by channel, see: http://www.pjb.com.au/muscript/gm.html for instruments "
                        "number correspondance:")
            for channel_nb in channels:
                logging.debug(f"\tchannel {channel_nb}: instrument {channels[channel_nb]['instrument']}")

        # One track, defaults to format 1 (tempo track automatically created)
        MyMIDI = MIDIFile(numTracks=1, adjust_origin=False)
        MyMIDI.addTempo(track, time, tempo)
        # add the channels (1 per instrument)
        for channel_nbr in channels:
            MyMIDI.addProgramChange(track, channel=channel_nbr, time=time, program=channels[channel_nbr]["instrument"])

        sequence_length = len(protein["seq"])

        durations_list = []

        for i in range(0, sequence_length):
            aa = protein["seq"][i]
            pitch_list = [midi_keys[aa]]

            if i == 0:
                prev_aa = protein["seq"][sequence_length - 1]
                next_aa = protein["seq"][i + 1]
            elif i == sequence_length - 1:
                prev_aa = protein["seq"][i - 1]
                next_aa = protein["seq"][0]
            else:
                prev_aa = protein["seq"][i - 1]
                next_aa = protein["seq"][i + 1]

            # set the duration of the key (current AA) depending on the number
            # of shared properties with the next AA
            if aa == "X" or next_aa == "X": # non determined AA
                shared_properties_current_next = 0
            else:
                shared_properties_current_next = len(set.intersection(aa_phy_chi[aa], aa_phy_chi[next_aa]))

            if shared_properties_current_next == 0:
                duration = 1
            elif shared_properties_current_next == 1:
                duration = 1.5
            elif shared_properties_current_next == 2:
                duration = 2
            else:
                duration = 4
            # add each duration adjusted with the tempo
            durations_list.append(float(duration) * (60 / tempo))

            # set the chords depending on the number of shared properties between the current AA and the previous AA
            if aa == "X" or prev_aa == "X":  # non determined AA
                shared_properties_current_previous = 0
            else:
                shared_properties_current_previous = len(set.intersection(aa_phy_chi[aa], aa_phy_chi[prev_aa]))

            if shared_properties_current_previous == 2:
                # 2 keys chord
                keys_in_chord = 2
                idx_key = keys_octave_alterations.index(midi_keys[aa])
                pitch_list = create_chord(pitch_list, keys_in_chord, idx_key, keys_octave_only)
            elif shared_properties_current_previous == 3:
                # 3 keys chord
                keys_in_chord = 3
                idx_key = keys_octave_alterations.index(midi_keys[aa])
                pitch_list = create_chord(pitch_list, keys_in_chord, idx_key, keys_octave_only)
            elif shared_properties_current_previous >= 4:
                # 4 keys chord
                keys_in_chord = 4
                idx_key = keys_octave_alterations.index(midi_keys[aa])
                pitch_list = create_chord(pitch_list, keys_in_chord, idx_key, keys_octave_only)

            # change the volume of each instrument depending on the structure
            if "structure" in protein.keys():
                if i in protein["structure"].keys():
                    logging.debug(f"{protein['structure'][i]}: {i}")
                    if protein["structure"][i] == "HELIX":
                        channels[0]["vol"] = 40
                        channels[1]["vol"] = 100
                        channels[2]["vol"] = 60
                    elif protein["structure"][i] == "STRAND":
                        channels[0]["vol"] = 60
                        channels[1]["vol"] = 40
                        channels[2]["vol"] = 100
                    elif protein["structure"][i] == "TURN":
                        channels[0]["vol"] = 100
                        channels[1]["vol"] = 40
                        channels[2]["vol"] = 60
                    else:
                        channels[0]["vol"] = 100
                        channels[1]["vol"] = 60
                        channels[2]["vol"] = 40

            if debug:
                logging.debug(f"position: {i}")
                logging.debug(f"AA: {aa}")
                logging.debug(f"pitch: {pitch_list}")
                logging.debug(f"time: {time}")
                logging.debug(f"duration: {duration}")
            for channel_nbr in channels:
                for pitch in pitch_list:
                    MyMIDI.addNote(track, channel=channel_nbr, pitch=pitch, time=time, duration=duration,
                                   volume=channels[channel_nbr]["vol"])

            time = time + duration
        MyMIDI.writeFile(midiFile)
    logging.info(f"\tMIDI file for {protein['entry_name']} {protein['organism']} created: {path_midi}")

    return durations_list
