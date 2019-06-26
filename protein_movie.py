#! /usr/bin/env python3

import os
from midi2audio import FluidSynth
import imageio



def create_movie(path_movie, dir_frames, duration_keys, logger):
    '''Create the movie of the protein.

    :param str path_movie: the path where the movie will be created.
    :param str dir_frames: the path of the pdb frames directory.
    :param list of floats duration_keys: the duration keys list.
    :param logger logger: the logger.
    :return: the movie path.
    :rtype: str
    '''

    # using the default sound font in 44100 Hz sample rate
    fs = FluidSynth()
    ######## TO CHANGE
    midi = '/home/nico/fablab/proteios_sounds/sounds_proteins/P10145_IL8_HUMAN_Homo_sapiens_100bpm_instrus-0-42-65.midi'
    flac = '/home/nico/fablab/proteios_sounds/proteios_sounds/P10145.flac'
    fs.midi_to_audio(midi, flac)

    # download ffmpeg if necessary
    imageio.plugins.ffmpeg.download()

    # create a dictionary of the frames directory with the index as keys
    frames_dict = {}
    for png in os.listdir(dir_frames):
        idx_frame_in_prot = os.path.splitext(png)[0].split('_')[1]
        frames_dict[idx_frame_in_prot] = os.path.join(dir_frames, png)

    # set the video parameters
    FPS = 24

    # set the video writer
    video_writer = imageio.get_writer(path_movie, fps=FPS)

    # create the frames per second
    print('Creating the movie file, please wait..')
    for idx, key_duration in enumerate(duration_keys):
        print('\tFrames for amino acid {}/{}'.format(idx + 1,
                                                     len(duration_keys)))
        nb_frames_on_key = int(FPS * key_duration)
        if str(idx) in frames_dict:
            for _ in range(nb_frames_on_key):
                video_writer.append_data(imageio.imread(frames_dict[str(idx)]))
        else:
            for _ in range(nb_frames_on_key):
                video_writer.append_data(imageio.imread(frames_dict['no-idx']))
    video_writer.close()
