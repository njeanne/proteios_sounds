#! /usr/bin/env python3

import os
from midi2audio import FluidSynth
import imageio
import subprocess



def create_movie(path_movie, dir_frames, duration_keys, midi_path, logger):
    '''Create the movie of the protein.

    :param str path_movie: the path where the movie will be created.
    :param str dir_frames: the path of the pdb frames directory.
    :param list of floats duration_keys: the duration keys list.
    :param logger logger: the logger.
    :return: the movie path.
    :rtype: str
    '''

    # convert the MIDI file to audio file
    # using the default sound font in 44100 Hz sample rate
    fs = FluidSynth(sound_font='FluidR3_GM.sf2')
    flac = '{}.flac'.format(os.path.splitext(midi_path)[0])
    fs.midi_to_audio(midi_path, flac)

    # create the movie
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

    cmd_audio_video = 'ffmpeg -i {0} -i {1} -c copy -map 0:v:0 -map 1:a:0 {0}'.format(path_movie,
                                                                                           flac)
    subprocess.run(cmd_audio_video, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
