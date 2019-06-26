#! /usr/bin/env python3

import os
import numpy as np
import cv2
from cv2 import VideoWriter, VideoWriter_fourcc

def create_movie(path_movie, dir_frames, duration_keys, logger):
    '''Create the movie of the protein.

    :param str path_movie: the path where the movie will be created.
    :param str dir_frames: the path of the pdb frames directory.
    :param list of floats duration_keys: the duration keys list.
    :param logger logger: the logger.
    :return: the movie path.
    :rtype: str
    '''
    # create a dictionary of the frames directory with the index as keys
    frames_dict = {}
    for png in os.listdir(dir_frames):
        idx_frame_in_prot = os.path.splitext(png)[0].split('_')[1]
        frames_dict[idx_frame_in_prot] = cv2.imread(os.path.join(dir_frames,
                                                                 png),
                                                    1)
    print('[protein_movie.py]: {}'.format(frames_dict))
    # set the video parameters
    width = 1280
    height = 720
    FPS = 24

    # set the video Codec
    fourcc = VideoWriter_fourcc(*'MP42')
    # set the video writer
    video = VideoWriter(path_movie, fourcc, float(FPS), (width, height))
    # create the frames per second
    for idx, key_duration in enumerate(duration_keys):
        nb_frames_on_key = int(FPS * key_duration)
        print('[protein_movie.py]: key duration: {}\t nb frames: {}'.format(key_duration,
                                                        nb_frames_on_key))
        if str(idx) in frames_dict:
            for _ in range(nb_frames_on_key):
                video.write(frames_dict[str(idx)])
        else:
            for _ in range(nb_frames_on_key):
                video.write(frames_dict['no-idx'])


    video.release()
