#! /usr/bin/env python3

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

    print('[protein_movie.py]: {}'.format(duration_keys))
    width = 1280
    height = 720
    FPS = 24

    # set the video Codec
    fourcc = VideoWriter_fourcc(*'MP42')
    # set the video writer
    video = VideoWriter(path_movie, fourcc, float(FPS), (width, height))
    # create the frames per second
    # for idx, key_duration in enumerate(duration_keys):
    #     video.write()
