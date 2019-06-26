#! /usr/bin/env python3

import cv2
from cv2 import VideoWriter, VideoWriter_fourcc


def create_movie(frames_dir, logger):
    '''Create the movie of the protein.

    :param str frames_dir: the path of the pdb frames directory.
    :param logger logger: the logger.
    :return: the movie path.
    :rtype: str
    '''

    
