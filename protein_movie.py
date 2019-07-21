#! /usr/bin/env python3

import os
import subprocess
from midi2audio import FluidSynth
import imageio


def create_movie(path_movie, dir_frames, duration_keys, midi_path, logger):
    '''Create the movie of the protein.

    :param str path_movie: the path where the movie will be created.
    :param str dir_frames: the path of the pdb frames directory.
    :param list of floats duration_keys: the duration keys list.
    :param logger logger: the logger.
    :return: the movie path.
    :rtype: str
    '''

    #TODO: Video creation with ffmpeg to reduce the number of libraries
    # convert the MIDI file to audio file
    # using the default sound font in 44100 Hz sample rate
    sound_font_path = os.path.join(os.path.dirname(__file__),
                                   'ressources', 'FluidR3_GM.sf2')
    fs = FluidSynth(sound_font=sound_font_path)
    flac = '{}.flac'.format(os.path.splitext(midi_path)[0])
    fs.midi_to_audio(midi_path, flac)
    logger.info('Audio file created: {}'.format(flac))

    #TODO: Video creation with ffmpeg to reduce the number of libraries
    # create the movie
    # create a dictionary of the frames directory with the index as keys
    frames_dict = {}
    for png in os.listdir(dir_frames):
        idx_frame_in_prot = os.path.splitext(png)[0].split('_')[1]
        frames_dict[idx_frame_in_prot] = os.path.join(dir_frames, png)

    # set the video writer
    FPS = 24
    path_tmp_movie = os.path.join(os.path.dirname(path_movie),
                                  'tmp_movie.avi')
    video_writer = imageio.get_writer(path_tmp_movie, fps=FPS)

    # create the frames per second
    msg = 'Movie creation'
    logger.info(msg)
    print('{}, please wait..'.format(msg))

    for idx, key_duration in enumerate(duration_keys):
        nb_frames_on_key = round(float(FPS) * key_duration)
        msg = 'Frames for amino acid {}/{}\n\tkey duration: {} seconds\n\tnb of frames: {}'.format(idx + 1,
                                                                                                   len(duration_keys),
                                                                                                   key_duration,
                                                                                                   nb_frames_on_key)
        logger.info(msg)
        print(msg)
        if str(idx) in frames_dict:
            for _ in range(nb_frames_on_key):
                video_writer.append_data(imageio.imread(frames_dict[str(idx)]))
        else:
            for _ in range(nb_frames_on_key):
                video_writer.append_data(imageio.imread(frames_dict['no-idx']))
    video_writer.close()

    cmd_audio_video = 'ffmpeg -y -i {} -i {} -c copy -map 0:v:0 -map 1:a:0 {}'.format(path_tmp_movie,
                                                                                      flac,
                                                                                      path_movie)
    try:
        logger.info('ffmpeg: add soundtrack to the movie.')
        ffmpeg_process = subprocess.run(cmd_audio_video,
                                        shell=True,
                                        stdout=subprocess.PIPE,
                                        stderr=subprocess.PIPE)
        if ffmpeg_process.stdout:
            logger.info(ffmpeg_process.stdout)
        if ffmpeg_process.stderr:
            logger.warn(ffmpeg_process.stderr.decode('utf-8'))
        # remove tmp movie file (file without sound)
        os.remove(path_tmp_movie)
        msg = 'Movie file created: {}'.format(path_movie)
        print(msg)
        logger.info(msg)
    except Exception as ex:
        logger.error(ex)
