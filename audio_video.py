#! /usr/bin/env python3

import os
import logging
import subprocess
from midi2audio import FluidSynth
import imageio

def convert_midi_to_flac(midi_path):
    '''Convert a MIDI file to a FLAC file.

    :param str midi_path: the path to the input MIDI file.
    :return: the path of the output FLAC file.
    :rtype: str
    '''
    # convert the MIDI file to audio file
    # using the default sound font in 44100 Hz sample rate

    # TODO: suppression du chemin FluidR3_GM.sf2 si pas de téléchargement dans le
    # projet à l'install dans ressources/
    sound_font_path = os.path.join(os.path.dirname(__file__),
                                   'ressources', 'FluidR3_GM.sf2')
    if not os.path.exists(sound_font_path):
        sound_font_path = '/usr/share/sounds/sf2/FluidR3_GM.sf2'
    logging.info('FluidSynth: soundfonts library path is {}'.format(sound_font_path))
    try:
        fs = FluidSynth(sound_font=sound_font_path)
        flac_path = '{}.flac'.format(os.path.splitext(midi_path)[0])
        fs.midi_to_audio(midi_path, flac_path)
    except Exception as ex:
        logging.error('FluidSynth: {}'.format(ex))
    logging.info('Audio file created: {}'.format(flac_path))
    return(flac_path)

def create_movie(path_movie, dir_frames, duration_keys, flac_path):
    '''Create the movie of the protein.

    :param str path_movie: the path where the movie will be created.
    :param str dir_frames: the path of the pdb frames directory.
    :param list of floats duration_keys: the duration keys list.
    :return: the movie path.
    :rtype: str
    '''

    # TODO: Video creation with ffmpeg to reduce the number of libraries
    # create the movie
    # create a dictionary of the frames directory with the index as keys
    frames_dict = {}
    for png in os.listdir(dir_frames):
        idx_frame_in_prot = os.path.splitext(png)[0].split('_')[1]
        frames_dict[idx_frame_in_prot] = os.path.join(dir_frames, png)

    # set the video writer
    FPS = 24
    path_tmp_movie = os.path.join(os.path.dirname(path_movie),
                                  'tmp_no_sound_movie.avi')

    video_writer = imageio.get_writer(path_tmp_movie, fps=FPS)

    # create the frames per second
    msg = 'Movie creation'
    logging.info(msg)
    print('{}, please wait..'.format(msg))

    for idx, key_duration in enumerate(duration_keys):
        nb_frames_on_key = round(float(FPS) * key_duration)
        msg = 'Frames for amino acid {}/{}\n\tkey duration: {} seconds\n\tnb of frames: {}'.format(idx + 1,
                                                                                                   len(duration_keys),
                                                                                                   key_duration,
                                                                                                   nb_frames_on_key)
        logging.info(msg)
        print(msg)
        if str(idx) in frames_dict:
            for _ in range(nb_frames_on_key):
                video_writer.append_data(imageio.imread(frames_dict[str(idx)]))
        else:
            for _ in range(nb_frames_on_key):
                video_writer.append_data(imageio.imread(frames_dict['no-idx']))
    video_writer.close()

    cmd_audio_video = 'ffmpeg -y -i {} -i {} -c copy -map 0:v:0 -map 1:a:0 {}'.format(path_tmp_movie,
                                                                                      flac_path,
                                                                                      path_movie)
    try:
        logging.info('ffmpeg: add soundtrack to the movie.')
        logging.info(cmd_audio_video)
        ffmpeg_process = subprocess.run(cmd_audio_video,
                                        shell=True,
                                        stdout=subprocess.PIPE,
                                        stderr=subprocess.PIPE)
        if ffmpeg_process.stdout:
            logging.info(ffmpeg_process.stdout)
        if ffmpeg_process.stderr:
            logging.warn(ffmpeg_process.stderr.decode('utf-8'))
        # remove tmp movie file (file without sound)
        os.remove(path_tmp_movie)
        msg = 'Movie file created: {}'.format(path_movie)
        print(msg)
        logging.info(msg)
    except Exception as ex:
        logging.error(ex)
