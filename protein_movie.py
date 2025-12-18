#! /usr/bin/env python3

import logging
import os
import subprocess
from midi2audio import FluidSynth
import imageio


def create_movie(path_movie, dir_frames, duration_keys, midi_path):
    """
    Create the movie of the protein.

    :param path_movie: the path where the movie will be created.
    :type path_movie: str
    :param dir_frames: the path of the pdb frames directory.
    :type dir_frames: str
    :param duration_keys: the duration keys list.
    :type duration_keys: list
    :param midi_path: the midi file path.
    :type midi_path: str
    :return: the movie path.
    :rtype: str
    """

    #TODO: Video creation with ffmpeg to reduce the number of libraries
    # convert the MIDI file to audio file
    # using the default sound font in 44100 Hz sample rate
    sound_font_path = os.path.join(os.path.dirname(__file__), "ressources", "FluidR3_GM.sf2")
    fs = FluidSynth(sound_font=sound_font_path)
    flac = f"{os.path.splitext(midi_path)[0]}.flac"
    fs.midi_to_audio(midi_path, flac)
    logging.info(f"Audio file created: {flac}")

    #TODO: Video creation with ffmpeg to reduce the number of libraries
    # create the movie
    # create a dictionary of the frames directory with the index as keys
    frames_dict = {}
    for png in os.listdir(dir_frames):
        idx_frame_in_prot = os.path.splitext(png)[0].split("_")[1]
        frames_dict[idx_frame_in_prot] = os.path.join(dir_frames, png)

    # set the video writer
    fps = 24
    path_tmp_movie = os.path.join(os.path.dirname(path_movie), "tmp_movie.avi")
    video_writer = imageio.get_writer(path_tmp_movie, fps=fps)

    # create the frames per second
    logging.info("Movie creation, please wait..")

    for idx, key_duration in enumerate(duration_keys):
        nb_frames_on_key = round(float(fps) * key_duration)
        logging.info(f"Frames for amino acid {idx + 1}/{len(duration_keys)}")
        logging.info(f"\tkey duration: {key_duration} seconds\n\tnb of frames: {nb_frames_on_key}")
        if str(idx) in frames_dict:
            for _ in range(nb_frames_on_key):
                video_writer.append_data(imageio.imread(frames_dict[str(idx)]))
        else:
            for _ in range(nb_frames_on_key):
                video_writer.append_data(imageio.imread(frames_dict["no-idx"]))
    video_writer.close()

    cmd_audio_video = f"ffmpeg -y -i {path_tmp_movie} -i {flac} -c copy -map 0:v:0 -map 1:a:0 {path_movie}"
    try:
        logging.info("ffmpeg: add soundtrack to the movie.")
        ffmpeg_process = subprocess.run(cmd_audio_video, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if ffmpeg_process.stdout:
            logging.info(ffmpeg_process.stdout)
        if ffmpeg_process.stderr:
            logging.warn(ffmpeg_process.stderr.decode("utf-8"))
        # remove tmp movie file (file without a sound)
        os.remove(path_tmp_movie)
        logging.info(f"Movie file created: {path_movie}")
    except Exception as ex:
        logging.error(ex)
