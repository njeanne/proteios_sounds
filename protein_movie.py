#! /usr/bin/env python3

import logging
import os
import subprocess
import urllib.request
import zipfile

from midi2audio import FluidSynth
import imageio


def convert_midi_to_aac(midi):
    """
    Convert the MIDI file to an AAC file.

    :param midi: the MIDI file path.
    :type midi: str
    :return: the AAC file path.
    :rtype: str
    """
    fluidr3_gm_dir = os.path.join(os.path.dirname(__file__), "resources", "FluidR3_GM")
    sound_font_path = os.path.join(fluidr3_gm_dir, "FluidR3_GM.sf2")

    # if the FluidR3_GM.sf2 does not exist, download it
    if not os.path.exists(sound_font_path):
        os.makedirs(fluidr3_gm_dir, exist_ok=True)
        logging.info("\tMissing FluidR3_GM.sf2 file in the resources directory, downloading it..")
        fluidr3_gm_zip = os.path.join(fluidr3_gm_dir, "FluidR3_GM.zip")
        filehandle, _ = urllib.request.urlretrieve("https://filedn.eu/lf4Nj1iOnB8JeCa95AXHCxy/FluidR3_GM.zip",
                                                   fluidr3_gm_zip)
        with zipfile.ZipFile(fluidr3_gm_zip, "r") as zip_ref:
            zip_ref.extractall(fluidr3_gm_dir)
        os.remove(fluidr3_gm_zip)

    fs = FluidSynth(sound_font=sound_font_path, sample_rate=44100)
    wav = f"{os.path.splitext(midi)[0]}.wav"
    fs.midi_to_audio(midi, wav)
    aac = f"{os.path.splitext(midi)[0]}.m4a"
    subprocess.run(["ffmpeg", "-y", "-i", wav, "-c:a", "aac", "-b:a", "192k", aac])
    os.remove(wav)
    logging.info(f"\tAudio file created: {aac}")
    return aac


def movie_creation(aac, dir_frames, duration_keys, path_movie):
    """
    Create the movie with the pictures without sound, and after add the soundtrack.

    :param aac: the AAC sound file path.
    :type aac: str
    :param dir_frames: the frames' directory.
    :type dir_frames: str
    :param duration_keys: the keys' durations.
    :type duration_keys: list
    :param path_movie: the movie's path.
    :type path_movie: str
    """
    frames_dict = {}
    for png in os.listdir(dir_frames):
        idx_frame_in_prot = os.path.splitext(png)[0].split("_")[1]
        frames_dict[idx_frame_in_prot] = os.path.join(dir_frames, png)
    # set the video writer
    fps = 24
    path_tmp_movie = os.path.join(os.path.dirname(path_movie), "tmp_movie.avi")
    video_writer = imageio.get_writer(path_tmp_movie, fps=fps, codec="libx264")
    # create the frames per second
    logging.info("Movie creation, please wait..")
    for idx, key_duration in enumerate(duration_keys):
        nb_frames_on_key = round(float(fps) * key_duration)
        logging.debug(f"\tFrames for amino acid {idx + 1}/{len(duration_keys)}")
        logging.debug(f"\t\tkey duration:\t{key_duration} seconds.")
        logging.debug(f"\t\tframes count:\t{nb_frames_on_key}")
        if str(idx) in frames_dict:
            for _ in range(nb_frames_on_key):
                video_writer.append_data(imageio.imread(frames_dict[str(idx)]))
        else:
            for _ in range(nb_frames_on_key):
                video_writer.append_data(imageio.imread(frames_dict["no-idx"]))
    video_writer.close()
    cmd_audio_video = (f"ffmpeg -y -fflags +genpts -i {path_tmp_movie} -i {aac} -c:v copy -c:a copy -shortest "
                       f"{path_movie}")
    try:
        logging.info("\tffmpeg: add soundtrack to the movie.")
        ffmpeg_process = subprocess.run(cmd_audio_video, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if ffmpeg_process.stdout:
            logging.info(ffmpeg_process.stdout)
        if ffmpeg_process.stderr:
            logging.warning(ffmpeg_process.stderr.decode("utf-8"))
        # remove tmp movie file (file without sounds)
        os.remove(path_tmp_movie)
        logging.info(f"\tMovie file created: {path_movie}")
    except Exception as ex:
        logging.error(ex)


def create_movie(path_of_the_movie, frames_directory, durations_of_the_keys, midi_path):
    """
    Create the movie of the protein.

    :param path_of_the_movie: the path where the movie will be created.
    :type path_of_the_movie: str
    :param frames_directory: the path of the pdb frames directory.
    :type frames_directory: str
    :param durations_of_the_keys: the duration keys list.
    :type durations_of_the_keys: list
    :param midi_path: the midi file path.
    :type midi_path: str
    """
    logging.info("Create the movie:")
    # sound conversion
    aac_path = convert_midi_to_aac(midi_path)
    # movie creation
    movie_creation(aac_path, frames_directory, durations_of_the_keys, path_of_the_movie)


