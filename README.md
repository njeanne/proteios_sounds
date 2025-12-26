# Proteios sounds project

Listen the sound of the proteins.

This project aims to transform the data collected by the biologists on proteins to a MIDI file, then play it.
3 instruments are chosen, from the [MIDI instruments numbering](https://en.wikipedia.org/wiki/General_MIDI), to play 
the score created from the protein amino acids sequence.
To operate, the script needs an internet connection.

## Data

The protein data are retrieved from the [UniProt](https://www.uniprot.org/) database.

## Installation with Conda
 
A [Conda](https://docs.conda.io/projects/conda/en/stable/index.html) environment is provided in the 
`conda_env/proteios_sonds_env.yml` file.
The file contains all the dependencies to run the script except MuseScore which installation procedure is explained with
the link provided on the next section.

The conda environment is generated using the command:
```shell script
# create the environment
conda env create -f conda_env/proteios_sonds_env.yml

# activate the environment
conda activate proteios_sonds
```

## Manual installation

The script was tested with [Python 3.12](https://www.python.org/downloads/release/python-31212/).

### Python libraries

All Python library can be installed with [pip](https://pypi.org/), using `pip install <LIBRARY>`, see:

- [Biopython](https://biopython.org/)
- [midiutil](https://pypi.org/project/MIDIUtil/)
- [midi2audio](https://pypi.org/project/midi2audio/)
- [imageio](https://pypi.org/project/imageio/)
- [imageio-ffmpeg](https://pypi.org/project/imageio-ffmpeg/)

### External softwares

- [pymol](https://github.com/schrodinger/pymol-open-source)
- [musescore](https://musescore.org/en/download)
- [FluidSynth](http://www.fluidsynth.org/)
- [ffmpeg](https://ffmpeg.org/)

## Automatic downloads

If the file `FluidR3_GM.sf2` is not present in the `resources/FluidR3_GM` directory at the first execution of the 
script, it will be automatically downloaded.

## Usage
```
proteios_sounds.py -o <results_directory> [-t <tempo>] [-i <INT INT INT>] [-m] [-d] <uniprot_accession_number>
```

- `-o --out <results_directory>`: the path of the output directory where the MIDI file is produced.
- `-t --tempo <tempo>`: optional, the tempo in BPM. Must be an integer between 60 and 150, default is 100.
- `-i --instruments <INT INT INT>`: optional, 3 integers separated by spaces between 0 and 127 to set the instruments on the 3 channels. Default are 0, 42 and 65. See the [General MIDI patch numbers](https://www.perfessorbill.com/lyrics/gsinst.htm) for the correspondances between the integers and the MIDI instruments.
- `-m --musescore`: optional, use musescore to create the score corresponding to the MIDI file.
- `-d --debug`: create a log file.
- `<uniprot_accession_number>`: the uniprot accession number of the protein to create the MIDI file from the Uniprot entry. i.e: human interleukine 8 accession number is [P10145](https://www.uniprot.org/uniprot/P10145).
