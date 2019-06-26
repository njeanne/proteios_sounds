# Proteios sounds project

Listen the sound of the proteins

This project aims to transform the data collected by the biologists on proteins to a MIDI file.
To operate, the script needs an internet connection.

## data

The proteins data are retrieved from the [UniProt](https://www.uniprot.org/) database.

## Python libraries

All Python library can be installed using `pip install <LIBRARY>`, see:

- [Biopython](https://biopython.org/)
- [midiutil](https://pypi.org/project/MIDIUtil/)
- [midi2audio](https://pypi.org/project/midi2audio/)
- [imageio](https://pypi.org/project/imageio/)

## External softwares

- [pymol](https://github.com/schrodinger/pymol-open-source)
- [musescore](https://musescore.org/en/download)
- [FluidSynth](http://www.fluidsynth.org/)

## Usage
```
proteios_sounds.py -o <results_directory> [-t <tempo>] [-i <INT INT INT>] [-m] [-d] <uniprot_accession_number>
```

- `-o --out <results_directory>`: the path of the output directory where the MIDI file is produced.
- `-t --tempo <tempo>`: optional, the tempo in BPM. Must be an integer between 60 and 150, default is 100.
- `-i --instruments <INT INT INT>`: optional, 3 integers separated by spaces between 0 and 127 to set the instruments on the 3 channels. Default are 0, 42 and 65. See the [General MIDI patch numbers](http://www.pjb.com.au/muscript/gm.html#patch) for the correspondances between the integers and the MIDI instruments.
- `-m --musescore`: optional, use musescore to create the score corresponding to the MIDI file.
- `-d --debug`: create a log file.
- `<uniprot_accession_number>`: the uniprot accession number of the protein to create the MIDI file from the Uniprot entry. i.e: human interleukine 8 accession number is [P10145](https://www.uniprot.org/uniprot/P10145).
