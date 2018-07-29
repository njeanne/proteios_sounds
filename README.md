# Proteios sounds project

Listen the sound of the proteins

This project aims to transform the data collected by the biologists on proteins to a MIDI file.
To operate, the script needs an internet connection.

## data

The proteins data are retrieved from the [UniProt](https://www.uniprot.org/) database.

## Python libraries

- [Biopython](https://biopython.org/)
- [midiutil](https://pypi.org/project/MIDIUtil/)

## Usage
```
proteios_sounds.py -o <results_directory> -m <mode> [-t <tempo>] <uniprot_accession_number>
```

- `-o --out <results_directory>`: the path of the output directory where the MIDI file is produced.
- `-m --mode <mode>`: the mode to apply (major, mixolydian, dorian or blues).
- `-t --tempo <tempo>`: optional, the tempo in BPM. Must be an integer between 60 and 150, default is 100.
- `-i --instruments <INT INT INT>`: otional, 3 integers separated by spaces between 0 and 127 to set the instruments on the 3 channels. Default are 0, 42 and 65. See the [General MIDI patch numbers](http://www.pjb.com.au/muscript/gm.html#patch) for the correspondances between the integers and the MIDI instruments.
- `-d --debug`: create a log file.
- `<uniprot_accession_number>`: the uniprot accession number of the protein to create the MIDI file from the Uniprot entry. i.e: human interleukine 8 accession number is [P10145](https://www.uniprot.org/uniprot/P10145).
