# Proteios sounds project

Listen the sound of the proteins

This project aims to transform the data collected by the biologists on proteins to a MIDI file.
To operate, the script needs an internet connection.

## data

The proteins data are retrieved from the [UniProt](https://www.uniprot.org/) database.

## Python libraries

- [Biopython](https://biopython.org/)
- [mudiutil](https://pypi.org/project/MIDIUtil/)

## Usage
```
proteios_sounds.py -o <results_directory> -m <mode> [-t <tempo>] <uniprot_accession_number>
```

- `-o --out <results_directory>`: the path of the output directory where the MIDI file is produced.
- `-m --mode <mode>`: the mode to apply (major, mixolydian, dorian or blues).
- `-t --tempo <tempo>`: optional, the tempo in BPM. Must be an integer between 60 and 150, default is 100.
- `<uniprot_accession_number>`: the uniprot accession number of the protein to transform to a MIDI file. i.e: human interleukine 8 accession number is [P10145](https://www.uniprot.org/uniprot/P10145).
