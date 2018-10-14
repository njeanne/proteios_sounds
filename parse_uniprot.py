#! /usr/bin/env python3

import re
from Bio import ExPASy
from Bio import SwissProt

def parse_entry(uniprot_accession_number):

    # Parse the Uniprot entry
    try:
        handle = ExPASy.get_sprot_raw(uniprot_accession_number)
    except urllib.error.HTTPError as err:
        msg = '{}{} accession number does not exists in UniProt database, check https://www.uniprot.org/ for the correct accession number.'.format(err, args.uniprot_accession_number)
        print('ERROR: {}'.format(msg))
        if args.debug:
            logger.error(msg)
        sys.exit(0)

    uniprot = SwissProt.read(handle)

    ### TOREMOVE
    # print(uniprot.sequence)
    # for item in uniprot.features:
    #     print(item)
    ##############

    # get the organism and entry name in UniProt
    organism = uniprot.organism
    if '(' in organism:
         organism = organism.split('(')[1].split(')')[0]
    organism = organism.rstrip('.,;:').replace(' ', '_')
    entry_name = uniprot.entry_name
    uniprot_sequence = uniprot.sequence

    # create a dictionary for the protein
    protein = {'seq': uniprot_sequence, 'organism': organism, 'entry_name': entry_name}


    # create a dictionary to get the structures positions => type
    structures = {}
    for v in uniprot.features:
        print('{}'.format(v))
    # parse the UniProt Features (see https://www.uniprot.org/help/?query=*&fil=category%3A%22PTM+processing%22+AND+section%3Amanual&columns=title)
    for feature in uniprot.features:
        # strucutre information
        if feature[0] == 'HELIX' or feature[0] == 'STRAND' or feature[0] == 'TURN':
            structures[feature[1]] = {'type': feature[0], 'end': feature[2]}
        # signal peptide
        if feature[0] == 'SIGNAL':
            if 'signal_peptide' in protein.keys():
                protein['signal_peptide'].append((feature[1], feature[2]))
            else:
                protein['signal_peptide']= [(feature[1], feature[2])]
        # Modified residue: Modified residues excluding lipids, glycans and protein crosslinks
        if feature[0] == 'MOD_RES':
            modif = re.split(' {', feature[3])[0].rstrip('.;,')
            if not 'modified_residue' in protein.keys():
                protein['modified_residue'] = {modif: [(feature[1], feature[2])]}
            else:
                if modif in protein['modified_residue'].keys():
                    protein['modified_residue'][modif].append((feature[1], feature[2]))
                else:
                    protein['modified_residue'][modif] = [(feature[1], feature[2])]
        # Glycosylation: Covalently attached glycan group(s)
        if feature[0] == 'CARBOHYD':
            glyco = re.split(' ', feature[3])[0].rstrip('.;,')
            if not 'glycosylation' in protein.keys():
                protein['glycosylation'] = {glyco: [(feature[1], feature[2])]}
            else:
                if glyco in protein['glycosylation'].keys():
                    protein['glycosylation'][glyco].append((feature[1], feature[2]))
                else:
                    protein['glycosylation'][glyco] = [(feature[1], feature[2])]
        # Site: Any interesting single amino acid site on the sequence
        if feature[0] == 'SITE':
            site = re.split(' {', feature[3])[0].rstrip('.;,')
            if site in protein['site'].keys():
                protein['site'][site].append((feature[1], feature[2]))
            else:
                protein['site'][site] = [(feature[1], feature[2])]
        # Propeptide: Part of a protein that is cleaved during maturation or activation
        if feature[0] == 'PROPEP':
            propep = re.split(' {', feature[3])[0].rstrip('.;,')
            if propep in protein['propeptide'].keys():
                protein['propeptide'][propep].append((feature[1], feature[2]))
            else:
                protein['propeptide'][propep] = [(feature[1], feature[2])]
        # Disulfid bond, get the two linked positions on the same chain (no interchain)
        if feature[0] == 'DISULFID':
            if not feature[3].startswith('Interchain'):
                if 'disulfid' in protein.keys():
                    protein['disulfid'].append((feature[1], feature[2]))
                else:
                    protein['disulfid'] = [(feature[1], feature[2])]


    # update protein with the structures
    sequence_length = len(protein['seq'])
    structure_last_position = 0
    for position in sorted(structures.keys()):
        if structure_last_position == 0 and position == 1:
            if 'structure' in protein.keys():
                protein['structure'][position] = structures[position]['type']
            else:
                protein['structure'] = {position: structures[position]['type']}
            structure_last_position = structures[position]['end']
        elif position >= structure_last_position + 1:
            if 'structure' in protein.keys():
                protein['structure'][structure_last_position + 1] = 'FREE'
            else:
                protein['structure'] = {structure_last_position + 1: 'FREE'}
            protein['structure'][position] = structures[position]['type']
            structure_last_position = structures[position]['end']
    if structure_last_position < sequence_length:
        if 'structure' in protein.keys():
            protein['structure'][structure_last_position + 1] = 'FREE'
        else:
            protein['structure'] = {structure_last_position + 1: 'FREE'}

    ### TOREMOVE
    for k,v in protein.items():
        if k == 'modified_residue':
            print(k)
            for k2,v2 in v.items():
                print('\t{}'.format(k2))
                for pos in v2:
                    print('\t\t{}'.format(pos))
        elif k == 'glycosylation':
            print(k)
            for k2,v2 in v.items():
                print('\t{}'.format(k2))
                for pos in v2:
                    print('\t\t{}'.format(pos))
        elif k == 'site':
            print(k)
            for k2,v2 in v.items():
                print('\t{}'.format(k2))
                for pos in v2:
                    print('\t\t{}'.format(pos))
        elif k == 'propeptide':
            print(k)
            for k2,v2 in v.items():
                print('\t{}'.format(k2))
                for pos in v2:
                    print('\t\t{}'.format(pos))
        elif k == 'disulfid':
            print(k)
            for bond in v:
                print('\t{} disulfid bond to {}'.format(bond[0], bond[1]))
        elif k == 'signal_peptide':
            print(k)
            for signal_peptide in v:
                print('\tsignal peptide from {} to {}'.format(signal_peptide[0], signal_peptide[1]))
    #################

    return protein
