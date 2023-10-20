import pandas as pd
import numpy as np
import helpers as hp
from decimal import Decimal

df_acids = pd.DataFrame()

def get_hydrophobicity(peptide):
    h = Decimal(0)
    for a in peptide:
        h += Decimal(float(df_acids.loc[df_acids['abbreviation'] == a, 'hydrophobicity'].values[0]))
    return h

def get_charge(peptide):
    c = Decimal(0)
    for a in peptide:
        c += Decimal(float(df_acids.loc[df_acids['abbreviation'] == a, 'charge'].values[0]))
    return c

def get_volume(peptide):
    v = Decimal(0)
    for a in peptide:
        v += Decimal(float(df_acids.loc[df_acids['abbreviation'] == a, 'volume'].values[0]))
    return v

def get_polarity(peptide):
    p = Decimal(0)
    for a in peptide:
        p += Decimal(float(df_acids.loc[df_acids['abbreviation'] == a, 'polarity'].values[0]))
    return p

def get_acid(peptide, i):
    if len(peptide) < i:
        return None
    else:
        return peptide[i-1]
    
def get_end_acid(peptide, i):
    return peptide[i*-1]
    
def get_PKY(peptide):
    if 'PKY' in peptide:
        return 1
    else:
        return 0
    
def get_YKV(peptide):
    if 'YKV' in peptide:
        return 1
    else:
        return 0
    
def get_KYV(peptide):
    if 'KYV' in peptide:
        return 1
    else:
        return 0
    
def get_QNT(peptide):
    if 'QNT' in peptide:
        return 1
    else:
        return 0
    
def get_NTL(peptide):
    if 'NTL' in peptide:
        return 1
    else:
        return 0

# get acid stats at each position

def main():
    global df_acids
    df_acids = pd.read_csv('amino_acid_properties.csv', sep=',', header=0, names=('aminoAcid', 'abbreviation', 'hydrophobicity', 'volume', 'chemical', 'physicochemical',  'charge', 'polarity', 'hydrogenDonor'))
    df = pd.read_csv('mhc2023.txt', delim_whitespace=True, names=['bind', 'peptide'])

    df['length'] = df['peptide'].apply(lambda x: len(x))
    df['hydrophobicity'] = df['peptide'].apply(get_hydrophobicity)
    df['charge'] = df['peptide'].apply(get_charge)
    df['volume'] = df['peptide'].apply(get_volume)
    df['polarity'] = df['peptide'].apply(get_polarity)

    df['PKY'] = df['peptide'].apply(get_PKY)
    df['YKV'] = df['peptide'].apply(get_YKV)
    df['KYV'] = df['peptide'].apply(get_KYV)
    df['QNT'] = df['peptide'].apply(get_QNT)
    df['NTL'] = df['peptide'].apply(get_NTL)

    df['active_aa_count'] = df['peptide'].apply(hp.get_active_aa_count)
    df['binding_aa_count'] = df['peptide'].apply(hp.get_binding_aa_count)
    df['protein_pK'] = df['peptide'].apply(lambda x: hp.calculate_protein_pK(x, 7))
    df['longest_repeat'] = df['peptide'].apply(hp.longest_repeat_length)
    df['hydrogen_donors'] = df['peptide'].apply(hp.count_hydrogen_donors)
    df['hydrogen_acceptors'] = df['peptide'].apply(hp.count_hydrogen_acceptors)

    for i in range(1, 31):
        df['aa_' + str(i)] = df['peptide'].apply(lambda x: get_acid(x, i))

    for i in range(1, 7):
        df['aa_end_' + str(i)] = df['peptide'].apply(lambda x: get_end_acid(x, i))

    df.to_csv('data.csv', index=False, sep=',')

if __name__ == "__main__":
    main()


