import pandas as pd
import sys

def to_numbers(x):
    result = []
    for part in x.split(';'):
        a, b = part.split('...')
        a, b = int(a), int(b)
        result.extend(range(a, b + 1))
    return result


def to_peri_numbers(x, y):
    result = []
    filtered = []
    for part in x.split(';'):
        a, b = part.split('...')
        a, b = int(a), int(b)
        result.extend(range(a - 8, a))
        result.extend(range(b + 1, b + 9))
    for i in result:
        if 0 < i <= len(y) and i not in to_numbers(x):
            filtered.append(i)
    return list(set(filtered))


def to_tele_numbers(x, y):
    result = []
    filtered = []
    result.extend(range(1, len(y) + 1))
    for i in result:
        if i not in to_numbers(x) and i not in to_peri_numbers(x, y):
            filtered.append(i)
    return list(set(filtered))


region = pd.read_table(sys.argv[1], names=['Uniprot_ID', 'Region'])
df = pd.read_table('../4.panMut/hs.protein.score.Nucleus', names=['Uniprot_ID', 'Gene_name', 'Sequence', 'Score']). \
    merge(region, how='inner', on='Uniprot_ID')

df['NLS'] = df.apply(lambda dataframe: to_numbers(dataframe['Region']), axis=1)
df['peri_NLS'] = df.apply(lambda dataframe: to_peri_numbers(dataframe['Region'], dataframe['Sequence']), axis=1)
df['tele_NLS'] = df.apply(lambda dataframe: to_tele_numbers(dataframe['Region'], dataframe['Sequence']), axis=1)

functional = pd.read_table('../0.dataset/pSTY.FuncScore.txt', names=['Uniprot_ID', 'Position', 'Functional_score']). \
    merge(df, how='inner', on='Uniprot_ID')

Fall_on = []
for i, row in functional.iterrows():
    if row['Position'] in row['NLS']:
        Fall_on.append('NLS')
    elif row['Position'] in row['peri_NLS']:
        Fall_on.append('peri_NLS')
    elif row['Position'] in row['tele_NLS']:
        Fall_on.append('tele_NLS')
    else:
        Fall_on.append('NaN')

functional = functional.assign(Fall_on=Fall_on)

functional.drop(columns=['NLS', 'peri_NLS', 'tele_NLS']).to_csv(sys.argv[2])
#functional.to_csv('position.csv')


