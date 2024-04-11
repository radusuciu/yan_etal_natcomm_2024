import itertools
import pathlib
import re
import shutil
import xml.etree.ElementTree as ET

import pandas as pd


AA_MONO_MASS = {
    'G': 57.0214636,
    'A': 71.0371136,
    'S': 87.0320282,
    'P': 97.0527636,
    'V': 99.0684136,
    'T': 101.0476782,
    # 'C' : 103.0091854,
    'C': 160.0306454,  # 103.0091854+57.02146,
    'L': 113.0840636,
    'I': 113.0840636,
    'X': 113.0840636,
    'N': 114.0429272,
    'O': 114.0793126,
    'B': 114.5349350,
    'D': 115.0269428,
    'Q': 128.0585772,
    'K': 128.0949626,
    'Z': 128.5505850,
    'E': 129.0425928,
    'M': 131.0404854,  # +15.99940,
    'H': 137.0589116,
    'F': 147.0684136,
    'R': 156.1011106,
    'Y': 163.0633282,
    'W': 186.0793126,
}

MONO_HPLUS = 1.0072765
MONO_H = 1.0078250
MONO_O = 15.9949146

# this matches customizations made to unimod.xml
# these are specific to this pipeline!
MONO_MASS_TO_MOD_NAME = {
    '15.994915': 'Oxidation',
    '57.021464': 'Carbamidomethyl',
    '28.031300': 'Dimethyl',
    '6.031817': 'Dimethyl Heavy Delta',
    # TMTPro
    '304.207146': '13C(7) 15N(2) C(8) H(25) N O(3)',
}

CIMAGE_MOD_REPLACEMENTS = {
    '.(Dimethyl Heavy Delta)': '',
    '.(Dimethyl)': '',
    '(Dimethyl Heavy Delta)': '',
    '(Dimethyl)': '',
    '(Carbamidomethyl)': '',
    '(Oxidation)': '*',
}


def add_mod_descriptions_to_comet_output(comet_output_path: pathlib.Path) -> None:
    """Adds description attribute to comet output

    Args:
        comet_output_path (pathlib.Path): Path to the comet .pep.xml file to be modified
    """

    ns = {'pepxml': 'http://regis-web.systemsbiology.net/pepXML'}
    ET.register_namespace('', ns['pepxml'])
    tree = ET.parse(comet_output_path)

    mods_aa = tree.findall('.//pepxml:search_summary/pepxml:aminoacid_modification', ns)
    mods_term = tree.findall('.//pepxml:search_summary/pepxml:terminal_modification', ns)

    for mod_el in itertools.chain(mods_aa, mods_term):
        # if it has an amino acid key, return unmodified. if it's a terminal mod,
        # we need to append -term to since that's what OpenMS expects
        amino_acid = mod_el.get('aminoacid', f'{mod_el.get("terminus")}-term')
        mono_mass = mod_el.get('massdiff')

        description = f'{MONO_MASS_TO_MOD_NAME.get(mono_mass)} ({amino_acid})'

        mod_el.set('description', description)
        mod_el.attrib.pop('symbol', None)

    tree.write(str(comet_output_path).replace('.pep.xml', '.fixed_mods.pep.xml'), encoding='UTF-8', xml_declaration=True)


def calculate_sequest_mono_mass(sequence):
    if not sequence:
        return 0

    m = 2 * MONO_H + MONO_O
    for aa in sequence:
        try:
            m += AA_MONO_MASS[aa]
        except KeyError:
            pass
    return m


def sequence_shorthand(sequence: str) -> str:
    for before, after in CIMAGE_MOD_REPLACEMENTS.items():
        sequence = sequence.replace(before, after)
    return sequence.replace('.', '')


def heavy_or_light(sequence: str) -> str:
    return 'heavy' if 'Heavy' in sequence else 'light'


def format_description(description):
    match = re.findall(r'(.+?)\|(.+?)\|(.+?)\s(.+?)(\s..=.+){1,10}', description)  # <5 params typically

    if not match:
        raise RuntimeError(f'Error parsing: {description}')

    source_db, accession, protein_name, description, params = match[0]

    match_gene_symbol = re.findall(r'GN=(.+?)\s|\Z', params)

    if match_gene_symbol:
        gene_symbol = match_gene_symbol[0]
    else:
        # assign protein_name to gene_symbol if it doesn't exist
        gene_symbol = protein_name

    # some headers are explicitly assigned dashes. yuck.
    if gene_symbol == '-' or gene_symbol == '=-' or gene_symbol == '':
        gene_symbol = protein_name

    # Handle prefixes
    if '_' in source_db:
        prefix = source_db.split('_')[0] + '_'
    else:
        prefix = ''

    # handle contaminants
    if 'CONTAMINANT' in description:
        gene_symbol = f'{gene_symbol}_CONTAMINANT'

    return f'{prefix}{gene_symbol} {protein_name} {description}'


def create_ipi_name_table(export_folder: pathlib.Path, dta_folder: pathlib.Path) -> None:
    """Creates ipi_name.table file needed by cimage

    Args:
        export_folder (pathlib.Path): Path containing OpenMS output files
        dta_folder (pathlib.Path): Path to Cimage dta folder where file will be output
    """
    proteins = pd.read_csv(export_folder.joinpath('proteins.csv'), skiprows=1)
    print('Imported proteins.csv.')

    # work on proteins
    proteins = proteins.loc[proteins['#PROTEIN'] == 'PROTEIN']
    proteins = proteins.iloc[:, 1:]

    # extract name
    proteins[['db', 'uniprot_accession', 'protein_name']] = proteins.accession.str.split('|', expand=True)
    proteins.accession = proteins.accession.str.strip()
    proteins.protein_description = proteins.protein_description.str.strip()
    proteins.sort_values('uniprot_accession', inplace=True)

    proteins['name'] = proteins.accession + ' ' + proteins.protein_description
    proteins.name = proteins.name.apply(format_description)
    proteins.name = proteins.uniprot_accession + '\t"' + proteins.name + '"'

    proteins = proteins.drop_duplicates('uniprot_accession')

    with dta_folder.joinpath('ipi_name.table').open('w') as f:
        f.write('name\n' + '\n'.join(proteins.name.values))

    print('Exported ipi_name.table')


def create_all_scan_table(
    df: pd.DataFrame,
    experiment_name: str,
    dta_folder: pathlib.Path,
) -> pd.DataFrame:
    """Create all_scan.table file needed by cimage

    Args:
        df (pd.DataFrame): Dataframe containing PSM information from OpenMS output
        experiment_name (str): The cimage experiment name

    Returns:
        pd.DataFrame: Processed PSM dataframe
    """
    # duplicate rows with non uniq proteins
    # such that we have 1 accession per row
    non_unique = df.loc[df.accessions.str.contains(';')].copy()
    non_unique.accessions = non_unique.accessions.apply(lambda s: s.split(';'))
    df = df.drop(non_unique.index)

    agg_non_uniq = []

    for i, row in non_unique.iterrows():
        accessions = row.accessions.copy()
        for accession in accessions:
            row = row.copy()
            row.accessions = accession
            agg_non_uniq.append(row)

    df = pd.concat((df, pd.DataFrame(agg_non_uniq)))

    # remove contaminants
    df = df.loc[~df.accessions.str.contains('CONTAMINANT')]
    df = df.reset_index()

    # parse annotation data
    df['run_num'] = df.file_origin.str.extract(r'.*_(\d+).idXML').fillna('01')
    df['scan'] = df.spectrum_reference.str.extract('.+scan=(\d+)').astype(int)
    df[['db', 'uniprot_accession', 'protein_name']] = df.accessions.str.split('|', expand=True)
    df['sequence_shorthand'] = df.sequence.apply(sequence_shorthand)
    df['HL'] = df.sequence.apply(heavy_or_light)
    df['key'] = (
        df.uniprot_accession
        + ':'
        + df.aa_before
        + '.'
        + df.sequence_shorthand
        + '.'
        + df.aa_after
        + ':'
        + df.charge.astype(str)
        + ':'
        + df.run_num
    )
    df['run'] = experiment_name

    df = df.sort_values(['HL', 'uniprot_accession', 'aa_before', 'sequence_shorthand'])

    # export all_scan table
    # contaminants may have same uniprot as true entries
    df[['key', 'run', 'scan', 'HL']].drop_duplicates().to_csv(
        dta_folder.joinpath('all_scan.table'), index=False, sep=' '
    )
    print('Exported all_scan.table')
    return df


def create_cross_scan_table(df: pd.DataFrame, experiment_name: str, dta_folder: pathlib.Path) -> None:
    """Create cross_scan.table file needed by cimage

    Args:
        df (pd.DataFrame): Partially processed dataframe containing PSM info from OpenMS output
        experiment_name (str): The cimage experiment name
    """
    df['mass'] = df.sequence_shorthand.apply(calculate_sequest_mono_mass)

    # sort dataframe by percolator score MS:1001492, bigger is better (HUPO-PSI)
    df = df.sort_values(['uniprot_accession', 'key', 'MS:1001492'], ascending=[True, True, False])

    # propagate the best ms2
    cross_scan = df[['key', 'mass', 'scan']].drop_duplicates('key').copy()
    cross_scan.columns = ['key', 'mass', experiment_name]
    cross_scan.to_csv(dta_folder.joinpath('cross_scan.table'), index=False, sep=' ')
    print('Exported cross_scan.table')


def create_input_files(experiment_name: str, cimage_folder: pathlib.Path) -> None:
    """Create cimage input files from OpenMS output.

    Args:
        experiment_name (str): The cimage experiment name
        cimage_folder (pathlib.Path): Path to the folder where cimage will operate
    """
    ratio_dir = 'LH'

    if experiment_name.endswith('_HL'):
        experiment_name = experiment_name[:-3]
        ratio_dir = 'HL'

    print(f'Running cimage on experiment {experiment_name} with ratios as {ratio_dir}.')
    print('Importing PSMs and protein information.')

    export_folder = cimage_folder.parent
    dta_folder = cimage_folder.joinpath('dta')

    df = pd.read_csv(export_folder.joinpath('psms.csv'))

    if df.empty:
        raise Exception('Cannot create Cimage input files due to empty PSM export')

    print('Imported psms.csv')

    df.columns = [c.replace('#', '') for c in df.columns]

    create_ipi_name_table(export_folder, dta_folder)

    # work on PSMs
    df = create_all_scan_table(df, experiment_name, dta_folder)
    create_cross_scan_table(df, experiment_name, dta_folder)


def setup_cimage_folder(base_path: pathlib.Path) -> None:
    cimage_path = base_path / 'cimage'
    dta_path = cimage_path / 'dta'
    output_path = dta_path / 'output'

    cimage_path.mkdir(exist_ok=True)
    dta_path.mkdir(exist_ok=True)
    output_path.mkdir(exist_ok=True)

    mzML_paths = list(base_path.glob('*.mzML'))

    if not mzML_paths:
        raise Exception('Cannot setup cimage folder due to missing mzML files')

    for p in mzML_paths:
        # cimage expects that files are named with suffixes that denote fractions
        # in the format _0N, so if there's no fraction specified in the .mzML
        # we add it here to the link name
        name = p.stem if re.search(r'_\d+$', p.name) else f'{p.stem}_01'
        link_path = cimage_path / f'{name}.mzML'
        link_path.unlink(missing_ok=True)

        # coping the file here because this pipeline is setup with snakemake
        # and some operations are done in a container whereas others are done in the host
        # a symlink would be fine in the host and snakemake but it would be broken inside
        # the container unless we exactly mirrored the directory structure
        # a hardlink also works, but it touches the original mzML and thus makes snakemake
        # think that the file has been updated, leading to the search workflow being re-run
        shutil.copy(p, link_path)
