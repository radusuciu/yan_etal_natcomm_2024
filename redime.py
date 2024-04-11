import pathlib
from functools import partial

import numpy as np
import pandas as pd
import polars as pl
from pyteomics import fasta


DATA_INPUT_PATH = pathlib.Path('data/input')
DATA_OUTPUT_PATH = pathlib.Path('data/output')
FASTA_PATH = pathlib.Path('data/fasta/uniprot_mouse_2020_10_10_long-ddhd1_contam_revcat.fasta')

R2_MIN = 0
MAX_RATIO = 32
BANNED_VALUES = (np.nan, np.inf, -np.inf)

sh_map_path = DATA_INPUT_PATH / 'serine_hydrolase_list.tsv'

# map of mouse serine hydrolase uniprot accession to symbol to use in plot
sh_map = dict(zip(*pl.read_csv(sh_map_path, separator='\t').to_dict(as_series=False).values()))


def apply_log2_correction(df_cimage: pd.DataFrame, log2_ratio_correction: float):
    # correct peak_area_ratios and cap
    df_cimage['log2_ratio'] = df_cimage.log2_ratio - log2_ratio_correction
    df_cimage.loc[df_cimage.log2_ratio > np.log2(MAX_RATIO), 'log2_ratio'] = np.log2(MAX_RATIO)
    df_cimage.loc[df_cimage.log2_ratio < -np.log2(MAX_RATIO), 'log2_ratio'] = -np.log2(MAX_RATIO)
    df_cimage.loc[df_cimage.peak_area_ratio == 0, 'log2_ratio'] = np.nan

    # CIMAGE caps values so applying corrections to these values will result
    # in even wonkier stats. As such, assign clipped values regardless of
    # correction factor.
    # use gte and lte comparisons because CIMAGE truncates decimals
    df_cimage.loc[df_cimage.peak_area_ratio >= MAX_RATIO, 'log2_ratio'] = np.log2(MAX_RATIO)
    df_cimage.loc[
        (df_cimage.peak_area_ratio <= 1 / MAX_RATIO) & ((df_cimage.peak_area_ratio != 0)), 'log2_ratio'
    ] = -np.log2(MAX_RATIO)

    df_cimage['peak_area_ratio_corrected'] = 2**df_cimage.log2_ratio.values

    return df_cimage


def read_cimage_output(data_path: pathlib.Path) -> pd.DataFrame:
    df_cimage = pd.read_csv(data_path, sep='\t')
    # drop "extra" column named index
    df_cimage = df_cimage.drop('index', axis=1)

    # column name clean-up
    df_cimage.columns = [col.split('.')[0] for col in df_cimage.columns]
    df_cimage = df_cimage.rename(
        columns={
            'ipi': 'uniprot_accession',
            'IR': 'peak_area_ratio',
            'INT': 'peak_intensity',
            'R2': 'coelution_rsquared',
        }
    )

    # fix min ratios of 0.03, which should be exactly 1/32
    df_cimage.peak_area_ratio = df_cimage.peak_area_ratio.replace(0.03, 1 / MAX_RATIO)

    # calculated fields
    # translation table to remove all non-alpha characters
    non_alpha_delete_translation = str.maketrans(
        '', '', ''.join(c for c in map(chr, range(256)) if not c.isalpha())
    )
    df_cimage['clean_sequence'] = (
        df_cimage.sequence.str.split('.')
        .str.get(1)
        .str.translate(non_alpha_delete_translation)
    )

    fasta_sequences = [
        sequence
        for description, sequence
        in fasta.read(str(FASTA_PATH), use_index=True)
        if not description.startswith('Reverse_')
    ]

    # checking for peptide uniqueness by counting occurences of substring in fasta file
    # this is a naive method, but using here because the actual method we use grabs this from
    # another OpenMS export, which we avoid bringing in here for simplicity
    df_cimage['sequence_counts'] = df_cimage.clean_sequence.apply(
        lambda x: sum(1 for _ in fasta_sequences if x in _)
    )

    # NB: is_unique clashes with pd.Series.is_unique
    df_cimage['is_unique_peptide'] = df_cimage.sequence_counts == 1

    # ignoring divide by zero warnings because we are going to filter out these values anyway
    with np.errstate(divide='ignore'):
        df_cimage['log2_ratio'] = np.log2(df_cimage.peak_area_ratio.values)

    df_cimage['is_valid'] = (df_cimage.coelution_rsquared >= R2_MIN) & (df_cimage.peak_area_ratio != 0)

    # calculate median protein ratio ("median-of-medians")
    good_ratios = df_cimage[(df_cimage.coelution_rsquared >= R2_MIN) & (np.isfinite(df_cimage.log2_ratio))]
    median_log2_protein_ratio = good_ratios.groupby('uniprot_accession').log2_ratio.median().median()
    # if the skew is too great, we do not want to automatically apply a correction
    is_skewed_and_correctable = median_log2_protein_ratio > -1 and median_log2_protein_ratio < 1

    if is_skewed_and_correctable:
        df_cimage = apply_log2_correction(
            log2_ratio_correction=median_log2_protein_ratio,
            df_cimage=df_cimage,
        )
    else:
        df_cimage['peak_area_ratio_corrected'] = df_cimage.peak_area_ratio

    return df_cimage


def get_cimage_protein_df(processed_df: pd.DataFrame) -> pd.DataFrame:
    peptide_df = processed_df.rename(
        columns={
            'is_unique_peptide': 'is_unique',
            'peak_area_ratio': 'original_ratio',
            'peak_area_ratio_corrected': 'ratio',
        }
    )

    protein_df = (
        peptide_df.assign(
            is_unique_quantified_peptide=peptide_df.is_unique & peptide_df.is_valid,
            # whether a peptide ratio is valid or not can depend on more than
            # there being a non-NaN ratio
            # eg. if rsquared cutoff is set to a non-zero value
            # so we set those to NaN here, and they are not included when aggregating
            # below
            log2_ratio=peptide_df.log2_ratio.where(peptide_df.is_valid, np.nan),
        )
        .groupby('uniprot_accession')
        .agg(
            log2_ratio=('log2_ratio', np.nanmedian),
            stdev_log2_ratio=('log2_ratio', partial(np.nanstd, ddof=1)),
            # we're adding spectral counts together even for invalid peptides
            # since spectral counts are already a very course metric and we'd
            # like to maximize discovery potential
            num_total_peptides=('sequence', 'count'),
            num_quantified_peptides=('is_valid', 'sum'),
            num_unique_peptides=('is_unique', 'sum'),
            num_unique_quantified_peptides=('is_unique_quantified_peptide', 'sum'),
            # a protein ratio is set to valid if any of its peptides are valid
            # this can be manually curated later
            is_valid=('is_valid', 'any'),
        )
        .assign(
            ratio=lambda x: np.exp2(x.log2_ratio),
            percent_inhibition=lambda x: (100 * (1 - 1 / x.ratio)).clip(upper=100),
            percent_activity=lambda x: 100 - x.percent_inhibition,
        )
        # we do not store NaN values in the database
        # mostly due to inconsistencies in how they are handled across the stack
        # eg. JSON serialization, differences in between IEEE 754 spec and Postgres implementation
        .replace(BANNED_VALUES, None)
    )

    return protein_df


def get_single_replicate_sh_df(data_path: pathlib.Path) -> pl.DataFrame:
    df_cimage = read_cimage_output(data_path)
    protein_df = get_cimage_protein_df(df_cimage)

    # switching to polars here since I started using it for new code and report generation
    # cimage related code is comparatively old and in "maintenance mode"
    return (
        pl
        .from_pandas(protein_df.reset_index())
        .select([
            'uniprot_accession',
            'percent_inhibition',
            'num_quantified_peptides',
            'num_unique_quantified_peptides',
        ])
        .rename({
            'percent_inhibition': 'percent_competition',
            'num_quantified_peptides': 'num_peptides',
            'num_unique_quantified_peptides': 'num_unique_peptides',
        })
        .filter(pl.col('uniprot_accession').is_in(sh_map.keys()))
        .with_columns(
            symbol=pl.col('uniprot_accession').replace(sh_map),
        )
    )


def get_sh_df() -> pl.DataFrame:
    dfs = []

    for data_path in DATA_INPUT_PATH.glob('*.to_excel'):
        # for simplicity, the tissue, concentration, and replicate are all
        # annotated in the file name, so we parse them back out here
        tissue, concentration_str, replicate_str = data_path.name.replace('.raw_output_rt_10_sn_2.5.to_excel', '').split('_')
        concentration = float(concentration_str.replace('conc', ''))
        replicate = int(replicate_str.replace('rep', ''))

        dfs.append(
            get_single_replicate_sh_df(data_path)
            .with_columns(
                tissue=pl.lit(tissue),
                concentration=pl.lit(concentration),
                replicate=pl.lit(replicate),
            )
        )

    return pl.concat(dfs)


def generate_tissue_report_df(sh_df: pl.DataFrame, tissue: str):
    df = (
        sh_df
        .filter(
            (pl.col('tissue') == tissue) &
            (
                (pl.col('num_unique_peptides') >= 1) &
                (pl.col('num_peptides') >= 3) |
                (pl.col('symbol') == 'AIG1')
            )
        )
        .to_pandas()
        .pivot_table(
            index='symbol',
            columns=['concentration', 'replicate'],
            values='percent_competition',
            sort=True,
        )
        .round(decimals=1)
        .fillna('-')
    )

    report_output_path = DATA_OUTPUT_PATH / f'{tissue}_redime.csv'

    df.to_csv(report_output_path)


def output_reports():
    sh_df = get_sh_df()

    for tissue in sh_df['tissue'].unique():
        generate_tissue_report_df(sh_df, tissue)


if __name__ == '__main__':
    output_reports()
