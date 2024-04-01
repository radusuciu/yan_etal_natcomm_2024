import pathlib

import polars as pl

DATA_INPUT_PATH = pathlib.Path('data/input')
DATA_OUTPUT_PATH = pathlib.Path('data/output')

sh_map_path = DATA_INPUT_PATH / 'serine_hydrolase_list.tsv'

# map of mouse serine hydrolase uniprot accession to symbol to use in plot
sh_map = dict(zip(*pl.read_csv(sh_map_path, separator='\t').to_dict(as_series=False).values()))


def get_single_replicate_sh_df(data_path: pathlib.Path) -> pl.DataFrame:
    # maps names of all columns we want to keep from the original PD names to
    # names that are easier to work with programatically
    desired_column_map = {
        'Accession': 'uniprot_accession',
        '# Unique Peptides': 'num_unique_peptides',
        '# Peptides': 'num_peptides',
        'Abundance Ratio: (Light) / (Heavy)': 'ratio',
    }

    return (
        pl
        .read_csv(
            source=data_path,
            separator='\t',
            columns=list(desired_column_map.keys()),
            # overriding dtype for ratio since it can be missing
            # and we want to make sure that it gets treated a float and not a string
            dtypes={
                'Abundance Ratio: (Light) / (Heavy)': pl.Float32,
            }
        )
        .rename(desired_column_map)
        .filter(pl.col('uniprot_accession').is_in(sh_map.keys()))
        .with_columns(
            symbol=pl.col('uniprot_accession').replace(sh_map)
        )
    )


def get_sh_df() -> pl.DataFrame:
    dfs = []

    for data_path in DATA_INPUT_PATH.glob('*_Proteins.txt'):
        # for simplicity, the tissue, concentration, and replicate are all
        # annotated in the file name, so we parse them back out here
        tissue, concentration_str, replicate_str = data_path.name.replace('_Proteins.txt', '').split('_')
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
        .with_columns(
            percent_competition=(100 * (1 - 1 / pl.col('ratio')))
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


if __name__ == "__main__":
    output_reports()
