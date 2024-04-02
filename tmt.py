"""
The code in this script was originally part of a larger codebase that processed
data from TMT experiments. For this reason more calculations and metadata are
kept than are strictly needed for the analysis carried out in this paper. The
code is nevertheless and tested to run as is to generate the data tables used
in the the supplementary figures of the paper this project references.

Please feel free to ask any questions in the GitHub repository for this project.
"""

import itertools
import pathlib
import textwrap
from copy import deepcopy
from functools import cached_property
from typing import Iterable, Literal, Optional, Union

import numpy as np
import pandas as pd
import polars as pl
import polars.selectors as cs
from typing_extensions import NotRequired, TypedDict


class ChannelMapDict(TypedDict):
    """Dict mapping condition name to list of channel indices."""

    channel_map: dict[str, list[int]]
    control_condition: NotRequired[str]
    condition_to_control_map: NotRequired[dict[str, str]]
    channel_to_sample_map: NotRequired[dict[Union[str, int], int]]


PlexOptions = Literal[6, 10, 11, 16, 18]

MIN_PER_CONTROL_CHANNEL_INTENSITY = 0
MIN_MEAN_SNR = 2
MAX_CONTROL_CV = 0.5
MAX_LOG2_COMPETITION_RATIO = 5
MIN_LOG2_COMPETITION_RATIO = -5
LOG2_MEDIAN_NORMALIZATION_THRESHOLD = 1


class ChannelMap:
    def __init__(
        self,
        quantification_params: ChannelMapDict,
        num_channels: Optional[int] = None,
        unmapped_condition_name='__unmapped__',
    ) -> None:
        self._quantification_params = quantification_params
        # if num_channels is not explicitly provided, then we just count the number of channels
        # specified in the map.
        self.num_channels = num_channels or sum(
            len(channels) for channels in quantification_params['channel_map'].values()
        )
        self.unmapped_condition_name = unmapped_condition_name
        self.condition_to_indices = self._get_condition_to_indices()

        self.channel_to_sample_map = quantification_params.get(
            'channel_to_sample_map', {channel: channel for channel in range(1, self.num_channels + 1)}
        )

        if 'control_condition' in quantification_params:
            self.condition_to_control_map = {
                condition: quantification_params['control_condition'] for condition in self.conditions
            }
        elif 'condition_to_control_map' in quantification_params:
            self.condition_to_control_map = quantification_params['condition_to_control_map']
        else:
            raise Exception(
                'Channel map must contain either a control condition or a condition to control map.'
            )

    def __repr__(self):
        return textwrap.dedent(
            f'''\
            {self.__class__.__name__}(
                quantification_params={self._quantification_params},
                num_channels={self.num_channels},
                unmapped_condition_name='{self.unmapped_condition_name}'
            )'''
        )

    @cached_property
    def index_to_condition(self) -> dict[int, str]:
        """Maps channel index to condition name.

        Returns:
            dict[int, str]: Dict mapping channel index to condition name.
        """
        return {
            channel_index: condition
            for condition, channel_indices in self.condition_to_indices.items()
            for channel_index in channel_indices
        }

    @cached_property
    def conditions(self) -> list[str]:
        """List of conditions in channel map.

        Returns:
            list[str]: The conditions names.
        """
        return list(self.condition_to_indices.keys())

    def _get_condition_to_indices(self) -> dict[str, list[int]]:
        user_channel_map = deepcopy(self._quantification_params['channel_map'])

        # NOTE: quantification params contain channel map with channel numbers (ie. 1 indexed)
        mapped_channel_numbers = [int(i) for i in itertools.chain.from_iterable(user_channel_map.values())]
        unique_mapped_channel_numbers = set(mapped_channel_numbers)

        all_channel_numbers = set(range(1, self.num_channels + 1))
        unmapped_channel_numbers = all_channel_numbers - unique_mapped_channel_numbers

        if unmapped_channel_numbers:
            user_channel_map[self.unmapped_condition_name] = list(unmapped_channel_numbers)

        return {condition: [num - 1 for num in numbers] for condition, numbers in user_channel_map.items()}


def get_normalized_channel_signals(df: pl.DataFrame, columns: Iterable[str], by: str) -> pl.DataFrame:
    # this is the total signal for each channel
    per_channel_total_signal = (
        pl.col(col).sum().over('channel_index').alias(f'channel_total_{col}') for col in columns
    )
    # this is the mean of the total signal for groups of channels, defined by the `by` column
    per_group_mean_total_signal = (
        pl.col(f'channel_total_{col}').mean().over(by).alias(f'group_mean_channel_total_{col}') for col in columns
    )
    # each individual signal is adjusted so that the total channel signal is equal to the mean of the total channel signal for its group
    normalized_signals = (
        (pl.col(col) * pl.col(f'group_mean_channel_total_{col}') / pl.col(f'channel_total_{col}')).alias(
            f'normalized_{col}'
        )
        for col in columns
    )

    return (
        df.with_columns(per_channel_total_signal)
        .with_columns(per_group_mean_total_signal)
        .with_columns(normalized_signals)
        .drop(cs.matches('^(group_mean_channel_total|channel_total)_.+$'))
        # shouldn't need the selector, see issue:
        # https://github.com/pola-rs/polars/issues/7753
        # .drop('^channel_total_.+$', '^group_mean_channel_total_.+$')
    )


def normalize_pipeline_output_df(
    pipeline_output_df: pl.DataFrame,
    channel_map: ChannelMap,
    normalize_by: Literal['control_condition', 'condition'] = 'control_condition',
) -> pl.DataFrame:
    num_psms = len(pipeline_output_df)
    num_channels = len(pipeline_output_df.get_column('ms3_intensities')[0])
    # condition and respective control condition for each index, in order
    # used to label each intensity value and channel index with its condition and control condition
    # which we need to do because these values determine how normalization is done.
    # for example, for screening experiments that include multiple probes, we want to normalize
    # between channels that use the same probe, which is the control condition.
    conditions = list(dict(sorted(channel_map.index_to_condition.items())).values())
    control_conditions = [channel_map.condition_to_control_map[cond] for cond in conditions]

    columns_to_explode = [
        'ms3_intensities',
        'ms3_signal_to_noise_ratios',
        'channel_index',
        'condition',
        'control_condition',
    ]

    return (
        pipeline_output_df.lazy()
        # this is basically the pk of the PSM before we insert it into the db
        # it's replaced with actual PKs at that time
        .with_row_index('psm_index')
        .with_columns(
            channel_index=pl.repeat([*range(0, num_channels)], num_psms, dtype=pl.List(pl.Int8)),
            condition=pl.repeat(conditions, num_psms, dtype=pl.List(pl.String)).alias('condition'),
            control_condition=pl.repeat(control_conditions, num_psms, dtype=pl.List(pl.String)).alias(
                'control_condition'
            ),
            original_ms3_intensities=pl.col('ms3_intensities'),
            original_ms3_signal_to_noise_ratios=pl.col('ms3_signal_to_noise_ratios'),
        )
        .explode(columns_to_explode)
        .pipe(get_normalized_channel_signals, ['ms3_intensities', 'ms3_signal_to_noise_ratios'], normalize_by)
        # grouping by all columns except the ones we exploded, and any normalized columns
        # also, excepting List[str] columns, because polars doesn't support grouping by them right now
        .group_by(
            ~cs.by_dtype(pl.List(str))
            & ~cs.by_name(
                [
                    *columns_to_explode,
                    'normalized_ms3_intensities',
                    'normalized_ms3_signal_to_noise_ratios',
                ]
            ),
            maintain_order=True,
        )
        .agg(
            # picking first values from List[str] columns because we couldn't group by them
            # even though we know they're all the same.
            cs.by_dtype(pl.List(str)).first(),
            # all of these we agg back to lists, which I chose to do because each row of this
            # df then corresponds to a row in the database.
            # TODO: consider moving this aggregation only where it's needed
            # TODO: consider adding channel_index column to the database and storing
            #       the long-form table instead of using ArrayFields.
            pl.col(
                'ms3_intensities',
                'ms3_signal_to_noise_ratios',
                'normalized_ms3_intensities',
                'normalized_ms3_signal_to_noise_ratios',
            ),
            pl.col('channel_index'),
        )
        # TODO: consider storing original ms3 intensities in database
        .drop(
            'original_ms3_intensities',
            'original_ms3_signal_to_noise_ratios',
            'channel_index',
            'condition',
            'control_condition',
            # dropping because we want to use the normalized columns
            'ms3_intensities',
            'ms3_signal_to_noise_ratios',
        )
        .rename(
            {
                'normalized_ms3_intensities': 'ms3_intensities',
                'normalized_ms3_signal_to_noise_ratios': 'ms3_signal_to_noise_ratios',
            }
        )
        .collect()
    )


def get_psm_df(
    pipeline_output_df: pl.DataFrame,
    channel_map: ChannelMap,
    normalize_intensities=True,
    normalize_by: Literal['control_condition', 'condition'] = 'control_condition',
) -> pl.DataFrame:
    if normalize_intensities:
        return normalize_pipeline_output_df(pipeline_output_df, channel_map, normalize_by=normalize_by)

    return pipeline_output_df.with_row_index('psm_index')


def base_filter(filter_name: str, filter_expr: pl.Expr, reason_template: str, *reason_args) -> pl.Expr:
    return (
        pl.when(filter_expr)
        .then(pl.format(reason_template, *reason_args))
        .alias(f'exclusion_reason_{filter_name}')
    )


def psm_control_cv_filter(max_control_cv: int) -> pl.Expr:
    return base_filter(
        'psm_control_cv_filter',
        pl.col('control_cv') > max_control_cv,
        'Control CV is {} which is greater than the maximum: {}',
        pl.col('control_cv'),
        max_control_cv,
    )


def psm_control_mean_snr_filter(min_mean_snr: int) -> pl.Expr:
    return base_filter(
        'psm_control_mean_snr_filter',
        pl.col('control_mean_snr') < min_mean_snr,
        'Control mean signal-to-noise ratio is {} which is less than the minimum: {}',
        pl.col('control_mean_snr'),
        min_mean_snr,
    )


def psm_control_sum_intensity_filter(min_per_control_channel_intensity: int) -> pl.Expr:
    return base_filter(
        'psm_control_sum_intensity_filter',
        pl.col('control_sum_intensity') < min_per_control_channel_intensity * pl.col('num_control_channels'),
        'Control sum intensity is {} which is less than the minimum: {}',
        pl.col('control_sum_intensity').round(2),
        min_per_control_channel_intensity * pl.col('num_control_channels'),
    )


def get_channel_median_log2_ratios(df: pl.DataFrame) -> pl.DataFrame:
    return (
        df.filter(pl.col('condition') != pl.col('control_condition'))
        .group_by('channel_index')
        .agg(
            pl.col('log2_competition_ratio')
            .filter((pl.col('log2_competition_ratio').is_finite()) & (pl.col('log2_competition_ratio') > 0))
            .median()
            .alias('median_log2_ratio')
        )
        .sort('channel_index')
    )


def median_normalize_log2_ratios(
    df: pl.DataFrame,
    normalization_threshold=LOG2_MEDIAN_NORMALIZATION_THRESHOLD,
    force_median_normalization=False,
) -> pl.DataFrame:
    is_within_norm_threshold = pl.col('median_log2_ratio').is_between(
        -normalization_threshold, normalization_threshold
    )

    normalization_factor_column = (
        pl.when(force_median_normalization | is_within_norm_threshold)
        .then(pl.col('median_log2_ratio'))
        .otherwise(pl.lit(0))
        .alias('normalization_factor')
    )

    normalization_factors_df = (
        get_channel_median_log2_ratios(df)
        .with_columns(normalization_factor_column)
        .select('channel_index', 'normalization_factor')
        .sort('channel_index')
    )

    return (
        df.join(normalization_factors_df, on='channel_index', how='left')
        .with_columns(
            log2_competition_ratio=(
                pl.col('log2_competition_ratio') - pl.col('normalization_factor').fill_nan(0).fill_null(0)
            ).fill_null(np.nan)
        )
        .drop('normalization_factor')
    )


@pl.StringCache()
def get_solo_quant_df(
    psm_df: pl.DataFrame,
    channel_map: ChannelMap,
    index_column: str,
    normalize_to_median=True,
    force_normalize_to_median=False,
    min_mean_snr=MIN_MEAN_SNR,
) -> pl.DataFrame:
    num_channels = channel_map.num_channels

    conditions = [cond for _, cond in sorted(channel_map.index_to_condition.items())]

    control_conditions = [channel_map.condition_to_control_map[cond] for cond in conditions]
    unique_control_conditions = set(control_conditions)

    sample_ids = [
        sample_id for _, sample_id in sorted(channel_map.channel_to_sample_map.items(), key=lambda x: int(x[0]))
    ]
    channel_indices = [*range(num_channels)]

    solo_quant_df = (
        psm_df.rename(
            {
                'ms3_intensities': 'intensity',
                'ms3_signal_to_noise_ratios': 'snr',
            }
        )
        .with_columns(
            pl.repeat(channel_indices, len(psm_df), dtype=pl.List(pl.Int8)).alias('channel_index'),
            pl.repeat(conditions, len(psm_df), dtype=pl.List(pl.String)).alias('condition'),
            pl.repeat(control_conditions, len(psm_df), dtype=pl.List(pl.String)).alias('control_condition'),
            pl.repeat(sample_ids, len(psm_df), dtype=pl.List(pl.Int32)).alias('sample_id'),
        )
        .explode('intensity', 'snr', 'channel_index', 'condition', 'control_condition', 'sample_id')
        .cast(
            {
                'condition': pl.Categorical,
                'control_condition': pl.Categorical,
                # we need this cast to keep things consistent when normalization is not applied
                # in which case polars casts to Int64.. which is fine except we have some code
                # that assumes it's a float that needs to be rounded and I'd rather not special
                # case that when most often we do apply normalization
                'intensity': pl.Float64,
                'snr': pl.Float64,
            }
        )
    )

    controls_grouped_stats = (
        solo_quant_df.filter(pl.col('condition').is_in(unique_control_conditions))
        .group_by(['psm_index', 'condition'])
        .agg(
            control_sum_intensity=pl.col('intensity').sum(),
            control_cv=pl.col('intensity').std() / pl.col('intensity').mean(),
            control_mean_snr=pl.col('snr').mean(),
            num_control_channels=pl.col('channel_index').n_unique(),
        )
        .sort(['psm_index', 'condition'])
        .rename({'condition': 'control_condition'})
    )

    solo_quant_df = (
        solo_quant_df.join(controls_grouped_stats, how='left', on=['psm_index', 'control_condition'])
        .with_columns(
            min_control_sum_intensity=MIN_PER_CONTROL_CHANNEL_INTENSITY * pl.col('num_control_channels'),
        )
        .with_columns(
            psm_control_mean_snr_filter(min_mean_snr),
            psm_control_sum_intensity_filter(MIN_PER_CONTROL_CHANNEL_INTENSITY),
            psm_control_cv_filter(MAX_CONTROL_CV),
        )
        .with_columns(exclusion_reason=pl.coalesce(pl.col('^exclusion_reason_.+$')))
        .with_columns(
            exclude=pl.col('exclusion_reason').is_not_null(),
            exclusion_reason=pl.col('exclusion_reason').fill_null(''),
        )
        .drop('min_control_sum_intensity', 'control_cv', 'control_mean_snr', 'num_control_channels')
        .drop(cs.matches('^exclusion_reason_.+$'))
    )

    # now we are grouping across PSMs, for each index_column (eg. uniprot_accession)
    # and condition among the control conditions, and then calculating a per channel mean
    mean_control_sum_snr = (
        solo_quant_df.filter(
            pl.col('condition').is_in(unique_control_conditions),
        )
        .group_by([index_column, 'condition', 'channel_index'])
        .agg(
            pl.col('snr')
            .filter(~(pl.col('exclude') & (pl.col('condition') == pl.col('control_condition'))))
            .sum()
            .alias('snr')
        )
        .group_by([index_column, 'condition'])
        .agg(pl.col('snr').mean().alias('mean_control_sum_snr'))
        .select([index_column, pl.col('condition').alias('control_condition'), 'mean_control_sum_snr'])
    )

    solo_quant_df = (
        solo_quant_df.join(mean_control_sum_snr, how='left', on=[index_column, 'control_condition'])
        .with_columns(
            snr=pl.when(pl.col('exclude')).then(pl.lit(None)).otherwise(pl.col('snr')),
            intensity=pl.when(pl.col('exclude')).then(pl.lit(None)).otherwise(pl.col('intensity')),
        )
        .group_by([index_column, 'channel_index', 'sample_id', 'condition', 'control_condition'])
        .agg(
            mean_control_sum_snr=pl.col('mean_control_sum_snr').first(),
            total_intensity=pl.col('intensity').sum(),
            total_snr=pl.col('snr').sum(),
            psm_index=pl.col('psm_index'),
            exclude=pl.col('exclude'),
            exclusion_reason=pl.col('exclusion_reason'),
        )
        .with_columns(
            log2_competition_ratio=(pl.col('mean_control_sum_snr') / pl.col('total_snr')).log(base=2),
        )
        .drop('mean_control_sum_snr', 'total_snr')
        .sort(index_column, 'channel_index')
    )

    if normalize_to_median:
        solo_quant_df = solo_quant_df.pipe(
            median_normalize_log2_ratios,
            normalization_threshold=LOG2_MEDIAN_NORMALIZATION_THRESHOLD,
            force_median_normalization=force_normalize_to_median,
        )

    clipped_solo_quant_df = solo_quant_df.with_columns(
        log2_competition_ratio=pl.col('log2_competition_ratio').clip(
            lower_bound=MIN_LOG2_COMPETITION_RATIO,
            upper_bound=MAX_LOG2_COMPETITION_RATIO,
        )
    ).with_columns(
        percent_competition=(100 * (1 - 1 / 2 ** pl.col('log2_competition_ratio'))).clip(lower_bound=0),
        percent_of_control=(100 - 100 * (1 - 1 / 2 ** pl.col('log2_competition_ratio'))).clip(lower_bound=0),
    )

    return clipped_solo_quant_df


def get_peptide_stats(psm_df: pl.DataFrame, solo_quant_df: pl.DataFrame) -> pl.DataFrame:
    return (
        solo_quant_df.explode('psm_index', 'exclude', 'exclusion_reason')
        .join(
            psm_df.select('psm_index', 'sequence', 'is_unique'),
            how='left',
            on='psm_index',
        )
        .with_columns(
            is_unique_quantified=pl.col('is_unique') & ~pl.col('exclude'),
            is_quantified=~pl.col('exclude'),
        )
        .group_by('uniprot_accession', 'sequence', 'condition')
        .agg(
            is_quantified=pl.col('is_quantified').any(),
            is_unique=pl.col('is_unique').any(),
            is_unique_quantified=pl.col('is_unique_quantified').any(),
        )
        .group_by('uniprot_accession', 'condition')
        .agg(
            num_total_peptides=pl.col('sequence').n_unique(),
            num_quantified_peptides=pl.col('is_quantified').sum(),
            num_unique_peptides=pl.col('is_unique').sum(),
            num_unique_quantified_peptides=pl.col('is_unique_quantified').sum(),
            is_valid=pl.col('is_quantified').any(),
        )
    )


def read_pd_output(pd_export_path: pathlib.Path) -> pl.DataFrame:
    # TODO: replace with pure polars approach
    df = pd.read_table(pd_export_path)

    desired_columns = [
        'clean_sequence',
        'species',
        'uniprot_accession',
        'charge',
        'sequence',
        'spectrum_reference',
        'ms2_mz',
        'ms2_rt',
        'ms3_intensities',
        'ms3_signal_to_noise_ratios',
        'score',
        'score_type',
        'is_unique',
    ]

    return pl.from_pandas(
        df.assign(
            clean_sequence=df['Annotated Sequence'].str.split('.').str[1].str.upper(),
            species='',
            uniprot_accession=df['Protein Accessions'].str.split('; '),
            charge=df['Charge'],
            sequence=lambda df: df['Annotated Sequence'].str.cat(
                df['Modifications'].fillna(''), sep='; Modifications: '
            ),
            original_sequence=df['Annotated Sequence'],
            spectrum_reference=df['First Scan'],
            ms2_mz=df['m/z [Da]'],
            ms2_rt=df['RT [min]'],
            ms3_intensities=df.filter(regex='^Abundance:').fillna(0).values.tolist(),
            ms3_signal_to_noise_ratios=df.filter(regex='^Abundance:').fillna(0).values.tolist(),
            score=df['Percolator q-Value'],
            score_type='Percolator q-Value',
            is_unique=df['# Protein Groups'] == 1,
        )[desired_columns]
        .explode('uniprot_accession')
        .reset_index(drop=True)
    )



def output_report():
    DATA_INPUT_PATH = pathlib.Path('data/input')
    DATA_OUTPUT_PATH = pathlib.Path('data/output')
    sh_map_path = DATA_INPUT_PATH / 'serine_hydrolase_list.tsv'
    pl.enable_string_cache()

    # map of mouse serine hydrolase uniprot accession to symbol to use in plot
    sh_map = dict(zip(*pl.read_csv(sh_map_path, separator='\t').to_dict(as_series=False).values()))

    pd_export_path = DATA_INPUT_PATH / 'wat_dose_response_PSMs.txt'

    quantification_params = {
        'channel_map': {
            'DMSO': [1, 3, 5],
            'ABD_1': [8, 10, 12],
            'ABD_10': [14, 16, 18],
            'ABD_0.1': [13, 15, 17],
            'ABD_0.01': [2, 4, 6],
            'ABD_0.001': [7, 9, 11]
        },
        'control_condition': 'DMSO'
    }

    channel_map = ChannelMap(quantification_params)

    output_df = read_pd_output(pd_export_path)
    psm_df = get_psm_df(output_df, channel_map)
    
    protein_quant_df = get_solo_quant_df(
        psm_df=psm_df,
        channel_map=channel_map,
        index_column='uniprot_accession',
        normalize_to_median=True,
        force_normalize_to_median=False,
        min_mean_snr=0,
    )

    peptide_stats = get_peptide_stats(psm_df, protein_quant_df)

    report_df = (
        protein_quant_df.join(
            peptide_stats,
            how='left',
            on=['uniprot_accession', 'condition'],
        )
        .filter(
            pl.col('uniprot_accession').is_in(sh_map.keys()),
            (pl.col('condition') != 'DMSO'),
        )
        .with_columns(
            symbol=pl.col('uniprot_accession').replace(sh_map),
            concentration=pl.col('condition').cast(str).str.split('_').list.last().cast(float),
        )
        .filter(
            (pl.col('num_unique_quantified_peptides') >= 1) &
            (pl.col('num_quantified_peptides') >= 3)
        )
        .to_pandas()
        .pivot_table(
            index='symbol',
            columns=['concentration', 'channel_index'],
            values='percent_competition',
            sort=True,
        )
        .round(decimals=1)
        .fillna('-')
    )

    report_output_path = DATA_OUTPUT_PATH / 'wat_dose_response_tmt.csv'
    report_df.to_csv(report_output_path)


if __name__ == '__main__':
    output_report()
