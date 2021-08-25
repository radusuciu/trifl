
from trifl.models import (
    database,
    create_tables,
)
from trifl.dataset import (
    TmtDataset,
    SilacDataset
)
import trifl.filters.isobaric as isobaric_filters
import trifl.filters.isotopic as isotopic_filters
import trifl.config as config
import trifl.models as m
import trifl.utils as utils

# installed modules
import pandas as pd

# base modules
from dataclasses import dataclass
from typing import Optional
from functools import partial, reduce
from collections import defaultdict
from collections.abc import Iterable
from abc import ABC, abstractmethod
import itertools
import dataclasses
import pathlib
import operator



@dataclass
class AnalysisParams:
    name: Optional[str] = None
    datasets: Optional[dict] = None
    filters: Optional[dict] = None
    fasta_db_path: Optional[str] = None
    datasets_base_path: Optional[str] = None
    output_folder: Optional[str] = None

    @classmethod
    def from_dict(cls, dct):
        return cls(**{f.name: dct[f.name] for f in dataclasses.fields(cls)})


@dataclass
class SilacAnalysisParams(AnalysisParams):
    parser: Optional[str] = config.DEFAULT_CIMAGE_PARSER
    combine: Optional[str] = config.DEFAULT_SILAC_COMBINE_LEVEL
    datasets_base_url: Optional[str] = ''
    dta_folder_name: Optional[str] = config.DEFAULT_DTA_FOLDER_NAME


class _Analysis(ABC):
    def __init__(self, params: AnalysisParams) -> None:
        # create all tables if they do not already exist
        create_tables()

        self.params = params
        self.name = params.name
        self.filters = self.init_filters(params.filters)
        self.fasta_db_path = pathlib.Path(params.fasta_db_path)
        self.output_path = pathlib.Path(params.output_folder)

        self.datasets = self.init_datasets(
            params.datasets,
            base_path=pathlib.Path(params.datasets_base_path)
        )

        self.fasta_db = m.FastaDatabaseIndex.get_or_create(self.fasta_db_path)
        self.results = None

    @property
    def filtered_out(self) -> set:
        filtered = set()

        for filtered_from_dataset in self._analysis.filters.values():
            for filters in filtered_from_dataset.values():
                for filtered_ids in filters.values():
                    filtered.update(filtered_ids)

        return filtered

    def _new_analysis(self):
        return m.AnalysisModel(
            name=self.name,
            datasets=[d.as_dict() for d in self.datasets]
        )

    def run(self):
        analysis = self._new_analysis()
        results, filtered_out = self.apply_filters()
        analysis.filters = filtered_out
        analysis.save()
        self._analysis = analysis
        database.close()

    def init_filters(self, filters_dict: dict):
        filter_objects = defaultdict(list)

        for level, fs in filters_dict.items():
            if fs is None:
                continue

            for f in fs:
                filter_objects[level].append(self.filter_module.Filter.from_dict(f, self.filter_module))

        return filter_objects

    def apply_filters(self):
        result = dict()

        for dataset in self.datasets:
            _, filtered_out = dataset.apply_filters()
            result[dataset.name] = filtered_out

        return (None, result)

    @property
    def experiments_included(self):
        return list(itertools.chain.from_iterable(d.experiments for d in self.datasets))

    @property
    def experiment_ids_included(self):
        return [e.id for e in self.experiments_included]

class TmtAnalysis(_Analysis):
    data_model = m.TmtData
    filter_module = isobaric_filters
    dataset_class = TmtDataset

    def init_datasets(self, datasets: dict, base_path=''):
        return [TmtDataset(
            name=dataset['name'],
            filters=self.filters,
            replicate_paths=[pathlib.Path(base_path).joinpath(r) for r in dataset['replicates']],
            channel_layout=dataset['channel_layout'],
            control_channels=dataset['control_channels'],
            fasta_db_path=self.fasta_db_path
        ) for dataset in datasets]

    def report(self, report_prefix=config.DEFAULT_REPORT_PREFIX, output_callback=None):
        m = self.data_model
        d = self.datasets[0]

        query = (m
            .select(
                m.id,
                m.experiment,
                m.filename,
                m.scan_num,
                m.uniprot,
                m.symbol,
                m.description,
                m.unique_sequence,
                m.clean_sequence,
                m.channel_intensities,
            )
            .where(
                (m.experiment_id.in_(self.experiment_ids_included))
            )
        )

        df = pd.DataFrame.from_records(list(query.dicts()))
        df = df.set_index('id')

        for channel in d.channel_layout:
            channel_num, channel_name = next(iter(channel.items()))
            new_col_name = f'{channel_name}.intensity'
            # append number to column to distinguish duplicates
            new_col_name = f'{new_col_name}_{len(df.filter(regex=new_col_name).columns)}'
            df[new_col_name] = df.channel_intensities.str.get(channel_num - 1)

        df = df.drop(columns=['channel_intensities'])

        control_cols = df.filter(regex=f'{d.control_channels}.intensity')
        data_cols = df.filter(regex=f'.+\.intensity')

        percentages = data_cols.div(control_cols.mean(axis=1), axis=0).mul(100, axis=0)
        percentages = percentages.rename(columns={
            c: c.replace('intensity', 'percent_of_control')
            for c in df.filter(regex='\.intensity').columns
        })

        df = df.join(percentages)

        report_output_name_template = '{}_{}_{}_{{}}.csv'.format(
            report_prefix,
            self.name,
            self._analysis.id
        )

        report_output_path = utils.get_timestamped_report_path(
            report_output_name_template,
            self.output_path
        )

        if output_callback:
            df = output_callback(df)

        df.to_csv(report_output_path, index=False, encoding='utf-8-sig')
        return df

    def filter_report(self, filename: str = ''):
        pass


class SilacAnalysis(_Analysis):
    data_model = m.SilacData
    filter_module = isotopic_filters
    dataset_class = SilacDataset

    def init_datasets(self, datasets: dict, base_path=''):
        return [self.dataset_class(
            name=dataset['name'],
            filters=self.filters,
            replicate_paths=[pathlib.Path(base_path).joinpath(r) for r in dataset['replicates']],
            fasta_db_path=self.fasta_db_path,
            parser=self.params.parser,
            base_url=self.params.datasets_base_url,
            combine_level=self.params.combine,
            dta_folder_name=dataset.get('dta_folder_name', self.params.dta_folder_name),
        ) for dataset in datasets]

    def report(self, report_prefix=config.DEFAULT_REPORT_PREFIX, whitelist=None, blacklist=None, output_callback=None):
        m = self.data_model

        expression_list = [
            m.experiment_id.in_(self.experiment_ids_included),
            m.id.not_in(self.filtered_out),
        ]

        if whitelist:
            expression_list.append(m.uniprot.in_(whitelist))

        if blacklist:
            expression_list.append(m.uniprot.not_in(blacklist))

        where_expression = reduce(operator.and_, expression_list)

        query = (m
            .select(
                m.id,
                m.experiment,
                m.uniprot,
                m.symbol,
                m.ratio,
            )
            .where(where_expression)
        )

        df = pd.DataFrame.from_records(list(query.dicts()))
        df = self.dataset_class.generate_id(df)

        for dataset in self.datasets:
            df.loc[df.experiment.isin(dataset.experiment_ids_included), 'condition'] = dataset.name

        df = df.set_index('id')

        result_df = (df
            .groupby(['uniprot', 'symbol', 'condition', 'experiment'])
            .agg(ratio=('ratio', 'median'), num_peptides=('ratio', len))
            .groupby(level=('uniprot', 'symbol', 'condition'))  
            .agg(ratio=('ratio', 'median'), num_peptides=('num_peptides', 'sum'), ratio_list=('ratio', list))
        )

        # result_df['num_peptides'] = result_df.num_peptides.astype(int)
        result_df = result_df.unstack(level='condition')
        result_df['num_peptides'] = result_df.num_peptides.fillna(0).astype(int)
        # result_df['ratio'] = result_df.ratio.fillna('-')

        def ratios_to_string(ratios, invert=True):
            if not isinstance(ratios, list):
                return
            if invert:
                ratios = [1/r for r in ratios]
            return ', '.join(map(str, ratios))

        result_df['ratio_list'] = result_df.ratio_list.transform(lambda x: x.apply(ratios_to_string))
        result_df['ratio'] = result_df.ratio.rdiv(1)
        
        result_df.columns = result_df.columns.swaplevel()

        sorted_col_multiindex, _ = result_df.columns.sortlevel()
        result_df = result_df[sorted_col_multiindex]

        # result_df.loc[:, pd.IndexSlice[:, 'ratio']] = (result_df
        #     .loc[:, pd.IndexSlice[:, 'ratio']]
        #     .rdiv(1)
        # )

        result_df = (result_df
            .reset_index()
            .sort_values(by=[(self.datasets[0].name, 'ratio'), 'symbol'], ascending=True)
            .fillna('-')
        )

        report_output_name_template = '{}_{}_{}_{{}}.csv'.format(
            report_prefix,
            self.name,
            self._analysis.id
        )

        report_output_path = utils.get_timestamped_report_path(
            report_output_name_template,
            self.output_path
        )

        if output_callback:
            result_df = output_callback(result_df)

        result_df.to_csv(report_output_path, index=False, encoding='utf-8-sig')
        return result_df

    def filter_report(self, merge_cells=True):
        datasets_to_filter = self.experiment_ids_included

        m = self.data_model

        query = (m
            .select(
                m.id,
                m.experiment,
                m.uniprot,
                m.symbol,
                m.description,
                m.sequence,
                m.clean_sequence,
                m.ratio,
                m.num_ms2,
                m.rsquared,
                m.charge,
                m.meta
            )
            .where(m.experiment_id.in_(datasets_to_filter))
        )

        df = pd.DataFrame.from_records(list(query.dicts()))
        df = self.dataset_class.generate_id(df)

        for dataset in self.datasets:
            df.loc[df.experiment.isin(dataset.experiment_ids_included), 'condition'] = dataset.name

        df = df.set_index('id')

        df = df[[
            '_id',
            'experiment',
            'condition',
            'uniprot',
            'symbol',
            'sequence',
            'meta',
            'num_ms2',
            'rsquared',
            'ratio',
        ]]

        for filtered_out in self._analysis.filters.values():
            for cat, f in filtered_out.items():
                for filter_name, filtered_ids in f.items():
                    df.loc[filtered_ids, f'{cat}.{filter_name}'] = False

        # this should probably be just done in SQL
        # experiment_ids = df.experiment.unique().tolist()
        # q2 = Experiment.select(Experiment.id, Experiment.source_url).where(Experiment.id.in_(experiment_ids))
        # experiments = dict(list(q2.tuples()))
        # df.link = df.apply(lambda x: experiments[x.experiment] + x.link.split('"')[1], axis=1)

        report_output_name_template = 'filter_report_{}_{}_{{}}.xlsx'.format(
            self.name,
            self._analysis.id
        )

        report_output_path = utils.get_timestamped_report_path(
            report_output_name_template,
            self.output_path
        )

        # # df.set_index(['seq_id', 'condition', 'experiment']).sort_index(level=0).to_excel(report_output_path)
        df.set_index(['_id', 'experiment', 'condition']).sort_index(level=0).to_excel(report_output_path, merge_cells=merge_cells)
        return df
    
    def unfiltered_report(self, merge_cells=True):
        m = self.data_model

        query = (m
            .select(
                m.id,
                m.experiment,
                m.uniprot,
                m.symbol,
                m.description,
                m.sequence,
                m.mass,
                m.charge,
                m.rsquared,
                m.ratio,
            )
            .where(
                m.experiment_id.in_(self.experiment_ids_included)
            )
        )

        df = pd.DataFrame.from_records(list(query.dicts()))

        for dataset in self.datasets:
            df.loc[df.experiment.isin(dataset.experiment_ids_included), 'condition'] = dataset.name

        df = df[~df.uniprot.str.startswith('Reverse_')]
        df['description'] = df.description.str.split().str[1:].str.join(' ')

        df = df.set_index(['uniprot', 'symbol', 'description', 'condition', 'experiment']).sort_index(level=0)

        report_output_name_template = 'unfiltered_{}_{}_{{}}.xlsx'.format(
            self.name,
            self._analysis.id
        )

        report_output_path = utils.get_timestamped_report_path(
            report_output_name_template,
            self.output_path
        )

        df.to_excel(report_output_path, encoding='utf-8-sig', merge_cells=merge_cells)
        return df


def analyze_tmt(analysis_params: dict):
    database.init(f'{analysis_params.get("user", config.DEFAULT_DB_NAME)}.db', pragmas=config.pragmas)
    params = AnalysisParams.from_dict(analysis_params)

    analysis = TmtAnalysis(params=params)
    analysis.run()
    return analysis

def analyze_silac(analysis_params: dict):
    database.init(f'{analysis_params.get("user", config.DEFAULT_DB_NAME)}.db', pragmas=config.pragmas)
    params = SilacAnalysisParams.from_dict(analysis_params)

    analysis = SilacAnalysis(params=params)
    analysis.run()
    return analysis
