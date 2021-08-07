from .shared import *
from .shared import filter_fun as f, Requirement as r

import numpy as np
from scipy import stats as stats
from typing import List

ChannelIntensities = List[int]


@f(requires=r('channel_intensities', 'control_intensities'), level='individual')
def control_sum(control_intensities: ChannelIntensities, min_intensity_per_control_channel: int = 5000) -> bool:
    min_control_sum = len(control_intensities) * min_intensity_per_control_channel
    control_sum = sum(control_intensities)
    return control_sum >= min_control_sum

@f(requires=r('channel_intensities', 'control_intensities'), level='individual')
def control_cv(control_intensities: ChannelIntensities, max_control_cv: float = 0.5) -> bool:
    control_cv = stats.variation(control_intensities)
    return control_cv <= max_control_cv

@f(requires=r('channel_intensities', 'channel_intensities'), level=['individual', 'experiment', 'aggregate'])
def max_cv(channel_intensities: ChannelIntensities, max_cv: float = 0.5) -> bool:
    cv = stats.variation(channel_intensities)
    return cv <= max_cv

@f(requires=r('channel_intensities', 'channel_intensities'), level=['individual'])
def max_condition_cv(channel_intensities: ChannelIntensities, channel_indices: List[int], max_cv: float = 0.5):
    cv = stats.variation(channel_intensities.T[:, channel_indices].T)
    return cv <= max_cv

@f(requires=r('channel_intensities', 'channel_intensities'), level=['individual'])
def channel_is_nonzero(channel_intensities: ChannelIntensities, channel_numbers: List[int]):
    return (channel_intensities.T[:, np.array(channel_numbers) - 1] != 0).all(axis=1)

@f(requires=r('channel_intensities', 'channel_intensities'), level=['individual'])
def total_intensity_sum(channel_intensities: ChannelIntensities, min_total_intensity_sum: int):
    intensity_sum = sum(channel_intensities)
    return intensity_sum >= min_total_intensity_sum

@f(requires=r('tmt_purity', 'tmt_purity'), level=['individual'])
def tmt_purity(tmt_purity: float, min_tmt_purity: float = 0.5):
    return tmt_purity > min_tmt_purity

# @f(requires=r('channel_intensities', 'channel_intensities'), level='aggregate')
# def min_number_of_datasets(channel_intensities: ChannelIntensities, max_cv: float = 0.5) -> bool:
#     cv = stats.variation(channel_intensities)
#     return cv <= max_cv
