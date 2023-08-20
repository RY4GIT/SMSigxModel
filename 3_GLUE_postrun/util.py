import numpy as np
import pandas as pd


# Global function
def weighted_quantile(
    values, quantiles, sample_weight=None, values_sorted=False, old_style=False
):
    # Code from https://stackoverflow.com/questions/21844024/weighted-percentile-using-numpy/32216049
    # Equation was taken from here "The weighted percentile method": https://en.wikipedia.org/wiki/Percentile#Definition_of_the_Weighted_Percentile_method

    """Very close to numpy.percentile, but supports weights.
    NOTE: quantiles should be within [0, 1]!
    :param values: numpy.array with data
    :param quantiles: array-like with many quantiles needed
    :param sample_weight: array-like of the same length as `array`
    :param values_sorted: bool, if True, then will avoid sorting of
        initial array
    :param old_style: if True, will correct output to be consistent
        with numpy.percentile.
    :return: numpy.array with computed quantiles.
    """

    # Check input variables
    values = np.array(values)

    quantiles = np.array(quantiles)
    assert np.all(quantiles >= 0) and np.all(
        quantiles <= 1
    ), "quantiles should be in [0, 1]"

    if sample_weight is None:
        sample_weight = np.ones(len(values))
    sample_weight = np.array(sample_weight)

    # Sort values
    if not values_sorted:
        sorter = np.argsort(values)
        values = values[sorter]
        sample_weight = sample_weight[sorter]

    # Calculate the percentiles of the weights
    weighted_quantiles = np.cumsum(sample_weight) - 0.5 * sample_weight
    if old_style:
        # To be convenient with numpy.percentile
        weighted_quantiles -= weighted_quantiles[0]
        weighted_quantiles /= weighted_quantiles[-1]
    else:
        weighted_quantiles /= np.sum(sample_weight)

    # Get the data corresponding to desired quantiles
    # Interpolate the values (y) x weighted quantiles (x1) relationship to quantiles (x2)
    return np.interp(quantiles, weighted_quantiles, values)


def triangle_weight(x, a, b, c):
    # a: lowerlim, c:upperlim, b:midpoint
    y = np.full((len(x),), np.nan)
    # for i in range(len(x)):
    y = np.where((x <= a) | (x >= c), 0, y)
    y = np.where((a <= x) & (x <= b), (x - a) / (b - a), y)
    y = np.where((b <= x) & (x <= c), (c - x) / (c - b), y)
    return y


def to_datetime(df, time_column, format):
    df = df.copy()
    df[time_column] = pd.to_datetime(df[time_column], format=format)
    return df.set_index(time_column)
