import numpy as np


def getWelchDF(var1, var2, N1, N2, df1, df2):
    """
    Calculate degree of freedom according to Welch-Satterthwaite equation.
    https://en.wikipedia.org/wiki/Welch's_t-test
    :param var1: float, variance of set1
    :param var2: float, variance of set2
    :param N1: float, size of set 1
    :param N2: float, size of set 2
    :param df1: degree of freedom of set1
    :param df2: degree of freedom of set2
    :return: int, calculated degrees of freedom
    """


    tmp1, tmp2 = var1 / N1, var2 / N2

    df = ((tmp1 + tmp2) ** 2) / (tmp1 ** 2 / (df1 - 1) + tmp2 ** 2 / (df2 - 1))

    if np.isnan(df):
        return df
    else:
        return int(round(df))