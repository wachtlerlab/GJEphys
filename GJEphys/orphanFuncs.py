import numpy as np
import rpy2
import rpy2.robjects

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

# shamelessly stolen from https://github.com/scipy/scipy/pull/4933#issuecomment-292314817
def r_mannwhitneyu(sample1, sample2, exact=True, alternative="g"):
    sample1 = "c({})".format(str(list(sample1))[1:-1])
    sample2 = "c({})".format(str(list(sample2))[1:-1])
    rpy2.robjects.R()("""wres <- wilcox.test({}, {}, alternative="{}"{});
                       rm(sample1);
                       rm(sample2);""".format(sample1, sample2, alternative,
                                              ", exact=TRUE" if exact else ""))
    wres = rpy2.robjects.r['wres']
    uval = wres[0][0]
    pval = wres[2][0]
    return uval, pval