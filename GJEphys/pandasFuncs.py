'''
This file contains functions for transforming pandas.DataFrames in certain specific ways.
'''

import pandas as pd


def dfCorr(data, x, ys, func, outLabels, hue=None):
    '''
    Applies func(y, x) for each y in ys. func must return an iterable of the same size as outLabels.
    Returns the results as a DataFrame with ys as indexes and outLabels as columns. If hue is not None, the analysis is
    done for each y grouped on the hue. In this case, columns would be multiindexed with the hue on the higher level and
    outLabels at the lower.
    :param data: a pandas DataFrame
    :param x: the column of data to be used as the independent variable
    :param ys: the columns of data to be used as
    :param func: the function to be applied
    :param outLabels: labels for the values returned by func
    :param hue: grouping key
    :return: a pandas DataFrame
    '''


    allDfs = []
    for y in ys:
        df = pd.DataFrame({'Overall (n={})'.format(len(data)): func(data[x], data[y])}, index=outLabels).T


        if hue is not None:
            groupedData = data.groupby(hue)

            for group in groupedData.groups:

                hueData = groupedData.get_group(group)

                dfT = pd.DataFrame({group + ' (n={})'.format(len(hueData)): func(hueData[x], hueData[y])},
                                   index=outLabels).T

                df = df.append(dfT)

        allDfs.append(df)

    return pd.concat(allDfs, axis=0, keys=ys)

def getLevelUniques(df, level):
    '''
    Returns the unique indexes of a multiindex at a level specified
    :param level: int, the level of the index values to be uniqued and returned
    :return: Index
    '''

    return df.index.get_level_values(level).unique()


def dfInterHueFunc(data, pars, hue, func, outLabels, kwargsDict={}):
    '''
    Splits each column named in pars into series' depending on hue values. Applies func with these series as argument.
    Compiles and returns the results as a dataframe with pars as index names and outLabels as column names.
    :param data: pandas dataframe
    :param pars: iterable strings, valid labels of columns of data
    :param hue: string, valid label of a column of data. The column with 'class labels'
    :param func: a function which takes data[hue].unique().size number of arguments
    :param outLabels: iterable of strings, corresponding to the return values of func
    :return: pandas dataframe
    '''

    hueGroupedData = data.groupby(hue)
    hueGroupedDFs = {hue: hueDF for hue, hueDF in hueGroupedData}

    allDFs = []

    for par in pars:

        df = pd.DataFrame({'{}'.format(par): func(*[hueDF[par] for hueDF in hueGroupedDFs.itervalues()], **kwargsDict)},
                          index=outLabels).T
        allDFs.append(df)

    return pd.concat(allDFs, axis=0)




