'''
Contains a convenience function for obtaining the labor state of DL-Int-1 experiments used in analyses.
'''

m = {
    '130313-4Rh':   'Forager',
    '130322-1LY':   'Forager',
    '130326-2Rh':   'Forager',
    '130408-1LY':   'Forager',
    '130425-1Al':   'Forager',
    '130501-2Rh':   'Forager',
    '130523-3LY':   'Newly Emerged',
    '130605-1LY':   'Newly Emerged',
    '130605-2LY':   'Newly Emerged',
    '130705-1LY':   'Forager',
    '140424-1LY':   'Forager',
    '140701-1Al':   'Newly Emerged',
    '140813-3Al':   'Newly Emerged',
    '140930-1Al':   'Newly Emerged',
    '140917-1Al':   'Newly Emerged',
    '141030-1Al':   'Newly Emerged',
    }

def expIDLaborStateMap(expName):
    '''
    Returns "Forager" or "Newly Emerged", the labor state of the experiment ID contained in expName
    :param expName: string, experiment ID
    :return:
    '''
    if expName in m:
        return m[expName]
    else:
        'I don\'t know the labor state of the experiment {}'.format(expName)


