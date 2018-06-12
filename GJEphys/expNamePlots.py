'''
This file contains convenience functions to call other python files with specified parameters. They were required to
avoid memory leakage when plotting and saving many figures in a loop.
'''

import subprocess
import os

def makeContImages(jsonParFile):
    '''
    Executes the file "makeContImages.py" in the same path with parameters 'jsonParFile' as a subprocess
    :param jsonParFile: string, path of a json par file
    :return:
    '''
    filePath = os.path.split(__file__)[0]
    subprocess.call(['python', os.path.join(filePath, 'makeContImages.py'), jsonParFile])

def makePulseImages(jsonParFile):
    '''
    Executes the file "makePulseImages.py" in the same path with parameters 'jsonParFile' as a subprocess
    :param jsonParFile: string, path of a json par file
    :return:
    '''
    filePath = os.path.split(__file__)[0]
    subprocess.call(['python', os.path.join(filePath, 'makePulseImages.py'), jsonParFile])

def makeSpikeTotalVsFreqImages(jsonParFile):
    '''
    Executes the file "makeSpikeTotalVsFreqImages.py" in the same path with parameters 'jsonParFile' as a subprocess
    :param jsonParFile: string, path of a json par file
    :return:
    '''
    filePath = os.path.split(__file__)[0]
    subprocess.call(['python', os.path.join(filePath, 'makeSpikeTotalVsFreqImages.py'), jsonParFile])
