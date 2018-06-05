import subprocess
import os

def makeContImages(jsonParFile):

    filePath = os.path.split(__file__)[0]
    subprocess.call(['python', os.path.join(filePath, 'makeContImages.py'), jsonParFile])

def makePulseImages(jsonParFile):

    filePath = os.path.split(__file__)[0]
    subprocess.call(['python', os.path.join(filePath, 'makePulseImages.py'), jsonParFile])

def makeSpikeTotalVsFreqImages(jsonParFile):

    filePath = os.path.split(__file__)[0]
    subprocess.call(['python', os.path.join(filePath, 'makeSpikeTotalVsFreqImages.py'), jsonParFile])
