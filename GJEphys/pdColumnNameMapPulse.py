'''
Contains dictionaries mapping certain feature name abbreviations to their longer forms.
'''

"""
Mappings for meta data
"""
mdFN = {
    "expID":                    "Experiment ID",
    "cat":                      "Neuron Category",
    "respType":                 "Response \nType (ON/OFF)",
    "pulseDur":                 "Stimulus\nPulse Duration (ms)",
    "pulseInt":                 "Stimulus\nPulse Interval (ms)",
    "pulseNumber":              "Pulse Number\nin order of\n application"
}

"""
Mappings for response features
"""
fFN = {
    "firstSpikeLat":            "First Spike\nLatency (ms)",
    "allSpikes":                "All Spike\nTimes (ms)",
    "totalSpikeNumber":         "Total Spike\nNumber",
}
