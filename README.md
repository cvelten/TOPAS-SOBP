# TOPAS-SOBP

The `SOBPSource` python script can be used to generate spread-out Bragg peaks for protons in [TOPAS](https://topasmc.org) through the use of time features. Script usage `python SOBPSource.py -h`:
```
usage: SOBPSource.py [-h] -e ENERGY -c CHI [-p POWERP] [-n NBEAMS] [--delta DELTA] [--recommended]

Generate TOPAS time features to generate spread-out Bragg peaks for protons in water.

optional arguments:
  -h, --help            show this help message and exit
  -e ENERGY, --energy ENERGY
                        maximum energy in MeV
  -c CHI, --chi CHI     SOBP width as fraction of max. range
  -p POWERP, --powerp POWERP
                        power law parameter p; if empty will interpolate from recommended values
  -n NBEAMS, --nbeams NBEAMS
                        minimum number of beamlets
  --delta DELTA         beamlet spacing in cm
  --recommended         print recommended p-value table
```

NB: The user needs to specify the source themselves and set the source's beam energy to
```
d:So/USERSOURCE/BeamEnergy = Tf/BeamEnergy/Value MeV
```

When using this code, please refer to and cite the [manuscript (PubMed)](https://pubmed.ncbi.nlm.nih.gov/33444283/):
> Velten C, Tomé WA. Simulation of spread-out bragg peaks in proton beams using Geant4/TOPAS. Biomed Phys Eng Express. 2020 May 14;6(4):047001. doi: 10.1088/2057-1976/ab8f6d. PMID: 33444283.
