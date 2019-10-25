# ------------------------------------------------------------------------------------------------ #
# electrochemicalconstants.py
# Buz Barstow
# Created: 2018-07-26
# ------------------------------------------------------------------------------------------------ #


from rewiredcarbon.physicalconstants import elementaryCharge
from rewiredcarbon.physicalconstants import avogadroNumber

# ------------------------------------------------------------------------------------------------ #

# Redox potential of H20/O2 couple vs. Standard Hydrogen Electrode
vH2O = 0.82

# Redox potential of H+/H2 couple vs. Standard Hydrogen Electrode
vH2 = -0.41

# Bandgap of Si photovoltaic
siBandGapEnergy = 1.1*elementaryCharge

# Redox potential of NADH/NAD+ couple vs. Standard Hydrogen Electrode
vNADH = -0.32

# Redox potential of reduced and oxidized Ferredoxin couple vs. Standard Hydrogen Electrode
vFd = -0.42


# Redox potential of menaquinone (MQ/MQH2) couple vs. Standard Hydrogen Electrode
vMenaquinone = -0.0885

# deltaG for the ADP/ATP reaction on a per molecule basis
deltaGADPATP = 54000./avogadroNumber


