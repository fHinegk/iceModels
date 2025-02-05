# iceModels

Matlab scripts for simple 1D mechanical and thermal expansion models of a purely elastic ice beam on hydraulic foundation.

The mechanical models return the lateral deflection, flexural moment and flexural stress within an ice beam strained by a drop in the water level supporting it, with different boundary conditions. The ends of the beam are hinged or fixed, to represent a more or less rigid connection of the ice sheet to the shores. An intermediate support is introduced for cases where the ice beam encounters a bathymetric obstacle (rock, island, sharp coastal feature) during the descent of water level.

The thermal expansion models return the thermal stress within a free floating ice sheet, approximated as an infinite beam in one case, and as a plate under plane stresss conditions in the other (see the work by Evans and Untersteiner: http://dx.doi.org/10.1029/JC076i003p00694). The gravitational terms in the model by Evans and Untersteiner were corrected for the missing gravitational acceleration in the attached scripts.
