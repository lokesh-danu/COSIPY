import math
import numpy as np
from constants import *
import sys

def penetrating_radiation(GRID, SWnet, dt):
    penetrating_allowed = ['Bintanja95']
    if penetrating_method == 'Bintanja95':
        subsurface_melt, Si = method_Bintanja(GRID, SWnet, dt)
    else:
        raise ValueError("Penetrating method = \"{:s}\" is not allowed, must be one of {:s}".format(penetrating_method, ", ".join(penetrating_allowed)))

    return subsurface_melt, Si


def method_Bintanja(GRID, SWnet, dt):

    # Total height of first layer
    total_height = 0.0

    # Store the total subsurface melt
    subsurface_melt = 0.0

    # Absorption of shortwave radiation
    if GRID.get_node_density(0) <= snow_ice_threshold:
        Si = float(SWnet)*0.1
        depth = np.asarray(GRID.get_depth())
        depth = np.insert(depth,0,0)
        decay = np.exp(17.1*-depth)
        E = Si*np.abs(np.diff(decay))
    else:
        Si = float(SWnet)*0.2
        depth = np.asarray(GRID.get_depth())
        depth = np.insert(depth,0,0)
        decay = np.exp(2.5*-depth)
        E = Si*np.abs(np.diff(decay))

    # List with layer numbers to be removed
    list_of_layers_to_remove = []

    for idxNode in range(0, GRID.number_nodes - 1):

        # New temperature due to penetrating shortwave radiation
        T_rad = float(GRID.get_node_temperature(idxNode) + (E[idxNode] /\
                    (GRID.get_node_density(idxNode) * spec_heat_ice)) * \
                    (dt / GRID.get_node_height(idxNode)))
        
        if (T_rad-zero_temperature>0.0):

            # Temperature difference between layer and freezing temperature
            dT = T_rad - zero_temperature

            # Compute conversion factor A
            A = (spec_heat_ice*ice_density)/(water_density*lat_heat_melting)

            # Changes in volumetric contents; dtheta_w change in water fraciont and dtheate_i change in ice fraction
            dtheta_w = A * dT * GRID.get_node_ice_fraction(idxNode)
            dtheta_i = (water_density/ice_density) * -dtheta_w
            dh = -dtheta_i/GRID.get_node_ice_fraction(idxNode)

            if (dh >= 1.0):
                list_of_layers_to_remove.append(idxNode)
            else:
                GRID.set_node_liquid_water_content(idxNode, \
                    GRID.get_node_liquid_water_content(idxNode)+dtheta_w)
                LWC_temp = GRID.get_node_liquid_water_content(idxNode) * GRID.get_node_height(idxNode)
                GRID.set_node_temperature(idxNode, zero_temperature)
                GRID.set_node_height(idxNode, (1-dh)*GRID.get_node_height(idxNode))
                GRID.set_node_liquid_water_content(idxNode, LWC_temp/GRID.get_node_height(idxNode))

            subsurface_melt = subsurface_melt + dtheta_w * GRID.get_node_height(idxNode)
        else:
            GRID.set_node_temperature(idxNode, T_rad)

    # Remove layers which have been melted
    GRID.remove_node(list_of_layers_to_remove)

    return subsurface_melt, Si
