# -*- coding: utf-8 -*-
"""
Created on Wed Jan 20 13:32:04 2016

Old functions created by Arjen Haag during his M.Sc. Thesis.
Not used and since moved to separate script for better oversight.

@author: J.M.Hoch, Department of Physical Geography, Utrecht University

"""
# ---------------------------------------------------------------------
# OLD FUNCTIONS, HAVE BEEN RE-WRITTEN, IMPROVED, OR ARE NO LONGER USED!
# ---------------------------------------------------------------------

def add_waterDepth(model, added_depth, river_cells, model_type, verbose, limit_waterDepth_to_zero = True):
    """
    Setting a user-specified initial water depth for a specified number of cells.
    Particularly useful for FM together with RFS to improve flow in the initial states
    of the simulation.
    Wehn coupling to LISFLOOD-FP, this function is not applicable. Better specify
    initial water level or depth with startfile.
    
    Input:
		- added_depth: user-specified depth to be added
		- river cells: list of cells to be updated; either all cells, or if RFS active, only river cells for instance	
    """
    
    # retrieve current water levels in FM
    if model_type == 'FM':
        updated_waterDepth = np.copy(model.get_var('s1'))
    else:
        os.sys.exit('Function add_waterDepth currently only applicable for Delft3D FM')
    
    # add specified additional water depth to each FM cell
    if (added_depth != 0.):      
        updated_waterDepth[np.array(river_cells)] += added_depth
        # set negative values to zero if specified
        if limit_waterDepth_to_zero == True:
            updated_waterDepth[np.where(updated_waterDepth<0.)] = 0.
        print 'initial water levels of', added_depth, 'added to river cells'  
        
    # overwrite current with updated water level
    model.set_var('s1', updated_waterDepth)
    
    return
    
# =============================================================================

def calculateDeltaVolumes_OLD(model_pcr, missing_value_pcr, secPerDay, CoupledPCRcellIndices, cellarea_data_pcr, use_RFS, model_type):
    """
    Calculating the delta volumes [m3/d] for all coupled PCR-cells.
    Delta volumes are based on discharge, surfaceRunoff, and topWaterLayer (only 2way-coupling) in PCR-GLOBWB.
    
    Input:
        - list with indeces pointing to coupled PCR-cells
        - PCR landmask data
        
    Output:
        - if RFS active, two arrays with delta volumes for river and floodplain cells, respectively
        - if RFS not active, one array with aggregated delta volumes
        - all outputs are in m3/day
    """
    #- update PCR for one day to get values
    model_pcr.update(1)

    #- retrieve data from PCR-GLOBWB
    current_discharge_pcr  = model_pcr.get_var('discharge')
    current_runoff_pcr     = model_pcr.get_var('landSurfaceRunoff')
    current_waterlayer_pcr = model_pcr.get_var('topWaterLayer')
    
    # 1a. Discharge
    
    # prepare empty array for all PCR-cells
    water_volume_PCR_rivers = np.zeros([len(current_discharge_pcr),len(current_discharge_pcr[0])])
    
    # loop over current discharge and convert to m3/d; missing values are replaced with zero
    for i in range(len(current_discharge_pcr)):
        for j in range(len(current_discharge_pcr[0])):
            if current_discharge_pcr[i][j] != missing_value_pcr:
                water_volume_PCR_rivers[i][j] = current_discharge_pcr[i][j] * secPerDay
            else:
                water_volume_PCR_rivers[i][j] = 0.
                
	# prepare empty array for coupled PCR-cells
    water_volume_PCR_rivers_coupled = np.zeros(len(CoupledPCRcellIndices))
    
    # get daily discharge volumes for all coupled PCR-cells [m3/day]    
    for i in range(len(CoupledPCRcellIndices)):
        water_volume_PCR_rivers_coupled[i] = water_volume_PCR_rivers[CoupledPCRcellIndices[i]]
        
    # 1b. Runoff and Waterlayer
        
    # prepare empty arrays for all PCR-cells
    water_volume_PCR_runoff      = np.zeros([len(current_runoff_pcr),len(current_runoff_pcr[0])])
    water_volume_PCR_waterlayer  = np.zeros([len(current_waterlayer_pcr),len(current_waterlayer_pcr[0])])
    
    # loop over current runoff and waterlayer and convert to m3/d; missing values are replaced with zero    
    for i in range(len(current_runoff_pcr)):
        for j in range(len(current_runoff_pcr[0])):
            if current_runoff_pcr[i][j] != missing_value_pcr:
                water_volume_PCR_runoff[i][j]     = current_runoff_pcr[i][j] * cellarea_data_pcr[i][j]
                water_volume_PCR_waterlayer[i][j] = current_waterlayer_pcr[i][j] * cellarea_data_pcr[i][j]
            else:
                water_volume_PCR_runoff[i][j]     = 0.
                water_volume_PCR_waterlayer[i][j] = 0.
                
    # prepare empty arrays for coupled PCR-cells
    water_volume_PCR_runoff_coupled      = np.zeros(len(CoupledPCRcellIndices))
    water_volume_PCR_waterlayer_coupled  = np.zeros(len(CoupledPCRcellIndices))
    water_volume_PCR_floodplains_coupled = np.zeros(len(CoupledPCRcellIndices))
    
    # prepare empty arrays for coupled PCR-cells if RFS is not active
    water_volume_PCR_coupled = np.zeros(len(CoupledPCRcellIndices))
    
    # get daily runoff and waterlayer volumes for all coupled PCR-cells [m3/day] 
    for i in range(len(CoupledPCRcellIndices)):
        water_volume_PCR_runoff_coupled[i]      = water_volume_PCR_runoff[CoupledPCRcellIndices[i]]
        water_volume_PCR_waterlayer_coupled[i]  = water_volume_PCR_waterlayer[CoupledPCRcellIndices[i]]
        water_volume_PCR_floodplains_coupled[i] = water_volume_PCR_runoff_coupled[i] + water_volume_PCR_waterlayer_coupled[i]
        # directly merge volumes from discharge, runoff, and waterlayer 
        water_volume_PCR_coupled[i] 			= water_volume_PCR_floodplains_coupled[i] + water_volume_PCR_rivers_coupled[i]
    
    # 2. calculate delta volumes for coupled cells [m3/day]
    
    if (use_RFS == True) and (model_type == 'FM'):
	
        delta_volume_PCR_rivers_coupled      	= water_volume_PCR_rivers_coupled - 0.
        delta_volume_PCR_floodplains_coupled 	= water_volume_PCR_floodplains_coupled - 0.
        delta_volume_PCR_coupled				= 0.
        
    if (use_RFS == False) and (model_type == 'FM'):
	
        delta_volume_PCR_rivers_coupled      	= 0.
        delta_volume_PCR_floodplains_coupled 	= 0.
        delta_volume_PCR_coupled 				= water_volume_PCR_coupled - 0.
        
    if model_type == 'FP':
	
        delta_volume_PCR_rivers_coupled      	= 0.
        delta_volume_PCR_floodplains_coupled 	= 0.
        delta_volume_PCR_coupled 				= water_volume_PCR_coupled - 0.
        
    return delta_volume_PCR_rivers_coupled, delta_volume_PCR_floodplains_coupled, delta_volume_PCR_coupled
 
# =============================================================================  

def calculateDeltaWater_OLD(model, CouplePCR2model, CoupleModel2PCR, delta_volume_PCR_rivers_coupled, delta_volume_PCR_floodplains_coupled, \
                                                    boolean_river_cell_in_coupled_PCR_cell, bottom_lvl, \
                                                    river_cells_per_coupled_PCR_cell, FM_floodplain_cells_per_coupled_PCR_cell, \
                                                    CellAreaSpherical, fraction_timestep, verbose_folder, model_type, useFluxes, verbose):
    """
   In this function the calculated daily delta volumes [m3/d] is translated to suitable units later to be used in the updating step.
    The input volumes of PCR-GLOBWB are here divided over the number of hydrodynamic cells within each PCR-cell.
    When fluxes are used, daily volumes are simply translated to averaged flux of m3/s per day.
    When states are used, daily volumes are first divided by the cell area of each coupled hydrodynamic cell, and then divided by the
    user-specified timestep, if provided.
    
    Input:
        - delta volumes
        - cell areas
        - model settings
    Output:
        - delta state [m/d] or [m/timestep]
        - delta flux [m3/s]
    """
    
    #- creating zero arrays to be populated with data
    additional_water_level = np.zeros(len(CoupleModel2PCR))
    additional_water_volume = np.zeros(len(CoupleModel2PCR))
    
    verbose_volume = np.zeros(len(CoupleModel2PCR))
    
    verbose_volume_rivers = np.zeros(len(CoupleModel2PCR))
    verbose_volume_floodplains = np.zeros(len(CoupleModel2PCR))
    verbose_volume_coupled = np.zeros(len(CoupleModel2PCR))
    
    # loop over all coupled PCR cells
    for i in range(len(CouplePCR2model)):
                   
        # check if there are any river cell within current PCR cell
        if boolean_river_cell_in_coupled_PCR_cell[i]:
            
            # 1. RIVER CELLS:
            
            # check if delta volume is positive
            if delta_volume_PCR_rivers_coupled[i] >= 0:
                
                # divided volume by number of hydrodynamic cells in each PCR-cell
                temp_water_volume_river_FM = delta_volume_PCR_rivers_coupled[i] / len(river_cells_per_coupled_PCR_cell[i])
                
                # at each PCR-cell...
                for j in range(len(river_cells_per_coupled_PCR_cell[i])):
					
                    # ...get current FM cell index
                    current_cell_index = river_cells_per_coupled_PCR_cell[i][j]
                    
                    if useFluxes == False:
						
						# ...calculate additional water levels [m/d] for current index
                        additional_water_level[current_cell_index] = temp_water_volume_river_FM / CellAreaSpherical[current_cell_index]
                        
                    elif useFluxes == True:
						
						# ...calculate additional water volume [m3/d] for current index
                        additional_water_volume[current_cell_index] = temp_water_volume_river_FM / 1.
                        
                    if verbose == True:
						
                        verbose_volume_rivers[current_cell_index] = temp_water_volume_river_FM / 1.
                        verbose_volume[current_cell_index] = temp_water_volume_river_FM / 1.
                        
            else:
                sys.exit('delta_volume_PCR_rivers_coupled is negative, should not be the case - check and fix this!')
            
            # 2. FLOODPLAIN CELLS:
            
            # check if delta volume is positive
            if delta_volume_PCR_floodplains_coupled[i] >= 0:
                
                # divided volume by number of hydrodynamic cells in each PCR-cell
                temp_water_volume_floodplains_FM = delta_volume_PCR_floodplains_coupled[i] / len(FM_floodplain_cells_per_coupled_PCR_cell[i])
                
                # at each PCR-cell...
                for j in range(len(FM_floodplain_cells_per_coupled_PCR_cell[i])):
                    
                    # ...get current FM cell index
                    current_cell_index = FM_floodplain_cells_per_coupled_PCR_cell[i][j]
                    
                    if useFluxes == False:
						
						# ...calculate additional water levels [m/d] for current index						
                        additional_water_level[current_cell_index] = temp_water_volume_floodplains_FM / CellAreaSpherical[current_cell_index]
                        
                    elif useFluxes == True:
						
						# ...calculate additional water volume [m3/d] for current index
                        additional_water_volume[current_cell_index] = temp_water_volume_floodplains_FM / 1.
                        
                    if verbose == True:

                        verbose_volume_floodplains[current_cell_index] = temp_water_volume_floodplains_FM / 1.						
                        verbose_volume[current_cell_index] = temp_water_volume_floodplains_FM / 1.

            else:
                sys.exit('delta_volume_PCR_floodplains_coupled is negative, should not be the case - check and fix this!')

        # if there is no river cell within the current PCR-cell
        else:
            
            # merge delta volumes from rivers and floodplains
            delta_volume_PCR_total_coupled = delta_volume_PCR_rivers_coupled[i] + delta_volume_PCR_floodplains_coupled[i]
            
            # check if delta volume is positive
            if delta_volume_PCR_total_coupled >= 0:
				
                # divided volume by number of hydrodynamic cells in each PCR-cell
                temp_water_volume = delta_volume_PCR_total_coupled / len(CouplePCR2model[i])
                
                # at each PCR-cell..
                for j in range(len(CouplePCR2model[i][1])):
                    
                    # ...get current FM cell index
                    current_cell_index = CouplePCR2model[i][1][j]
                    
                    if useFluxes == False:
						
                        # ...calculate additional water levels [m/d] for current index
                        additional_water_level[current_cell_index] = temp_water_volume / CellAreaSpherical[current_cell_index]
                        
                    elif useFluxes == True:
						
						# ...calculate additional water volume [m3/d] for current index
                        additional_water_volume[current_cell_index] = temp_water_volume / 1.
                        
                    if verbose == True:
						
                        verbose_volume_coupled[current_cell_index] = temp_water_volume / 1.                        
                        verbose_volume[current_cell_index] = temp_water_volume / 1.

            else:
                sys.exit('delta_volume_PCR_total_coupled is negative, should not be the case - check and fix this!')
    
    # calculate additional water levels or fluxes based on chosen settings
    if (useFluxes == False) and (model_type == 'FM'):
        delta_water = additional_water_level / fraction_timestep        # [m/timestep]
        
    elif (useFluxes == False) and (model_type == 'FP'):
        delta_water = np.copy(additional_water_level)                   # [m/day]
   
    elif (useFluxes == True) and (model_type == 'FM'):              
        delta_water = additional_water_level * 1000.                    # [mm/day] 
        
    elif (useFluxes == True) and (model_type == 'FP'):              
        delta_water = additional_water_volume / 86400.                  # [m3/s]
    
    return delta_water, verbose_volume, verbose_volume_rivers, verbose_volume_floodplains, verbose_volume_coupled

# =============================================================================
    
def extractAndConvertVolumes_OLD(CouplePCR2model, model, boolean_river_cell_in_coupled_PCR_cell, river_cells_per_coupled_PCR_cell, \
                                                             FM_floodplain_cells_per_coupled_PCR_cell, inundated_area_floodplains_model_2_PCR_coupled):
    """
    Only required for 2way-coupling!    
    extracting FM water volume for each coupled PCR cell and converting these to water depths to add back to PCR
    Exists, but not actually used for 1way-coupling.
    """
    current_volume = model.get_var('vol1')
    
    water_volume_rivers_FM_2_PCR_coupled = np.zeros(len(CouplePCR2model))
    water_volume_floodplains_FM_2_PCR_coupled = np.zeros(len(CouplePCR2model))
    water_depths_floodplains_FM_2_PCR_coupled = np.zeros(len(CouplePCR2model))
    
    # loop over all coupled PCR cells
    for i in range(len(CouplePCR2model)):
        # check if there are any river cells:
        if boolean_river_cell_in_coupled_PCR_cell[i]:
            
            # 1. RIVER CELLS
            water_volume_rivers_FM_2_PCR_coupled[i] = np.sum(current_volume[river_cells_per_coupled_PCR_cell[i]])
            # 2. FLOODPLAIN CELLS
            water_volume_floodplains_FM_2_PCR_coupled[i] = np.sum(current_volume[FM_floodplain_cells_per_coupled_PCR_cell[i]])
            # divide total volume by inundated FM area to obtain water depth
            if inundated_area_floodplains_model_2_PCR_coupled[i] > 0. :
                water_depths_floodplains_FM_2_PCR_coupled[i] = water_volume_floodplains_FM_2_PCR_coupled[i] / inundated_area_floodplains_model_2_PCR_coupled[i]
            else:
                water_depths_floodplains_FM_2_PCR_coupled[i] = 0.    
                
        # if there are no river cells:
        else:
            # extracting total water volume from all FM cell coupled to current PCR cel
            water_volume_floodplains_FM_2_PCR_coupled[i] = np.sum(current_volume[CouplePCR2model[i][1]])
            # divide total volume by inundated FM area to obtain water depth
            if inundated_area_floodplains_model_2_PCR_coupled[i] > 0. :
                water_depths_floodplains_FM_2_PCR_coupled[i] = water_volume_floodplains_FM_2_PCR_coupled[i] / inundated_area_floodplains_model_2_PCR_coupled[i]
            else:
                water_depths_floodplains_FM_2_PCR_coupled[i] = 0.
    
    return water_volume_rivers_FM_2_PCR_coupled, water_volume_floodplains_FM_2_PCR_coupled, water_depths_floodplains_FM_2_PCR_coupled

# =============================================================================

def convertMeterToMeter3(inputArray, cellAreaMap, missingValues=-999):
    """
    Converts units of a numpy array created from a PCRaster map from meters [m] to cubic meters [m3] by multiplying with the cell area [m2]
    
    Input:  numpy array with input values, cell area map, missing values (optional, default = -999)
    
    Output: numpy array with output values
    
    NOTE: currently requires a cell area map of the same area as the landmask map of the model!
	      (could also be done with data read from model initialized with BMI)
    """
    # get cell area
    try:
        # if cell area input is path to PCR map
        cell_area_map = pcr.readmap(cellAreaMap)
    except:
        # if cell area input is already loaded PCR map
        cell_area_map = cellAreaMap
    # convert input array to PCRaster map (for easier multiplication)
    input_values_map  = pcr.numpy2pcr(pcr.Scalar, inputArray, missingValues)
    # multiply maps ([m]*[m2]=[m3])
    output_values_map = input_values_map * cell_area_map
    # convert resulting map back to numpy array
    output_values_np  = pcr.pcr2numpy(output_values_map, missingValues)
    
    return output_values_np

def convertPCRindicesOLD(indexList, PCRmap, missingValues=255):
    """
    Converts a list with single indices referring to a point on a PCR map to double indices (array,column).
    This is required to call PCR cells from maps read during the coupling, since these work with double (array,column) indices,
    and this basically reverses the single index that was created in the coupling functions to find coupled cells more easily.
    
    Input:  - list with single indices (output of coupleAllCells[3] or [4])
            - reference PCR map (e.g. landmask)
            - missing values of map (optional, default = 255 for landmask)
    
    Output: list with double indices (array,column)
    """
    # get reference map
    try:
        # if map is given as path to PCR map
        reference_map = pcr.readmap(PCRmap)
    except:
        # if map is an already loaded PCR map
        reference_map = PCRmap
    # create zero arrays for filling in values
    indices_array  = np.zeros(len(indexList))
    indices_column = np.zeros(len(indexList))
    # loop over all single indices
    for i in range(len(indexList)):
        # reset break
        break_loop = False
        # get the current single index
        try:
            # if single index list is given as list of lists (from function coupleAllCells[3])
            current_index = indexList[i][0]
        except:
            # if single index list is given as list (from function coupleAllCells[4])
            current_index = indexList[i]
        # start counting towards single index at zero
        count = 0
        # loop over all arrays
        for p in range(len(pcr.pcr2numpy(reference_map, missingValues))):
            # break out of this loop if index was found in next loop
            if break_loop:
                break
            # loop over all columns
            for q in range(len(pcr.pcr2numpy(reference_map, missingValues)[0])):
                # check for missing values
                if pcr.pcr2numpy(reference_map, missingValues)[p][q] != missingValues:
                    # count corresponding single index
                    count += 1
                    # check if count matches single index
                    if count == current_index+1:
                        # assign array and column index
                        indices_array[i]  = p
                        indices_column[i] = q
                        # break out of the array and column loops
                        break_loop = True
                        break
    # create a single list with combined (array,column) indices
    indexListDouble = zip(indices_array,indices_column)
    
    return indexListDouble

def coupleCellsFMtoPCRold(FMcoords, FMcoordsCentroids, PCRcoords, PCRcoordsMinMax, no_couple_value=None):
    """
    Creates a list of coupled cells by finding a coupled PCR cell for each FM cell
    
    Input:  - list of (x,y) coordinates of each FM polygon      (output of getFMcoords)
            - tuple of all coordinates of FM polygon centroids  (output of getFMcoordsCentroids)
            - list of (x,y) coordinates of each PCR polygon     (output of getPCRcoords)
            - tuple of all min/max coordinates of PCR polygons  (output of getPCRcoordsMinMax)
            - value that is used to indicate no coupled cell was found
            
    Output: - list of (x,y) coordinates of each FM polygon that is coupled
            - list of (x,y) coordinates of each PCR polygon that is coupled
            - list of all FM cells and its coupled PCR cells (if any)
            - array of PCR cells that were coupled to FM cells       (used later on as input for coupleCellsPCRtoFMold)
            - list of unique PCR cells that were coupled to FM cells (used later on as input for coupleCellsPCRtoFMold)
    """
    # splitting tuples with data into seperate arrays (so it is easier to see what is done in the code)
    centroids_x_fm            = FMcoordsCentroids[0]
    centroids_y_fm            = FMcoordsCentroids[1]
    x_min_all_cell_coords_pcr = PCRcoordsMinMax[0]
    x_max_all_cell_coords_pcr = PCRcoordsMinMax[1]
    y_min_all_cell_coords_pcr = PCRcoordsMinMax[2]
    y_max_all_cell_coords_pcr = PCRcoordsMinMax[3]
    
    # List with all FM cells and coupled PCR cells (if any):
    
    # initial list with FM cell numbers and zero values for their coupled PCR cells
    coupling_fm_pcr = []
    
    for i in range(len(FMcoords)):
        coupling_fm_pcr.append([(i),(0)])
        
    # empty array for assigning only PCR cells (or 'no couple value')
    coupling_pcr_only = np.zeros(len(FMcoords))
        
    # fill in zero array of coupled PCR cells:
    # looping over all FM cells
    for i in range(len(FMcoords)):
        
        # finding the to-be-coupled PCR cell
        temp = np.where((centroids_x_fm[i] >= x_min_all_cell_coords_pcr) & \
                        (centroids_x_fm[i] <  x_max_all_cell_coords_pcr) & \
                        (centroids_y_fm[i] >= y_min_all_cell_coords_pcr) & \
                        (centroids_y_fm[i] <  y_max_all_cell_coords_pcr))[0]
        
        # check if a cell was found
        if temp.size != 0:
            # if so, assign this value to the coupling list and PCR array
            coupling_fm_pcr[i][1] = temp[0]
            coupling_pcr_only[i]  = temp[0]
        else:
            # if not, assign the chosen 'no couple value' to the coupling list and PCR array
            coupling_fm_pcr[i][1] = no_couple_value
            coupling_pcr_only[i]  = no_couple_value
    
    # Lists with coordinates of FM and PCR cells that are coupled:
    
    # Create lists with coupled cells only:
    # empty lists
    FM_coupled_cells = []
    PCR_coupled_cells = []
    
    # loop over all indices of coupled cell list
    for i in range(len(coupling_fm_pcr)):
        
        # check if a PCR cell was found for coupling
        if coupling_fm_pcr[i][1] != no_couple_value:
            
            # if so, assign the current cell to the FM list
            FM_coupled_cells.append((coupling_fm_pcr[i][0]))
            
            # and the found PCR cell to the PCR list
            PCR_coupled_cells.append((coupling_fm_pcr[i][1]))
            
    # get only the unique PCR cells 
    # (there will be many duplicate cells, since the code tried to find a PCR cell for every single FM cell)
    PCR_coupled_cells_unique = list(set(PCR_coupled_cells))
            
    # Create lists with coordinates of coupled cells only:
    # empty lists
    all_cell_coords_FM_coupled = []
    all_cell_coords_PCR_coupled = []
    
    # loop over all coupled FM cells and get their coordinates
    for i in range(len(FM_coupled_cells)):
        all_cell_coords_FM_coupled.append(FMcoords[int(FM_coupled_cells[i])])
    
    # loop over all coupled PCR cells and get their coordinates    
    for i in range(len(PCR_coupled_cells_unique)):
        all_cell_coords_PCR_coupled.append(PCRcoords[int(PCR_coupled_cells_unique[i])])
        
    return all_cell_coords_FM_coupled, all_cell_coords_PCR_coupled, \
           coupling_fm_pcr, coupling_pcr_only, PCR_coupled_cells_unique
           
def coupleCellsPCRtoFMold(allCoupledPCRcells, uniqueCoupledPCRcells):
    """
    Creates a list of coupled cells by finding a coupled FM cell for each PCR cell
    
    Input:  - array of PCR cells that were coupled to FM cells       (output of coupleCellsFMtoPCRold)
            - list of unique PCR cells that were coupled to FM cells (output of coupleCellsFMtoPCRold)
            
    Output: list of each coupled PCR cell and all FM cell coupled to them
    """    
    # initial list with zero values for the FM cells
    coupling_pcr_fm = []
    
    for i in range(len(uniqueCoupledPCRcells)):
        coupling_pcr_fm.append([(uniqueCoupledPCRcells[i]),(0)])
    
    # loop over all unique coupled PCR cells
    for i in range(len(coupling_pcr_fm)):
        
        # find all FM cells that are coupled to the current PCR cell, and assign these to the list
        coupling_pcr_fm[i][1] = np.where(uniqueCoupledPCRcells[i]==allCoupledPCRcells)[0]
        
    return coupling_pcr_fm

def getFMcellsRiverFloodplainOLD(coupling_pcr_fm, FMcoords, FMareas, area_threshold='all', threshold_factor_3vertices=1.5, threshold_factor_4vertices=3.0):
    """
    Creates lists of FM river and floodplain cells
    
    Input:  - list of all coupled PCR cells and their coupled FM cells [3rd output of coupleAllCells]
            - list of FM cell coordinates                              [output of getFMcoords]
            - list of FM cell areas                                    [output of calculateAllCellAreaSpherical or obtained through BMI get_var('ba')]
            - option to use the minimum area of all cells or only current PCR cells in loop as threshold base-level
            - multiplication factors for threshold base-level for cells with 3 and 4 vertices
    
    Output: - list of FM river cells
            - list of FM floodplain cells
            - list of coordinates of FM river cells (useful for plotting/checking results)
            - list of coordinates of FM floodplain cells (useful for plotting/checking results)
    """
    
    # create empty list for filling in river cells
    FM_cells_river = []
    
    # create empty list for filling in floodplain cells
    FM_cells_floodplain = []
    
    for i in range(len(coupling_pcr_fm)):
    
        # -------------------------------------------------------------------------------------------------------
        # Create a separate array for FM cells of this coupled PCR cell only
        # -------------------------------------------------------------------------------------------------------
    
        # empty array
        FM_cells = np.zeros(len(coupling_pcr_fm[i][1]), dtype=np.int64)
    
        # loop over all FM cells and assign to empty array
        for j in range(len(coupling_pcr_fm[i][1])):
            FM_cells[j] = coupling_pcr_fm[i][1][j]
    
        # get coordinates of these cells
        FM_cells_coords = []
        for j in range(len(FM_cells)):
            FM_cells_coords.append(FMcoords[FM_cells[j]])
    
        # get areas of these cells
        FM_cells_areas = []
        for j in range(len(FM_cells)):
            FM_cells_areas.append(FMareas[FM_cells[j]])
    
        # -------------------------------------------------------------------------------------------------------
        # Get threshold
        # -------------------------------------------------------------------------------------------------------
    
        if area_threshold == 'all':
    
            # set thresholds
            threshold_4vertices = min(FMareas)*threshold_factor_4vertices
            threshold_3vertices = min(FMareas)*threshold_factor_3vertices
    
        else:
    
            # set thresholds
            threshold_4vertices = min(FM_cells_areas)*threshold_factor_4vertices
            threshold_3vertices = min(FM_cells_areas)*threshold_factor_3vertices
    
        # -------------------------------------------------------------------------------------------------------
        # River cells
        # -------------------------------------------------------------------------------------------------------
    
        # assign to this list the cells that have 4 vertices and that fall below the relevant threshold
        FM_cells_river_temp = FM_cells[np.where(FM_cells_areas<=threshold_4vertices)[0]]
        for j in range(len(FM_cells_river_temp)):
            if len(FMcoords[FM_cells_river_temp[j]]) == 4:
                FM_cells_river.append(FM_cells_river_temp[j])
    
        # and all cells that fall below the threshold for polygons with 3 vertices
        FM_cells_river_temp = FM_cells[np.where(FM_cells_areas<=threshold_3vertices)[0]]
        for j in range(len(FM_cells_river_temp)):
            FM_cells_river.append(FM_cells_river_temp[j])
    
        # get coords of these cells
        FM_cells_river_coords = []
        for j in range(len(FM_cells_river)):
            FM_cells_river_coords.append(FMcoords[FM_cells_river[j]])
    
        # -------------------------------------------------------------------------------------------------------
        # Floodplain cells
        # -------------------------------------------------------------------------------------------------------
    
        # assign to this list all cells that are not river cells
        for j in range(len(FM_cells)):
            if FM_cells[j] not in FM_cells_river:
                FM_cells_floodplain.append(FM_cells[j])
    
        # get coords of these cells
        FM_cells_floodplain_coords = []
        for j in range(len(FM_cells_floodplain)):
            FM_cells_floodplain_coords.append(FMcoords[FM_cells_floodplain[j]])
    
    return FM_cells_river, FM_cells_floodplain, FM_cells_river_coords, FM_cells_floodplain_coords

def getFMcoordsCentroidsOLD(cell_pointsFM, FMcoords):
    """
    Get the centroids x and y coordinates of each FM polygon
    
    Input:  flow element nodes, list of (x,y) coordinates of each polygon (output of getFMcoords)
    
    Output: seperate arrays of x and y centroid coordinates
    """
    cell_sizes_fm = (cell_pointsFM > 0).sum(1)
    cell_count_fm = len(cell_pointsFM)
    
    # zero arrays
    centroids_x_fm = np.zeros(cell_count_fm)
    centroids_y_fm = np.zeros(cell_count_fm)
    
    # loop over all cells
    for i in range(cell_count_fm):
        
        # current cell size
        cell_size = cell_sizes_fm[i]
        
        # reset temporary values
        temp_value_for_A = 0
        temp_value_for_x = 0
        temp_value_for_y = 0
        
        # Calculate area of current cell:
        # loop over all vertices of current cell
        for j in range(cell_size):
            # checking if not at last vertex
            if j < (cell_size-1):
                # if not, calculation based on current and next vertex
                temp_value_for_A = temp_value_for_A + ((FMcoords[i][j][0] * FMcoords[i][j+1][1]) - \
                                                       (FMcoords[i][j+1][0] * FMcoords[i][j][1]))
            else:
                # otherwise, calculation based on current and first vertex (since this is the 'next' vertex in this case)
                temp_value_for_A = temp_value_for_A + ((FMcoords[i][j][0] * FMcoords[i][0][1]) - \
                                                       (FMcoords[i][0][0] * FMcoords[i][j][1]))
            
        cell_area = 0.5 * temp_value_for_A
        
        # Calculate x and y coordinate of cell centroid:
        # loop over all vertices of current cell
        for j in range(cell_size):
            # checking if not at last vertex
            if j < (cell_size-1):
                # if not, calculation based on current and next vertex
                temp_value_for_x = temp_value_for_x + ((FMcoords[i][j][0] + FMcoords[i][j+1][0]) * \
                                                       ((FMcoords[i][j][0] * FMcoords[i][j+1][1]) - \
                                                        (FMcoords[i][j+1][0] * FMcoords[i][j][1])))
                
                temp_value_for_y = temp_value_for_y + ((FMcoords[i][j][1] + FMcoords[i][j+1][1]) * \
                                                       ((FMcoords[i][j][0] * FMcoords[i][j+1][1]) - \
                                                        (FMcoords[i][j+1][0] * FMcoords[i][j][1])))
            else:
                # otherwise, calculation based on current and first vertex (since this is the 'next' vertex in this case)
                temp_value_for_x = temp_value_for_x + ((FMcoords[i][j][0] + FMcoords[i][0][0]) * \
                                                       ((FMcoords[i][j][0] * FMcoords[i][0][1]) - \
                                                        (FMcoords[i][0][0] * FMcoords[i][j][1])))
                
                temp_value_for_y = temp_value_for_y + ((FMcoords[i][j][1] + FMcoords[i][0][1]) * \
                                                       ((FMcoords[i][j][0] * FMcoords[i][0][1]) - \
                                                        (FMcoords[i][0][0] * FMcoords[i][j][1])))
            
        # use cell area and calculated temporary values to find the x and y coordinate of the cell centroids
        centroids_x_fm[i] = (1/(6*cell_area))*temp_value_for_x
        centroids_y_fm[i] = (1/(6*cell_area))*temp_value_for_y
        
    return centroids_x_fm, centroids_y_fm

def outputPCRtoFM(coupledCells, coupledCellsIndices, valueArray):
    """
    Creates an array with values, with output of 1 PCR cell distributed over its coupled FM cells. The indices of the array correspond to those of the input list coupledCells (which is output of the function coupleCellsPCRtoFM)
    
    Input:  - list with PCR cells and their coupled FM cells (output of coupleCellsPCRtoFMold)
            - double indices (array,column) of coupled cells (output of convertPCRindices)
            - array of values (can be obtained with BMI function get_var)
            
    Output: array with values that can be used as input to FM cells
    
    NOTE: currently the distribution over the FM cells is only based on the number of cells, this can be improved! (e.g. by taking into account cell size, elevation, bathymetry)
    """
    # create zero array
    input_values_fm = np.zeros(len(coupledCells))
    #  loop over all coupled PCR cells
    for i in range(len(coupledCells)):
        # find the number of coupled FM cells
        nr_coupled_FM_cells = len(coupledCells[i][1])
        # find PCR output value
        output_value_pcr = valueArray[coupledCellsIndices[i]]
        # calculate the value that each FM cell should receive as input
        input_values_fm[i] = output_value_pcr/nr_coupled_FM_cells
        
    return input_values_fm
	
# =============================================================================
    
def checkLocationCoupledPCRcells(model_pcr, CoupledPCRcellIndices):
    """
    This is mainly a debug-function to check whether PCR-cells are correctly coupled.
    
    Output:
		- boolean map that allows to assess the locations of the coupled PCR-GLOBWB cells.    
    """              
    
    # retrieving array representing PCR-GLOBWB data extent
    shapemap = np.copy(model_pcr.get_var(('routing', 'lddMap')))
    
    # creating empty array with PCR-GLOBWB extent 
    check_map = np.zeros([len(shapemap), len(shapemap[0])])

    # assigning True value for all coupled PCR-GLOBWB cells
    for i in range(len(CoupledPCRcellIndices)):
        check_map[CoupledPCRcellIndices[i]] = 1.    

    return check_map
	
# =============================================================================

def removeSmallWaterDepths(array_in, timestep, d_thres):
    """
    removing small delta depths as defined by water depth below user-specified threshold.
    returns array with same dimension as input that has adapted water depths as content.
    
    Input:
		- list with values
	Output:
		- list with values, values below threshold are set to zero
    """
    
    array_out = np.zeros(len(array_in))
    
    for i in xrange(len(array_in)):
        if array_in[i] > (d_thres / timestep):
            array_out[i] = array_in[i]
        else:
            pass
    
    return array_out

def getDFMcellsRiverFloodplain(coupling_pcr_DFM, hydroCoords, DFMareas, area_threshold='all', threshold_factor_3vertices=1.5, threshold_factor_4vertices=3.0):
    """
    Creates lists of DFM river and floodplain cells
    
    Input:  - list of all coupled PCR cells and their coupled DFM cells [3rd output of coupleAllCells]
            - list of DFM cell coordinates                              [output of gethydroCoords]
            - list of DFM cell areas                                    [output of calculateAllCellAreaSpherical or obtained through BMI get_var('ba')]
            - option to use the minimum area of all cells or only current PCR cells in loop as threshold base-level
            - multiplication factors for threshold base-level for cells with 3 and 4 vertices
    
    Output: - list of DFM river cells
            - list of DFM floodplain cells
            - list of coordinates of DFM river cells (useful for plotting/checking results)
            - list of coordinates of DFM floodplain cells (useful for plotting/checking results)
            - list of list of DFM river cells per coupled PCR cell
            - boolean array indicating whether coupled PCR cells have DFM river cells (True) or not (False)
            - list of list of DFM floodplain cells per coupled PCR cell
    """
    
    # create empty lists for filling in river cell indices and cpordinates
    DFM_cells_river = []
    DFM_cells_river_coords = []
    
    # create empty lists for filling in floodplain cell indices and cpordinates
    DFM_cells_floodplain = []
    DFM_cells_floodplain_coords = []
    
    # create empty list for filling in all coupled DFM cells
    DFM_cells_total = []
    
    # create False array for filling in True values if a river cell is found for a coupled PCR cell
    PCR_coupled_cell_has_river_boolean = np.zeros(len(coupling_pcr_DFM),dtype=bool)
    
    # create list filled with 'None' for filling in the river cells found for that coupled PCR cell
    PCR_coupled_cells_DFM_river_cells      = [[None]] * len(coupling_pcr_DFM)
    PCR_coupled_cells_DFM_floodplain_cells = [[None]] * len(coupling_pcr_DFM)
    
    for i in range(len(coupling_pcr_DFM)):
    
        # -------------------------------------------------------------------------------------------------------
        # Create a separate array for DFM cells of this coupled PCR cell only
        # -------------------------------------------------------------------------------------------------------
    
        # empty array
        DFM_cells = np.zeros(len(coupling_pcr_DFM[i][1]), dtype=np.int64)
    
        # loop over all DFM cells and assign to empty array
        for j in range(len(coupling_pcr_DFM[i][1])):
            DFM_cells[j] = coupling_pcr_DFM[i][1][j]
            DFM_cells_total.append(DFM_cells[j])
    
        # get areas of these cells
        DFM_cells_areas = []
        for j in range(len(DFM_cells)):
            DFM_cells_areas.append(DFMareas[DFM_cells[j]])
    
        # -------------------------------------------------------------------------------------------------------
        # Get threshold
        # -------------------------------------------------------------------------------------------------------
    
        if area_threshold == 'all':
    
            # set thresholds
            threshold_4vertices = min(DFMareas)*threshold_factor_4vertices
            threshold_3vertices = min(DFMareas)*threshold_factor_3vertices
    
        else:
    
            # set thresholds
            threshold_4vertices = min(DFM_cells_areas)*threshold_factor_4vertices
            threshold_3vertices = min(DFM_cells_areas)*threshold_factor_3vertices
    
        # -------------------------------------------------------------------------------------------------------
        # River cells
        # -------------------------------------------------------------------------------------------------------
        
        # create empty list to fill in all river cells for current coupled PCR cell
        DFM_cells_river_temp = []
    
        # assign to the list all cells that have 4 vertices and that fall below the relevant threshold
        DFM_cells_river_temp_1 = DFM_cells[np.where(DFM_cells_areas<=threshold_4vertices)[0]]
        for j in range(len(DFM_cells_river_temp_1)):
            if len(hydroCoords[DFM_cells_river_temp_1[j]]) == 4:
                DFM_cells_river_temp.append(DFM_cells_river_temp_1[j])
    
        # and all cells that fall below the threshold for polygons with 3 vertices
        DFM_cells_river_temp_2 = DFM_cells[np.where(DFM_cells_areas<=threshold_3vertices)[0]]
        for j in range(len(DFM_cells_river_temp_2)):
            if DFM_cells_river_temp_2[j] not in DFM_cells_river_temp:
                DFM_cells_river_temp.append(DFM_cells_river_temp_2[j])
        
        # assign this to list of DFM river cells
        for j in range(len(DFM_cells_river_temp)):
            DFM_cells_river.append(DFM_cells_river_temp[j])
    
        # if a river cell was found, assign True value (1) to list and add the cell indices to the relevant list
        if len(DFM_cells_river_temp) >= 1:
            PCR_coupled_cell_has_river_boolean[i] = 1
            PCR_coupled_cells_DFM_river_cells[i]   = DFM_cells_river_temp
        
        # -------------------------------------------------------------------------------------------------------
        # Floodplain cells
        # -------------------------------------------------------------------------------------------------------
        
        # if a river cell was found for this coupled PCR cell, it will probably also have floodplain cells
        if PCR_coupled_cell_has_river_boolean[i] == 1:
            
            # create empty list to fill in all floodplain cells for current coupled PCR cell
            DFM_cells_floodplain_temp = []
            
            # loop over all current DFM cells and assign to the list all cells that are not river cells
            for k in range(len(DFM_cells)):
                if DFM_cells[k] not in DFM_cells_river_temp:
                    DFM_cells_floodplain_temp.append(DFM_cells[k])
            
            # and assign these to the list of floodplain cells for each coupled PCR cell
            PCR_coupled_cells_DFM_floodplain_cells[i] = DFM_cells_floodplain_temp
        
        # if not river cell was found, all current DFM cells are assigned as floodplain cells
        else:
            PCR_coupled_cells_DFM_floodplain_cells[i] = DFM_cells
    
    # -------------------------------------------------------------------------------------------------------
    # River and floodplain cell coordinates
    # -------------------------------------------------------------------------------------------------------
    
    # get coords of river cells
    for i in range(len(DFM_cells_river)):
        DFM_cells_river_coords.append(hydroCoords[DFM_cells_river[i]])
    
    # assign to the list all cells that are not river cells
    for i in range(len(DFM_cells_total)):
        if DFM_cells_total[i] not in DFM_cells_river:
            DFM_cells_floodplain.append(DFM_cells_total[i])
    
    # get coordinates of these cells
    for i in range(len(DFM_cells_floodplain)):
        DFM_cells_floodplain_coords.append(hydroCoords[DFM_cells_floodplain[i]])
    
    return DFM_cells_river, DFM_cells_floodplain, DFM_cells_river_coords, DFM_cells_floodplain_coords, \
	       PCR_coupled_cells_DFM_river_cells, PCR_coupled_cell_has_river_boolean, PCR_coupled_cells_DFM_floodplain_cells

def plotGridValues(ax, cell_coords, cell_values, model_name, limits=None, mask_values=True, cell_cmap='jet', cell_aLFPha=1.0, \
                        cell_color='black', cells_colorbar=True, nr_ticks=10, padding_colorbar=None, label_colorbar=None, \
                        padding_label_colorbar=None, cell_edgecolor='none', cell_linewidth=1, extend_colorbar='neither', \
                        x1=None, x2=None, y1=None, y2=None):
    """
    Creates a 2D colour plot of cell values on the grid with specified coordinates.
    The function is created specifically to plot water depths, but can also plot other values on the grid. 
    It automatically adjusts the input array values by removing negative values and masking zero values,
    to allow for a clear distinction between wet and dry cells.
    
    Input:  - empty plot call 		   	[ax]
            - cell coordinates of grid 	[output of gethydroCoords or PCRcoords]
            - cell values              	[e.g. output of BMI function get_var('s1')]
            - model name               	[e.g. 'DFM' or 'PCR']
            - limits for plotted values	[min. and max. value to be plotted, colorbar is adjusted accordingly]
            - there are many other optional input variables that control the outcome of the plot
    
    Output: plot of model grid with coloured cell values
    
    NOTE: adjusted from Gena's 'plot_grid' function.
    """
    # create a collection of all polygons and set the colormap
    polygons = matplotlib.collections.PolyCollection(cell_coords, cmap=cell_cmap)
    polygons.set_aLFPha(cell_aLFPha)
    
    # check if model name suggests data from PCR-GLOBWB (which has to be adjusted)
    if model_name not in 'PCR-GLOBWB':
        pass
    else:	
        temp = np.zeros(len(cell_coords))
        count = 0
        for i in range(len(cell_values)):
            for j in range(len(cell_values[0])):
                if cell_values[i][j] != -999:
                    temp[count] = cell_values[i][j]
                    count += 1
        cell_values = temp
    # set all negative values to zero (array can contain very small negative values due to numerics involved)
    cell_values[np.where(cell_values<0)] = 0
    if mask_values:
        # mask all zero values (to be able to plot these using a specified color)
        cell_values = np.ma.masked_values(cell_values, 0)
    # pass the values to the collection of polygons
    polygons.set_array(cell_values)
    if limits != None:
        polygons.set_clim(vmin=limits[0], vmax=limits[1])
    
    if cells_colorbar:
        if padding_colorbar == None:
            cb = plt.colorbar(polygons, ax=ax, extend=extend_colorbar, pad=0.01)
        else:
            cb = plt.colorbar(polygons, ax=ax, extend=extend_colorbar, pad=padding_colorbar)
        if limits != None:
            cb.set_ticks(np.linspace(limits[0], limits[1], nr_ticks+1, endpoint=True), update_ticks=True)
        if label_colorbar == None:
            if padding_label_colorbar == None:
                cb.set_label('cell values')
            else:
                cb.set_label('cell values', labeLFPad=padding_label_colorbar)
        else:
            if padding_label_colorbar == None:
                cb.set_label(label_colorbar)
            else:
                cb.set_label(label_colorbar, labeLFPad=padding_label_colorbar)
    
    polygons.set_edgecolor(cell_edgecolor)
    polygons.set_linewidth(cell_linewidth)
    ax.add_collection(polygons)
    
    # set axis
    if (x1==None) and (x2==None) and (y1==None) and (y2==None):
        ax.autoscale()
    else:
        ax.axis((x1,x2,y1,y2))

def zeroMapArray(PCRmap, missingValuesOutput=-999, missingValuesMap=255):
    """
    Creates a numpy array with zero values at every map location with no missing value. This can be used to fill in data that will be extracted from a model.
    
    Input:  - PCR map showing location of (non-)missing values (e.g. landmask)
            - missing values for output numpy array (optional, default at -999)
            - missing values of map (optional, default = 255 for landmask)
            
    Output: numpy array with zero values at non-missing values of map and output missing values at missing values of map
    """
    reference_map    = pcr.readmap(PCRmap)
    reference_map_np = pcr.pcr2numpy(reference_map, missingValuesMap)
    output_map_np    = np.zeros([len(reference_map_np), len(reference_map_np[0])])
    
    for i in range(len(reference_map_np)):
        for j in range(len(reference_map_np[0])):
            if reference_map_np[i][j] == missingValuesMap:
                output_map_np[i][j] = missingValuesOutput
            
    return output_map_np
	
def gethydroCoordsZ(cell_pointsDFM, zDFM):
    """
    Get all z-coordinates of a DFM grid.
    
    Input:  flow element nodes, z-values of an DFM grid
    
    Output: list of z-coordinates of each polygon
    """
    # get all cell sizes and total cell count
    cell_sizes_DFM = (cell_pointsDFM > 0).sum(1)
    cell_count_DFM = len(cell_pointsDFM)
    # Create list of list with all DFM cell coordinates
    # initial empty list
    z_cell_coords_DFM = []
    # loop over all cells
    for i in range(cell_count_DFM):
        # get current cell size
        cell_size = cell_sizes_DFM[i]
        # empty list
        cell_coords = []
        # loop over all cell vertices
        for j in range(cell_size):
            # get coordinates of current cell vertex and append to list
            cell_coords.append((zDFM[cell_pointsDFM[i, j] - 1]))
        # add list of coordinates of current cell to total list
        z_cell_coords_DFM.append(cell_coords)
        
    return z_cell_coords_DFM

def calculateSingleCellAreaSpherical(cellCoords):
    """
    Calculates the area of a single cell, based on its spherical coordinates (latitude, longitude).
    
    Input:  spherical coordinate of a cell (list with (x,y) coordinates)
    
    Output: area of cells
    
    References: - http://gis.stackexchange.com/questions/8495/converting-longitude-and-latitude-coordinates-to-square-miles
                - http://dev.openlayers.org/releases/OpenLayers-2.10/lib/OpenLayers/Geometry/LinearRing.js
                - Robert. G. Chamberlain and William H. Duquette, "Some Algorithms for Polygons on a Sphere",
                  JPL Publication 07-03, Jet Propulsion Laboratory, Pasadena, CA, June 2007
                  http://trs-new.jpl.nasa.gov/dspace/handle/2014/40409
    
    NOTE: This is an approximation and not an exact solution!
    """
    # temporary value used in calculation
    temp_value_for_A = 0
    # cell size
    cell_size = len(cellCoords)
    # loop over all vertices
    for i in range(cell_size):
        # check if not at last vertex
        if i < (cell_size-1):
            # if not, calculation based on current and next vertex
            temp_value_for_A = temp_value_for_A + (np.radians(cellCoords[i+1][0] - cellCoords[i][0]) * \
                                                   (2 + (np.sin(np.radians(cellCoords[i][1]))) + (np.sin(np.radians(cellCoords[i+1][1])))))
        else:
            # otherwise, based on current and first vertex (since this is 'next' vertex)
            temp_value_for_A = temp_value_for_A + (np.radians(cellCoords[0][0] - cellCoords[i][0]) * \
                                                   (2 + (np.sin(np.radians(cellCoords[i][1]))) + (np.sin(np.radians(cellCoords[0][1])))))
    # calculate area of cell
    cell_area = abs(temp_value_for_A * 6378137 * 6378137/2)
    
    return cell_area
	
def calculateSingleCellArea(cellCoords):
    """
    Calculates the area of a single cell, based on its coordinates
    
    Input:  coordinates of one cell (list with (x,y) coordinates)
    
    Output: area of cell
    
    NOTE: area depends on the coordinate system used!
    """
    # temporary value used in calculation
    temp_value_for_A = 0
    # cell size
    cell_size = len(cellCoords)
    # loop over all vertices
    for i in range(cell_size):
        # check if not at last vertex
        if i < (cell_size-1):
            # if not, calculation based on current and next vertex
            temp_value_for_A = temp_value_for_A + ((cellCoords[i][0] * cellCoords[i+1][1]) -\
                                                   (cellCoords[i+1][0] * cellCoords[i][1]))
        else:
            # otherwise, based on current and first vertex (since this is 'next' vertex)
            temp_value_for_A = temp_value_for_A + ((cellCoords[i][0] * cellCoords[0][1]) -\
                                                   (cellCoords[0][0] * cellCoords[i][1]))
    cell_area = 0.5 * temp_value_for_A
    
    return cell_area
	
def convertPCRindices(indexList, PCRmap, missingValues=255):
    """
    Converts a list with single indices referring to a point on a PCR map to double indices (array,column).
    This is required to call PCR cells from maps read during the coupling, since these work with double (array,column) indices,
    and this basically reverses the single index that was created in the coupling functions to find coupled cells more easily.
    
    Input:  - list with single indices (output of coupleAllCells[3] or [4])
            - reference PCR map (e.g. landmask)
            - missing values of map (optional, default = 255 for landmask)
    
    Output: list with double indices (array,column)
    """
    # get reference map
    try:
        # if map is given as path to PCR map
        reference_map = pcr.readmap(PCRmap)
    except:
        # if map is an already loaded PCR map
        reference_map = PCRmap
    reference_map= pcr.cover(reference_map,0)
    #pdb.set_trace()
    # convert map to numpy array
    reference_map_np = pcr.pcr2numpy(reference_map, missingValues)
    # find the number (and cumulative number) of non-missing values per row
    number_non_missing_each_row            = np.zeros(len(reference_map_np))
    cumulative_number_non_missing_each_row = np.zeros(len(reference_map_np))
    for i in range(len(reference_map_np)):
        number_non_missing_each_row[i]            = (len(np.where(reference_map_np[i] != missingValues)[0]))
        cumulative_number_non_missing_each_row[i] = (np.sum(number_non_missing_each_row))
        #print reference_map_np[i], cumulative_number_non_missing_each_row[i], number_non_missing_each_row[i]
    
    # create zero arrays for filling in values
    indices_array  = np.zeros(len(indexList))
    indices_column = np.zeros(len(indexList))
    # loop over all single indices
    for i in range(len(indexList)):
        # get the current single index
        try:
            # if single index list is given as list of lists (from function coupleAllCells[3])
            current_index = indexList[i][0]
            print 'option 1 chosen'
        except:
            # if single index list is given as list (from function coupleAllCells[4])
            current_index = indexList[i]
            #print 'option 2 chosen'
        # find the array number of the current index
        temp_array = np.where((cumulative_number_non_missing_each_row - current_index) > 0)[0][0]
        # find the number of non-missing values that have to be looped over to find the right column index
        temp_remaining = (current_index - cumulative_number_non_missing_each_row[temp_array-1]) + 1 # old version, works for Amazon, but not for Elbe
        #temp_remaining = (current_index + cumulative_number_non_missing_each_row[temp_array]) # new version, works for Elbe, but not Amazon and not for Niger
        #print i, cumulative_number_non_missing_each_row, current_index, temp_remaining, temp_array
        # find the column number of the current index
        temp_count = 0
        for j in range(len(reference_map_np[temp_array])):
            if reference_map_np[temp_array][j] != missingValues:
                #print 'temp_count: ',temp_count
                #print 'temp_remaining:', temp_remaining
                temp_count += 1
                if temp_count == temp_remaining:
                    #print 'temp_count == temp_remaining'
                    temp_column = j
                    break
        # assign the found values to the total arrays
        indices_array[i]  = temp_array
        indices_column[i] = temp_column
        
    #pdb.set_trace()
    # create a single list with combined (array,column) indices
    indexListDouble = zip(indices_array,indices_column)
    
    return indexListDouble
	
def calculateAllCellEmptyVolumes(cellCoordsZ, cellAreas, BedlevType=None, option='adjusted'):
    """
    Calculates the volume of all cells up until the level that the elevation values of all vertices are equal,
    (so the cell is considered "full") based on their (spherical) coordinates.
    
    Input:  - z-coordinates of all cells
            - cell areas of all cells
            - BedlevType (default=None). If not specified, z-values will be taken as mean of closest min or max z-value.
              If specified, it follows DFM .mdu file numbers, so 3 calculated mean values of all z-values and 4 uses min of all values.
              Currently only working for BedlevType 3 or 4 (or default=None of course)!
			- option keyword, either 'original' or 'adjusted' (default)
              (this specifies which z-values will be used for the calculation of cell volumes)
            
    Output: - (empty) volume of all cells
            - shape identifier for all cells
            - list with adjusted z-values
    
    NOTE: These are not the exact volumes! A number of simplifications/approximations are used!
    """
    
    # calculate number of cells
    cell_count = len(cellCoordsZ)
    
    # zero array for filling in cell volumes
    cell_volumes = np.zeros(cell_count)
    
    # zero array for filling in cell structure identifier
    cell_shapes = np.zeros(cell_count)
    
    # empty list for setting adjusted z-values
    adjusted_z_values_list = []
    
    # loop over all cells
    for i in range(cell_count):
        
        # current cell
        current_cell = cellCoordsZ[i]
        
        # current cell length
        current_cell_length = len(current_cell)
        
        # zero array and filling in values (and to more easily access values, using an array instead of list)
        temp_z_values = np.zeros(current_cell_length)
        for p in range(current_cell_length):
            temp_z_values[p] = current_cell[p]
        
        # mean, min and max value of current cell
        temp_mean = np.mean(current_cell)
        temp_min  = np.min(current_cell)
        temp_max  = np.max(current_cell)
        
        # perform all calculations (if min != max, since that means there are elevation differences in the cell)
        if temp_min != temp_max:
            
            # index of min/max value of current cell
            temp_min_index = np.where(current_cell==min(current_cell))[0][0]
            temp_max_index = np.where(current_cell==max(current_cell))[0][0]
            
            # empty arrays for filling in differences
            temp_diff_min = np.zeros(len(current_cell))
            temp_diff_max = np.zeros(len(current_cell))
            
            # empty array stating where point belongs to (0/False = min, 1/True = max)
            temp_diff_array = np.zeros(len(current_cell),dtype=bool)
            
            # loop over all cell points to find differences
            for j in range(current_cell_length):
                
                # calculate z-values difference between min and other points
                temp_diff_min[j] = current_cell[j] - temp_min
                
                # calculate z-values difference between max and other points
                temp_diff_max[j] = temp_max - current_cell[j]
                
                # assign every point to closest match (default is min, so only checking for max)
                if temp_diff_max[j] < temp_diff_min[j]:
                    temp_diff_array[j] = 1
            
            # adjusting z-values
            if BedlevType == None:
                temp_z_values[~temp_diff_array] = np.mean(temp_z_values[~temp_diff_array])
                temp_z_values[temp_diff_array] = np.mean(temp_z_values[temp_diff_array])
            elif BedlevType == 3:
                temp_z_values[~temp_diff_array] = np.mean(temp_z_values)
                temp_z_values[temp_diff_array] = np.mean(temp_z_values[temp_diff_array])
            elif BedlevType == 4:
                temp_z_values[~temp_diff_array] = np.min(temp_z_values)
                temp_z_values[temp_diff_array] = np.mean(temp_z_values[temp_diff_array])
            
            adjusted_z_values_list.append(temp_z_values)
            
            # Calculating volume
            
            # check specified option to see which z-values to use
            if option=='original':
                
                pass
            
            elif option=='adjusted':
                
                # calculate new min/max values, based on adjusted z-values
                temp_min = np.min(temp_z_values)
                temp_max = np.max(temp_z_values)
            
            # checking cell shape (3 or 4 vertices)
            if current_cell_length == 4:
                
                # cell with 1 high vertex, 3 low
                if np.sum(temp_diff_array) == 1:
                    cell_shapes[i] = 1
                    cell_volumes[i] = (2./3) * cellAreas[i] * (temp_max - temp_min)
                    
                # cell with 3 high vertices, 1 low
                if np.sum(temp_diff_array) == 3:
                    cell_shapes[i] = 2
                    cell_volumes[i] = (1./3) * cellAreas[i] * (temp_max - temp_min)
                    
                # cell with 2 high vertices, 2 low
                if np.sum(temp_diff_array) == 2:
                    
                    # Check if similar points lie next to each other (or not)
                    temp_check_value = 0
                    # loop over all cell vertices
                    for k in range(current_cell_length):
                        # check if not at last index
                        if k != (current_cell_length-1):
                            # if current point equals next point, similar points lie next to each other
                            if temp_diff_array[k] == temp_diff_array[k+1]:
                                temp_check_value = 1
                                break
                        # if at last index, use first index as next point
                        else:
                            if temp_diff_array[k] == temp_diff_array[0]:
                                temp_check_value = 1
                    
                    # cells with similar points next to each other
                    if temp_check_value == 1:
                        cell_shapes[i] = 3
                        cell_volumes[i] = (1./2) * cellAreas[i] * (temp_max - temp_min)
                    # or similar points opposite to each other (diagonally for 4 vertices)
                    else:
                        cell_shapes[i] = 4
                        cell_volumes[i] = (2./3) * cellAreas[i] * (temp_max - temp_min)
                        
            elif current_cell_length == 3:
                
                # cell with 1 high vertex, 2 low
                if np.sum(temp_diff_array) == 1:
                    cell_shapes[i] = 5
                    cell_volumes[i] = (2./3) * cellAreas[i] * (temp_max - temp_min)
                    
                # cell with 2 high vertices, 1 low
                if np.sum(temp_diff_array) == 2:
                    cell_shapes[i] = 6
                    cell_volumes[i] = (1./3) * cellAreas[i] * (temp_max - temp_min)
                    
        # but if min and max are equal, all cell vertices have the same z-value and there is no "empty volume"
        else:
            
            # add current z-values to adjusted z-values list
            adjusted_z_values_list.append(current_cell)
            
            # check cell length to assign cell shape identifier
            if current_cell_length == 4:
                cell_shapes[i] = 7
            else:
                cell_shapes[i] = 8
                
    return cell_volumes, cell_shapes, adjusted_z_values_list
	
def calculateAllCellAreaSpherical(cellCoords):
    """
    Calculates the area of all cells, based on their spherical coordinates (latitude, longitude).
    
    Input:  spherical coordinates of all cells (list of list with (x,y) coordinates)
    
    Output: area of all cells (array)
    
    References: - http://gis.stackexchange.com/questions/8495/converting-longitude-and-latitude-coordinates-to-square-miles
                - http://dev.openlayers.org/releases/OpenLayers-2.10/lib/OpenLayers/Geometry/LinearRing.js
                - Robert. G. Chamberlain and William H. Duquette, "Some Algorithms for Polygons on a Sphere",
                  JPL Publication 07-03, Jet Propulsion Laboratory, Pasadena, CA, June 2007
                  http://trs-new.jpl.nasa.gov/dspace/handle/2014/40409
    
    NOTE: This is an approximation and not an exact solution!
    """
    # calculate number of cells
    cell_count = len(cellCoords)
    # zero array for filling in cell areas
    cell_areas = np.zeros(cell_count)
    # loop over all cells
    for i in range(cell_count):
        # current cell size
        cell_size = len(cellCoords[i])
        # temporary value used in calculation
        temp_value_for_A = 0
        # loop over all vertices
        for j in range(cell_size):
            # check if not at last vertex
            if j < (cell_size-1):
                # if not, calculation based on current and next vertex
                temp_value_for_A = temp_value_for_A + (np.radians(cellCoords[i][j+1][0] - cellCoords[i][j][0]) * \
                                                       (2 + (np.sin(np.radians(cellCoords[i][j][1]))) + (np.sin(np.radians(cellCoords[i][j+1][1])))))
            else:
                # otherwise, based on current and first vertex (since this is 'next' vertex)
                temp_value_for_A = temp_value_for_A + (np.radians(cellCoords[i][0][0] - cellCoords[i][j][0]) * \
                                                       (2 + (np.sin(np.radians(cellCoords[i][j][1]))) + (np.sin(np.radians(cellCoords[i][0][1])))))
        # calculate area of cell
        cell_area = abs(temp_value_for_A * 6378137 * 6378137/2)
        # assign area of current cell to array
        cell_areas[i] = cell_area
    
    return cell_areas
	
def plotCoupledPCRCell(PCRcellIndex, coupling_pcr_DFM, PCRcoords, hydroCoords, bufferPlotAxis=0.1):
    """
    Creates 2D a plot of one coupled PCR cell and its coupled DFM cells
    
    Input:  - index of PCR cell (index in input list)
            - list of all coupled PCR cells										(output of coupleAllCells[3])
            - list of (x,y) coordinates of each DFM polygon						(output of gethydroCoords)
            - list of (x,y) coordinates of each PCR polygon						(output of getPCRcoords)
            - buffer around axis limit (to be able to see all vertices of the cell, optional, default is 0.1)
            
    Output: plot of coupled PCR cell with its coupled DFM cells
    """
    # Get data from list:
    # get cells based on chosen index
    if len(coupling_pcr_DFM) > 5:
        # input is only the relevant list with coupled PCR cells
        PCRcell = coupling_pcr_DFM[PCRcellIndex][0]
        DFMcells = coupling_pcr_DFM[PCRcellIndex][1]
    else:
        # input is full output of function coupleAllCells (list of 5 lists)
        PCRcell = coupling_pcr_DFM[3][PCRcellIndex][0]
        DFMcells = coupling_pcr_DFM[3][PCRcellIndex][1]
    
    # get coords of these cells
    PCRcellCoords = PCRcoords[PCRcell]
    DFMcellCoords  = [hydroCoords[i] for i in DFMcells]
    
    # find plot limits
    xMinPlot = min(PCRcellCoords)[0]
    xMaxPlot = max(PCRcellCoords)[0]
    yMinPlot = min(PCRcellCoords)[1]
    yMaxPlot = max(PCRcellCoords)[1]
    
    # Plot:
    # calling empty plot
    f, ax = plt.subplots(1, 1)
    
    # set axis
    ax.set_xlim(xMinPlot - bufferPlotAxis, xMaxPlot + bufferPlotAxis)
    ax.set_ylim(yMinPlot - bufferPlotAxis, yMaxPlot + bufferPlotAxis)
    
    # create polygons
    test_polygons_PCR_coupled = matplotlib.collections.PolyCollection([PCRcellCoords])
    test_polygons_DFM_coupled  = matplotlib.collections.PolyCollection(DFMcellCoords)
    
    # set polygon appearance
    test_polygons_PCR_coupled.set_aLFPha(1)
    test_polygons_PCR_coupled.set_facecolor('white')
    test_polygons_PCR_coupled.set_edgecolor('black')
    
    test_polygons_DFM_coupled.set_aLFPha(0.2)
    test_polygons_DFM_coupled.set_facecolor('black')
    test_polygons_DFM_coupled.set_edgecolor('black')
    
    # add PCR polygons to plot
    ax.add_collection(test_polygons_PCR_coupled)
    ax.add_collection(test_polygons_DFM_coupled)
	
def gethydroCoords(cell_pointsDFM, xDFM, yDFM):
    """
    Get all relevant coordinates of a DFM grid.
    
    Input:  flow element nodes, x, y values of an DFM grid
    
    Output: list of list of (x,y) coordinates of each polygon
    """
    # get all cell sizes and total cell count
    cell_sizes_DFM = (cell_pointsDFM > 0).sum(1)
    cell_count_DFM = len(cell_pointsDFM)
    # Create list of list with all DFM cell coordinates
    # initial empty list
    all_cell_coords_DFM = []
    # loop over all cells
    for i in range(cell_count_DFM):
        # get current cell size
        cell_size = cell_sizes_DFM[i]
        # empty list
        cell_coords = []
        # loop over all cell vertices
        for j in range(cell_size):
            # get coordinates of current cell vertex and append to list
            cell_coords.append((xDFM[cell_pointsDFM[i, j] - 1], yDFM[cell_pointsDFM[i, j] - 1]))
        # add list of coordinates of current cell to total list
        all_cell_coords_DFM.append(cell_coords)
        
    return all_cell_coords_DFM
