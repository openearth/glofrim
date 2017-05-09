

def inundatedArea_noRFS(model, bottom_lvl, CouplePCR2model, threshold_inundated_depth, CellAreaSpherical, model_type, waterDepth=0.):
    """
    Only required for 2way-coupling!    
    calculating fraction and area of (inundated) water per coupled PCR cell 
    """
    if model_type == 'DFM':
        current_water_depth = model.get_var('s1') - bottom_lvl
    else:
        current_water_depth = np.copy(waterDepth)
    
    inundated_area_model_2_PCR_coupled = np.zeros(len(CouplePCR2model))
    
    # loop over all coupoled PCR cells and fill in and/or calculate values
    for i in range(len(CouplePCR2model)):
        
        # get total area of all flooded DFM cells coupled to current PCR cell
        
        # temporary varialbe for storing total ara of flooded DFM cells coupled to current PCR cell
        temp_inundated_area = 0.
        
        # loop over all coupled DFM cells
        for j in range(len(CouplePCR2model[i][1])):
            
            # get current DFM cell index
            current_cell_index = CouplePCR2model[i][1][j]
            
            # check if water depth of current DFM cell is above chosen threshold
            if current_water_depth[current_cell_index] > threshold_inundated_depth:
                
                # if so, add cell are to temp var
                temp_inundated_area += CellAreaSpherical[current_cell_index]
                
        # at end of loop, assign temporary variable to array storing inundated areas
        inundated_area_model_2_PCR_coupled[i] = temp_inundated_area
        
    return inundated_area_model_2_PCR_coupled

# =============================================================================

def inundatedArea(model, bottom_lvl, CouplePCR2model, boolean_river_cell_in_coupled_PCR_cell, river_cells_per_coupled_PCR_cell, \
                    DFM_floodplain_cells_per_coupled_PCR_cell, threshold_inundated_depth_rivers, threshold_inundated_depth_floodplains, CellAreaSpherical):
    """
    Only required for 2way-coupling!     
    calculating fraction and area of (inundated) water per coupled PCR cell 
    """
    current_water_depth_DFM = model.get_var('s1') - bottom_lvl
    
    inundated_area_rivers_model_2_PCR_coupled          = np.zeros(len(CouplePCR2model))
    inundated_area_floodplains_model_2_PCR_coupled     = np.zeros(len(CouplePCR2model))
    
    # loop over all coupled PCR cells
    for i in range(len(CouplePCR2model)):
        
        # check if there are any river cells:
        if boolean_river_cell_in_coupled_PCR_cell[i]:
            # calculating fraction and area of wet DFM river cells per coupled PCR cell
            temp_inundated_area = 0.
            # 1. RIVER CELLS
            for j in range(len(river_cells_per_coupled_PCR_cell[i])):
                # get current DFM cell index
                current_cell_index = river_cells_per_coupled_PCR_cell[i][j]
                # check if water depth of current DFM cell is above chosen threshold
                if current_water_depth_DFM[current_cell_index] > threshold_inundated_depth_rivers:
                    # if so, add cell area to temporary variable
                    temp_inundated_area += CellAreaSpherical[current_cell_index]
            # at end of loop, assign temporary variable to array storing inundated areas
            inundated_area_rivers_model_2_PCR_coupled[i] = temp_inundated_area
            # calculating fraction and area of inundated DFM floodplain cells per coupled PCR cell
            temp_inundated_area = 0.
            # 2. FLOODPLAIN CELLS
            for j in range(len(DFM_floodplain_cells_per_coupled_PCR_cell[i])):
                # get current DFM cell index
                current_cell_index = DFM_floodplain_cells_per_coupled_PCR_cell[i][j]
                # check if water depth of current DFM cell is above chosen threshold
                if current_water_depth_DFM[current_cell_index] > threshold_inundated_depth_floodplains:
                    # if so, add cell area to temporary variable
                    temp_inundated_area += CellAreaSpherical[current_cell_index]
            # at end of loop, assign temporary variable to array storing inundated areas
            inundated_area_floodplains_model_2_PCR_coupled[i] = temp_inundated_area
            
        # if there are no river cells:
        else:
            temp_inundated_area = 0.
            # loop over all coupled DFM cells
            for j in range(len(CouplePCR2model[i][1])):
                # get current DFM cell index
                current_cell_index = CouplePCR2model[i][1][j]
                # check if water depth of current DFM cell is above chosen threshold
                if current_water_depth_DFM[current_cell_index] > threshold_inundated_depth_floodplains:
                    # if so, add cell area to temporary variable
                    temp_inundated_area += CellAreaSpherical[current_cell_index]
            # at end of loop, assign temporary variable to array storing inundated areas
            inundated_area_floodplains_model_2_PCR_coupled[i] = temp_inundated_area
    
    return inundated_area_rivers_model_2_PCR_coupled, inundated_area_floodplains_model_2_PCR_coupled

# =============================================================================

def extractAndConvertVolumes(model, CouplePCR2model, inundated_area_model_2_PCR_coupled, CellAreaSpherical, model_type, waterDepth=0.):
    """
    Only required for 2way-coupling!
    extracting DFM water volume for each coupled PCR cell and converting these to water depths to add back to PCR
    Exists, but not actually used for 1way-coupling.    
    """
    
    if model_type == 'DFM':
        current_volume = model.get_var('vol1')
    else:
        current_volume = waterDepth.ravel() * CellAreaSpherical

    water_volume_model_2_PCR_coupled = np.zeros(len(CouplePCR2model))
    water_depths_model_2_PCR_coupled = np.zeros(len(CouplePCR2model))
    
    # loop over all coupled PCR cells and fill in and/or calculate values
    for i in range(len(CouplePCR2model)):
        
        # get total volume of all DFM cells coupled to current PCR cell
        water_volume_model_2_PCR_coupled[i] = np.sum(current_volume[CouplePCR2model[i][1]])
        
        # divide total volume by inundated DFM area to obtain water depth
        if inundated_area_model_2_PCR_coupled[i] > 0.:
            water_depths_model_2_PCR_coupled[i] = water_volume_model_2_PCR_coupled[i] / inundated_area_model_2_PCR_coupled[i]
        else:
            water_depths_model_2_PCR_coupled[i] = 0.
    
    return water_volume_model_2_PCR_coupled, water_depths_model_2_PCR_coupled
       
# =============================================================================