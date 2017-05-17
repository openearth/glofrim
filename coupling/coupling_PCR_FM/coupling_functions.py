# -*- coding: utf-8 -*-
"""
Set of functions related to the coupling of PCR-GLOBWB with Delft 3D Flexible Mesh and LISFLOOD-FP.
Unlike the "model_functions" where primarily the translation of water volumes is computed, the subsequent set
of functions is used to spatially couple the grids of the models.

@author: Jannis Hoch, Department of Physical Geography, Geosciences, Utrecht University (j.m.hoch@uu.nl)
@date: 20-03-2017
"""

import pcraster as pcr
import numpy as np
import os
import pdb
import spatialDataSet2PCR
import pcrGlobalGeometry
import pyproj as pyproj

import matplotlib.collections
import matplotlib.pyplot as plt

# =============================================================================

def getFMcoords(cell_pointsFM, xFM, yFM):
    """
    Get all relevant coordinates of a FM grid.
    
    Input:  flow element nodes, x, y values of an FM grid
    
    Output: list of list of (x,y) coordinates of each polygon
    """
    # get all cell sizes and total cell count
    cell_sizes_fm = (cell_pointsFM > 0).sum(1)
    cell_count_fm = len(cell_pointsFM)
    # Create list of list with all FM cell coordinates
    # initial empty list
    all_cell_coords_fm = []
    # loop over all cells
    for i in range(cell_count_fm):
        # get current cell size
        cell_size = cell_sizes_fm[i]
        # empty list
        cell_coords = []
        # loop over all cell vertices
        for j in range(cell_size):
            # get coordinates of current cell vertex and append to list
            cell_coords.append((xFM[cell_pointsFM[i, j] - 1], yFM[cell_pointsFM[i, j] - 1]))
        # add list of coordinates of current cell to total list
        all_cell_coords_fm.append(cell_coords)
        
    return all_cell_coords_fm
	
def getVerticesFromMidPoints(xCoord, yCoord, xResolution, yResolution, verbose = False, proj0 = pyproj.Proj(init='epsg:4326'), proj1 = pyproj.Proj(init='epsg:4326')):
    
    """
     For rectangular, regular grids, this function returns the x- and y-coordinates of the
     four corners as a list, similar to function gethydroCoords.
	 Furthermore reprojecting coordinates from projected to spherical coordinate system if required.
     
     Input:
     -----
     array or list with xCoords of centre points
     array of list with yCoords of centre points
     map resolution in x- and y-direction
     
     Output:
     ------
     list with length of cell number with four tupels containing corner coordinates
	
	"""
    if isinstance(xCoord, np.ndarray):
		xCoord= xCoord.tolist()
    if isinstance(yCoord, np.ndarray):
		yCoord= yCoord.tolist()
	
    coordVertices= []
    for lCnt in xrange(len(xCoord)):
        x_ll= xCoord[lCnt]-0.5*xResolution 
        y_ll= yCoord[lCnt]-0.5*yResolution
        x_ll_1, y_ll_1 = pyproj.transform(proj0, proj1, x_ll, y_ll)
        ll = (x_ll_1, y_ll_1)
              
        x_ul= xCoord[lCnt]-0.5*xResolution
        y_ul= yCoord[lCnt]+0.5*yResolution
        x_ul_1, y_ul_1 = pyproj.transform(proj0, proj1, x_ul, y_ul)
        ul = (x_ul_1, y_ul_1)
              
        x_ur= xCoord[lCnt]+0.5*xResolution
        y_ur= yCoord[lCnt]+0.5*yResolution
        x_ur_1, y_ur_1 = pyproj.transform(proj0, proj1, x_ur, y_ur)
        ur = (x_ur_1, y_ur_1)
              
        x_lr= xCoord[lCnt]+0.5*xResolution
        y_lr= yCoord[lCnt]-0.5*yResolution
        x_lr_1, y_lr_1 = pyproj.transform(proj0, proj1, x_lr, y_lr)
        lr = (x_lr_1, y_lr_1)
              
              
        coordVertices.append([ll, ul, ur, lr])
           
    if verbose == True:
        print 'length x coords', len(xCoord) 
        print 'length y coords', len(yCoord)
        print 'length model coords', len(coordVertices)
    
    return coordVertices
	
# =============================================================================
	
def getPCRcoords(PCRmap, missing_value_pcr=-999):
    """
    Get all vertices coordinates of a PCRaster map.

    Input:
	-----
	pcraster map (preferrably landmask)
	value for MV (optional, default at -999)

    Output:
	------
	list of (x,y) coordinates of each polygon

    """
    # Get coordinates as numpy array:
    # upper left coordinates
    pcr.setglobaloption("coorul")

    xcoord_pcr_ul_map           = pcr.xcoordinate(PCRmap)
    xcoord_pcr_ul_np            = pcr.pcr2numpy(xcoord_pcr_ul_map,missing_value_pcr)

    ycoord_pcr_ul_map           = pcr.ycoordinate(PCRmap)
    ycoord_pcr_ul_np            = pcr.pcr2numpy(ycoord_pcr_ul_map,missing_value_pcr)

    # lower right coordinates
    pcr.setglobaloption("coorlr")

    xcoord_pcr_lr_map           = pcr.xcoordinate(PCRmap)
    xcoord_pcr_lr_np            = pcr.pcr2numpy(xcoord_pcr_lr_map,missing_value_pcr)

    ycoord_pcr_lr_map           = pcr.ycoordinate(PCRmap)
    ycoord_pcr_lr_np            = pcr.pcr2numpy(ycoord_pcr_lr_map,missing_value_pcr)

    # centroid coordinates
    pcr.setglobaloption("coorcentre")

    xcoord_pcr_centr_map        = pcr.xcoordinate(PCRmap)
    xcoord_pcr_centr_np         = pcr.pcr2numpy(xcoord_pcr_centr_map,missing_value_pcr)

    ycoord_pcr_centr_map        = pcr.ycoordinate(PCRmap)
    ycoord_pcr_centr_np         = pcr.pcr2numpy(ycoord_pcr_centr_map,missing_value_pcr)

    # Construct collection of polygon vertices:
    # number of arrays/elements to loop over and/or construct new arrays
    array_count_pcr             = len(ycoord_pcr_lr_np)
    elements_per_array_pcr      = np.size(ycoord_pcr_lr_np)/array_count_pcr
    nonmiss_val_per_array_pcr   = np.sum(ycoord_pcr_lr_np!=missing_value_pcr)


    # filling empty arrays while looping over data
    i, j = np.where(xcoord_pcr_lr_np != missing_value_pcr)
    xcoord_pcr_lr_np_nonmiss = xcoord_pcr_lr_np[i, j]
    xcoord_pcr_ul_np_nonmiss = xcoord_pcr_ul_np[i, j]
    xcoord_pcr_ll_np_nonmiss = xcoord_pcr_ul_np[i, j]
    xcoord_pcr_ur_np_nonmiss = xcoord_pcr_lr_np[i, j]

    ycoord_pcr_lr_np_nonmiss = ycoord_pcr_lr_np[i, j]
    ycoord_pcr_ul_np_nonmiss = ycoord_pcr_ul_np[i, j]
    ycoord_pcr_ll_np_nonmiss = ycoord_pcr_lr_np[i, j]
    ycoord_pcr_ur_np_nonmiss = ycoord_pcr_ul_np[i, j]

    xcoord_pcr_centr_np_nonmiss = xcoord_pcr_centr_np[i, j]
    ycoord_pcr_centr_np_nonmiss = ycoord_pcr_centr_np[i, j]

    # empty collection for polygons
    ll = zip(xcoord_pcr_ll_np_nonmiss, ycoord_pcr_ll_np_nonmiss)
    lr = zip(xcoord_pcr_lr_np_nonmiss, ycoord_pcr_lr_np_nonmiss)
    ur = zip(xcoord_pcr_ur_np_nonmiss, ycoord_pcr_ur_np_nonmiss)
    ul = zip(xcoord_pcr_ul_np_nonmiss, ycoord_pcr_ul_np_nonmiss)
    # wrap all cell coordinates into a list of lists (one list per cell, with multiple tuples per cell corner)
    all_cell_coords_pcr = [[ll[i], lr[i], ur[i], ul[i]] for i in range(len(ll))]

    return all_cell_coords_pcr
	
# =============================================================================

def coupleAllCells(hydroCoords, PCRcoords, no_couple_value=None):
    """
    Creates a list of coupled cells by finding a coupled PCR cell for each DFM cell.
	This function is only required to compute array that can be used as input for
	plotting of spatially coupled grids.
	Actual coupling performed in function "assignPCR2cells".
    
    Input:
	-----
	list of (x,y) coordinates of each hydrodynamic cell
    list of (x,y) coordinates of each PCR polygon 
    value that is used to indicate no coupled cell was found (optional, default = None)
            
    Output:
	------
	list of (x,y) coordinates of each DFM polygon that is coupled
    list of (x,y) coordinates of each PCR polygon that is coupled
    list of list of all DFM cells and its coupled PCR cells (if any)
    list of list of all unique coupled PCR cells including an array of all their coupled DFM cells
    list of list of all unique coupled PCR cells only
	
    """
    # calculate DFM cell centroid coordinates
    hydroCoordsCentroids = getHydroCoordsCentroids(hydroCoords)
    # calculate PCR cell min/max values
    PCRcoordsMinMax = getPCRcoordsMinMax(PCRcoords)
    # splitting this data into seperate arrays (so it is easier to see what is done in the code)
    centroids_x_DFM            = hydroCoordsCentroids[0]
    centroids_y_DFM            = hydroCoordsCentroids[1]
    x_min_all_cell_coords_pcr = PCRcoordsMinMax[0]
    x_max_all_cell_coords_pcr = PCRcoordsMinMax[1]
    y_min_all_cell_coords_pcr = PCRcoordsMinMax[2]
    y_max_all_cell_coords_pcr = PCRcoordsMinMax[3]
    
    # List with all DFM cells and coupled PCR cells (if any):
    # initial list with DFM cell numbers and zero values for their coupled PCR cells
    coupling_DFM_pcr = []
    for i in range(len(hydroCoords)):
        coupling_DFM_pcr.append([(i),(0)])
    # empty array for assigning only PCR cells (or 'no couple value')
    coupling_pcr_only = np.zeros(len(hydroCoords))
    # fill in zero array of coupled PCR cells:
    # looping over all DFM cells
    for i in range(len(hydroCoords)):
        # finding the to-be-coupled PCR cell
        temp = np.where((centroids_x_DFM[i] >= x_min_all_cell_coords_pcr) & \
                        (centroids_x_DFM[i] <  x_max_all_cell_coords_pcr) & \
                        (centroids_y_DFM[i] >= y_min_all_cell_coords_pcr) & \
                        (centroids_y_DFM[i] <  y_max_all_cell_coords_pcr))[0]
        # check if a cell was found
        if temp.size != 0:
            # if so, assign this value to the coupling list and PCR array
            coupling_DFM_pcr[i][1] = temp[0]
            coupling_pcr_only[i]  = temp[0]
        else:
            # if not, assign the chosen 'no couple value' to the coupling list and PCR array
            coupling_DFM_pcr[i][1] = no_couple_value
            coupling_pcr_only[i]  = no_couple_value
    
    # List with unique coupled PCR cells and another list that also includes their coupled DFM cells
    # Create lists with coupled cells only:
    # empty lists
    DFM_coupled_cells = []
    PCR_coupled_cells = []
    # loop over all indices of coupled cell list
    for i in range(len(coupling_DFM_pcr)):
        # check if a PCR cell was found for coupling
        if coupling_DFM_pcr[i][1] != no_couple_value:
            # if so, assign the current cell to the DFM list
            DFM_coupled_cells.append((coupling_DFM_pcr[i][0]))
            # and the found PCR cell to the PCR list
            PCR_coupled_cells.append((coupling_DFM_pcr[i][1]))
    # get only the unique coupled PCR cells 
    # (there will be many duplicate cells, since the code tried to find a PCR cell for every single DFM cell)
    PCR_coupled_cells_unique = list(set(PCR_coupled_cells))
    
    # Also include DFM cells coupled to unique coupled PCR cells:
    # initial list with zero values for the DFM cells
    coupling_pcr_DFM = []
    for i in range(len(PCR_coupled_cells_unique)):
        coupling_pcr_DFM.append([(PCR_coupled_cells_unique[i]),(0)])
    # loop over all unique coupled PCR cells
    for i in range(len(coupling_pcr_DFM)):
        # find all DFM cells that are coupled to the current PCR cell, and assign these to the list
        coupling_pcr_DFM[i][1] = np.where(PCR_coupled_cells_unique[i]==coupling_pcr_only)[0]
            
    # Create lists with coordinates of coupled cells only:
    # empty lists
    all_cell_coords_DFM_coupled = []
    all_cell_coords_PCR_coupled = []
    # loop over all coupled DFM cells and get their coordinates
    for i in range(len(DFM_coupled_cells)):
        all_cell_coords_DFM_coupled.append(hydroCoords[int(DFM_coupled_cells[i])])
    # loop over all coupled PCR cells and get their coordinates    
    for i in range(len(PCR_coupled_cells_unique)):
        all_cell_coords_PCR_coupled.append(PCRcoords[int(PCR_coupled_cells_unique[i])])
        
    return all_cell_coords_DFM_coupled, all_cell_coords_PCR_coupled, \
           coupling_DFM_pcr, coupling_pcr_DFM, PCR_coupled_cells_unique
	
# =============================================================================

def getHydroCoordsCentroids(hydroCoords):
    """
    Computing the centroids of each hydrodynamic cell
    
    Input:
	-----
	list of (x,y) coordinates of each hydrodynamic cell
    
    Output: 
	------
	separate arrays of x and y centroid coordinates
	
    """
    # total number of cells
    cell_count_DFM = len(hydroCoords)
    # zero arrays
    centroids_x_DFM = np.zeros(cell_count_DFM)
    centroids_y_DFM = np.zeros(cell_count_DFM)
    # calculate cell areas
    cell_areas = calculateAllCellArea(hydroCoords)
    
    # loop over all cells
    for i in range(cell_count_DFM):
        # current cell size
        cell_size = len(hydroCoords[i])
        # current cell area
        cell_area = cell_areas[i]
        # Calculate x and y coordinate of cell centroid:
        # reset temporary values
        temp_value_for_x = 0
        temp_value_for_y = 0
        # loop over all vertices of current cell
        for j in range(cell_size):
            # checking if not at last vertex
            if j < (cell_size-1):
                # if not, calculation based on current and next vertex
                temp_value_for_x = temp_value_for_x + ((hydroCoords[i][j][0] + hydroCoords[i][j+1][0]) * \
                                                       ((hydroCoords[i][j][0] * hydroCoords[i][j+1][1]) - \
                                                        (hydroCoords[i][j+1][0] * hydroCoords[i][j][1])))
                
                temp_value_for_y = temp_value_for_y + ((hydroCoords[i][j][1] + hydroCoords[i][j+1][1]) * \
                                                       ((hydroCoords[i][j][0] * hydroCoords[i][j+1][1]) - \
                                                        (hydroCoords[i][j+1][0] * hydroCoords[i][j][1])))
            else:
                # otherwise, calculation based on current and first vertex (since this is the 'next' vertex in this case)
                temp_value_for_x = temp_value_for_x + ((hydroCoords[i][j][0] + hydroCoords[i][0][0]) * \
                                                       ((hydroCoords[i][j][0] * hydroCoords[i][0][1]) - \
                                                        (hydroCoords[i][0][0] * hydroCoords[i][j][1])))
                
                temp_value_for_y = temp_value_for_y + ((hydroCoords[i][j][1] + hydroCoords[i][0][1]) * \
                                                       ((hydroCoords[i][j][0] * hydroCoords[i][0][1]) - \
                                                        (hydroCoords[i][0][0] * hydroCoords[i][j][1])))
            
        # use cell area and calculated temporary values to find the x and y coordinate of the cell centroids
        centroids_x_DFM[i] = (1/(6*cell_area))*temp_value_for_x
        centroids_y_DFM[i] = (1/(6*cell_area))*temp_value_for_y
        
    return centroids_x_DFM, centroids_y_DFM

# =============================================================================

def calculateAllCellArea(cellCoords):
    """
    Calculates the area of all cells, based on their coordinates
    
    Input:  
	-----
	list of (x,y) coordinates of each hydrodynamic cell
    
    Output: 
	------
	area of all cells (array)
   
    """
    # calculate number of cells
    cell_count = len(cellCoords)
    # zero array for filling in cell areas
    cell_areas = np.zeros(cell_count)
    # loop over all cells
    for i in range(cell_count):
        # current cell size
        cell_size = len(cellCoords[i])
        # reset temporary value
        temp_value_for_A = 0
        # Calculate area of current cell:
        # loop over all vertices
        for j in range(cell_size):
            # check if not at last vertex
            if j < (cell_size-1):
                # if not, calculation based on current and next vertex
                temp_value_for_A = temp_value_for_A + ((cellCoords[i][j][0] * cellCoords[i][j+1][1]) -\
                                                       (cellCoords[i][j+1][0] * cellCoords[i][j][1]))
            else:
                # otherwise, based on current and first vertex (since this is 'next' vertex)
                temp_value_for_A = temp_value_for_A + ((cellCoords[i][j][0] * cellCoords[i][0][1]) -\
                                                       (cellCoords[i][0][0] * cellCoords[i][j][1]))
        cell_area = 0.5 * temp_value_for_A
        # assign area of current cell to array
        cell_areas[i] = cell_area
    
    return cell_areas

# =============================================================================

def getPCRcoordsMinMax(PCRcoords):
    """
    Get minimum and maximum of x and y coordinates for each PCR polygon
    
    Input:  
	-----
	list of (x,y) coordinates of each polygon (output of getPCRcoords)
    
    Output: 
	------
	separate arrays with xMin, xMax, ymin and yMax coordinates
    
    NOTE: the used indices (i.e. 0,1,2) depend on the way the coordinates are ordered in the function getPCRcoords (currently ll, lr, ur, ul)
    """
	
    cell_count_pcr = len(PCRcoords)
    
    x_min_all_cell_coords_pcr = np.zeros(cell_count_pcr)
    x_max_all_cell_coords_pcr = np.zeros(cell_count_pcr)
    y_min_all_cell_coords_pcr = np.zeros(cell_count_pcr)
    y_max_all_cell_coords_pcr = np.zeros(cell_count_pcr)
    
    for i in range(cell_count_pcr):
        x_min_all_cell_coords_pcr[i] = min(PCRcoords[i][0][0],PCRcoords[i][1][0])
        x_max_all_cell_coords_pcr[i] = max(PCRcoords[i][0][0],PCRcoords[i][1][0])
        y_min_all_cell_coords_pcr[i] = min(PCRcoords[i][1][1],PCRcoords[i][2][1])
        y_max_all_cell_coords_pcr[i] = max(PCRcoords[i][1][1],PCRcoords[i][2][1])
    
    return x_min_all_cell_coords_pcr, x_max_all_cell_coords_pcr, \
           y_min_all_cell_coords_pcr, y_max_all_cell_coords_pcr

# =============================================================================

def assignPCR2cells(landmask_pcr, hydroCoords, verbose):
    """
	Function converting single indices of coupled PCR cells to double (array,column) indices.
	This is the key function coupling the grids of PCR-GLOBWB and the hydrodynamic models.
	
	Input:  
	-----
    :param landmask_pcr: land mask of hydrological model
    :param hydroCoords: list of (x,y) coordinates of each hydrodynamic cell
    :param verbose: print information yes/no

	Output:
	------
	coupledHydro2PCR: double containing coupled PCR cell per cell of hydrodynamic model
	coupledPCR2hydro: double containing all hydrodynamic cells per coupled PCR cell
	zipPCRIndices: indices pointing to all coupled PCR cells in an array
    """

    #-read in x and y coordinates set clone, get its attributes and read in landmask
    pcr.setclone(landmask_pcr)
    landMask= pcr.readmap(landmask_pcr)
    cloneAttributes= spatialDataSet2PCR.spatialAttributes(landmask_pcr)
    #-call function to get cell centres
    lon_hydro, lat_hydro= getMidPointsFromVertices(hydroCoords)
	#-print extent of coupled cells of hydrodynamic models if specified
    if verbose == True:
	    print ' - coordinates of cell centres contain %d entries with a latitudes spanning %f-%f and longitudes %f-%f' %\
		    (len(lat_hydro), min(lat_hydro), max(lat_hydro), min(lon_hydro), max(lon_hydro))
    #-read in lat and lon for DFM
    index_hydro = np.arange(len(lon_hydro))
    #-get corresponding PCRaster indices
    cellIndexPCR, rowIndexPCR, colIndexPCR= \
        getPCRIndices(index_hydro, lat_hydro, lon_hydro, landMask, cloneAttributes, verbose= False)
    #-check on assignment
    numberDFMCells= np.zeros((cloneAttributes.numberRows, cloneAttributes.numberCols))
    pcrCellID= np.zeros((cloneAttributes.numberRows, cloneAttributes.numberCols))
    for iCnt in xrange(len(cellIndexPCR)):
        rowCnt= rowIndexPCR[iCnt]
        colCnt= colIndexPCR[iCnt]
        numberDFMCells[rowCnt, colCnt]+= 1
        pcrCellID[rowCnt, colCnt]= cellIndexPCR[iCnt]
    #-return information for coupling
    coupledHydro2PCR, coupledPCR2hydro, zipPCRIndices= \
		returnCoupledPCRIndices(index_hydro, lat_hydro, lon_hydro, landMask, cloneAttributes, verbose= False)
		
    #-check on assignment
    # check for PCRaster first
    numberDFMCells= np.zeros((cloneAttributes.numberRows, cloneAttributes.numberCols))
    pcrCellID= np.zeros((cloneAttributes.numberRows, cloneAttributes.numberCols))
    sortedRows, sortedCols= zip(*zipPCRIndices)
    for iCnt in xrange(len(coupledPCR2hydro)):
        cellID= coupledPCR2hydro[iCnt][0]
        rowCnt= sortedRows[iCnt]
        colCnt= sortedCols[iCnt]
        numberEntries= np.size(coupledPCR2hydro[iCnt][1])
   
    return coupledHydro2PCR, coupledPCR2hydro, zipPCRIndices

# =============================================================================

def getMidPointsFromVertices(coordVertices):
	
	"""
	Reads in a compound list of coordinates of the vertices of each cell on a cartesian grid
	in x, y coordinates and returns two arrays with the x, y coordinate of the mid point
	(geometric centre) of the corresponding cell
	
	Input:
	------
	coordVertices:	the coordinates of the vertices consisting of N(x, y pairs) with
									the entries of each cell given per row; it is assumed coordinate pairs
									are entered as tupples and grouped into a list of N entries
									
	Output:
	-------
	xCoord:					list of N entries containing the x coordinate of the centre point of each cell
	yCoord:					list of N entries containing the y coordinate of the centre point of each cell
	"""
	
	#-create empty lists with floats of length specified by the list of vertices
	xCoord= [0.0 for iCnt in xrange(len(coordVertices))]
	yCoord= [0.0 for iCnt in xrange(len(coordVertices))]
	#-get vertices and pop entries as pairs
	for iCnt in xrange(len(coordVertices)):
		vertices= coordVertices[iCnt]
		x= 0.0
		y= 0.0
		vCnt= 0
		for x_v, y_v in vertices:
			x+= x_v
			y+= y_v
			vCnt+= 1
		#-compute average and insert in lists
		x/= vCnt
		y/= vCnt
		xCoord[iCnt]= x
		yCoord[iCnt]= y
	#-return lists
	return xCoord, yCoord
 
# ============================================================================= 

def getPCRIndices(index_hydro, lat_hydro, lon_hydro, pcrLandMask, mapAttributes, verbose= False):
	"""

	Returns for each entry of the flexible mesh, specified as ID, latitude and longitude of its cell centre,
	the corresponding cell number and its location as row and column number for a PCRaster map that identifies
	the landmask (1, True; 0: False) on the basis	of the nearest distance between cell centres.
	
	Input:
	------
	index_hydro:				sequence of indices to identify DFM cells
	lat_hydro:					latitude  of cell centre (dec. degrees)
	lon_hydro:					longitude of cell centre (dec. degrees)
	pcrLandMask:		PCRaster field of land mask (boolean, 1 denoting land mask)
	mapAttributes:	instance of map attributes (may be replaced by direct access to PCRaster map attributes)
	verbose:		test on output, if True, assignment is printed to screen (optional, default is False)

	Output:
	-------
	cellIndexPCR:		sequence (array) of same size as index_hydro identifying the assigned cells
									on the PCRaster landmask; cells are numbered from upper-left to bottom-right
									and start at 0
	rowIndexPCR:		sequence (array) of same size as index_hydro listing the row corresponding to the assigned cell
	colIndexPCR:		sequence (array) of same size as index_hydro listing the column corresponding to the assigned cell
	
	Row and column indices are numbered according to python convention, i.e., starting at zero.

	"""
	
	#-get np indices referring to the rows and columns of the PCRaster file
	rowsPCR= np.arange(mapAttributes.numberRows)
	colsPCR= np.arange(mapAttributes.numberCols)
	#-get unique latitudes and longitudes
	latPCR= mapAttributes.yUR-\
		(rowsPCR+0.5)*np.abs(mapAttributes.yResolution)
	lonPCR= mapAttributes.xLL+(colsPCR+0.5)*mapAttributes.xResolution
	#-get landmask as numpy array, MVs set to zero
	npLandMask= pcr.pcr2numpy(pcrLandMask,0)
	#-create raveled arrays of selected coordinates of land mask
	selectedLatPCR= np.zeros(npLandMask.shape)
	selectedLonPCR= np.zeros(npLandMask.shape)
	for colCnt in colsPCR:
		selectedLatPCR[:, colCnt]= latPCR.copy()
	for rowCnt in rowsPCR:
		selectedLonPCR[rowCnt, :]= lonPCR.copy()
	selectedLatPCR= selectedLatPCR[npLandMask == 1]
	selectedLonPCR= selectedLonPCR[npLandMask == 1]
	#-create arrays to hold information derived from landmask
	minimumDistance= np.zeros((index_hydro.size))-1.0
	rowIndexPCR= np.zeros((index_hydro.size))-1
	colIndexPCR= np.zeros((index_hydro.size), dtype= int)-1
	cellIndexPCR= np.zeros((index_hydro.size), dtype= int)-1
	#-iterate over all DFM indices and find minimum distance and corresponding
	# entries in the PCRaster map
	for iCnt in xrange(index_hydro.size):
		#-for current lat and lon, find the arc distance and select the row, column entry corresponding
		# to the minimum distance; if two or more cells are the closest, the first is chosen by default
		arcDistance= pcrGlobalGeometry.getArcDistance(lat_hydro[iCnt], lon_hydro[iCnt], selectedLatPCR, selectedLonPCR)
		minimumDistance[iCnt]= np.min(arcDistance)
		selectionMask= arcDistance == minimumDistance[iCnt]
		#-check on double entries
		if np.sum(selectionMask) > 1:
			validIndex= np.arange(selectionMask.size)[selectionMask][0]
			selectionMask[:]= False
			selectionMask[validIndex]= True
		#-get masks for latitude and longitude
		latMask= latPCR == selectedLatPCR[selectionMask]
		lonMask= lonPCR == selectedLonPCR[selectionMask]
		#-get corresponding row, column and compute cell index
		rowIndexPCR[iCnt]= rowsPCR[latMask]
		colIndexPCR[iCnt]= colsPCR[lonMask]
		cellIndexPCR[iCnt]= colIndexPCR[iCnt]+\
			rowIndexPCR[iCnt]*colsPCR.size
		#-echo to check assignment if required
		if verbose:
			print ' - DFM entry %5d at lat, lon %7.1f, %7.1f is matched to PCRaster cell %5d at %7.1f, %7.1f at %5.1f km located at row %5d, col %5d' %\
				(iCnt, lat_hydro[iCnt], lon_hydro[iCnt], cellIndexPCR[iCnt],\
				selectedLatPCR[selectionMask], selectedLonPCR[selectionMask], 1.e-3*minimumDistance[iCnt],\
					rowIndexPCR[iCnt], colIndexPCR[iCnt])
	#-all entries processed, return cell index and corresponding row and column numbers
	return cellIndexPCR, rowIndexPCR, colIndexPCR
 
# ============================================================================= 
   
def returnCoupledPCRIndices(index_hydro, lat_hydro, lon_hydro, pcrLandMask, mapAttributes, verbose):
	"""

	Returns for each entry of the hydrodynamic model, specified as ID, latitude and longitude of its cell centre,
	the information that is required for its coupling with cells in the PCR fields that match the spatial
	attributes of the specified land mask. The function calls getPCRIndices() first, and then post-processes 
	the results to obtain the required output, containing a nested list of hydrodynamic model IDs coupled one-to-one
	to PCR IDs, a compound list of sorted PCR IDs linked to an array of hydrodynamic model IDs (one-to-many) and
	a zipped list of row and column indices of all sorted sorted PCR IDs.
	
	Input:
	------
	index_hydro:				sequence of indices to identify hydrodynamic cells
	lat_hydro:					latitude  of cell centre (dec. degrees)
	lon_hydro:					longitude of cell centre (dec. degrees)
	pcrLandMask:		        PCRaster field of land mask (boolean, 1 denoting land mask)
	mapAttributes:	            instance of map attributes (may be replaced by direct access to PCRaster map attributes)
	verbose:		            test on output, if True, assignment is printed to screen (optional, default is False)

	Output:
	-------
	coupledHydro2PCR:	nested list of hydrodynamic model IDs and PCR IDs that are paired on a one-to-one basis;
									the list has as many entries as hydrodynamic model IDs where each entry consists of a list with
									two entries, detailing the hydrodynamic model ID and the corresponding PCRaster ID. Flexible mesh
									IDs are sorted in ascending order
	coupledPCR2hydro:	compound list of PCRaster IDs that are paired on a one-to-many basis: the list contains as many
									entries as there are unique PCRaster IDs that are sorted in ascending order. Every entry
									contains a list with two entries, the PCRaster ID first, followed by a 1D numpy array of var-
									iable size holding all hydrodynamic model IDs that are matched by the PCRaster ID under consideration 
	zipPCRIndices:	zipped list holding row and column indices per unique PCRaster ID -sorted in ascending order
	
	Row and column indices are numbered according to python convention, i.e., starting at zero.

	"""
	
	#-call function and get the PCRaster information		
	cellIndexPCR, rowIndexPCR, colIndexPCR= getPCRIndices(index_hydro, lat_hydro, lon_hydro, pcrLandMask, mapAttributes, verbose)
	
	#-get the unique IDs from the PCRaster IDs and sort both these and the index_hydro
	uniqueIndexPCR= np.unique(cellIndexPCR)
	uniqueIndexPCR.sort()
	orderindex_hydro= np.argsort(index_hydro)
	#-sort DFM indices and corresponding PCRaster IDs
	index_hydro= index_hydro[orderindex_hydro]
	sortedIndices= cellIndexPCR[orderindex_hydro]
	coupledHydro2PCR= [[-1, -1] for iCnt in xrange(index_hydro.size)]
	for iCnt in xrange(index_hydro.size):
		coupledHydro2PCR[iCnt]= [index_hydro[iCnt], sortedIndices[iCnt]]
	#-write data for the unique and sorted PCRaster indices
	# as well as the corresponding rows and columns 
	coupledPCR2hydro= [[-1, -1] for iCnt in xrange(uniqueIndexPCR.size)]
	sortedRows= [0 for iCnt in xrange(uniqueIndexPCR.size)]
	sortedCols= [0 for iCnt in xrange(uniqueIndexPCR.size)]
	# info on min and max length of array with DFM IDs
	minSize= index_hydro.size
	minSizeIndex= 0
	maxSize= 0
	maxSizeIndex= 0
	for iCnt in xrange(uniqueIndexPCR.size):
		sortedRows[iCnt]= rowIndexPCR[cellIndexPCR == uniqueIndexPCR[iCnt]][0]
		sortedCols[iCnt]= colIndexPCR[cellIndexPCR == uniqueIndexPCR[iCnt]][0]
		coupledPCR2hydro[iCnt]= [uniqueIndexPCR[iCnt], index_hydro[cellIndexPCR == uniqueIndexPCR[iCnt]]]
		if np.size(coupledPCR2hydro[iCnt][1]) > maxSize:
			maxSize= np.size(coupledPCR2hydro[iCnt][1])
			maxSizeIndex= iCnt
		if np.size(coupledPCR2hydro[iCnt][1]) < minSize:
			minSize= np.size(coupledPCR2hydro[iCnt][1])
			minSizeIndex= iCnt
	#-zip row and column arrays
	zipPCRIndices= zip(sortedRows, sortedCols)
	#-test
	if verbose == True:
		print '\n - %d unique IDs out of %d PCRaster IDs are found starting at %d and ending at %d' %\
			(uniqueIndexPCR.size, cellIndexPCR.size, uniqueIndexPCR[0], uniqueIndexPCR[-1])
		print '   minimum array of matched flexible mesh IDs of length %d found for %d' %\
			(minSize, uniqueIndexPCR[minSizeIndex])
		print '   ', coupledPCR2hydro[minSizeIndex][1]
		print '   maximum array of matched flexible mesh IDs of length %d found for %d' %\
			(maxSize, uniqueIndexPCR[maxSizeIndex])
		print '   ', coupledPCR2hydro[maxSizeIndex][1]
		print ' - %d IDs are provided for the flexible mesh' % index_hydro.size
		if np.all(orderindex_hydro == np.arange(index_hydro.size)):
			print '   no sorting is required for the flexible mesh IDs'
		print '   coupled flexible mesh to PCRaster IDs:', coupledHydro2PCR[0], coupledHydro2PCR[-1]

	#-return values
	return coupledHydro2PCR, coupledPCR2hydro, zipPCRIndices

# =============================================================================

def plotGridfromCoords(coordsGrid1, coordsGrid2=None, linewidthGrid1=1, linewidthGrid2=1):
    """
    Creates a 2D plot from polygon coordinates (output of coupleAllCells)
    
    Input:
	-----
	list of (x,y) coordinates of each PCR polygon that is coupled
    list of (x,y) coordinates of each hydrodynamic polygon that is coupled (optional)
            
    Output: 
	------
	plot of model grids
    
    NOTE: to get the plots for 2 grids right, it is important that:
            - 1st grid = PCR grid
            - 2nd grid = hydrodynamic grid
			- not yet possible to visualize 1-D channels in DFM as they are vectors not grids
    """
    # calling empty plot
    f, ax = plt.subplots(1, 1)
    
    # create polygons
    polygonsGrid1 = matplotlib.collections.PolyCollection(coordsGrid1)
    
    # set polygon appearance
    if coordsGrid2 == None:
        polygonsGrid1.set_aLFPha(0.2)
        polygonsGrid1.set_facecolor('black')
        polygonsGrid1.set_edgecolor('black')
        polygonsGrid1.set_linewidth(linewidthGrid1)
    else:
        polygonsGrid1.set_aLFPha(1)
        polygonsGrid1.set_facecolor('white')
        polygonsGrid1.set_edgecolor('black')
        polygonsGrid1.set_linewidth(linewidthGrid1)
    
    # add PCR polygons to plot
    ax.add_collection(polygonsGrid1)
    
    # checking for second grid
    if coordsGrid2 != None:
        polygonsGrid2  = matplotlib.collections.PolyCollection(coordsGrid2)
    
        polygonsGrid2.set_aLFPha(0.2)
        polygonsGrid2.set_facecolor('black')
        polygonsGrid2.set_edgecolor('black')
        polygonsGrid2.set_linewidth(linewidthGrid2)
    
        ax.add_collection(polygonsGrid2)
        
    # adjust scale
    ax.autoscale()
	


