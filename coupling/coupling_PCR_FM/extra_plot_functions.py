# -*- coding: utf-8 -*-
"""
Created on Mon Dec 15 11:26:01 2014

@author: haag

Some extra functions for plotting geographical data with my model grids.

"""

from coupling_PCR_FM import coupling_functions

import shapefile

import numpy as np

import matplotlib.patches as patches
from matplotlib.lines import Line2D
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection

def plotCities(ax, path, plot_labels=None, plot_markersize=50, plot_markercolor='black', plot_markertype='o', plot_fontsize=14, \
               lon_offset=0.1, lat_offset=0.0, plot_horizontalalignment='right', plot_verticalalignment='bottom', \
               box_facecolor='white', box_edgecolor='none', box_alpha=1.0, x1=None, x2=None, y1=None, y2=None):
    """
    Creates a plot of cities (optionally with their names) from shapefile input data.
    
    TO-DO: allow for manual changing of location of specified country names
    """
    
    # read shapefile
    sf     = shapefile.Reader(path)
    
    # get shapefile data
    recs   = sf.records()
    shapes = sf.shapes()
    
    # plot cities using specified parameters
    for i in range(len(recs)):
        lon_loc_station  = recs[i][0]
        lat_loc_station  = recs[i][1]
        name_loc_station = recs[i][2]
        # only plot if lat/lon falls within specified extent
        if (x1 != None) and (x2 != None) and (y1 != None) and (y2 != None):
            if (lon_loc_station >= x1) and (lon_loc_station <= x2) and (lat_loc_station >= y1) and (lat_loc_station <= y2):
                ax.scatter(lon_loc_station, lat_loc_station, s=plot_markersize, c=plot_markercolor, marker=plot_markertype)
                if plot_labels != None:
                    ax.text(lon_loc_station + lon_offset, lat_loc_station + lat_offset, name_loc_station, fontsize=plot_fontsize, \
                            horizontalalignment=plot_horizontalalignment, verticalalignment=plot_verticalalignment,\
                            bbox={'facecolor':box_facecolor, 'edgecolor':box_edgecolor, 'alpha':box_alpha})
        else:
            ax.scatter(lon_loc_station, lat_loc_station, s=plot_markersize, c=plot_markercolor, marker=plot_markertype)
            if plot_labels != None:
                ax.text(lon_loc_station + lon_offset, lat_loc_station + lat_offset, name_loc_station, fontsize=plot_fontsize, \
                        horizontalalignment=plot_horizontalalignment, verticalalignment=plot_verticalalignment,\
                        bbox={'facecolor':box_facecolor, 'edgecolor':box_edgecolor, 'alpha':box_alpha})

def plotCountries(ax, path, return_data='no', plot_labels=None, label_loc_edit=None, plot_facecolor='none', plot_edgecolor='k', plot_linewidths=0.2, \
                  plot_fontsize=14, plot_horizontalalignment='center', plot_verticalalignment='bottom', \
                  x1=None, x2=None, y1=None, y2=None):
    """
    Creates a plot of countries (optionally with their names) from shapefile input data.
    
    It is also be possible to return the data used to create this plot instead of plotting it, to manually edit this further if desired.
    
    TO-DO: allow for manual changing of location of specified country names
    """
    
    # read shapefile
    sf     = shapefile.Reader(path)
    
    # get shapefile data
    recs   = sf.records()
    shapes = sf.shapes()
    
    # create list of country polygons
    country_polygons = []
    Nshp             = len(shapes)
    for nshp in xrange(Nshp):
        ptchs   = []
        pts     = np.array(shapes[nshp].points)
        prt     = shapes[nshp].parts
        par     = list(prt) + [pts.shape[0]]
        for pij in xrange(len(prt)):
            ptchs.append(Polygon(pts[par[pij]:par[pij+1]]))
        country_polygons.append(ptchs)
        
    # create list of country names
    country_names = []
    for i in range(len(recs)):
        country_names.append(recs[i][3])
        
    # create list of country polygon extents
    country_polygons_extent = []
    for i in range(len(country_polygons)):
        temp_list = []
        ll_temp   = tuple(np.array(country_polygons[i][0].get_extents())[0])
        ur_temp   = tuple(np.array(country_polygons[i][0].get_extents())[1])
        temp_list.append(ll_temp)
        temp_list.append((ur_temp[0],ll_temp[1]))
        temp_list.append(ur_temp)
        temp_list.append((ll_temp[0],ur_temp[1]))
        country_polygons_extent.append(temp_list)
        
    # calculate country polygon centroids (used to display country names)
    country_polygons_centroid = coupling_functions.getFMcoordsCentroids(country_polygons_extent)
    
    # edit these so they better fit on the map (when displayed at full extent)
    country_polygons_centroid_edit = coupling_functions.getFMcoordsCentroids(country_polygons_extent)
    for i in range(len(country_polygons_centroid_edit[0])):
        if country_names[i] == 'Central African Republic':
            country_polygons_centroid_edit[0][i] += 0.2
            country_polygons_centroid_edit[1][i] -= 0.8
        elif country_names[i] == 'Chad':
            country_polygons_centroid_edit[0][i] += 1.0
            country_polygons_centroid_edit[1][i] -= 0.5
        elif country_names[i] == 'Ivory Coast':
            country_polygons_centroid_edit[0][i] -= 2.5
            country_polygons_centroid_edit[1][i] += 2.2
        elif country_names[i] == 'Cameroon':
            country_polygons_centroid_edit[1][i] -= 3.5
        elif country_names[i] == 'Guinea':
            country_polygons_centroid_edit[1][i] += 0.8
        elif country_names[i] == 'Nigeria':
            country_polygons_centroid_edit[1][i] += 5.0
        elif country_names[i] == 'Senegal':
            country_polygons_centroid_edit[0][i] -= 0.6
            country_polygons_centroid_edit[1][i] -= 0.2
        elif country_names[i] == 'Togo':
            country_polygons_centroid_edit[0][i] += 0.2
    
    if label_loc_edit == None:
        pass
    elif label_loc_edit == 'pcr_full':
        for i in range(len(country_polygons_centroid_edit[0])):
            if country_names[i] == 'Burkina Faso':
                country_polygons_centroid_edit[0][i] -= 1.0
                country_polygons_centroid_edit[1][i] -= 0.4
            elif country_names[i] == 'Guinea':
                country_polygons_centroid_edit[0][i] -= 1.0
                country_polygons_centroid_edit[1][i] += 0.2
            elif country_names[i] == 'Mali':
                country_polygons_centroid_edit[0][i] += 1.5
                country_polygons_centroid_edit[1][i] += 2.5
            elif country_names[i] == 'Niger':
                country_polygons_centroid_edit[0][i] += 2.0
                country_polygons_centroid_edit[1][i] += 1.0
            elif country_names[i] == 'Nigeria':
                country_polygons_centroid_edit[0][i] += 2.0
                country_polygons_centroid_edit[1][i] += 2.0
    elif len(label_loc_edit) >= 2:
        for i in range(len(label_loc_edit)/2):
            temp_index = country_names.index(label_loc_edit[(i*2)])
            country_polygons_centroid_edit[0][temp_index] = label_loc_edit[((i*2)+1)][0]
            country_polygons_centroid_edit[1][temp_index] = label_loc_edit[((i*2)+1)][1]
    
    # check whether the function is called to return data or to plot values directly
    if return_data == 'no':
            
        # plot the countries using specified parameters
        for i in range(len(country_polygons)):
            ax.add_collection(PatchCollection(country_polygons[i], facecolor=plot_facecolor, edgecolor=plot_edgecolor, linewidths=plot_linewidths))
            if plot_labels != None:
                if i not in [5,9,10,12,13,18]:
                    #if country_names[i] == 'Niger':
                    #   ax.annotate(country_names[i], xy=(1.8, 14.7), fontsize=14, horizontalalignment='center', verticalalignment='bottom')
                    #elif country_names[i] == 'Mali':
                    #    ax.annotate(country_names[i], xy=(-6.5, 14.5), fontsize=14, horizontalalignment='center', verticalalignment='bottom')
                    #elif country_names[i] == 'Mauritania':
                    #    ax.annotate(country_names[i], xy=(-7.5, 16.5), fontsize=14, horizontalalignment='center', verticalalignment='bottom')
                    #else:
                        # only plot if location of text falls within specified extent
                    #    if (country_polygons_centroid_edit_2[0][i] >= x1) and (country_polygons_centroid_edit_2[0][i] <= x2) and \
                    #       (country_polygons_centroid_edit_2[1][i] >= y1) and (country_polygons_centroid_edit_2[1][i] <= y2):
                    #        ax.annotate(country_names[i], xy=(country_polygons_centroid_edit_2[0][i], country_polygons_centroid_edit_2[1][i]), \
                    #                    fontsize=14, horizontalalignment='center', verticalalignment='bottom')
                    # only plot if location of text falls within specified extent
                    if (x1 != None) and (x2 != None) and (y1 != None) and (y2 != None):
                        if (country_polygons_centroid_edit[0][i] >= x1) and (country_polygons_centroid_edit[0][i] <= x2) and \
                           (country_polygons_centroid_edit[1][i] >= y1) and (country_polygons_centroid_edit[1][i] <= y2):
                            ax.annotate(country_names[i], xy=(country_polygons_centroid_edit[0][i], country_polygons_centroid_edit[1][i]), \
                                        fontsize=plot_fontsize, horizontalalignment=plot_horizontalalignment, verticalalignment=plot_verticalalignment)
                    else:
                        ax.annotate(country_names[i], xy=(country_polygons_centroid_edit[0][i], country_polygons_centroid_edit[1][i]), \
                                    fontsize=plot_fontsize, horizontalalignment=plot_horizontalalignment, verticalalignment=plot_verticalalignment)
                        
    else:
        
        return country_polygons, country_names, country_polygons_centroid, country_polygons_centroid_edit

def plotRivers(ax, path, return_data='no', plot_color='blue', plot_linewidth=1.0):
    """
    Creates a plot of rivers from shapefile input data.
    
    It is also be possible to return the data used to create this plot instead of plotting it, to manually edit this further if desired.
    """
    
    # read shapefile
    sf     = shapefile.Reader(path)
    
    # get shapefile data
    recs   = sf.records()
    shapes = sf.shapes()
    
    # create a list of (x,y) points of rivers
    river_points = []
    Nshp = len(shapes)
    for nshp in xrange(Nshp):
        pts     = np.array(shapes[nshp].points)
        prt     = shapes[nshp].parts
        par     = list(prt) + [pts.shape[0]]
        for pij in xrange(len(prt)):
            x_pts = zip(*pts[par[pij]:par[pij+1]])[0]
            y_pts = zip(*pts[par[pij]:par[pij+1]])[1]
            river_points.append((x_pts,y_pts))

    # check whether the function is called to return data or to plot values directly
    if return_data == 'no':
        
        # plot rivers using specified parameters
        for i in range(len(river_points)):
            ax.add_line(Line2D(river_points[i][0], river_points[i][1], color=plot_color, linewidth=plot_linewidth))
    
    else:
        
        return river_points

