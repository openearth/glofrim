import os, sys
import numpy as np
import pcraster as pcr
from math import pi as m_pi
from types import NoneType

pcr.setglobaloption('radians')
deg2Rad= m_pi/180.

def getArcDistance(latA, lonA, latB, lonB, radius= 6371221.3, testVerbose= False):
  '''Computes the distance between two points, positioned by \
their geographic coordinates along the surface of a perfect sphere \
in units given by the radius used. Input variables include:\n
- latA, latB: the latitude of the points considered in decimal degrees,\n
- lonA, lonB: the longitude of the points considered in decimal degrees,\n
- radius: the radius of the sphere in metres, set by default to \
that of Earth (6371221.3 m).'''
	#-make arrays if needed
  if isinstance(latA,float):
		latA= np.array(latA)
  if isinstance(lonA,float):
		lonA= np.array(lonA)
  if isinstance(latB,float):
		latB= np.array(latB)
  if isinstance(lonB,float):
		lonB= np.array(lonB)
  #-pad latitudes, longitudes
  if latA.size <> latB.size:
    latA= np.ones(latB.shape)*latA
  if lonA.size <> lonB.size:
    lonA= np.ones(lonB.shape)*lonA
  #-convert all coordinates to radians
  pA= latA*deg2Rad; pB= latB*deg2Rad
  lA= lonA*deg2Rad; lB= lonB*deg2Rad
  dL= lB-lA; dP= pB-pA
  #-set arcDist default to zero and create mask of cells to be processed
  a= np.sin(0.5*dP)*np.sin(0.5*dP)+\
    np.cos(pA)*np.cos(pB)*\
    np.sin(0.5*dL)*np.sin(0.5*dL)
  arcDist= 2*np.arctan2(a**0.5,(1-a)**0.5)
  arcDist*= radius
  if testVerbose:
    print ' * along an ideal sphere of radius %f, the distance between the points at lat, lon \
  %f, %f and %f, %f respectively amounts to %f' %\
      (radius, latA[arcDist == arcDist.max()][0].astype(float), lonA[arcDist == arcDist.max()][0].astype(float),\
       latB[arcDist == arcDist.max()][0].astype(float),\
       lonB[arcDist == arcDist.max()][0].astype(float), arcDist[arcDist == arcDist.max()][0].astype(float))
  #-return arcDist
  return arcDist

def getArcDistance2(latA, lonA, latB, lonB, radius= 6371221.3, testVerbose= False):
  '''Computes the distance between two points, positioned by \
their geographic coordinates along the surface of a perfect sphere \
in units given by the radius used. Input variables include:\n
- latA, latB: the latitude of the points considered in decimal degrees,\n
- lonA, lonB: the longitude of the points considered in decimal degrees,\n
- radius: the radius of the sphere in metres, set by default to \
that of Earth (6371221.3 m).
'''
	#-make arrays if needed
  if isinstance(latA,float):
		latA= np.array(latA)
  if isinstance(lonA,float):
		lonA= np.array(lonA)
  if isinstance(latB,float):
		latB= np.array(latB)
  if isinstance(lonB,float):
		lonB= np.array(lonB)
  #-pad latitudes, longitudes
  if latA.size <> latB.size:
    latA= np.ones(latB.shape)*latA
  if lonA.size <> lonB.size:
    lonA= np.ones(lonB.shape)*lonA
  #-convert all coordinates to radians
  pA= latA*deg2Rad; pB= latB*deg2Rad
  lA= lonA*deg2Rad; lB= lonB*deg2Rad
  Z= np.sin(pA)*np.sin(pB)+\
    np.cos(pA)*np.cos(pB)*\
    np.cos(lA-lB)
  #-set arcDist default to zero and create mask of cells to be processed
  if Z.shape == ():
    Z= np.array([Z])
  mask= np.zeros(Z.shape)
  mask= (Z < 1) & (Z > -1)
  arcDist= np.zeros(Z.shape)
  arcDist[Z <= -1]= m_pi
  arcDist[mask]= np.arccos(Z[mask])
  arcDist*= radius
  if testVerbose:
    print ' * along an ideal sphere of radius %f, the distance between the points at lat, lon \
  %f, %f and %f, %f respectively amounts to %f' %\
      (radius, latA[arcDist == arcDist.max()][0].astype(float), lonA[arcDist == arcDist.max()][0].astype(float),\
       latB[arcDist == arcDist.max()][0].astype(float),\
       lonB[arcDist == arcDist.max()][0].astype(float), arcDist[arcDist == arcDist.max()][0].astype(float))
  #-return arcDist
  return arcDist

def getArcDistancePCR(latA, lonA, latB, lonB, radius= 6371221.3):
  '''Computes the distance using PCRaster between two points, positioned by \
their geographic coordinates along the surface of a perfect sphere \
in units given by the radius used. Input variables include:\n
- latA, latB: the latitude of the points considered in decimal degrees,\n
- lonA, lonB: the longitude of the points considered in decimal degrees,\n
- radius: the radius of the sphere in metres, set by default to \
that of Earth (6371221.3 m).
'''
  #-pad latitudes, longitudes
  #-convert all coordinates to radians
  pA= latA*deg2Rad; pB= latB*deg2Rad
  lA= lonA*deg2Rad; lB= lonB*deg2Rad
  dL= lB-lA; dP= pB-pA
  #-set arcDist default to zero and create mask of cells to be processed
  a= pcr.sin(0.5*dP)*pcr.sin(0.5*dP)+\
    pcr.cos(pA)*pcr.cos(pB)*\
    pcr.sin(0.5*dL)*pcr.sin(0.5*dL)
  arcDist= 2*pcr.scalar(pcr.atan(a**0.5/(1-a)**0.5))
  arcDist*= radius
  #-return arcDist
  return arcDist

def getAzimuth(latA, lonA, latB, lonB, radius= 6371221.3, testVerbose= False):
  '''Returns the array of the azimuth between two points, positioned by \
their geographic coordinates along the surface of a perfect sphere \
in units given by the radius used. Input variables include:\n
- latA, latB: the latitude of the points considered in decimal degrees,\n
- lonA, lonB: the longitude of the points considered in decimal degrees,\n
- radius: the radius of the sphere in metres, set by default to \
that of Earth (6371221.3 m).
NOTE: azimuth is computed relative to the point specified by A and positive CW from N'''
	#-make arrays if needed
  if isinstance(latA,float):
		latA= np.array(latA)
  if isinstance(lonA,float):
		lonA= np.array(lonA)
  if isinstance(latB,float):
		latB= np.array(latB)
  if isinstance(lonB,float):
		lonB= np.array(lonB)
  #-pad latitudes, longitudes
  if latA.size <> latB.size:
    latA= np.ones(latB.shape)*latA
  if lonA.size <> lonB.size:
    lonA= np.ones(lonB.shape)*lonA
  pA= latA*deg2Rad; pB= latB*deg2Rad
  lA= lonA*deg2Rad; lB= lonB*deg2Rad
  dL= lB-lA; dP= pB-pA
  #-azimuth
  x= np.sin(dL)*np.cos(pB)
  y= np.cos(pA)*np.sin(pB)-np.sin(pA)*np.cos(pB)*np.cos(dL)
  azimuth= deg2Rad**-1*np.arctan2(x, y)+360.
  azimuth= azimuth % 360.
  if testVerbose:
		for f_a in [0,0.25,0.5,0.75,1.0,1.25,1.5,1.75]:
			angle= deg2Rad**-1*f_a*m_pi
			mask= selectNearestAngle(azimuth*deg2Rad,angle*deg2Rad,None)
			wStr= ' * angle %f matched to ' % angle
			for a in azimuth[mask]:
				wStr+= '%f' % a
			print wStr
  #-return azimuth
  return azimuth

def getIntermediatePoint(latA, lonA, latB, lonB, radius= 6371221.3, fractionPath= 0.5):
	'''Computes the location between two points'''
	#-make arrays if needed
	if isinstance(latA,float):
		latA= np.array(latA)
	if isinstance(lonA,float):
		lonA= np.array(lonA)
	if isinstance(latB,float):
		latB= np.array(latB)
	if isinstance(lonB,float):
		lonB= np.array(lonB)
	#-pad latitudes, longitudes
	if latA.size <> latB.size:
		latA= np.ones(latB.shape)*latA
	if lonA.size <> lonB.size:
		lonA= np.ones(lonB.shape)*lonA
	#-get latitude and longitude in radians
	pA= latA*deg2Rad; pB= latB*deg2Rad
	lA= lonA*deg2Rad; lB= lonB*deg2Rad
	#-get distance 
	d= getArcDistance(latA, lonA, latB, lonB)/radius
	#-compute arguments a, b, x, y, z
	a= np.sin((1.-fractionPath)*d)/np.sin(d)
	b= np.sin(fractionPath*d)/np.sin(d)
	x= a*np.cos(pA)*np.cos(lA)+b*np.cos(pB)*np.cos(lB)
	y= a*np.cos(pA)*np.sin(lA)+b*np.cos(pB)*np.sin(lB)
	z= a*np.sin(pA)+b*np.sin(pB)
	pC= np.arctan2(z,(x*x+y*y)**0.5)
	lC= np.arctan2(y, x)
	return pC*deg2Rad**-1, lC*deg2Rad**-1

def getDestinationPoint(latA, lonA, arcDistance, azimuth, radius= 6371221.3):
	'''Computes the destination point given the location of latA, lonA, the azimuth and arc distance'''
	#-make arrays if needed
	if isinstance(latA,float):
		latA= np.array(latA)
	if isinstance(lonA,float):
		lonA= np.array(lonA)
	if isinstance(arcDistance,float):
		arcDistance= np.array(arcDistance)
	if isinstance(azimuth,float):
		azimuth= np.array(azimuth)
	#-pad latitudes, longitudes
	maxSize= max(latA.size, lonA.size, arcDistance.size, azimuth.size)
	if latA.size <> maxSize:
		latA= np.ones((maxSize))*latA.ravel()
	if lonA.size <> maxSize:
		lonA= np.ones((maxSize))*lonA.ravel()
	if arcDistance.size <> maxSize:
		arcDistance= np.ones((maxSize))*arcDistance.ravel()
	if azimuth.size <> maxSize:
		azimuth= np.ones((maxSize))*azimuth.ravel()
	#-get latitude and longitude in radians
	pA= latA*deg2Rad; lA= lonA*deg2Rad
	d= arcDistance/radius
	theta= azimuth*deg2Rad
	#-compute pB, lB
	pB= np.arcsin(np.sin(pA)*np.cos(d)+np.cos(pA)*np.sin(d)*np.cos(theta))
	lB= lA+np.arctan2(np.sin(theta)*np.sin(d)*np.cos(pA),\
		np.cos(d)-np.sin(pA)*np.sin(pB))
	latB= pB*deg2Rad**-1
	lonB= lB*deg2Rad**-1
	lonB+= 540.
	lonB= (lonB % 360.)-180.
	#-return latitude, longitude
	return latB, lonB

def getEquivalentCellArea(arcResolution, radius= 6371221.3,\
	latOrigin= 90.0, lonOrigin= -180., latExtent= 180., lonExtent= 360., testVerbose= False):
	'''Computes an array with the equivalent cell area in metres squared \
for a perfect sphere of given radius subdivided into geographic degrees. \
Input variables include:\n
- arcResolution: the resolution of the map in decimal degrees,\n
- radius: the radius of the sphere in metres, set by default to \
that of Earth (6371221.3 m).'''
	#-set resolution and derive number of rows
	arcResolution= np.array([arcResolution], dtype= np.float64)
	#-set number of rows and columns
	nrRows= int(round(latExtent/arcResolution))
	nrCols= int(round(lonExtent/arcResolution))
	#-echo:
	if testVerbose:
		print ' * generating the cell area in the units of the radius specified over a geographic grid with \
%d, %d rows and columns with a resolution of %f decimal degrees' % (nrRows, nrCols, arcResolution)
	#-generate arrays with latitudes and longitudes at cell-centres
	latitude= np.ones((nrRows, nrCols))*latOrigin
	longitude= np.ones((nrRows, nrCols))*lonOrigin
	for colCnt in xrange(nrCols):
		latitude[:,colCnt]-= (np.arange(nrRows)+0.5)*arcResolution
	for rowCnt in xrange(nrRows):
		longitude[rowCnt,:]+= (np.arange(nrCols)+0.5)*arcResolution
	#-check on generated extent
	deltaLat= (latitude[-1,-1]/(latOrigin-(latExtent-0.5*arcResolution))-1.)*100.
	deltaLon= (longitude[-1,-1]/(lonOrigin+(lonExtent-0.5*arcResolution))-1.)*100.
	if testVerbose:
		print '   the differences between the generated latitude, longitude amount to \
 %.4f, %.4f %s' % (deltaLat, deltaLon,'%')
	#-compute cell area in units of radius
	cellArea= np.abs(np.sin((latitude+0.5*arcResolution)*deg2Rad)-np.sin((latitude-0.5*arcResolution)*deg2Rad))
	cellArea*= arcResolution*deg2Rad*radius**2
	if testVerbose:
			print '   the average, min and max cell areas in units of the radius of length %f are \
%.f, %.f, and %.f' % (radius,cellArea.mean(), cellArea.min(), cellArea.max()),
	totalArea= cellArea.sum()
	sphereArea= 4*m_pi*radius**2 
	deltaArea= (totalArea/sphereArea-1.)*100.
	if testVerbose:
			print 'equivalent to a total area of %f compared to %f for an ideal sphere, \
with a deviation of %f %s' % (totalArea, sphereArea, deltaArea, '%')	
	#-return cellArea and coordinates
	return cellArea, latitude, longitude
	
def findDistance(distance,azimuth,func,searchDirection= None, devAzimuth= None):
	'''Finds the distance and corresponding azimuth specified by func \
from arrays of distance and azimuth within a predefined sector'''
	if isinstance(searchDirection, NoneType):
		mask= np.ones(distance.shape, dtype= bool)
		selectedDistance= func(distance[mask])
		selectedAzimuth= azimuth[distance[mask] == selectedDistance]
	else:
		mask= selectNearestAngle(azimuth,searchDirection, devAzimuth)
		selectedDistance= func(distance[mask])
		selectedAzimuth= azimuth[distance == selectedDistance][0]
	return float(selectedDistance), float(selectedAzimuth)

def selectNearestAngle(a,t_a, d_a):
        '''returns from an array of the angle, a, the mask where a is nearest to t_a'''
        if isinstance(d_a, NoneType):
                absDev= np.abs(np.sin(a)-np.sin(t_a))+\
                        np.abs(np.cos(a)-np.cos(t_a))
                return absDev == absDev.min()
        else:
                return np.abs(a-t_a) < d_a

def mapEllipse(r_a, r_b, azimuth, latC, lonC, lat, lon, numberStepsArc= 360):
	'''Maps an ellipse for the given spatial attributes using a predefined number \
of steps along the arc:
 r_a: length of semi_major axis
 r_b: length of semi_major axis
 azimuth: angle of semi_major axis with the centroid to north
 latC: latitude of centroid
 lonC: longitude of centroid
 lat: map of latitude
 lon: map of longitude
 numberStepsArc= 360 (default)

where r= r_a*r_b/(r_a**2*sin(theta)**2+r_b**2*cos(theta)**2)**0.5 \n'''
	#-compute theta as the azimuth and rotated towards the semi-major axis
	theta= np.linspace(0.,360.-360./numberStepsArc,numberStepsArc)
	radius= (r_a*r_b)/(r_a**2*(np.sin(theta*deg2Rad))**2+r_b**2*(np.cos(theta*deg2Rad))**2)**0.5
	theta= (theta+azimuth) % 360.
	lat1, lon1= getDestinationPoint(latC, lonC, radius, theta)
	ellipse= pcr.boolean(0)
	for iCnt in xrange(lat1.size):
		ellipse= ellipse | \
			(pcr.abs(lat-lat1[iCnt]) == pcr.mapminimum(pcr.abs(lat-lat1[iCnt]))) & \
			(pcr.abs(lon-lon1[iCnt]) == pcr.mapminimum(pcr.abs(lon-lon1[iCnt])))
	ellipse= pcr.ifthen(ellipse, ellipse)
	#-return ellipse
	return ellipse
  
def main():
  latA= np.array([35.0])
  lonA= np.array([45.0])
  latB= np.array([35.0])
  lonB= np.array([135.0])
  print getAzimuth(latA,lonA,latB,lonB, testVerbose= True)
  print getAzimuth(latB,lonB,latA,lonA, testVerbose= True)
  print getDestinationPoint(53.33, -1.75, 124800., 96.)

if __name__ == '__main__':
  sys.exit(main())

