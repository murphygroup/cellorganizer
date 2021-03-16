# -*- coding: utf-8 -*-

import math
import numpy as np
import sys
import statsmodels.api as sm
from scipy import stats
from scipy import spatial as sps
import os
import re

#############################################################
# Classes

## this is used to check if optional arguments are supplied
# DEFAULT = object()
class DEFAULT:
	def __init__(self):
		self.self = 'default'

class PointPattern:
	"""

	Point pattern object, holds all necessary components of a point pattern.

	dimension : dimension of point pattern (int)
	points : array of point coordinates (nparray)
	npoints : number of observed points (int)
	marked : flag indicating whether a pattern is marked (bool)
	multimarked : flag indicating whether a pattern has multiple marks (bool)
	marks : a 1 by npoints list of strings, tuples, ints, or floats corresponding to the mark/label of each point. ORDER MATTERS HERE (list)
	umarks : unique marks in the point pattern (list)
	cellpat : array of cell membrane coordinates (nparray)
	cellflag : flag, true if cellpat is provided (bool)
	cellarea : area contained within cell (float)
	nucpat : array of nuclear membrane coordinates (nparray)
	nucflag : flag, true if nucpat is provided (bool)
	box : minimum box containing all points (np.array, [xmin, xmax; ymin, ymax; zmin, zmax])


	"""

	## main function to create a point pattern from data
	def __init__(self,pos,dimension=DEFAULT()):
		## get correct dimension
		if not isinstance(dimension,DEFAULT):
			self.dimension = dimension
		else:
			if np.shape(pos)[1] != dimension:
				print 'Specified dimension does not match point dimension, reverting to point dimension'
				self.dimension = np.shape(pos)[1]
			else:
				self.dimension = dimension

		## store the number of points in the pattern
		self.npoints = np.shape(pos)[0]
		
		## sort the points by x values in ascending order (this will speed calculations of nearest neighbors, covar dist, etc later)
		self.points = pos[pos[:,0].argsort()]

##### delete many helper functions of PointPattern

class QuadScheme:
	"""
	Quadrature scheme for Berman Turner

	dimension : dimension of quad points (int)
	points : array of all point coordinates, real followed by dummy (nparray)
	ptlabs : array of point identities, 0 if dummy, 1 if real (nparray)
	quadweights : array of quadrature weights corresponding to each point in both the pattern and quadscheme, the first npoints are real, all following are dummy points (nparray)
	responses : array of calculated responses, equal to 1/quadweight for each real point, 0 for each dummy point (nparray)
	ptsonedge : array of points that happen to fall on the edge of a quad tile (unlikely), nparray
	nrp : number of real data points (int)
	numdummy : number of dummy points per dimension (int)
	pat : point pattern that the quadscheme was generated from/for (PointPattern)

	"""
	
	## initialize the quadrature scheme
		# numdummy : desired number of points on any given dimension (actually the cube (or square) root of the desired total number of dummy points)
		# creates a quadrature scheme based on equal sized tiling
	def __init__(self , pat , numdummy):

		## store the quad dimension
		self.dimension = pat.dimension

		## calculate the tile sizes
		# tilesizex = (pat.box[0,1]-pat.box[0,0])/numdummy
		# tilesizey = (pat.box[1,1]-pat.box[1,0])/numdummy
		tilesizex = 1.0/numdummy
		tilesizey = 1.0/numdummy

		## determine the dummy point coordinates by spacing evenly over each dimension
		# xvals = np.linspace(pat.box[0,0]+tilesizex/2,pat.box[0,1]-tilesizex/2,numdummy)
		# yvals = np.linspace(pat.box[1,0]+tilesizey/2,pat.box[1,1]-tilesizey/2,numdummy)
		xvals = np.linspace(tilesizex/2,1-tilesizex/2,numdummy)
		yvals = np.linspace(tilesizey/2,1-tilesizey/2,numdummy)
		
		## if 3D pattern, calculated z dimension for points
		if pat.dimension == 3:
			tilesizez = 1.0/numdummy

			# zvals = np.linspace(pat.box[2,0]+tilesizez/2,pat.box[2,1]-tilesizez/2,numdummy)
			zvals = np.linspace(0+tilesizez/2,1-tilesizez/2,numdummy)

			## manipulate the coordinate arrays to match point pattern formatting
			xv, yv, zv = np.meshgrid(xvals, yvals, zvals)
			xv = np.ndarray.flatten(xv)
			yv = np.ndarray.flatten(yv)
			zv = np.ndarray.flatten(zv)

			dummypts_all = np.column_stack((xv,yv,zv))

		else:
			## manipulate the coordinate arrays to match point pattern formatting
			xv, yv = np.meshgrid(xvals, yvals)
			xv = np.ndarray.flatten(xv)
			yv = np.ndarray.flatten(yv)
			dummypts_all = np.column_stack((xv,yv))

		# c_hull = sps.ConvexHull(pat.cellpat.data)
		# n_hull = sps.ConvexHull(pat.nucpat.data)
		# c_region = sps.Delaunay(pat.cellpat.data[c_hull.vertices])
		# n_region = sps.Delaunay(pat.nucpat.data[n_hull.vertices])
		
		## sort the dummy points for faster calculation
		rdps = []
		dummypts_all = dummypts_all[dummypts_all[:,0].argsort()]
		for pind in xrange(0,np.shape(dummypts_all)[0]):
			# if in_hull(dummypts_all[pind,:],c_region) and not in_hull(dummypts_all[pind,:],n_region):
			# 	rdps.append(dummypts_all[pind,:])
			rdps.append(dummypts_all[pind,:])
		dummypts = np.array(rdps)

		## pt2tile dict to store which point is in which tile
		pt2tile = {}
		tilecounter = 0

		## calculate the number of points in each tile for each point in both dummy and observed pattern
		if pat.dimension == 3:
			## set quadweights and quad responses
			self.quadweights = np.zeros([np.shape(dummypts)[0]+pat.npoints,1])
			self.responses = np.zeros([np.shape(dummypts)[0]+pat.npoints,1])
			tilecounts = np.zeros(np.shape(dummypts)[0])
			pt2tile = {}
			tilecounter = 0
			rps = sps.cKDTree(pat.points)

			for dpind in xrange(0,np.shape(dummypts)[0]):
				curdp = dummypts[dpind,:]
				tx = (curdp[0]-tilesizex/2,curdp[0]+tilesizex/2)
				ty = (curdp[1]-tilesizey/2,curdp[1]+tilesizey/2)
				tz = (curdp[2]-tilesizez/2,curdp[2]+tilesizez/2)
				cand_intile = rps.query_ball_point(curdp,round(max(tilesizex,tilesizey,tilesizez)))
				in_interval = lambda x,txl,tyl,tzl: True if txl[0] < x[0] < txl[1] and tyl[0] < x[1] < tyl[1] and tzl[0] < x[2] < tzl[1] else False
				intile = [pt for pt in cand_intile if in_interval(pat.points[pt,:],tx,ty,tz)]
				for val in intile:
					pt2tile[val] = tilecounter
				tilecounts[tilecounter] = len(intile)+1
				tilecounter += 1
		# else:
		# 	self.quadweights = np.zeros([np.shape(dummypts)[0]+pat.npoints,1])
		# 	self.responses = np.zeros([np.shape(dummypts)[0]+pat.npoints,1])
		# 	tilecounts = np.zeros(np.shape(dummypts)[0])
		# 	pt2tile = {}
		# 	tilecounter = 0
		# 	rps = sps.cKDTree(pat.points)

		# 	for dpind in xrange(0,np.shape(dummypts)[0]):
		# 		curdp = dummypts[dpind,:]
		# 		tx = (curdp[0]-tilesizex/2,curdp[0]+tilesizex/2)
		# 		ty = (curdp[1]-tilesizey/2,curdp[1]+tilesizey/2)
		# 		cand_intile = rps.query_ball_point(curdp,round(max(tilesizex,tilesizey)))
		# 		in_interval = lambda x,txl,tyl: True if txl[0] < x[0] < txl[1] and tyl[0] < x[1] < tyl[1] else False
		# 		intile = [pt for pt in cand_intile if in_interval(pat.points[pt,:],tx,ty)]
		# 		for val in intile:
		# 			pt2tile[val] = tilecounter
		# 		tilecounts[tilecounter] = len(intile)+1
		# 		tilecounter += 1
			
		## if any points fall exactly on a tile edge (unlikely), record them
		
		## calculcate quadweights (1/number of points in tile) and quadresponses (0 if dummy, 1/weight if observed)
		for pt in xrange(0,pat.npoints):
			if pt in pt2tile:
				curtile = pt2tile[pt]
				if self.dimension == 3:
					self.quadweights[pt] = (tilesizex*tilesizey*tilesizez)/tilecounts[curtile]
					self.responses[pt] = 1/self.quadweights[pt]
				else:
					self.quadweights[pt] = (tilesizex*tilesizey)/tilecounts[curtile]
					self.responses[pt] = 1/self.quadweights[pt]
			# else:
				# ptsonedge.append(pt)

		for dpt in xrange(0,np.shape(dummypts)[0]):
			if self.dimension == 3:
				self.quadweights[dpt+pat.npoints] = (tilesizex*tilesizey*tilesizez)/tilecounts[dpt]
			else:
				self.quadweights[dpt+pat.npoints] = (tilesizex*tilesizey)/tilecounts[dpt]
			self.responses[dpt+pat.npoints] = 0

		for ind in xrange(0,pat.npoints+np.shape(dummypts)[0]):
			if self.quadweights[ind] == 0.0:
				if self.dimension == 3:
					self.quadweights[ind] = (tilesizex*tilesizey*tilesizez)/2
				# else:
				# 	self.quadweights[ind] = (tilesizex*tilesizey)/2
				self.responses[ind] = 1/self.quadweights[ind]
				# print self.quadweights[ind]

		self.points = np.concatenate((pat.points,dummypts)) ### !!!!
		self.numdummy=np.shape(dummypts)[0]

		## store the pattern for future calculation
		self.pat = pat


#############################################################
# Model fitting:
# main model fitting function, uses weighted general linear models to optimize parameters uses beta version of weighted glm from statsmodels

def FitModel_all(quads,temp_folder_name,img_names):
# def FitModel( quad,temp_folder_name,img_name,covar):
	responses_all=[]
	quadweights_all=[]
	covars_all=[]

	betas = []
	failed = []
	for i in xrange(0,len(quads)):
		img_name=img_names[i]
		quad=quads[i]
		# try:
			covar = np.loadtxt(temp_folder_name+'others/'+img_names[i]+'__factor_matrix.txt')
			zmax=find_zmax(temp_folder_name,img_name)
			dummy_pts=quad.points[-quad.numdummy:]
			dummy_pts=np.column_stack((dummy_pts[:,0]*2048,dummy_pts[:,1]*2048,dummy_pts[:,2]*zmax)) #param
			dummy_pts=dummy_pts.astype(int)
			covar2=points2factor(dummy_pts,temp_folder_name,img_name)
			cols=()
			for ii in xrange(0,np.shape(covar2)[0]):
			    cols=cols+(covar2[ii],)
			covars=np.concatenate((covar,np.column_stack(cols)))

			mod_glm = sm.GLM(quad.responses[:,0],sm.tools.tools.add_constant(covars),family=sm.families.Poisson(sm.families.links.log), freq_weights = quad.quadweights[:,0])
			res_glm = mod_glm.fit(disp=0)
			beta = res_glm.params

			betas.append(beta)
			covars_all.append(covars)
			if i==0:
				responses_all=quad.responses[:,0]
				quadweights_all=quad.quadweights[:,0]
			else:
				responses_all=np.concatenate((responses_all,quad.responses[:,0]))
				quadweights_all=np.concatenate((quadweights_all,quad.quadweights[:,0]))
		# except:
		# 	failed.append(i)
		# 	print 'failed',i
		# 	pass

	covars_all2=covars_all[0]
	for i in xrange(0,len(quads)-1):
		covars_all2=np.concatenate((covars_all2,covars_all[i+1]))

	mod_glm = sm.GLM(responses_all,sm.tools.tools.add_constant(covars_all2),family=sm.families.Poisson(sm.families.links.log), freq_weights = quadweights_all)
	res_glm = mod_glm.fit(disp=0)
	full_model = res_glm.params
	return full_model,betas,failed

#############################################################
# cv likelihood tests
def inhomoLogLike(pat,coefs,temp_folder_name,img_name):
	num_rand = 70000
	# generate a uniform sample from the cell region
	rp = genRandPPBox(pat,num_rand)
	# calculate the covariate*coefficient value at all points
	trends = trendValue(coefs,temp_folder_name,img_name,rp)
	norm = np.mean(trends)
	logdens = sum(np.log(trendValue(coefs,temp_folder_name,img_name)/norm))
	return (logdens/pat.npoints + np.log(.001))

def genRandPPBox(pat,n):
	if pat.dimension == 3:
		x=stats.uniform.rvs(0,1,n)
		y=stats.uniform.rvs(0,1,n)
		z=stats.uniform.rvs(0,1,n)
		pos = np.column_stack((x,y,z))
	else:
		x=stats.uniform.rvs(0,1,n)
		y=stats.uniform.rvs(0,1,n)
		pos = np.column_stack((x,y))
	return PointPattern(pos,dimension=pat.dimension)

def trendValue(coefs,temp_folder_name,img_name,pat=DEFAULT()):
	covars = []
	if isinstance(pat,DEFAULT):
		factor=np.loadtxt(temp_folder_name+'others/'+img_name+'__factor_matrix.txt',dtype=float)
		covars.append(np.ones((np.shape(factor)[0],1))) #be param
		for ii in xrange(0,np.shape(factor)[1]):
			covars.append(factor[:,ii])
	else:
		covars.append(np.ones((np.shape(pat.points)[0],1))) #be param
		zmax=find_zmax(temp_folder_name,img_name)
		points=np.column_stack((pat.points[:,0]*2047,pat.points[:,1]*2047,pat.points[:,2]*(zmax-1))) #param
		points=points.astype(int)
		factor=points2factor(points,temp_folder_name,img_name)
		covars=covars+factor
	return np.exp([sum((coefs[i]*covars[i][pind] for i in xrange(0,len(coefs)))) for pind in xrange(0,len(covars[0]))])

def points2factor(points,temp_folder_name,img_name):
	factor=[]
	files=os.listdir(temp_folder_name+'dist_trans_matrix/')
	for fi in xrange(1,3): ###3 be param
		# fi=1
		regex = re.compile(img_name+'__dist_trans_matrix_'+str(fi))
		cur_img_dm_names=filter(regex.search, files)
		cur_img_dm_names=sorted_nicely(cur_img_dm_names)

		dist_matrix=[]
		for cur_img_dm_name in cur_img_dm_names:
		    dm_z=np.loadtxt(temp_folder_name+'dist_trans_matrix/'+cur_img_dm_name,dtype=float)
		    dist_matrix.append(dm_z)
		    print cur_img_dm_name
		
		factor_col=np.zeros((np.shape(points)[0],1)) ###2 be param
		for ci in xrange(0,np.shape(points)[0]):
		    x=int(points[ci][0])-1
		    y=int(points[ci][1])-1
		    z=int(points[ci][2])-1
		    factor_col[ci]=dist_matrix[z][x][y]
		factor.append(factor_col)
	return factor

def sorted_nicely( l ):
	    """ Sorts the given iterable in the way that is expected.
	 
	    Required arguments:
	    l -- The iterable to be sorted.
	 
	    """
	    convert = lambda text: int(text) if text.isdigit() else text
	    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
	    return sorted(l, key = alphanum_key)

def getAbsoluteFilePathsList(directory):
    files=[]
    for dirpath,_,filenames in os.walk(directory):
        for f in filenames:
            files.append(os.path.abspath(os.path.join(dirpath, f)))
    return files

def find_zmax(temp_folder_name,img_name):
	files=getAbsoluteFilePathsList(temp_folder_name+'others/')
	regex = re.compile(img_name+'__dist_trans_matrix_1')
	z_files=filter(regex.search, files)
	return len(z_files)
