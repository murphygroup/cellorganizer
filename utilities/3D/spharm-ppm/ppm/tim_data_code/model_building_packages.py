# -*- coding: utf-8 -*-

########## !!!! for now, QuadScheme only for 3D case, not 2D

import math
import numpy as np
import sys
import statsmodels.api as sm
from scipy import stats
from scipy import spatial as sps
import os
import scipy.io as sio
import re
import random

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
                                print ('Specified dimension does not match point dimension, reverting to point dimension')
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
                xvals = np.linspace(tilesizex/2.0,1-tilesizex/2.0,numdummy)
                yvals = np.linspace(tilesizey/2.0,1-tilesizey/2.0,numdummy)
                
                ## if 3D pattern, calculated z dimension for points
                if pat.dimension == 3:
                        tilesizez = 1.0/numdummy

                        # zvals = np.linspace(pat.box[2,0]+tilesizez/2,pat.box[2,1]-tilesizez/2,numdummy)
                        zvals = np.linspace(tilesizez/2.0,1-tilesizez/2.0,numdummy)

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
                for pind in range(0,np.shape(dummypts_all)[0]):
                        # if in_hull(dummypts_all[pind,:],c_region) and not in_hull(dummypts_all[pind,:],n_region):
                        #       rdps.append(dummypts_all[pind,:])
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

                        for dpind in range(0,np.shape(dummypts)[0]):
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
                #       self.quadweights = np.zeros([np.shape(dummypts)[0]+pat.npoints,1])
                #       self.responses = np.zeros([np.shape(dummypts)[0]+pat.npoints,1])
                #       tilecounts = np.zeros(np.shape(dummypts)[0])
                #       pt2tile = {}
                #       tilecounter = 0
                #       rps = sps.cKDTree(pat.points)

                #       for dpind in range(0,np.shape(dummypts)[0]):
                #               curdp = dummypts[dpind,:]
                #               tx = (curdp[0]-tilesizex/2,curdp[0]+tilesizex/2)
                #               ty = (curdp[1]-tilesizey/2,curdp[1]+tilesizey/2)
                #               cand_intile = rps.query_ball_point(curdp,round(max(tilesizex,tilesizey)))
                #               in_interval = lambda x,txl,tyl: True if txl[0] < x[0] < txl[1] and tyl[0] < x[1] < tyl[1] else False
                #               intile = [pt for pt in cand_intile if in_interval(pat.points[pt,:],tx,ty)]
                #               for val in intile:
                #                       pt2tile[val] = tilecounter
                #               tilecounts[tilecounter] = len(intile)+1
                #               tilecounter += 1
                        
                ## if any points fall exactly on a tile edge (unlikely), record them
                
                ## calculcate quadweights (1/number of points in tile) and quadresponses (0 if dummy, 1/weight if observed)
                for pt in range(0,pat.npoints):
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

                for dpt in range(0,np.shape(dummypts)[0]):
                        if self.dimension == 3:
                                self.quadweights[dpt+pat.npoints] = (tilesizex*tilesizey*tilesizez)/tilecounts[dpt]
                        else:
                                self.quadweights[dpt+pat.npoints] = (tilesizex*tilesizey)/tilecounts[dpt]
                        self.responses[dpt+pat.npoints] = 0

                for ind in range(0,pat.npoints+np.shape(dummypts)[0]):
                        if self.quadweights[ind] == 0.0:
                                if self.dimension == 3:
                                        self.quadweights[ind] = (tilesizex*tilesizey*tilesizez)/2
                                # else:
                                #       self.quadweights[ind] = (tilesizex*tilesizey)/2
                                self.responses[ind] = 1/self.quadweights[ind]
                                # print self.quadweights[ind]

                self.points = np.concatenate((pat.points,dummypts)) ### !!!!
                # self.numdummy=np.shape(dummypts)[0] ### origin
                self.numdummy = np.shape(self.points)[0]-pat.npoints ### new 20181210

                ## store the pattern for future calculation
                self.pat = pat


#############################################################
# Model fitting:
# main model fitting function, uses weighted general linear models to optimize parameters uses beta version of weighted glm from statsmodels
def FitModel_all(quads,intermediate_output_path,img_roi_names):
        responses_all=[]
        quadweights_all=[]
        covars_all=[]

        betas = []
        failed = []
        for i in range(0,len(quads)):
                img_roi_name=img_roi_names[i]
                quad=quads[i]
                # try:
                covar=np.loadtxt(os.path.join(intermediate_output_path,'others','factor_matrix__'+img_roi_name+'.txt'))
                try:
                        factor_channel_num=get_factor_channel_num()
                        for fi in range(0,len(factor_channel_num)):
                                if fi==0:
                                        covar_=covar[:,fi]
                                else:
                                        covar_=np.column_stack((covar_,covar[:,fi]))
                        covar=covar_
                except IndexError:
                        pass
                xmax,ymax,zmax=find_xmax_ymax_zmax(intermediate_output_path,img_roi_name)
                dummy_pts=quad.points[-quad.numdummy:]
                dummy_pts=np.column_stack((dummy_pts[:,0]*xmax,dummy_pts[:,1]*ymax,dummy_pts[:,2]*zmax)) #param
                dummy_pts=dummy_pts.astype(int)
                print ('generating factor matrix for '+img_roi_name+' dummy points...')
                covar2=points2factor(dummy_pts,intermediate_output_path,img_roi_name)
                cols=()
                for ii in range(0,np.shape(covar2)[0]):
                    cols=cols+(covar2[ii],)
                if len(np.shape(covar))==1:
                        covar=np.asarray([np.asarray([cr]) for cr in covar])
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
                #       failed.append(i)
                #       print 'failed',i
                #       pass

        covars_all2=covars_all[0]
        for i in range(0,len(quads)-1):
                covars_all2=np.concatenate((covars_all2,covars_all[i+1]))

        # print np.mean(responses_all)
        # print np.mean(covars_all2)
        # print np.mean(quadweights_all)

        mod_glm = sm.GLM(responses_all,sm.tools.tools.add_constant(covars_all2),family=sm.families.Poisson(sm.families.links.log), freq_weights = quadweights_all)
        res_glm = mod_glm.fit(disp=0)
        full_model = res_glm.params
        return full_model,betas,failed

#############################################################
# cv likelihood tests
def inhomoLogLike(test_img_roi_names,coefs,intermediate_output_path,num_rand):
        # num_rand = 70000
        # generate a uniform sample from the cell region
        # rp = genRandPPBox(num_rand) ##### share rd ?
        # calculate the covariate*coefficient value at all points
        print ('generating factor matrix for random points...')
        pre_logls=[]
        for test_img_roi_name in test_img_roi_names:
                print (test_img_roi_name)
                rp = genRandPPBox(num_rand)
                trends = trendValue(coefs,intermediate_output_path,test_img_roi_name,rp)
                norm = sum(trends)*(len(trends)**(-2.0/3)) 
                trends2 = trendValue(coefs,intermediate_output_path,test_img_roi_name)
                pre_logls.append(np.mean(np.log(trends2/norm)))
        return np.mean(pre_logls)

def genRandPPBox(n):
        x=stats.uniform.rvs(0,1,n)
        y=stats.uniform.rvs(0,1,n)
        z=stats.uniform.rvs(0,1,n)
        pos = np.column_stack((x,y,z))
        return PointPattern(pos,dimension=3)

def genRandPP(pat,n,win,c_region,n_region):
        # win is an np.array where win=[xmin xmax; ymin ymax; zmin zmax]
        if np.shape(win)[0] == 3:
                x=stats.uniform.rvs(win[0,0],win[0,1],n)
                y=stats.uniform.rvs(win[1,0],win[1,1],n)
                z=stats.uniform.rvs(win[2,0],win[2,1],n)
                pos = np.column_stack((x,y,z))
        else:
                x=stats.uniform.rvs(win[0,0],win[0,1],n)
                y=stats.uniform.rvs(win[1,0],win[1,1],n)
                pos = np.column_stack((x,y))

        for p in range(0,np.shape(pos)[0]):
                if in_hull(pos[p,:],c_region) and not in_hull(pos[p,:],n_region):
                        continue
                else: 
                        if np.shape(win)[0] == 3:
                                inflag = False
                                while inflag is False:
                                        x=stats.uniform.rvs(win[0,0],win[0,1],1)
                                        y=stats.uniform.rvs(win[1,0],win[1,1],1)
                                        z=stats.uniform.rvs(win[2,0],win[2,1],1)

                                        if in_hull(np.column_stack((x,y,z)),c_region) and not in_hull(np.column_stack((x,y,z)),n_region):
                                                pos[p,0] = x
                                                pos[p,1] = y
                                                pos[p,2] = z
                                                inflag = True
                        else:
                                inflag = False
                                while inflag is False:
                                        x=stats.uniform.rvs(win[0,0],win[0,1],1)
                                        y=stats.uniform.rvs(win[1,0],win[1,1],1)

                                        if in_hull(np.column_stack((x,y)),c_region) and not in_hull(np.column_stack((x,y)),n_region):
                                                pos[p,0] = x
                                                pos[p,1] = y
                                                inflag = True
        return PointPattern(pos,mark=pat.marks,area=pat.cellarea,cellpat=pat.cellpat.data,nucpat=pat.nucpat.data,dimension=pat.dimension,box=pat.box)

def trendValue(coefs,intermediate_output_path,img_roi_name,pat=DEFAULT()):
        covars=[]
        if isinstance(pat,DEFAULT):
                factor_channel_num=get_factor_channel_num()
                # cov_all=[]
                # for img_roi_name in img_roi_names:
                covars_current=[]
                factor=np.loadtxt(os.path.join(intermediate_output_path,'others','factor_matrix__'+img_roi_name+'.txt'),dtype=float)
                covars_current.append(np.ones((np.shape(factor)[0],1))) #be param
                try:
                        for fi in range(0,len(factor_channel_num)):
                                covars_current.append(factor[:,fi])
                except IndexError:
                        covars_current.append(factor) 
                # cov_all.append(covars_current)
                # covars=concatenate(cov_all) #concatenate all test factor_matrices
                covars=covars_current
        else:
                covars.append(np.ones((np.shape(pat.points)[0],1))) #be param
                xmax,ymax,zmax=find_xmax_ymax_zmax(intermediate_output_path,img_roi_name)
                # points=np.column_stack((pat.points[:,0]*(xmax-1),pat.points[:,1]*(ymax-1),pat.points[:,2]*(zmax-1))) #param
                points=np.column_stack((pat.points[:,0]*xmax,pat.points[:,1]*ymax,pat.points[:,2]*zmax)) #param
                points=points.astype(int)
                factor=points2factor(points,intermediate_output_path,img_roi_name)
                covars=covars+factor
        trends=np.exp([sum((coefs[i]*covars[i][pind] for i in range(0,len(coefs)))) for pind in range(0,len(covars[0]))])
        return trends

def points2factor(points,intermediate_output_path,img_roi_name):
        factor=[]
        factor_channel_num=get_factor_channel_num()
        for fi in range(0,len(factor_channel_num)): ###3 be param
                mat_contents = sio.loadmat(os.path.join(intermediate_output_path,'dist_trans_matrix','dist_trans_matrix__'+img_roi_name+'-'+str(factor_channel_num[fi])+'.mat'))
                dist_matrix=mat_contents['m']
                factor_col=np.zeros((np.shape(points)[0],1)) ###2 be param
                # print temp_folder_name+'dist_trans_matrix/'+img_name+'__dist_trans_matrix_'+str(factor_channel_num[fi])+'.mat'
                for ci in range(0,np.shape(points)[0]):
                    x=int(points[ci][0])-1
                    y=int(points[ci][1])-1
                    z=int(points[ci][2])-1
                    # try:
                    factor_col[ci]=dist_matrix[x][y][z]
                    # except:
                    #   print x
                    #   print y
                    #   return 
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

def find_xmax_ymax_zmax(intermediate_output_path,img_ROI_name):
        # factor_channel_num = []
        # with open('factor_channel_num.txt') as f:
        #       for line in f:
        #               factor_channel_num=[int(x) for x in line.split(',')]
        path=os.path.join(intermediate_output_path,'dist_trans_matrix','dist_trans_matrix__'+img_ROI_name+'-1.mat')
        mat_contents = sio.loadmat(path)
        dist_matrix=mat_contents['m']
        # factor_col[ci]=dist_matrix[x][y][z]   
        xmax=np.shape(dist_matrix)[0]
        ymax=np.shape(dist_matrix)[1]
        zmax=np.shape(dist_matrix)[2]
        # #for zmax
        # files=getAbsoluteFilePathsList(intermediate_output_path+'others/')
        # regex = re.compile(img_name+'__dist_trans_matrix_1')
        # z_files=filter(regex.search, files)
        # return xmax,ymax,len(z_files)
        return xmax,ymax,zmax

def get_files_by_regex(folder,regex_in):
        files=os.listdir(folder)
        regex=re.compile(regex_in)
        searched_list=sorted_nicely(filter(regex.search,files))
        return searched_list

def get_substring_between(str1,str2,str_orign):
        result=re.search(str1+'(.*)'+str2,str_orign)
        return result.group(1)

def concatenate(cov_all):
        covars=[]
        for cov_ind in range(0,len(cov_all)):
            col=cov_all[cov_ind]
            for col_ind in range(0,len(col)):
                if cov_ind==0:
                    covars.append([col[col_ind]])
                else:
                    covars[col_ind]=np.concatenate((covars[col_ind][0],col[col_ind]))
        return covars

def get_factor_channel_num():
        factor_channel_num = []
        with open('factor_channel_num.txt') as f:
                for line in f:
                        factor_channel_num=[int(x) for x in line.split(',')]
        return factor_channel_num

def get_cv_assignment_real(cv_assignment,img_names,img_roi_names):
    cv_assignment_real=[]
    for group in cv_assignment:
        cur_group=[]
        for element in group:
                new_group_str=get_strs_by_regex(img_names[element],img_roi_names)
                cur_group=cur_group+[img_roi_names.index(cv_str) for cv_str in new_group_str]
                cv_assignment_real.append(cur_group)
    return cv_assignment_real

def get_strs_by_regex(str_key,str_list):
        regex=re.compile('.*'+str_key+'.*')
        searched_list=sorted_nicely(filter(regex.search,str_list))
        return searched_list

def set2list(seq): 
   # order preserving
   checked = []
   for e in seq:
       if e not in checked:
           checked.append(e)
   return checked