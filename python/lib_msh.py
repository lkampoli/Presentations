"""
Author: Mathieu Deuse

This module contains classes and functions to manipulate 2D grids and curves.

Classes:

 - CurveFunc: describes a curve function, i.e. a function of a single parameter
   which associate a point [x(s),y(s)] to any value 's' of this parameter. This
   class is designed to manipulate a curve geometry without needing to discretize it.

 - Curve: describes a curve given by a list of point coordinates. This class is
   mainly a wrapper for 2D Numpy array (self.coord) of shape [N,2] where N is the
   number of points in the curve and the second dimension gives the x and y
   coordinates of each point.

 - Block: describes a single structured grid. This class is mainly a wrapper for 3D
   Numpy array (self.coord) of shape [Ni,Nj,2] where Ni and Nj are the grid dimensions
   and the third dimension gives the x and y coordinates of each point.

Functions:

 - read_plot3d and save_plot3d for the conversion between Plot3D files and Numpy arrays.

 - read_mask and save_mask to read/write masks files of overset grids
"""

import os
import numpy as np
import struct as st
from pandas import read_csv
from matplotlib import pyplot as plt
from scipy.interpolate import interp1d
from scipy.optimize import root,minimize


def read_plot3d(filename,has_header=True):

    s_int = st.calcsize('1i')
    s_float = st.calcsize('1f')

    f = open(filename,'rb')

    packed = f.read(3*s_int)
    (nxp,nyp,nzp) = st.unpack('3i',packed)
    n_pts = nxp*nyp*nzp

    if has_header:
        packed = f.read(4*s_float)
        header = st.unpack('4f',packed)
    else:
        header = None

    data = f.read()

    f.close()

    arr = np.frombuffer(data,dtype=np.float32).reshape(nxp,nyp,nzp,-1,order='F').astype(float)

    return (arr,header)


def save_plot3d(arr,filename,header=(0.,0.,0.,0.)):

    s_int = st.calcsize('1i')
    s_float = st.calcsize('1f')

    (nxp,nyp,nzp,nvar) = arr.shape

    f = open(filename,'wb+')

    f.write(st.pack('3i',nxp,nyp,nzp))

    if header is not None:
        f.write(st.pack('4f',*header))

    data = np.getbuffer(arr.flatten(order='F').astype(np.float32))

    f.write(data)

    f.close()


class CurveFunc:
    def __init__(self,func,bounds=(0.0,1.0)):
        """
        "func" is a function of one single parameter, taking values in bounds (default is [0,1]) an return [x,y] for this parameter.
        "func" must be such that it can be called with a 1D Numpy array 's' and return two 1D Numpy arrays 'x' and 'y' (same size)
        """
        self.func = func
        self.bounds = bounds

    def compute(self,s):
        """
        Return x and y arrays of the curve function for parameter array s
        """
        # make it work whether s is a float or a numpy array
        if not isinstance(s,np.ndarray):
            s = np.array([s])

        return self.func(s)
    
    def __call__(self,s):
        """
        Return a single array which concatenates x and y for parameter array s
        """
        (x,y) = self.compute(s)
        return np.array([x,y]).T

    def draw(self,fig=None,color='k',fill_color=None,n_points=1000):
        s = np.linspace(self.bounds[0],self.bounds[1],n_points)
        (x,y) = self.compute(s)

        if fill_color is None:
            plt.plot(x,y,color,figure=fig)
        else:   
            plt.fill(x,y,edgecolor=color,facecolor=fill_color)
    
    def show(self,color='k',fill_color=None,n_points=1000):
        fig = plt.figure()
        self.draw(fig=fig,color=color,fill_color=fill_color,n_points=n_points)
        plt.axis('equal')
        plt.show()

    def curvature(self,n_discretization=1000):
        """
        Returns a function that gives the curvature (inverse of the radius of curvature) of the curve
        """
        
        s = np.linspace(self.bounds[0],self.bounds[1],n_discretization)
        (x,y) = self.compute(s)
        
        dx = (x[2:] - x[:-2])/2.0
        dy = (y[2:] - y[:-2])/2.0
        
        ddx = x[2:] - 2.0*x[1:-1] + x[:-2]
        ddy = y[2:] - 2.0*y[1:-1] + y[:-2]
        
        kappa = np.abs((dx*ddy-dy*ddx)/(dx**2+dy**2)**(3.0/2.0))
        
        # constant extrapolation at extremities
        kappa = np.insert(kappa,0,kappa[0])
        kappa = np.insert(kappa,-1,kappa[-1])
        
        curvature = interp1d(s,kappa,kind='linear',copy=False,assume_sorted=True)
        
        return lambda t : curvature(t)

    def length(self,n_discretization=10000):
        """
        Return the length of the curve
        'n_discretization' corresponds to the number of points used to compute the length of the curve
        """

        s = np.linspace(self.bounds[0],self.bounds[1],n_discretization)
        (x,y) = self.compute(s)

        return np.sum(np.sqrt((x[1:]-x[:-1])**2 + (y[1:]-y[:-1])**2))

    def scale_param(self,new_bounds=[0.0,1.0]):
        
        def func(s):
            # change for parameter in [0,1]
            s = (s-new_bounds[0])/(new_bounds[1]-new_bounds[0])

            # change to the old bounds
            s = s*(self.bounds[1]-self.bounds[0])+self.bounds[0]

            return self.func(s)

        return CurveFunc(func,new_bounds)
    
    def crop(self,new_bounds):
        """
        Return a new curve function corresponding to the initial curve cropped in the given interval
        """

        return CurveFunc(self.func,new_bounds)

    def merge(self,other):
        """
        Merge two curve functions: the parameter interval of the second curved is moved after the first one
        """
        
        def func(s):
            x = np.empty(s.shape)
            y = np.empty(s.shape)

            ind = s <= self.bounds[1]
            (x[ind],y[ind]) = self.compute(s[ind])

            ind = s > self.bounds[1]
            (x[ind],y[ind]) = other.compute(s[ind]-self.bounds[1]+other.bounds[0])

            return(x,y)

        return CurveFunc(func,[self.bounds[0],self.bounds[1]+other.bounds[1]-other.bounds[0]])

    def natural(self,n_discretization=10000,kind='linear'):
        """
        Returns a new curve function with natural (arc-length) parametrization
        'n_discretization' corresponds to the number of points used to compute the length of the curve
        """
        
        s = np.linspace(self.bounds[0],self.bounds[1],n_discretization)
        (x,y) = self.compute(s)
        
        dl = np.sqrt((x[1:]-x[:-1])**2 + (y[1:]-y[:-1])**2)
        l = np.zeros(s.shape)
        l[1:] = np.cumsum(dl) # l[-1] is the length of the curve

        natural_param = interp1d(l,s,kind='linear',copy=False,assume_sorted=True,fill_value='extrapolate')

        func_nat = lambda t : self.func(natural_param(t))
        
        return CurveFunc(func_nat,[0.0,l[-1]])

    def intersect_axis(self,direction,loc,guess=None,dist_check=None,check_bounds=True):
        """
        Return the parameter value and the location of the intersection of the curve with
        an horizontal axis (direction is 'x' or 0) or a vertical axis (direction is 'y' or 1)
        For horizontal (vertical) axis, 'loc' is the y coordinate (x coordinate) of the axis
        An initial guess can be defined for the parameter value to improve convergence
        If a value is given for 'dist_check', raise an exception if an intersection is not found for that tolerance
        If 'check_bounds' is True, raise an exception if the intersection is found outside of the curve bounds
        """

        if direction=='x': direction = 0
        if direction=='y': direction = 1

        if direction==0:
            def func(s):
                (x,y) = self.compute(s)
                return y-loc
        elif direction==1:
            def func(s):
                (x,y) = self.compute(s)
                return x-loc
        else:
            raise ValueError("Incorrect direction")

        if guess is None:
            guess = 0.5*(self.bounds[0]+self.bounds[1])

        s = root(func,guess).x[0]
        p = self(s)[0,:]

        if check_bounds:
            if s<self.bounds[0] or s>self.bounds[1]:
                raise RuntimeError("Unable to find an intersection in curve bounds")

        if dist_check is not None:
            if direction==0:
                dist = abs(p[1]-loc)
            elif direction==1:
                dist = abs(p[0]-loc)

            if dist > dist_check:
                raise RuntimeError("Unable to find an intersection")

        return (s,p)

    def intersect(self,other,guess=None,dist_check=None):
        """
        Return the parameter value in the two curves, and the location, of the intersection of two curves functions
        An initial guess can be defined for [s1,s2] to improve convergence
        If a value is given for 'dist_check', raise an exception if an intersection is not found for that tolerance
        """
        
        def dist_func(s):
            c1 = self(s[0])
            c2 = other(s[1])
            return (c1[0,0]-c2[0,0])**2 + (c1[0,1]-c2[0,1])**2

        if guess is None:
            guess1 = 0.5*(self.bounds[0]+self.bounds[1])
            guess2 = 0.5*(other.bounds[0]+other.bounds[1])
            guess = [guess1,guess2]
        res = minimize(dist_func,guess,bounds=[self.bounds,other.bounds],tol=1e-20)

        s1 = res.x[0]    # the value at which distance is minimum for airfoil curve
        s2 = res.x[1]    # the value at which distance is minimum for vertical line
        p1 = self(s1)[0,:]   # the coordinates on airfoil curve from optimization
        p2 = other(s2)[0,:]  # The coordinates on vertical line from optimization

        if dist_check is not None:
            dist = np.sqrt((p1[0]-p2[0])**2 + (p1[1]-p2[1])**2)
            if dist > dist_check:
                raise RuntimeError("Unable to find an intersection")

        p = (p1+p2)/2.

        return (s1,s2,p)

    def tangent(self,s,ds=1e-4):
        """
        Return normalized tangent vector at parameter s
        Computed using 4th order central difference with step ds
        """

        tang = (self(s-2.0*ds) - 8.0*self(s-ds) + 8.0*self(s+ds) - self(s+2.0*ds))/(12.0*ds)

        norm = np.sqrt(tang[:,0]**2 + tang[:,1]**2)
        tang[:,0] = tang[:,0]/norm
        tang[:,1] = tang[:,1]/norm

        # eliminate useless dimension
        if tang.shape[0] == 1:
            tang = tang[0,:]

        return tang

    def normal(self,s,ds=1e-4):
        """
        Return normalized normal vector at parameter s
        Computed using 4th order central difference with step ds
        """

        tang = self.tangent(s,ds=ds)

        n = np.empty(tang.shape)

        if len(n.shape) == 2:
            n[:,0] =  tang[:,1]
            n[:,1] = -tang[:,0]
        elif len(n.shape) == 1:
            n[0] =  tang[1]
            n[1] = -tang[0]

        return n

    def translate(self,vec):
        """
        Return a new curve function corresponding to the initial curve function
        translated by a given vector
        """

        v_x = vec[0]
        v_y = vec[1]

        def func(s):
            (x,y) = self.compute(s)
            return (x+v_x,y+v_y)

        return CurveFunc(func,bounds=self.bounds)

    def rotate(self,angle,center=[0.,0.]):
        """
        Return a new curve function corresponding to the initial curve function
        rotated by a specified angle (in degrees) around a specified center
        """

        theta = angle*np.pi/180.0
        c_x = center[0]
        c_y = center[1]

        def func(s):

            (x,y) = self.compute(s)

            x -= c_x
            y -= c_y

            x_new = x*np.cos(theta) - y*np.sin(theta)
            y_new = x*np.sin(theta) + y*np.cos(theta)

            x_new += c_x
            y_new += c_y

            return (x_new,y_new)

        return CurveFunc(func,bounds=self.bounds)

    def scale(self,fact,center=[0.,0.]):
        """
        Return a new curve function corresponding to the initial curve function
        scaled by some factor around a specified center
        If 'fact' is a single real number then same scaling in all directions
        Different scaling can be used in x and y direction with fact=[fact_x,fact_y]
        """

        c_x = center[0]
        c_y = center[1]

        if isinstance(fact,list) or isinstance(fact,tuple):
            fact_x = fact[0]
            fact_y = fact[1]
        else:
            fact_x = fact
            fact_y = fact

        def func(s):

            (x,y) = self.compute(s)

            x_new = (x - c_x)*fact_x + c_x
            y_new = (y - c_y)*fact_y + c_y

            return (x_new,y_new)

        return CurveFunc(func,bounds=self.bounds)

    @staticmethod
    def line(p_start,p_end):
        """
        Curve function of line between p_start and p_end
        """
        
        def func(s):
            x = (1.0-s)*p_start[0] + s*p_end[0]
            y = (1.0-s)*p_start[1] + s*p_end[1]
            return (x,y)
        
        return CurveFunc(func)
 
    @staticmethod
    def circle_arc(center,start_point,ang):
        """
        Curve function of a circle arc (parametrized by s in [0,1]), first point is "start_point", and the opening angle is "ang" (in radians, positive means counter-clockwise)
        """
        
        alpha = np.arctan2(start_point[1]-center[1],start_point[0]-center[0])
        r = np.sqrt((start_point[0]-center[0])**2+(start_point[1]-center[1])**2)

        def func(s):
            x = r*np.cos(alpha+ang*s) + center[0]
            y = r*np.sin(alpha+ang*s) + center[1]
            return (x,y)

        return CurveFunc(func)

    @staticmethod
    def bezier(points):
        """
        Bezier curve with a given array of control points
        """

        if len(points) == 1:
            def func(s):
                return (points[0][0],points[0][1])
        else:
            def func(s):
                (x1,y1) = CurveFunc.bezier(points[:-1]).compute(s)
                (x2,y2) = CurveFunc.bezier(points[1:]).compute(s)
                x = (1.0-s)*x1 + s*x2
                y = (1.0-s)*y1 + s*y2
                return (x,y)

        return CurveFunc(func)

    @staticmethod
    def naca(t=0.12,m=0.0,p=0.4):
        """
        Curve function (parametrized by s in [0,1]) of a naca 4-digit airfoil of thickness 'r', camber 'm', and location of max camber 'p'
        """
        
        def func(s):
            x = (np.cos(2.0*np.pi*s)+1.0)/2.0
            
            yt = t/0.2*(0.2969*np.sqrt(x)-0.1260*x-0.3516*x**2+0.2843*x**3-0.1036*x**4)
            
            if m == 0.0:
                y = np.empty(s.size)
                y[s<0.5] = -yt[s<0.5]
                y[s>=0.5] = yt[s>=0.5]
            
            else:
                yc = np.empty(s.size)
                yc[x<p] = m*x[x<p]/p**2*(2.0*p-x[x<p])
                yc[x>=p] = m*(1.0-x[x>=p])/(1.0-p)**2*(1.0+x[x>=p]-2.0*p)
                
                theta = np.empty(s.size)
                theta[x<p] = 2.0*m/p**2*(p-x[x<p])
                theta[x>=p] = 2.0*m/(1.0-p)**2*(p-x[x>=p])
                
                x[s<0.5] = x[s<0.5] + yt[s<0.5]*np.sin(theta[s<0.5])
                x[s>=0.5] = x[s>=0.5] - yt[s>=0.5]*np.sin(theta[s>=0.5])
                
                y = np.empty(s.size)
                y[s<0.5] = yc[s<0.5] - yt[s<0.5]*np.cos(theta[s<0.5])
                y[s>=0.5] = yc[s>=0.5] + yt[s>=0.5]*np.cos(theta[s>=0.5])
            
            return (x,y)

        return CurveFunc(func)

    @staticmethod
    def cd_airfoil():
        """
        Curve function of the controlled-diffusion airfoil at design angle of attack (8 deg)
        """

        # find the directory of 'lib_msh.py', which should also contain 'CD_coord.npz'
        path = os.path.realpath(__file__)
        i = path.rfind('/')
        path = path[:i]

        cd_npz = np.load(path+'/CD_coord.npz')
        coord = cd_npz['coord']

        return Curve(coord=coord).curve_func(kind='cubic')


class Curve:
    def __init__(self,coord=None,curve_func=None,n_points=None,density=None):
        """
        Either the coordinates are given directly, or the curve is created from a curve function object, with a given number of points distributed according to a density function defined in [0,1] (if no density is given, constant density is assumed)
        """
    
        if coord is not None:
            self.coord = coord
            
        elif (curve_func is not None and n_points is not None):
            if (density is None):
                # constant density is assumed
                s = np.linspace(0.0,1.0,n_points)
            else:
                # generate a repartition from the given density
                t = np.linspace(0.0,1.0,n_points-1)
                d = density(t)
                s = np.cumsum(1.0/d)
                s = np.insert(s,0,0.0)
                s = s/s[-1]
            
            # scale parameter in [0,1]
            self.coord = curve_func.scale_param()(s)
    
    def __getitem__(self,key):
        return self.coord[key]

    def __setitem__(self,key,value):
        self.coord[key] = value

    def size(self):
        return self.coord.shape[0]
    
    def flip(self):
        self.coord = np.flip(self.coord,0)
        return self
    
    def draw(self,fig=None,color='.-k'):
        plt.plot(self[:,0],self[:,1],color,figure=fig)
    
    def show(self,color='.-k'):
        fig = plt.figure()
        self.draw(fig=fig,color=color)
        plt.axis('equal')
        plt.show()

    def saveas(self,filename):

        f = open(filename+'.dat','w+')

        f.write(' %i'%self.size())

        f.writelines(['\n  %.18f %.18f'%(p[0],p[1]) for p in self.coord])

        f.close()

    def split(self,indexes):
        curves = []
        
        ind = list(indexes)
        
        ind.insert(0,0)
        ind.append(self.size()-1)
        
        for i in range(len(ind)-1):
            coord = self[ind[i]:(ind[i+1]+1),:]
            curves.append(Curve(coord=coord))
   
        return curves

    def tangent(self,i=None):
        """
        Return normalized tangent vector at index i, or at each point if i is not specified
        Computed using 4th order central difference
        Periodicity is assumed at the extremities, with overlapping point
        """

        # -1 not to use the overlapping point twice
        n = self.size() - 1

        if i is None:
            tang = np.empty(self.coord.shape)

            for i in range(n):
                imm = (i-2)%n
                im  = (i-1)%n
                ip  = (i+1)%n
                ipp = (i+2)%n
                tang[i,:] = self[imm,:] - 8.0*self[im,:] + 8.0*self[ip,:] - self[ipp,:]

            # copy overlapping point
            tang[-1,:] = tang[0,:]

            # normalize
            norm = np.sqrt(tang[:,0]**2 + tang[:,1]**2)
            tang[:,0] = tang[:,0]/norm
            tang[:,1] = tang[:,1]/norm

        else:
            imm = (i-2)%n
            im  = (i-1)%n
            ip  = (i+1)%n
            ipp = (i+2)%n
            tang = self[imm,:] - 8.0*self[im,:] + 8.0*self[ip,:] - self[ipp,:]
            norm = np.sqrt(tang[0]**2 + tang[1]**2)
            tang = tang/norm

        return tang

    def normal(self,i=None):
        """
        Return normalized normal vector at index i, or at each point if i is not specified
        Computed using 4th order central difference
        """

        tang = self.tangent(i)
        n = np.empty(tang.shape)

        if i is None:
            n[:,0] =  tang[:,1]
            n[:,1] = -tang[:,0]
        else:
            n[0] =  tang[1]
            n[1] = -tang[0]

        return n

    def arc_length(self):
        """
        Compute arc length coordinate to be used to integrate along the curve
        """

        s = np.empty(self.size())
        s[0] = 0.0

        for i in range(1,self.size()):
            l = np.sqrt((self[i,0]-self[i-1,0])**2 + (self[i,1]-self[i-1,1])**2)
            s[i] = s[i-1] + l

        return s

    def length(self):
        """
        Return the length of the curve
        """

        return self.arc_length()[-1]

    def smooth_repartition(self,n=1):
        """
        Smooth points repartition by updating point location with the average of their neighbours 'n' times
        """

        for _ in range(n):
            self[1:-1,:] = 0.5*(self[:-2,:]+self[2:,:])

        return self

    def curve_func(self,kind='cubic'):
        """
        Return a curve function (parameter in [0,1]) that interpolates this curve
        """
        
        s = np.linspace(0.0,1.0,self.size())

        interp_x = interp1d(s,self.coord[:,0],kind=kind,copy=False,assume_sorted=True,fill_value='extrapolate')
        interp_y = interp1d(s,self.coord[:,1],kind=kind,copy=False,assume_sorted=True,fill_value='extrapolate')

        func = lambda s : (interp_x(s),interp_y(s))

        return CurveFunc(func,[0.0,1.0])

    def density_func(self):
        """
        Return a function (parameter in [0,1]) that interpolates the point density of the curve
        """

        d = 1.0/np.sqrt((self[1:,0]-self[:-1,0])**2 + (self[1:,1]-self[:-1,1])**2)

        s = np.linspace(0,1,d.size)

        d_func = interp1d(s,d,kind='cubic',copy=False,assume_sorted=True,fill_value='extrapolate')

        return d_func

    @staticmethod
    def merge(curves,ind=False):
        """
        For each pair of curves, the linking point is supposed to be present in both curves
        If 'ind' is true, a list of index that can be used again to split is returned in addition
        """
        
        new_curve = Curve()
        ind_list = []

        # the linking points are not repeated
        
        new_curve.coord = np.reshape(curves[0].coord[0,:],(1,2))
        
        for c in curves:
            new_curve.coord = np.concatenate([new_curve.coord,c.coord[1:,:]],axis=0)
            ind_list.append(new_curve.size() - 1)
        
        # remove the last index as it is meaningless
        ind_list.pop()

        if ind:
            return (new_curve,ind_list)
        else:
            return new_curve
    
    def copy(self):
        c = Curve()
        c.coord = self.coord.copy()
        return c

    @staticmethod
    def filter_array(arr,amp,n_steps,peri=True,order=2,overlap=True):
        """
        Filter a 1D Numpy array to improve metrics continuity. Modify the array directly, no copy
        IMPORTANT: if periodic filtering is used, must specify if the first and last point overlap or not
        """

        if order not in [2,4]:
            raise ValueError("Only order 2 and 4 implemented")

        if peri and overlap:
            # remove the last point for the filtering
            arr_ = arr[:-1].copy()
        else:
            arr_ = arr.copy()

        tmp = arr_.copy()

        if peri:
            if order==2:
                for _ in range(n_steps):
                    tmp[0,   :] = arr_[0   ] + amp*(arr_[1   ] - 2.0*arr_[0   ] + arr_[  -1])
                    tmp[1:-1,:] = arr_[1:-1] + amp*(arr_[ :-2] - 2.0*arr_[1:-1] + arr_[2:  ])
                    tmp[  -1,:] = arr_[  -1] + amp*(arr_[  -2] - 2.0*arr_[  -1] + arr_[0   ])
                    arr_[:] = tmp
            elif order==4:
                for _ in range(n_steps):
                    tmp[0   ,:] = arr_[0   ] - amp*(arr_[  -2] - 4.0*arr_[  -1] + 6.0*arr_[0   ] - 4.0*arr_[1   ] + arr_[2   ])
                    tmp[1   ,:] = arr_[1   ] - amp*(arr_[  -1] - 4.0*arr_[0   ] + 6.0*arr_[1   ] - 4.0*arr_[2   ] + arr_[3   ])
                    tmp[2:-2,:] = arr_[2:-2] - amp*(arr_[ :-4] - 4.0*arr_[1:-3] + 6.0*arr_[2:-2] - 4.0*arr_[3:-1] + arr_[4:  ])
                    tmp[  -2,:] = arr_[  -2] - amp*(arr_[  -4] - 4.0*arr_[  -3] + 6.0*arr_[  -2] - 4.0*arr_[  -1] + arr_[0   ])
                    tmp[  -1,:] = arr_[  -1] - amp*(arr_[  -3] - 4.0*arr_[  -2] + 6.0*arr_[  -1] - 4.0*arr_[0   ] + arr_[1   ])
                    arr_[:] = tmp
        else:
            if order==2:
                for _ in range(n_steps):
                    arr_[1:-1] = arr_[1:-1] + amp*(arr_[:-2] - 2.0*arr_[1:-1] + arr_[2:])
            else:
                raise Exception("Non periodic filter only implemented for order 2")

        if peri and overlap:
            # add the last point that has been removed
            arr[:-1] = arr_
            arr[ -1] = arr_[0]
        else:
            arr[:] = arr_

    def filter(self,amp,n_steps,peri=True,order=2,overlap=True):
        """
        Filter the curve to improve metrics continuity
        IMPORTANT: if periodic filtering is used, must specify if the first and last point overlap or not
        """

        Curve.filter_array(self[:,0],amp,n_steps,peri=peri,order=order,overlap=overlap)
        Curve.filter_array(self[:,1],amp,n_steps,peri=peri,order=order,overlap=overlap)
    
        return self

    def rotate(self,angle,center=[0.,0.]):
        """
        Rotation of a specified angle (in degrees) around a specified center
        """

        theta = angle*np.pi/180.0

        for i in range(self.size()):
            x = self[i,0] - center[0]
            y = self[i,1] - center[1]

            x_new = x*np.cos(theta) - y*np.sin(theta)
            y_new = x*np.sin(theta) + y*np.cos(theta)

            self[i,0] = x_new + center[0]
            self[i,1] = y_new + center[1]

        return self


class Block:
    def __init__(self,arg=None,fmt='ascii'):
        """
        Initialises a block
        If 'arg' is not provided then the coordinates are not initialised
        If 'arg' is a string then read the grid from a file, 'fmt' (format) can be 'ascii' (default) or 'plot3d'
        If 'arg' is a numpy array then this array correspond to the block coordinates
        If 'arg' is a list of four curves then the block is generated as a Coons patch
        If 'arg' is a list containing x and y arrays then the block is generated as a meshgrid
        """
        
        if arg is not None:
            if isinstance(arg,str):

                if fmt == 'ascii':

                    f = open(arg+'.dat','r')

                    tmp = f.readline().split()
                    try:
                        tmp.remove(',')
                    except ValueError:
                        pass
                    (n_i,n_j) = [int(n) for n in tmp[:2]]

                    self.coord = read_csv(f,delim_whitespace=True,header=None).values.reshape(n_i,n_j,2,order='F')

                    f.close()

                elif fmt == 'plot3d':

                    (self.coord,_) = read_plot3d(arg+'.xyz',has_header=False)

                    (n_i,n_j,_,_) = self.coord.shape

                    self.coord = self.coord[:,:,0,:2]

                else:
                    raise ValueError("Incorrect format")

            elif isinstance(arg,np.ndarray):
                self.coord = np.asfortranarray(arg)
                    
            elif isinstance(arg,list):

                if len(arg) == 2:
                    x = arg[0]
                    y = arg[1]
                    self.coord = np.empty((len(x),len(y),2),order='F')
                    (self.coord[:,:,0],self.coord[:,:,1]) = np.meshgrid(x,y,indexing='ij')

                elif len(arg) == 4:
                    self.coord = Block.coons_patch(arg[0],arg[1],arg[2],arg[3])

    def __getitem__(self,key):
        return self.coord[key]

    def __setitem__(self,key,value):
        self.coord[key] = value

    @staticmethod
    def coons_patch(xm,xp,ym,yp):
        """
        xm,xp,xm,yp are the boundary curves ('Curve' objects or numpy arrays) as named in hipstar
        """

        if isinstance(ym,Curve):
            n_i = ym.size()
        else:
            n_i = ym.shape[0]

        if isinstance(xm,Curve):
            n_j = xm.size()
        else:
            n_j = xm.shape[0]
        
        coord = np.empty((n_i,n_j,2),order='F')
        
        for i in range(n_i):
            for j in range(n_j):
                s = float(i)/(n_i-1)
                t = float(j)/(n_j-1)
                
                lc = (1.0-t)*ym[i,:] + t*yp[i,:]
                ld = (1.0-s)*xm[j,:] + s*xp[j,:]
                
                b = ym[0,:]*(1.0-s)*(1.0-t) + ym[-1,:]*s*(1.0-t) + yp[0,:]*(1.0-s)*t + yp[-1,:]*s*t
                
                coord[i,j,:] = lc + ld - b
        
        return coord
    
    def size(self):
        return (self.coord.shape[0],self.coord.shape[1])

    def num_points(self):
        (n_i,n_j) = self.size()
        return n_i*n_j

    def saveas(self,filename,fmt='ascii'):
        """
        'fmt' (format) can be 'ascii' (default) or 'plot3d'
        """
        (n_i,n_j) = self.size()
        
        if fmt == 'ascii':

            f = open(filename+'.dat','w+')

            f.write(' %i %i'%(n_i,n_j))

            f.writelines(['\n  %.18f %.18f'%(c[0],c[1]) for c in np.reshape(self.coord,[-1,2],order='F')])

            f.close()

        elif fmt == 'plot3d':

            # zeros for third coordinate
            z = np.full((n_i,n_j,1,1),0.0,order='F')

            arr = np.concatenate((self.coord.reshape(n_i,n_j,1,2,order='F'),z),axis=3)

            save_plot3d(arr,filename+'.xyz',header=None)

        else:
            raise ValueError("Incorrect format")
    
    def draw(self,fig=None,color='k',option='wireframe',stride=1):
        if option == 'wireframe':
            plt.pcolormesh(self[::stride,::stride,0],
                           self[::stride,::stride,1],
                           np.zeros(self[::stride,::stride,0].shape),
                           edgecolor=color,
                           cmap='binary',
                           linewidth=0.4,
                           antialiased=True,
                           figure=fig)
        else:
            sides = [self.extract_line(ind) for ind in ['ym','xp','yp','xm']]
            if option == 'contour':
                for c in sides:
                    c.draw(fig=fig,color=color)
            elif option == 'filled':
                sides[2].flip()
                sides[3].flip()
                x = np.concatenate([s[:,0] for s in sides])
                y = np.concatenate([s[:,1] for s in sides])
                plt.fill(x,y,color,figure=fig)

    def show(self,color='k',option='wireframe',stride=1):
        fig = plt.figure()
        self.draw(fig=fig,color=color,option=option,stride=stride)
        plt.axis('equal')
        plt.show()

    def show_boundaries(self):

        print self.size()

        fig = plt.figure(figsize=(9.0,6.5))

        xm = self.extract_line('xm')
        plt.plot(xm[:,0],xm[:,1],label='xm')

        xp = self.extract_line('xp')
        plt.plot(xp[:,0],xp[:,1],label='xp')

        ym = self.extract_line('ym')
        plt.plot(ym[:,0],ym[:,1],label='ym')

        yp = self.extract_line('yp')
        plt.plot(yp[:,0],yp[:,1],label='yp')

        plt.axis('equal')
        plt.tight_layout()
        plt.legend(loc='best')
        plt.show()

    def extrapolate(self,bound,n_points):
        """
        Use second order extrapolation to add 'n_points' points next to boundary 'bound'
        """

        (n_i,n_j) = self.size()

        if bound == 'xm':
            ext = np.empty((n_points,n_j,2))
            self.coord = np.concatenate([ext,self.coord],axis=0)
            for i in range(n_points-1,-1,-1):
                self[i,:,:] = 3.0*self[i+1,:,:] - 3.0*self[i+2,:,:] + self[i+3,:,:]
        elif bound == 'xp':
            ext = np.empty((n_points,n_j,2))
            self.coord = np.concatenate([self.coord,ext],axis=0)
            for i in range(n_i,n_i+n_points):
                self[i,:,:] = 3.0*self[i-1,:,:] - 3.0*self[i-2,:,:] + self[i-3,:,:]
        elif bound == 'ym':
            ext = np.empty((n_i,n_points,2))
            self.coord = np.concatenate([ext,self.coord],axis=1)
            for j in range(n_points-1,-1,-1):
                self[:,j,:] = 3.0*self[:,j+1,:] - 3.0*self[:,j+2,:] + self[:,j+3,:]
        elif bound == 'yp':
            ext = np.empty((n_i,n_points,2))
            self.coord = np.concatenate([self.coord,ext],axis=1)
            for j in range(n_j,n_j+n_points):
                self[:,j,:] = 3.0*self[:,j-1,:] - 3.0*self[:,j-2,:] + self[:,j-3,:]

        return self
                
 
    @staticmethod
    def merge(blocks,axis=0,ind=False):
        """
        The boundary along which two blocks are merge is supposed to be identical in both blocks
        If 'ind' is true, a list of index that can be used again to split is returned in addition
        """
        
        new_block = Block()
        ind_list = []

        if axis == 0:
            new_block.coord = np.reshape(blocks[0].coord[0,:,:],(1,-1,2))
            
            # To avoid duplicate the boundaries, the first row of the first block has been appended separately, then each block is appended without its first row
            for b in blocks:
                new_block.coord = np.concatenate([new_block.coord,b.coord[1:,:,:]],axis=axis)
                ind_list.append(new_block.coord.shape[0] - 1)
                
        elif axis == 1:
            new_block.coord = np.reshape(blocks[0].coord[:,0,:],(-1,1,2))
            
            # To avoid duplicate the boundaries, the first row of the first block has been appended separately, then each block is appended without its first row
            for b in blocks:
                new_block.coord = np.concatenate([new_block.coord,b.coord[:,1:,:]],axis=axis)
                ind_list.append(new_block.coord.shape[1] - 1)
        
        # remove the last index as it is meaningless
        ind_list.pop()

        if ind:
            return (new_block,ind_list)
        else:
            return new_block
    
    def split(self,indexes,axis=0):
        blocks = []
        
        ind = list(indexes)
        
        ind.insert(0,0)
        ind.append(self.coord.shape[axis]-1)
        
        if axis == 0:
            for i in range(len(ind)-1):
                b = Block()
                b.coord = self.coord[ind[i]:(ind[i+1]+1),:,:]
                blocks.append(b)
        elif axis == 1:
            for i in range(len(ind)-1):
                b = Block()
                b.coord = self.coord[:,ind[i]:(ind[i+1]+1),:]
                blocks.append(b)
        
        return blocks
    
    def flip(self,axis=0):
        """
        Flip the current block, not a copy
        """
        self.coord = np.flip(self.coord,axis)
        return self

    def swap_xy(self):
        """
        Swap the current block, not a copy
        """
        self.coord = np.swapaxes(self.coord,0,1)
        return self
    
    def copy(self):
        return Block(self.coord.copy(order='F'))
    
    def extract_line(self,index,axis=0):
        """
        Index can be an integer or a string ('xm','xp','ym','yp') 
        """
        
        if isinstance(index,str):
            if index == 'xm':
                c = Curve(self[0,:,:])
            if index == 'xp':
                c = Curve(self[-1,:,:])
            if index == 'ym':
                c = Curve(self[:,0,:])
            if index == 'yp':
                c = Curve(self[:,-1,:])
            
        else:
            if axis == 0:
                c = Curve(self[index,:,:])
            elif axis == 1:
                c = Curve(self[:,index,:])
        
        return c

    def find_nearest(self,p):
        """
        Returns the (i,j) indexes of the nearest point to point 'p' in the block
        """

        dist2 = (self[:,:,0]-p[0])**2 + (self[:,:,1]-p[1])**2

        ind_flat = np.argmin(dist2.flatten(order='F'))

        (i,j) = np.unravel_index(ind_flat,dist2.shape,order='F')

        return (i,j)

    def filter(self,amp,n_steps,axis=0,peri=True,order=2,overlap=True):
        """
        Filter the grid to improve metrics continuity
        IMPORTANT: if periodic filtering is used, must specify if the first and last point overlap or not
        """

        if axis is not None and axis not in [0,1]:
            raise ValueError("Incorrect axis")

        if order not in [2,4]:
            raise ValueError("Only order 2 and 4 implemented")

        if peri and overlap:
            # remove the last point for the filtering
            if axis==0:
                self.coord = self[:-1,:,:]
            elif axis==1:
                self.coord = self[:,:-1,:]

        tmp = np.copy(self.coord)

        if peri:
            if order==2:
                if axis==0:
                    for _ in range(n_steps):
                        tmp[0,   :,:] = self[0   ,:,:] + amp*(self[1   ,:,:] - 2.0*self[0   ,:,:] + self[  -1,:,:])
                        tmp[1:-1,:,:] = self[1:-1,:,:] + amp*(self[ :-2,:,:] - 2.0*self[1:-1,:,:] + self[2:  ,:,:])
                        tmp[  -1,:,:] = self[  -1,:,:] + amp*(self[  -2,:,:] - 2.0*self[  -1,:,:] + self[0   ,:,:])
                        self[:,:,:] = tmp
                elif axis==1:
                    for _ in range(n_steps):
                        tmp[:,0,   :] = self[:,0   ,:] + amp*(self[:,1   ,:] - 2.0*self[:,0   ,:] + self[:,  -1,:])
                        tmp[:,1:-1,:] = self[:,1:-1,:] + amp*(self[:, :-2,:] - 2.0*self[:,1:-1,:] + self[:,2:  ,:])
                        tmp[:,  -1,:] = self[:,  -1,:] + amp*(self[:,  -2,:] - 2.0*self[:,  -1,:] + self[:,0   ,:])
                        self[:,:,:] = tmp
            elif order==4:
                if axis==0:
                    for _ in range(n_steps):
                        tmp[0   ,:,:] = self[0   ,:,:] - amp*(self[  -2,:,:] - 4.0*self[  -1,:,:] + 6.0*self[0   ,:,:] - 4.0*self[1   ,:,:] + self[2  ,:,:])
                        tmp[1   ,:,:] = self[1   ,:,:] - amp*(self[  -1,:,:] - 4.0*self[0   ,:,:] + 6.0*self[1   ,:,:] - 4.0*self[2   ,:,:] + self[3  ,:,:])
                        tmp[2:-2,:,:] = self[2:-2,:,:] - amp*(self[ :-4,:,:] - 4.0*self[1:-3,:,:] + 6.0*self[2:-2,:,:] - 4.0*self[3:-1,:,:] + self[4: ,:,:])
                        tmp[  -2,:,:] = self[  -2,:,:] - amp*(self[  -4,:,:] - 4.0*self[  -3,:,:] + 6.0*self[  -2,:,:] - 4.0*self[  -1,:,:] + self[0  ,:,:])
                        tmp[  -1,:,:] = self[  -1,:,:] - amp*(self[  -3,:,:] - 4.0*self[  -2,:,:] + 6.0*self[  -1,:,:] - 4.0*self[0   ,:,:] + self[1  ,:,:])
                        self[:,:,:] = tmp
                elif axis==1:
                    for _ in range(n_steps):
                        tmp[:,0   ,:] = self[:,0   ,:] - amp*(self[:,  -2,:] - 4.0*self[:,  -1,:] + 6.0*self[:,0   ,:] - 4.0*self[:,1   ,:] + self[:,2  ,:])
                        tmp[:,1   ,:] = self[:,1   ,:] - amp*(self[:,  -1,:] - 4.0*self[:,0   ,:] + 6.0*self[:,1   ,:] - 4.0*self[:,2   ,:] + self[:,3  ,:])
                        tmp[:,2:-2,:] = self[:,2:-2,:] - amp*(self[:, :-4,:] - 4.0*self[:,1:-3,:] + 6.0*self[:,2:-2,:] - 4.0*self[:,3:-1,:] + self[:,4: ,:])
                        tmp[:,  -2,:] = self[:,  -2,:] - amp*(self[:,  -4,:] - 4.0*self[:,  -3,:] + 6.0*self[:,  -2,:] - 4.0*self[:,  -1,:] + self[:,0  ,:])
                        tmp[:,  -1,:] = self[:,  -1,:] - amp*(self[:,  -3,:] - 4.0*self[:,  -2,:] + 6.0*self[:,  -1,:] - 4.0*self[:,0   ,:] + self[:,1  ,:])
                        self[:,:,:] = tmp
        else:
            if order==2:
                if axis==0:
                    for _ in range(n_steps):
                        self[1:-1,:,:] = self[1:-1,:,:] + amp*(self[ :-2,:,:] - 2.0*self[1:-1,:,:] + self[2:  ,:,:])
                elif axis==1:
                    for _ in range(n_steps):
                        self[:,1:-1,:] = self[:,1:-1,:] + amp*(self[:, :-2,:] - 2.0*self[:,1:-1,:] + self[:,2:  ,:])
            else:
                raise Exception("Non periodic filter only implemented for order 2")

        if peri and overlap:
            # add the last point that has been removed
            if axis==0:
                self.coord = np.concatenate([self.coord,self[0:1,:,:]],axis=0)
            elif axis==1:
                self.coord = np.concatenate([self.coord,self[:,0:1,:]],axis=1)

        return self


class Mesh:
    def __init__(self,blocks,fmt='ascii'):
        """
        Initialization either by a given string corresponding to the mesh file name, or a list of either block objects or arguments to initialize new blocks
        If initialised from files, 'fmt' (format) can be 'ascii' (default) or 'plot3d'
        """

        if isinstance(blocks,str):
            # the argument is the prefix of block file names

            if fmt == 'ascii':
                ext = '.dat'
            elif fmt == 'plot3d':
                ext = '.xyz'
            else:
                raise ValueError("Incorrect format")

            # get the directory to be searched
            i = blocks.rfind('/') + 1
            if i == 0:
                directory = './'
                name = blocks
            else:
                directory = blocks[:i]
                name = blocks[i:]
            filenames = [f[:-4] for f in os.listdir(directory) if len(f) > len(name)
                                                      and f[:len(name)] == name
                                                      and f[-len(ext):] == ext       ]
            # remove invalid names:
            for s in filenames:
                try:
                    int(s[(len(name)+1):])
                except ValueError:
                    filenames.remove(s)

            filenames.sort(key = lambda s : int(s[(len(name)+1):]))

            self.blocks = [Block(directory+f,fmt=fmt) for f in filenames]

        elif isinstance(blocks,list):

            self.blocks = []
            for b in blocks:
                if isinstance(b,Block):
                    self.blocks.append(b)
                else:
                    self.blocks.append(Block(b,fmt=fmt))

        else:
            raise ValueError("Incorrect argument, mesh cannot be initialised")

    def __getitem__(self,key):
        return self.blocks[key]

    def __setitem__(self,key,value):
        self.blocks[key] = value

    def num_blocks(self):
        return len(self.blocks)

    def num_points(self):
        n_points = 0
        for b in self.blocks:
            n_points += b.num_points()
        return n_points

    def draw(self,fig=None,colors=None,option='wireframe',stride=1):
        if colors is None:
            colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
        
        for i,block in enumerate(self.blocks):
            block.draw(fig=fig,color=colors[i%len(colors)],option=option,stride=stride)

    def show(self,colors=None,option='wireframe',legend=True,stride=1):
        if colors is None:
            colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

        fig = plt.figure(figsize=(9.0,6.5))
        
        if legend:
            # dummy plot to set legend labels
            for i in range(len(self.blocks)):
                plt.plot(np.nan,np.nan,'s',color=colors[i%len(colors)],figure=fig,label='Block '+str(i+1))
        
        self.draw(fig=fig,colors=colors,option=option,stride=stride)
        
        plt.axis('equal')
        plt.tight_layout()
        if legend:
            plt.legend(loc='best')
        plt.show()
    
    def saveas(self,filename,fmt='ascii'):
        for i in range(len(self.blocks)):
            self.blocks[i].saveas(filename+'_'+str(i+1),fmt=fmt)


def read_mask(filename):
    """
    Reads the mask of an overset grid (ascii format)
    """

    f = open(filename,'r')

    tmp = f.readline().split()
    try:
        tmp.remove(',')
    except ValueError:
        pass
    (n_i,n_j) = [int(n) for n in tmp[:2]]

    m = read_csv(f,delim_whitespace=True,header=None).values.reshape(n_i,n_j,order='F')

    f.close()

    return m

def save_mask(m,filename):
    """
    Save the mask of an overset grid (ascii format)
    """

    (n_i,n_j) = m.shape

    f = open(filename,'w+')

    f.write(' %i %i'%(n_i,n_j))

    f.writelines(['\n  %i'%(flag) for flag in np.reshape(m,[-1],order='F')])

    f.close()
