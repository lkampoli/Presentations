#!/usr/bin/python

import os
import sys
import numpy as np
from matplotlib import pyplot as plt
from scipy.interpolate import lagrange
from lib_msh import *
from lib_stretch_func import *
from hyperbolic_grid import *
from elliptic_grid import *

############################## Parameters ##############################

print "-------- Parameters --------"

# Airfoil parameters
x_blunt    = 0.95    # X coordinate where the airfoil TE is cut
n_foil     =  501    # Number of points along the airfoil surface
ds_0       = 1e-3    # Thickness of the first cells on the airfoil surface
ds_max_og  = 6e-3    # Thickness of the cell at the outer boundary of the O-grid
l_og       = 0.25    # Thickness of the BL block around the airfoil (approximation)

# Density functions at airfoil LE and TE
density_le = lambda s : 6.0*np.exp(-(s-0.5)**2*200.0)
density_te = lambda s : 40.0*np.exp(-s**2*250.0)

# Parameters for hyperbolic meshing
alpha    = 10.0
eps_e    = 0.1
v_a      = 0.2

# Parameters of the Cartesian grid
x_front    = -8.0
x_back     = 7.0
y_top      = 7.0
y_bot      = -y_top
y_wake     = 0.2
x_b_corner = 0.3
fact       = 1.5
n_crop     = 4
d_max_wake = 0.075
d_max_out  = 0.2

wake_stretch = 'geom'
out_stretch  = 'tanh'

n_sponge = 20
n_zcbc   = 100

############################## Airfoil surface ##############################

print "-------- Airfoil surface --------"

# force to have an odd number of points along the airfoil (so that there is a point on the LE if symmetric)
if n_foil % 2 == 0:
    n_foil = n_foil + 1

# compute naca function
naca = CurveFunc.naca().natural()
l = naca.length()

# compute point and parameter corresponding to x_blunt
tmp1 = CurveFunc.line([x_blunt,0.0],[x_blunt,0.1])
(s_u,_,p) = naca.intersect(tmp1,guess=[0.99*l,0.5])
tmp1 = CurveFunc.line([x_blunt,0.0],[x_blunt,-0.1])
(s_l,_,_) = naca.intersect(tmp1,guess=[0.01*l,0.5])

n = naca.normal(s_u)

# compute center of the circle arc of the blunt TE
tmp1 = CurveFunc.line(p,p+n)
tmp2 = CurveFunc.line([0.0,0.0],[1.0,0.0])
(_,_,center) = tmp1.intersect(tmp2)

# compute angle of the circle arc
tmp1 = p - center
ang = -np.angle(tmp1[0]+1j*tmp1[1])

# circle arcs of the blunt TE, un upper and lower surface
arc_u = CurveFunc.circle_arc(center,p,ang).natural()
(x,y) = arc_u.compute(arc_u.bounds[1])
arc_l = CurveFunc.circle_arc(center,[x,y],ang).natural()

# create curve function of the airfoil with blunt TE
foil_func = naca.crop([s_l,s_u])
foil_func = arc_l.merge(foil_func)
foil_func = foil_func.merge(arc_u)

# create the actual curve with given density
density = lambda s : 1.0 + density_le(s) + density_te(s) + density_te(1.0-s)
foil = Curve(curve_func=foil_func,n_points=n_foil,density=density)

# correction for unit cord length, and place the TE at x=0
c = foil[0,0]
foil.coord = foil.coord / c
foil[:,0] -= 1.0

#foil.show()


############################## O-grid ##############################

print "-------- O-grid --------"

delta_s = np.diff(stretch('tanh',0.0,l_og,d_min=ds_0,d_max=ds_max_og))

ogrid = Block(hyperbolic_grid(foil.coord,delta_s,bc='periodic',
                              alpha=alpha,eps_e=eps_e,v_a=v_a,
                              dissip_type='linear',
                              alpha_mode='linear',
                              overlap=True))

# Filter to improve metrics continuity (very small not to modify the airfoil geometry)
ogrid.filter(0.02,40,peri=True,axis=0,order=2,overlap=True)

# remove overlap
ogrid.coord = ogrid[:-1,:,:]

#ogrid.show()


############################## Wake block ##############################

print "-------- Wake block --------"

# Split the O-grid in two blocks

c = ogrid.extract_line('yp')
(nx,_) = ogrid.size()

for i1 in range(nx):
    if c[i1,1] < -y_wake: break
for i2 in range(nx-1,-1,-1):
    if c[i2,1] > y_wake: break

i1 += 1

ogrid_2 = Block(np.concatenate((ogrid[i2:,:,:],ogrid[:i1,:,:]),axis=0))

ogrid.coord = ogrid[i1:i2,:,:]

#Mesh([ogrid,ogrid_2]).show()

# Create the initial wake block

c = ogrid_2.extract_line('yp')
(nx,_) = ogrid_2.size()

nj = 8
n_blend = 8

n_bezier_cp = 7
x_bezier_cp = np.linspace(x_b_corner,x_back,n_bezier_cp)

xm = x_back

for i in range(nx):

    x_i = ogrid_2[i,-nj:,0]
    y_i = ogrid_2[i,-nj:,1]

    s = np.arange(1,nj+1)

    ###

    fx = lagrange(s,x_i)
    fy = lagrange(s,y_i)

    s = np.arange(nj,nj+n_blend)
    coord = np.empty((len(s),2))
    coord[:,0] = fx(s)
    coord[:,1] = fy(s)

    c_i = Curve(coord=coord)

    ###

    p0 = ogrid_2[i,-2,:]
    p1 = ogrid_2[i,-1,:]

    n = p1 - p0
    d = np.linalg.norm(n)
    n /= d

    a = (x_b_corner-p1[0])/n[0]
    y_cp = p1[1] + a*n[1]

    p_list = [p1]
    for x_cp in x_bezier_cp:
        p_list.append([x_cp,y_cp])
    cf = CurveFunc.bezier(p_list).natural()
    l = cf.length()

    if i == 0:
        s = stretch(wake_stretch,0,l,d_min=d,d_max=d_max_wake)
        n_wake = len(s)
        b_coord = np.empty((nx,n_wake,2))
    else:
        s1 = stretch(wake_stretch,0,l,d_min=d,n=n_wake)
        s2 = stretch(wake_stretch,0,l,n=n_wake,d_max=d_max_wake)
        s = poly_blend(s1=s1,s2=s2)

    c = Curve(coord=cf(s))

    ###

    s = poly_blend(n_blend)

    for k in range(n_blend):
        c[k,:] = s[k]*c[k,:] + (1-s[k])*c_i[k,:]

    b_coord[i,:,:] = c.coord

wake = Block(b_coord[:,1:,:])

#Mesh([ogrid,ogrid_2,wake]).show()

# Improve orthogonality

coord_list = [b.coord for b in [ogrid,ogrid_2,wake]]

intfs = [(0,'xm',1,'xp'),
         (0,'xp',1,'xm'),
         (1,'yp',2,'ym')]

(new_ogrid,new_ogrid_2,wake) = Mesh(elliptic_grid(coord_list,intfs,n_iter=1)).blocks

#Mesh([new_ogrid,new_ogrid_2,wake]).show()

for i in range(ogrid.size()[0]):
    ogrid[i,:,0] = poly_blend(s1=ogrid[i,:,0],s2=new_ogrid[i,:,0])
    ogrid[i,:,1] = poly_blend(s1=ogrid[i,:,1],s2=new_ogrid[i,:,1])

for i in range(ogrid_2.size()[0]):
    ogrid_2[i,:,0] = poly_blend(s1=ogrid_2[i,:,0],s2=new_ogrid_2[i,:,0])
    ogrid_2[i,:,1] = poly_blend(s1=ogrid_2[i,:,1],s2=new_ogrid_2[i,:,1])

#Mesh([ogrid,ogrid_2,wake]).show()

# Remove some grid lines to avoid bad grid quality near the corner

n_crop = 4

new_ogrid = np.concatenate(( ogrid_2[-n_crop:,:-n_crop,:],
                             ogrid[:,:-n_crop,:],
                             ogrid_2[:n_crop,:-n_crop,:] ),axis=0)

new_ogrid_2 = ogrid_2[n_crop:-n_crop,:-n_crop,:]

new_wake = np.concatenate(( ogrid_2[n_crop:-n_crop,-n_crop:,:],
                            wake[n_crop:-n_crop,:,:] ),axis=1)

ogrid.coord = new_ogrid
ogrid_2.coord = new_ogrid_2
wake.coord = new_wake

#Mesh([ogrid,ogrid_2,wake]).show()

# Add overlapping points for overset grids method

n_ov = 8

ogrid.coord = np.concatenate(( ogrid_2[-n_ov:,:],
                               ogrid[:,:,:],
                               ogrid_2[:n_ov,:,:] ),axis=0)

wake.coord = np.concatenate((ogrid_2.coord,wake.coord),axis=1)

#Mesh([wake,ogrid]).show()


############################## Cartesian block ##############################

print "-------- Cartesian grid --------"

ds = fact*np.amin(np.sqrt((ogrid[:,-1,0]-ogrid[:,-2,0])**2 + (ogrid[:,-1,1]-ogrid[:,-2,1])**2))

x_min = np.amin(ogrid[:,-1,0])
x_max = np.amax(ogrid[:,-1,0])

nx = int((x_max-x_min)/ds)
x = np.linspace(x_min,x_max,nx)
dx = (x_max-x_min)/float(nx-1)

xm = stretch(out_stretch ,x_min,x_front,d_min=dx,d_max=d_max_out)
xp = stretch(out_stretch,x_max,x_back ,d_min=dx,d_max=d_max_wake)

x = np.concatenate((np.flip(xm[1:],0),x,xp[1:]))

###

y = np.flip(wake[:,-1,1])

ym = stretch(out_stretch,y[0],y_bot,d_min=y[1]-y[0],d_max=d_max_out)
yp = stretch(out_stretch,y[-1],y_top,d_min=y[-1]-y[-2],d_max=d_max_out)

n_blend = 8
s = np.arange(n_blend)
fym = lagrange(s,np.flip(y[:n_blend]))
fyp = lagrange(s,y[-n_blend:])

s = np.arange(n_blend,2*n_blend) - 1
ymm = fym(s)
ypp = fyp(s)

ym[:n_blend] = poly_blend(s1=ymm,s2=ym[:n_blend])
yp[:n_blend] = poly_blend(s1=ypp,s2=yp[:n_blend])

y = np.concatenate((np.flip(ym[1:],0),y,yp[1:]))

# Add sponge layer and ZCBC

sponge = np.arange(1,n_sponge+1)*d_max_out
zcbc = np.arange(1,n_zcbc+1)*d_max_wake

x = np.concatenate((x_front-np.flip(sponge,0),x,x_back+zcbc))
y = np.concatenate((y_bot-np.flip(sponge,0),y,y_top+sponge))

# Improve metrics continuity

amp_avg = 0.2
step_avg = 200
for i in range(step_avg):
    x[1:-1] = (1.0-amp_avg)*x[1:-1] + amp_avg*0.5*(x[2:]+x[:-2])
    #y[1:-1] = (1.0-amp_avg)*y[1:-1] + amp_avg*0.5*(y[2:]+y[:-2])

#show_metrics(x)
#show_metrics(y)

# Create block

c_block = Block([x,y])


############################## Save grids ##############################

print "-------- Save grids --------"

mesh = Mesh([c_block,ogrid,wake])
mesh.saveas('z_r_grid')

os.system("show_grid.py z_r_grid")

