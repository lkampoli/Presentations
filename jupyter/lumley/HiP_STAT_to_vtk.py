infile = './STAT/STAT.xmf'

out_dir = './vtk/'
out_file = out_dir + 'LS5TI5.vtk'
out_header = out_dir + 'header.txt'

var_mean = ['$ \\overline{\\rho}$']
var_fluc = ['$ \\tau_{11}$','$ \\tau_{12}$','$ \\tau_{22}$','$ \\tau_{33}$','$ \\tau_{13}$','$ \\tau_{23}$']
var_list=var_mean + var_fluc
print var_list
fhead=open(out_header,'w')
for var_str in var_list:
	var_tmp = var_str+'\n'
	fhead.write(var_tmp)
fhead.close()

scale_x = 1.0  #/1.011
scale_y = 1.0  #-1.0/1.011
delta = 1.57322  #/1.011
print scale_x,scale_y

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'XDMF Reader'
sTATxmf = XDMFReader(FileNames=[infile])

# Properties modified on sTATxmf
sTATxmf.PointArrayStatus = var_list
sTATxmf.GridStatus = ['Time 0', 'Time 0[1]', 'Time 0[2]', 'Time 0[3]', 'Time 0[4]', 'Time 0[5]', 'Time 0[6]', 'Time 0[7]', 'Time 0[8]', 'Time 0[9]', 'Time 0[10]']

# create a new 'Transform'
transform1 = Transform(Input=sTATxmf)
# Properties modified on transform1.Transform
transform1.Transform.Scale = [scale_x, scale_y, 1.0]

# create a new 'Transform'
transform2 = Transform(Input=transform1)
# Properties modified on transform1.Transform
transform2.Transform.Translate = [0.0, delta, 0.0]

# create a new 'Transform'
transform3 = Transform(Input=transform1)
# Properties modified on transform1.Transform
transform3.Transform.Translate = [0.0, -delta, 0.0]

# create a new 'Group Datasets'
groupDatasets1 = GroupDatasets(Input=[transform1, transform2])
# create a new 'Group Datasets'
groupDatasets2 = GroupDatasets(Input=[groupDatasets1, transform3])
# create a new 'Merge Blocks'
mergeBlocks1 = MergeBlocks(Input=groupDatasets2)

# save data
SaveData(out_file, proxy=mergeBlocks1)

