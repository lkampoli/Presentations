#### import the simple module from the paraview
from paraview.simple import *

#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'Legacy VTK Reader'
field = LegacyVTKReader(FileNames=['/home/unimelb.edu.au/lcampoli/UoM/Testcases/45_degree_Yaw/VT_NASA_BeVERLI_3D_Hill_Baseline_RANS_k_Omega_SST_Mes_Lev0_64mm_MeshB_run2_SELECTED/postProcessing/surfaces/120000/wallShearStress_Hill_wall.vtk'])

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [1544, 802]

# show data in view
fieldDisplay = Show(field, renderView1)
# trace defaults for the display properties.
fieldDisplay.Representation = 'Surface'
fieldDisplay.ColorArrayName = [None, '']
fieldDisplay.OSPRayScaleArray = 'wallShearStress'
fieldDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
fieldDisplay.SelectOrientationVectors = 'None'
fieldDisplay.ScaleFactor = 0.625
fieldDisplay.SelectScaleArray = 'None'
fieldDisplay.GlyphType = 'Arrow'
fieldDisplay.GlyphTableIndexArray = 'None'
fieldDisplay.DataAxesGrid = 'GridAxesRepresentation'
fieldDisplay.PolarAxes = 'PolarAxesRepresentation'

# update the view to ensure updated data information
renderView1.Update()

# set scalar coloring
ColorBy(fieldDisplay, ('POINTS', 'wallShearStress', 'Magnitude'))

# rescale color and/or opacity maps used to include current data range
fieldDisplay.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
fieldDisplay.SetScalarBarVisibility(renderView1, True)

# hide data in view
Hide(field, renderView1)

# reset view to fit data
renderView1.ResetCamera()

# set active source
SetActiveSource(field)

# show data in view
fieldDisplay = Show(field, renderView1)

# show color bar/color legend
fieldDisplay.SetScalarBarVisibility(renderView1, True)

# destroy renderView1
Delete(renderView1)
del renderView1

# Create a new 'SpreadSheet View'
spreadSheetView1 = CreateView('SpreadSheetView')
spreadSheetView1.ColumnToSort = ''
spreadSheetView1.BlockSize = 1024L
# uncomment following to set a specific view size
# spreadSheetView1.ViewSize = [400, 400]

# get layout
layout1 = GetLayout()

# place view in the layout
layout1.AssignView(0, spreadSheetView1)

# show data in view
fieldDisplay = Show(field, spreadSheetView1)

# export view
ExportView('field.csv', view=spreadSheetView1, FilterColumnsByVisibility=1)
