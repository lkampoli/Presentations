#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'Legacy VTK Reader'
k_Hill_wallvtk = LegacyVTKReader(FileNames=['/home/unimelb.edu.au/lcampoli/UoM/Testcases/45_degree_Yaw/VT_NASA_BeVERLI_3D_Hill_Baseline_RANS_k_Omega_SST_Mes_Lev0_64mm_MeshB_run2_SELECTED/postProcessing/surfaces/120000/k_Hill_wall.vtk'])

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [1544, 802]

# show data in view
k_Hill_wallvtkDisplay = Show(k_Hill_wallvtk, renderView1)
# trace defaults for the display properties.
k_Hill_wallvtkDisplay.Representation = 'Surface'
k_Hill_wallvtkDisplay.ColorArrayName = [None, '']
k_Hill_wallvtkDisplay.OSPRayScaleArray = 'k'
k_Hill_wallvtkDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
k_Hill_wallvtkDisplay.SelectOrientationVectors = 'None'
k_Hill_wallvtkDisplay.ScaleFactor = 0.625
k_Hill_wallvtkDisplay.SelectScaleArray = 'None'
k_Hill_wallvtkDisplay.GlyphType = 'Arrow'
k_Hill_wallvtkDisplay.GlyphTableIndexArray = 'None'
k_Hill_wallvtkDisplay.DataAxesGrid = 'GridAxesRepresentation'
k_Hill_wallvtkDisplay.PolarAxes = 'PolarAxesRepresentation'

# reset view to fit data
renderView1.ResetCamera()

# update the view to ensure updated data information
renderView1.Update()

# set scalar coloring
ColorBy(k_Hill_wallvtkDisplay, ('POINTS', 'k'))

# rescale color and/or opacity maps used to include current data range
k_Hill_wallvtkDisplay.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
k_Hill_wallvtkDisplay.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'k'
kLUT = GetColorTransferFunction('k')

# reset view to fit data
renderView1.ResetCamera()

# reset view to fit data
renderView1.ResetCamera()

# reset view to fit data
renderView1.ResetCamera()

# reset view to fit data
renderView1.ResetCamera()

# reset view to fit data
renderView1.ResetCamera()

# reset view to fit data
renderView1.ResetCamera()

# reset view to fit data
renderView1.ResetCamera()

# reset view to fit data
renderView1.ResetCamera()

# Properties modified on k_Hill_wallvtkDisplay.DataAxesGrid
k_Hill_wallvtkDisplay.DataAxesGrid.GridAxesVisibility = 1

# Properties modified on k_Hill_wallvtkDisplay.DataAxesGrid
k_Hill_wallvtkDisplay.DataAxesGrid.GridAxesVisibility = 0

# reset view to fit data
renderView1.ResetCamera()

# reset view to fit data
renderView1.ResetCamera()

# reset view to fit data
renderView1.ResetCamera()

# reset view to fit data
renderView1.ResetCamera()

# reset view to fit data
renderView1.ResetCamera()

# reset view to fit data
renderView1.ResetCamera()

# create a new 'Legacy VTK Reader'
r_Hill_wallvtk = LegacyVTKReader(FileNames=['/home/unimelb.edu.au/lcampoli/UoM/Testcases/45_degree_Yaw/VT_NASA_BeVERLI_3D_Hill_Baseline_RANS_k_Omega_SST_Mes_Lev0_64mm_MeshB_run2_SELECTED/postProcessing/surfaces/120000/R_Hill_wall.vtk'])

# show data in view
r_Hill_wallvtkDisplay = Show(r_Hill_wallvtk, renderView1)
# trace defaults for the display properties.
r_Hill_wallvtkDisplay.Representation = 'Surface'
r_Hill_wallvtkDisplay.ColorArrayName = [None, '']
r_Hill_wallvtkDisplay.OSPRayScaleArray = 'R'
r_Hill_wallvtkDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
r_Hill_wallvtkDisplay.SelectOrientationVectors = 'None'
r_Hill_wallvtkDisplay.ScaleFactor = 0.625
r_Hill_wallvtkDisplay.SelectScaleArray = 'None'
r_Hill_wallvtkDisplay.GlyphType = 'Arrow'
r_Hill_wallvtkDisplay.GlyphTableIndexArray = 'None'
r_Hill_wallvtkDisplay.DataAxesGrid = 'GridAxesRepresentation'
r_Hill_wallvtkDisplay.PolarAxes = 'PolarAxesRepresentation'

# update the view to ensure updated data information
renderView1.Update()

# set scalar coloring
ColorBy(r_Hill_wallvtkDisplay, ('POINTS', 'R', 'Magnitude'))

# rescale color and/or opacity maps used to include current data range
r_Hill_wallvtkDisplay.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
r_Hill_wallvtkDisplay.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'R'
rLUT = GetColorTransferFunction('R')

# hide data in view
Hide(k_Hill_wallvtk, renderView1)

# Rescale transfer function
rLUT.RescaleTransferFunction(6.05641332641e-08, 0.00143673620187)

# get opacity transfer function/opacity map for 'R'
rPWF = GetOpacityTransferFunction('R')

# Rescale transfer function
rPWF.RescaleTransferFunction(6.05641332641e-08, 0.00143673620187)

# set active source
SetActiveSource(k_Hill_wallvtk)

# show data in view
k_Hill_wallvtkDisplay = Show(k_Hill_wallvtk, renderView1)

# show color bar/color legend
k_Hill_wallvtkDisplay.SetScalarBarVisibility(renderView1, True)

# hide data in view
Hide(r_Hill_wallvtk, renderView1)

# set active source
SetActiveSource(r_Hill_wallvtk)

# show data in view
r_Hill_wallvtkDisplay = Show(r_Hill_wallvtk, renderView1)

# show color bar/color legend
r_Hill_wallvtkDisplay.SetScalarBarVisibility(renderView1, True)

# create a new 'Legacy VTK Reader'
u_Hill_wallvtk = LegacyVTKReader(FileNames=['/home/unimelb.edu.au/lcampoli/UoM/Testcases/45_degree_Yaw/VT_NASA_BeVERLI_3D_Hill_Baseline_RANS_k_Omega_SST_Mes_Lev0_64mm_MeshB_run2_SELECTED/postProcessing/surfaces/120000/U_Hill_wall.vtk'])

# show data in view
u_Hill_wallvtkDisplay = Show(u_Hill_wallvtk, renderView1)
# trace defaults for the display properties.
u_Hill_wallvtkDisplay.Representation = 'Surface'
u_Hill_wallvtkDisplay.ColorArrayName = [None, '']
u_Hill_wallvtkDisplay.OSPRayScaleArray = 'U'
u_Hill_wallvtkDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
u_Hill_wallvtkDisplay.SelectOrientationVectors = 'None'
u_Hill_wallvtkDisplay.ScaleFactor = 0.625
u_Hill_wallvtkDisplay.SelectScaleArray = 'None'
u_Hill_wallvtkDisplay.GlyphType = 'Arrow'
u_Hill_wallvtkDisplay.GlyphTableIndexArray = 'None'
u_Hill_wallvtkDisplay.DataAxesGrid = 'GridAxesRepresentation'
u_Hill_wallvtkDisplay.PolarAxes = 'PolarAxesRepresentation'

# update the view to ensure updated data information
renderView1.Update()

# Rescale transfer function
rLUT.RescaleTransferFunction(6.05641332641e-08, 0.00248850009871)

# Rescale transfer function
rPWF.RescaleTransferFunction(6.05641332641e-08, 0.00248850009871)

# set scalar coloring
ColorBy(u_Hill_wallvtkDisplay, ('POINTS', 'U', 'Magnitude'))

# rescale color and/or opacity maps used to include current data range
u_Hill_wallvtkDisplay.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
u_Hill_wallvtkDisplay.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'U'
uLUT = GetColorTransferFunction('U')

# hide data in view
Hide(r_Hill_wallvtk, renderView1)

# hide data in view
Hide(k_Hill_wallvtk, renderView1)

# reset view to fit data
renderView1.ResetCamera()

# reset view to fit data
renderView1.ResetCamera()

# reset view to fit data
renderView1.ResetCamera()

# reset view to fit data
renderView1.ResetCamera()

# reset view to fit data
renderView1.ResetCamera()

# reset view to fit data
renderView1.ResetCamera()

# rescale color and/or opacity maps used to exactly fit the current data range
u_Hill_wallvtkDisplay.RescaleTransferFunctionToDataRange(False, True)

# Rescale transfer function
uLUT.RescaleTransferFunction(0.0, 1.17578e-38)

# get opacity transfer function/opacity map for 'U'
uPWF = GetOpacityTransferFunction('U')

# Rescale transfer function
uPWF.RescaleTransferFunction(0.0, 1.17578e-38)

# Rescale transfer function
uLUT.RescaleTransferFunction(0.0, 1.17578133675e-38)

# Rescale transfer function
uPWF.RescaleTransferFunction(0.0, 1.17578133675e-38)

# create a new 'Legacy VTK Reader'
wallGradU_Hill_wallvtk = LegacyVTKReader(FileNames=['/home/unimelb.edu.au/lcampoli/UoM/Testcases/45_degree_Yaw/VT_NASA_BeVERLI_3D_Hill_Baseline_RANS_k_Omega_SST_Mes_Lev0_64mm_MeshB_run2_SELECTED/postProcessing/surfaces/120000/wallGradU_Hill_wall.vtk'])

# show data in view
wallGradU_Hill_wallvtkDisplay = Show(wallGradU_Hill_wallvtk, renderView1)
# trace defaults for the display properties.
wallGradU_Hill_wallvtkDisplay.Representation = 'Surface'
wallGradU_Hill_wallvtkDisplay.ColorArrayName = [None, '']
wallGradU_Hill_wallvtkDisplay.OSPRayScaleArray = 'wallGradU'
wallGradU_Hill_wallvtkDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
wallGradU_Hill_wallvtkDisplay.SelectOrientationVectors = 'None'
wallGradU_Hill_wallvtkDisplay.ScaleFactor = 0.625
wallGradU_Hill_wallvtkDisplay.SelectScaleArray = 'None'
wallGradU_Hill_wallvtkDisplay.GlyphType = 'Arrow'
wallGradU_Hill_wallvtkDisplay.GlyphTableIndexArray = 'None'
wallGradU_Hill_wallvtkDisplay.DataAxesGrid = 'GridAxesRepresentation'
wallGradU_Hill_wallvtkDisplay.PolarAxes = 'PolarAxesRepresentation'

# update the view to ensure updated data information
renderView1.Update()

# set scalar coloring
ColorBy(wallGradU_Hill_wallvtkDisplay, ('POINTS', 'wallGradU', 'Magnitude'))

# rescale color and/or opacity maps used to include current data range
wallGradU_Hill_wallvtkDisplay.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
wallGradU_Hill_wallvtkDisplay.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'wallGradU'
wallGradULUT = GetColorTransferFunction('wallGradU')

# create a new 'Legacy VTK Reader'
wallShearStress_Hill_wallvtk = LegacyVTKReader(FileNames=['/home/unimelb.edu.au/lcampoli/UoM/Testcases/45_degree_Yaw/VT_NASA_BeVERLI_3D_Hill_Baseline_RANS_k_Omega_SST_Mes_Lev0_64mm_MeshB_run2_SELECTED/postProcessing/surfaces/120000/wallShearStress_Hill_wall.vtk'])

# show data in view
wallShearStress_Hill_wallvtkDisplay = Show(wallShearStress_Hill_wallvtk, renderView1)
# trace defaults for the display properties.
wallShearStress_Hill_wallvtkDisplay.Representation = 'Surface'
wallShearStress_Hill_wallvtkDisplay.ColorArrayName = [None, '']
wallShearStress_Hill_wallvtkDisplay.OSPRayScaleArray = 'wallShearStress'
wallShearStress_Hill_wallvtkDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
wallShearStress_Hill_wallvtkDisplay.SelectOrientationVectors = 'None'
wallShearStress_Hill_wallvtkDisplay.ScaleFactor = 0.625
wallShearStress_Hill_wallvtkDisplay.SelectScaleArray = 'None'
wallShearStress_Hill_wallvtkDisplay.GlyphType = 'Arrow'
wallShearStress_Hill_wallvtkDisplay.GlyphTableIndexArray = 'None'
wallShearStress_Hill_wallvtkDisplay.DataAxesGrid = 'GridAxesRepresentation'
wallShearStress_Hill_wallvtkDisplay.PolarAxes = 'PolarAxesRepresentation'

# update the view to ensure updated data information
renderView1.Update()

# set scalar coloring
ColorBy(wallShearStress_Hill_wallvtkDisplay, ('POINTS', 'wallShearStress', 'Magnitude'))

# rescale color and/or opacity maps used to include current data range
wallShearStress_Hill_wallvtkDisplay.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
wallShearStress_Hill_wallvtkDisplay.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'wallShearStress'
wallShearStressLUT = GetColorTransferFunction('wallShearStress')

# hide data in view
Hide(wallShearStress_Hill_wallvtk, renderView1)

# hide data in view
Hide(wallGradU_Hill_wallvtk, renderView1)

# hide data in view
Hide(u_Hill_wallvtk, renderView1)

# set active source
SetActiveSource(k_Hill_wallvtk)

# show data in view
k_Hill_wallvtkDisplay = Show(k_Hill_wallvtk, renderView1)

# show color bar/color legend
k_Hill_wallvtkDisplay.SetScalarBarVisibility(renderView1, True)

# reset view to fit data
renderView1.ResetCamera()

# set active source
SetActiveSource(r_Hill_wallvtk)

# show data in view
r_Hill_wallvtkDisplay = Show(r_Hill_wallvtk, renderView1)

# show color bar/color legend
r_Hill_wallvtkDisplay.SetScalarBarVisibility(renderView1, True)

# set active source
SetActiveSource(u_Hill_wallvtk)

# show data in view
u_Hill_wallvtkDisplay = Show(u_Hill_wallvtk, renderView1)

# show color bar/color legend
u_Hill_wallvtkDisplay.SetScalarBarVisibility(renderView1, True)

# set active source
SetActiveSource(wallGradU_Hill_wallvtk)

# show data in view
wallGradU_Hill_wallvtkDisplay = Show(wallGradU_Hill_wallvtk, renderView1)

# show color bar/color legend
wallGradU_Hill_wallvtkDisplay.SetScalarBarVisibility(renderView1, True)

# set active source
SetActiveSource(wallShearStress_Hill_wallvtk)

# show data in view
wallShearStress_Hill_wallvtkDisplay = Show(wallShearStress_Hill_wallvtk, renderView1)

# show color bar/color legend
wallShearStress_Hill_wallvtkDisplay.SetScalarBarVisibility(renderView1, True)

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
wallShearStress_Hill_wallvtkDisplay = Show(wallShearStress_Hill_wallvtk, spreadSheetView1)

# export view
ExportView('/home/unimelb.edu.au/lcampoli/UoM/Testcases/45_degree_Yaw/VT_NASA_BeVERLI_3D_Hill_Baseline_RANS_k_Omega_SST_Mes_Lev0_64mm_MeshB_run2_SELECTED/postProcessing/surfaces/120000/wallShearStress.csv', view=spreadSheetView1, FilterColumnsByVisibility=1)

# create a new 'Legacy VTK Reader'
wallShearStress_X_0_Planevtk = LegacyVTKReader(FileNames=['/home/unimelb.edu.au/lcampoli/UoM/Testcases/45_degree_Yaw/VT_NASA_BeVERLI_3D_Hill_Baseline_RANS_k_Omega_SST_Mes_Lev0_64mm_MeshB_run2_SELECTED/postProcessing/surfaces/120000/wallShearStress_X_0_Plane.vtk'])

# export view
ExportView('/home/unimelb.edu.au/lcampoli/UoM/Testcases/45_degree_Yaw/VT_NASA_BeVERLI_3D_Hill_Baseline_RANS_k_Omega_SST_Mes_Lev0_64mm_MeshB_run2_SELECTED/postProcessing/surfaces/120000/wallShearStressX.csv', view=spreadSheetView1, FilterColumnsByVisibility=1)

# destroy spreadSheetView1
Delete(spreadSheetView1)
del spreadSheetView1

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [1544, 802]
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.StereoType = 0
renderView1.Background = [0.32, 0.34, 0.43]

# place view in the layout
layout1.AssignView(0, renderView1)

# set active source
SetActiveSource(wallShearStress_X_0_Planevtk)

# show data in view
wallShearStress_X_0_PlanevtkDisplay = Show(wallShearStress_X_0_Planevtk, renderView1)
# trace defaults for the display properties.
wallShearStress_X_0_PlanevtkDisplay.Representation = 'Surface'
wallShearStress_X_0_PlanevtkDisplay.ColorArrayName = [None, '']
wallShearStress_X_0_PlanevtkDisplay.OSPRayScaleArray = 'wallShearStress'
wallShearStress_X_0_PlanevtkDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
wallShearStress_X_0_PlanevtkDisplay.SelectOrientationVectors = 'None'
wallShearStress_X_0_PlanevtkDisplay.ScaleFactor = 0.9789028167724609
wallShearStress_X_0_PlanevtkDisplay.SelectScaleArray = 'None'
wallShearStress_X_0_PlanevtkDisplay.GlyphType = 'Arrow'
wallShearStress_X_0_PlanevtkDisplay.GlyphTableIndexArray = 'None'
wallShearStress_X_0_PlanevtkDisplay.DataAxesGrid = 'GridAxesRepresentation'
wallShearStress_X_0_PlanevtkDisplay.PolarAxes = 'PolarAxesRepresentation'

# reset view to fit data
renderView1.ResetCamera()

# show data in view
wallShearStress_X_0_PlanevtkDisplay = Show(wallShearStress_X_0_Planevtk, renderView1)

# reset view to fit data
renderView1.ResetCamera()

# update the view to ensure updated data information
renderView1.Update()

# reset view to fit data
renderView1.ResetCamera()

# reset view to fit data
renderView1.ResetCamera()

# set scalar coloring
ColorBy(wallShearStress_X_0_PlanevtkDisplay, ('POINTS', 'wallShearStress', 'Magnitude'))

# rescale color and/or opacity maps used to include current data range
wallShearStress_X_0_PlanevtkDisplay.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
wallShearStress_X_0_PlanevtkDisplay.SetScalarBarVisibility(renderView1, True)

# Rescale transfer function
wallShearStressLUT.RescaleTransferFunction(0.0, 0.00466002258136)

# get opacity transfer function/opacity map for 'wallShearStress'
wallShearStressPWF = GetOpacityTransferFunction('wallShearStress')

# Rescale transfer function
wallShearStressPWF.RescaleTransferFunction(0.0, 0.00466002258136)

# reset view to fit data
renderView1.ResetCamera()

# reset view to fit data
renderView1.ResetCamera()

# reset view to fit data
renderView1.ResetCamera()

# reset view to fit data
renderView1.ResetCamera()

# reset view to fit data
renderView1.ResetCamera()

# reset view to fit data
renderView1.ResetCamera()

# reset view to fit data
renderView1.ResetCamera()

# current camera placement for renderView1
renderView1.CameraPosition = [26.74412231094651, 4.894514083862305, 0.0]
renderView1.CameraFocalPoint = [-2.6183244237508517e-19, 4.894514083862305, 0.0]
renderView1.CameraViewUp = [0.0, 0.0, 1.0]
renderView1.CameraParallelScale = 6.921888198624195

# save screenshot
SaveScreenshot('/home/unimelb.edu.au/lcampoli/UoM/Testcases/45_degree_Yaw/VT_NASA_BeVERLI_3D_Hill_Baseline_RANS_k_Omega_SST_Mes_Lev0_64mm_MeshB_run2_SELECTED/postProcessing/surfaces/120000/pii.png', renderView1, ImageResolution=[1544, 802])

#### saving camera placements for all active views

# current camera placement for renderView1
renderView1.CameraPosition = [26.74412231094651, 4.894514083862305, 0.0]
renderView1.CameraFocalPoint = [-2.6183244237508517e-19, 4.894514083862305, 0.0]
renderView1.CameraViewUp = [0.0, 0.0, 1.0]
renderView1.CameraParallelScale = 6.921888198624195

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).