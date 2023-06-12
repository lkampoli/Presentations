#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# find source
frozen_sq_cyl_U_k_bij_hifi_Solve_only_Omega_with_Rterm_NEW_ProdLimiter2_Rtermclipped_1_Selectedfoam = FindSource('frozen_sq_cyl_U_k_bij_hifi_Solve_only_Omega_with_Rterm_NEW_ProdLimiter2_Rtermclipped_1_Selected.foam')

# create a new 'Slice'
slice1 = Slice(Input=frozen_sq_cyl_U_k_bij_hifi_Solve_only_Omega_with_Rterm_NEW_ProdLimiter2_Rtermclipped_1_Selectedfoam)
slice1.SliceType = 'Plane'
slice1.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice1.SliceType.Origin = [4.0, 0.0, 6.0]

# find source
a14 = FindSource('14')

# find source
a4 = FindSource('4')

# find source
a2 = FindSource('2')

# find source
a8 = FindSource('8')

# find source
a1 = FindSource('1')

# find source
a01 = FindSource('0.1')

# find source
a6 = FindSource('6')

# find source
a10 = FindSource('10')

# find source
a12 = FindSource('12')

# find source
a16 = FindSource('16')

# Properties modified on slice1.SliceType
slice1.SliceType.Origin = [-4.0, 0.0, 6.0]

# Properties modified on slice1.SliceType
slice1.SliceType.Origin = [-4.0, 0.0, 6.0]

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [1484, 802]

# get color transfer function/color map for 'k'
kLUT = GetColorTransferFunction('k')

# show data in view
slice1Display = Show(slice1, renderView1)
# trace defaults for the display properties.
slice1Display.Representation = 'Surface'
slice1Display.ColorArrayName = ['POINTS', 'k']
slice1Display.LookupTable = kLUT
slice1Display.OSPRayScaleArray = 'k'
slice1Display.OSPRayScaleFunction = 'PiecewiseFunction'
slice1Display.SelectOrientationVectors = 'None'
slice1Display.ScaleFactor = 1.6
slice1Display.SelectScaleArray = 'None'
slice1Display.GlyphType = 'Arrow'
slice1Display.GlyphTableIndexArray = 'None'
slice1Display.DataAxesGrid = 'GridAxesRepresentation'
slice1Display.PolarAxes = 'PolarAxesRepresentation'

# hide data in view
Hide(frozen_sq_cyl_U_k_bij_hifi_Solve_only_Omega_with_Rterm_NEW_ProdLimiter2_Rtermclipped_1_Selectedfoam, renderView1)

# show color bar/color legend
slice1Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# rename source object
RenameSource('-4', slice1)

# hide data in view
Hide(a16, renderView1)

# hide data in view
Hide(a14, renderView1)

# hide data in view
Hide(a12, renderView1)

# hide data in view
Hide(a10, renderView1)

# hide data in view
Hide(a6, renderView1)

# hide data in view
Hide(a2, renderView1)

# hide data in view
Hide(a01, renderView1)

# hide data in view
Hide(a1, renderView1)

# hide data in view
Hide(a8, renderView1)

# hide data in view
Hide(a4, renderView1)

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
slice1Display = Show(slice1, spreadSheetView1)

# export view
ExportView('/home/unimelb.edu.au/lcampoli/UoM/Testcases/frozen_sq_cyl_U_k_bij_hifi_Solve_only_Omega_with_Rterm_NEW_ProdLimiter2_Rtermclipped_1_Selected/frozen_sq_cyl_U_k_bij_hifi_Solve_only_Omega_with_Rterm_NEW_ProdLimiter2_Rtermclipped_1_Selected/k-4.csv', view=spreadSheetView1, FilterColumnsByVisibility=1)

# show data in view
a4Display = Show(a4, spreadSheetView1)

# export view
ExportView('/home/unimelb.edu.au/lcampoli/UoM/Testcases/frozen_sq_cyl_U_k_bij_hifi_Solve_only_Omega_with_Rterm_NEW_ProdLimiter2_Rtermclipped_1_Selected/frozen_sq_cyl_U_k_bij_hifi_Solve_only_Omega_with_Rterm_NEW_ProdLimiter2_Rtermclipped_1_Selected/k+4.csv', view=spreadSheetView1, FilterColumnsByVisibility=1)

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).