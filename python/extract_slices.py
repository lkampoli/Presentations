#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'OpenFOAMReader'
frozen_sq_cyl_U_k_bij_hifi_Solve_only_Omega_with_Rterm_NEW_ProdLimiter2_Rtermclipped_1_Selectedfoam = OpenFOAMReader(FileName='/home/unimelb.edu.au/lcampoli/UoM/Testcases/frozen_sq_cyl_U_k_bij_hifi_Solve_only_Omega_with_Rterm_NEW_ProdLimiter2_Rtermclipped_1_Selected/frozen_sq_cyl_U_k_bij_hifi_Solve_only_Omega_with_Rterm_NEW_ProdLimiter2_Rtermclipped_1_Selected/frozen_sq_cyl_U_k_bij_hifi_Solve_only_Omega_with_Rterm_NEW_ProdLimiter2_Rtermclipped_1_Selected.foam')
frozen_sq_cyl_U_k_bij_hifi_Solve_only_Omega_with_Rterm_NEW_ProdLimiter2_Rtermclipped_1_Selectedfoam.MeshRegions = ['internalMesh']
frozen_sq_cyl_U_k_bij_hifi_Solve_only_Omega_with_Rterm_NEW_ProdLimiter2_Rtermclipped_1_Selectedfoam.CellArrays = ['I1', 'I2', 'I3', 'I4', 'I5', 'Omegaij', 'R_res', 'R_sgs', 'R_tot', 'Rall', 'Rterm', 'RtermNorm', 'Sij', 'T1', 'T10', 'T2', 'T3', 'T4', 'T5', 'T6', 'T7', 'T8', 'T9', 'U', 'aijx', 'bij', 'gradU', 'k', 'k_res', 'nut', 'omega', 'xswtch']

# Properties modified on frozen_sq_cyl_U_k_bij_hifi_Solve_only_Omega_with_Rterm_NEW_ProdLimiter2_Rtermclipped_1_Selectedfoam
frozen_sq_cyl_U_k_bij_hifi_Solve_only_Omega_with_Rterm_NEW_ProdLimiter2_Rtermclipped_1_Selectedfoam.ScalarSize = '32-bit (SP)'
frozen_sq_cyl_U_k_bij_hifi_Solve_only_Omega_with_Rterm_NEW_ProdLimiter2_Rtermclipped_1_Selectedfoam.CellArrays = ['k']

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [1484, 802]

# show data in view
frozen_sq_cyl_U_k_bij_hifi_Solve_only_Omega_with_Rterm_NEW_ProdLimiter2_Rtermclipped_1_SelectedfoamDisplay = Show(frozen_sq_cyl_U_k_bij_hifi_Solve_only_Omega_with_Rterm_NEW_ProdLimiter2_Rtermclipped_1_Selectedfoam, renderView1)
# trace defaults for the display properties.
frozen_sq_cyl_U_k_bij_hifi_Solve_only_Omega_with_Rterm_NEW_ProdLimiter2_Rtermclipped_1_SelectedfoamDisplay.Representation = 'Surface'
frozen_sq_cyl_U_k_bij_hifi_Solve_only_Omega_with_Rterm_NEW_ProdLimiter2_Rtermclipped_1_SelectedfoamDisplay.ColorArrayName = [None, '']
frozen_sq_cyl_U_k_bij_hifi_Solve_only_Omega_with_Rterm_NEW_ProdLimiter2_Rtermclipped_1_SelectedfoamDisplay.OSPRayScaleArray = 'k'
frozen_sq_cyl_U_k_bij_hifi_Solve_only_Omega_with_Rterm_NEW_ProdLimiter2_Rtermclipped_1_SelectedfoamDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
frozen_sq_cyl_U_k_bij_hifi_Solve_only_Omega_with_Rterm_NEW_ProdLimiter2_Rtermclipped_1_SelectedfoamDisplay.SelectOrientationVectors = 'None'
frozen_sq_cyl_U_k_bij_hifi_Solve_only_Omega_with_Rterm_NEW_ProdLimiter2_Rtermclipped_1_SelectedfoamDisplay.ScaleFactor = 2.4000000000000004
frozen_sq_cyl_U_k_bij_hifi_Solve_only_Omega_with_Rterm_NEW_ProdLimiter2_Rtermclipped_1_SelectedfoamDisplay.SelectScaleArray = 'None'
frozen_sq_cyl_U_k_bij_hifi_Solve_only_Omega_with_Rterm_NEW_ProdLimiter2_Rtermclipped_1_SelectedfoamDisplay.GlyphType = 'Arrow'
frozen_sq_cyl_U_k_bij_hifi_Solve_only_Omega_with_Rterm_NEW_ProdLimiter2_Rtermclipped_1_SelectedfoamDisplay.GlyphTableIndexArray = 'None'
frozen_sq_cyl_U_k_bij_hifi_Solve_only_Omega_with_Rterm_NEW_ProdLimiter2_Rtermclipped_1_SelectedfoamDisplay.DataAxesGrid = 'GridAxesRepresentation'
frozen_sq_cyl_U_k_bij_hifi_Solve_only_Omega_with_Rterm_NEW_ProdLimiter2_Rtermclipped_1_SelectedfoamDisplay.PolarAxes = 'PolarAxesRepresentation'
frozen_sq_cyl_U_k_bij_hifi_Solve_only_Omega_with_Rterm_NEW_ProdLimiter2_Rtermclipped_1_SelectedfoamDisplay.ScalarOpacityUnitDistance = 0.11999587205850736

# reset view to fit data
renderView1.ResetCamera()

# update the view to ensure updated data information
renderView1.Update()

# set scalar coloring
ColorBy(frozen_sq_cyl_U_k_bij_hifi_Solve_only_Omega_with_Rterm_NEW_ProdLimiter2_Rtermclipped_1_SelectedfoamDisplay, ('POINTS', 'k'))

# rescale color and/or opacity maps used to include current data range
frozen_sq_cyl_U_k_bij_hifi_Solve_only_Omega_with_Rterm_NEW_ProdLimiter2_Rtermclipped_1_SelectedfoamDisplay.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
frozen_sq_cyl_U_k_bij_hifi_Solve_only_Omega_with_Rterm_NEW_ProdLimiter2_Rtermclipped_1_SelectedfoamDisplay.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'k'
kLUT = GetColorTransferFunction('k')

# create a new 'Slice'
slice1 = Slice(Input=frozen_sq_cyl_U_k_bij_hifi_Solve_only_Omega_with_Rterm_NEW_ProdLimiter2_Rtermclipped_1_Selectedfoam)
slice1.SliceType = 'Plane'
slice1.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice1.SliceType.Origin = [4.0, 0.0, 6.0]

# rename source object
RenameSource('x=0', slice1)

# Properties modified on slice1.SliceType
slice1.SliceType.Origin = [0.0, 0.0, 6.0]

# Properties modified on slice1.SliceType
slice1.SliceType.Origin = [0.0, 0.0, 6.0]

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

# set active source
SetActiveSource(frozen_sq_cyl_U_k_bij_hifi_Solve_only_Omega_with_Rterm_NEW_ProdLimiter2_Rtermclipped_1_Selectedfoam)

# create a new 'Slice'
slice1_1 = Slice(Input=frozen_sq_cyl_U_k_bij_hifi_Solve_only_Omega_with_Rterm_NEW_ProdLimiter2_Rtermclipped_1_Selectedfoam)
slice1_1.SliceType = 'Plane'
slice1_1.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice1_1.SliceType.Origin = [4.0, 0.0, 6.0]

# rename source object
RenameSource('x=1', slice1_1)

# Properties modified on slice1_1.SliceType
slice1_1.SliceType.Origin = [1.0, 0.0, 6.0]

# Properties modified on slice1_1.SliceType
slice1_1.SliceType.Origin = [1.0, 0.0, 6.0]

# show data in view
slice1_1Display = Show(slice1_1, renderView1)
# trace defaults for the display properties.
slice1_1Display.Representation = 'Surface'
slice1_1Display.ColorArrayName = ['POINTS', 'k']
slice1_1Display.LookupTable = kLUT
slice1_1Display.OSPRayScaleArray = 'k'
slice1_1Display.OSPRayScaleFunction = 'PiecewiseFunction'
slice1_1Display.SelectOrientationVectors = 'None'
slice1_1Display.ScaleFactor = 1.6
slice1_1Display.SelectScaleArray = 'None'
slice1_1Display.GlyphType = 'Arrow'
slice1_1Display.GlyphTableIndexArray = 'None'
slice1_1Display.DataAxesGrid = 'GridAxesRepresentation'
slice1_1Display.PolarAxes = 'PolarAxesRepresentation'

# hide data in view
Hide(frozen_sq_cyl_U_k_bij_hifi_Solve_only_Omega_with_Rterm_NEW_ProdLimiter2_Rtermclipped_1_Selectedfoam, renderView1)

# show color bar/color legend
slice1_1Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# set active source
SetActiveSource(frozen_sq_cyl_U_k_bij_hifi_Solve_only_Omega_with_Rterm_NEW_ProdLimiter2_Rtermclipped_1_Selectedfoam)

# create a new 'Slice'
slice1_2 = Slice(Input=frozen_sq_cyl_U_k_bij_hifi_Solve_only_Omega_with_Rterm_NEW_ProdLimiter2_Rtermclipped_1_Selectedfoam)
slice1_2.SliceType = 'Plane'
slice1_2.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice1_2.SliceType.Origin = [4.0, 0.0, 6.0]

# rename source object
RenameSource('x=2', slice1_2)

# Properties modified on slice1_2.SliceType
slice1_2.SliceType.Origin = [2.0, 0.0, 6.0]

# Properties modified on slice1_2.SliceType
slice1_2.SliceType.Origin = [2.0, 0.0, 6.0]

# show data in view
slice1_2Display = Show(slice1_2, renderView1)
# trace defaults for the display properties.
slice1_2Display.Representation = 'Surface'
slice1_2Display.ColorArrayName = ['POINTS', 'k']
slice1_2Display.LookupTable = kLUT
slice1_2Display.OSPRayScaleArray = 'k'
slice1_2Display.OSPRayScaleFunction = 'PiecewiseFunction'
slice1_2Display.SelectOrientationVectors = 'None'
slice1_2Display.ScaleFactor = 1.6
slice1_2Display.SelectScaleArray = 'None'
slice1_2Display.GlyphType = 'Arrow'
slice1_2Display.GlyphTableIndexArray = 'None'
slice1_2Display.DataAxesGrid = 'GridAxesRepresentation'
slice1_2Display.PolarAxes = 'PolarAxesRepresentation'

# hide data in view
Hide(frozen_sq_cyl_U_k_bij_hifi_Solve_only_Omega_with_Rterm_NEW_ProdLimiter2_Rtermclipped_1_Selectedfoam, renderView1)

# show color bar/color legend
slice1_2Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# set active source
SetActiveSource(frozen_sq_cyl_U_k_bij_hifi_Solve_only_Omega_with_Rterm_NEW_ProdLimiter2_Rtermclipped_1_Selectedfoam)

# create a new 'Slice'
slice1_3 = Slice(Input=frozen_sq_cyl_U_k_bij_hifi_Solve_only_Omega_with_Rterm_NEW_ProdLimiter2_Rtermclipped_1_Selectedfoam)
slice1_3.SliceType = 'Plane'
slice1_3.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice1_3.SliceType.Origin = [4.0, 0.0, 6.0]

# rename source object
RenameSource('x=3', slice1_3)

# Properties modified on slice1_3.SliceType
slice1_3.SliceType.Origin = [3.0, 0.0, 6.0]

# Properties modified on slice1_3.SliceType
slice1_3.SliceType.Origin = [3.0, 0.0, 6.0]

# show data in view
slice1_3Display = Show(slice1_3, renderView1)
# trace defaults for the display properties.
slice1_3Display.Representation = 'Surface'
slice1_3Display.ColorArrayName = ['POINTS', 'k']
slice1_3Display.LookupTable = kLUT
slice1_3Display.OSPRayScaleArray = 'k'
slice1_3Display.OSPRayScaleFunction = 'PiecewiseFunction'
slice1_3Display.SelectOrientationVectors = 'None'
slice1_3Display.ScaleFactor = 1.6
slice1_3Display.SelectScaleArray = 'None'
slice1_3Display.GlyphType = 'Arrow'
slice1_3Display.GlyphTableIndexArray = 'None'
slice1_3Display.DataAxesGrid = 'GridAxesRepresentation'
slice1_3Display.PolarAxes = 'PolarAxesRepresentation'

# hide data in view
Hide(frozen_sq_cyl_U_k_bij_hifi_Solve_only_Omega_with_Rterm_NEW_ProdLimiter2_Rtermclipped_1_Selectedfoam, renderView1)

# show color bar/color legend
slice1_3Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# set active source
SetActiveSource(frozen_sq_cyl_U_k_bij_hifi_Solve_only_Omega_with_Rterm_NEW_ProdLimiter2_Rtermclipped_1_Selectedfoam)

# create a new 'Slice'
slice1_4 = Slice(Input=frozen_sq_cyl_U_k_bij_hifi_Solve_only_Omega_with_Rterm_NEW_ProdLimiter2_Rtermclipped_1_Selectedfoam)
slice1_4.SliceType = 'Plane'
slice1_4.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice1_4.SliceType.Origin = [4.0, 0.0, 6.0]

# rename source object
RenameSource('x=4', slice1_4)

# show data in view
slice1_4Display = Show(slice1_4, renderView1)
# trace defaults for the display properties.
slice1_4Display.Representation = 'Surface'
slice1_4Display.ColorArrayName = ['POINTS', 'k']
slice1_4Display.LookupTable = kLUT
slice1_4Display.OSPRayScaleArray = 'k'
slice1_4Display.OSPRayScaleFunction = 'PiecewiseFunction'
slice1_4Display.SelectOrientationVectors = 'None'
slice1_4Display.ScaleFactor = 1.6
slice1_4Display.SelectScaleArray = 'None'
slice1_4Display.GlyphType = 'Arrow'
slice1_4Display.GlyphTableIndexArray = 'None'
slice1_4Display.DataAxesGrid = 'GridAxesRepresentation'
slice1_4Display.PolarAxes = 'PolarAxesRepresentation'

# hide data in view
Hide(frozen_sq_cyl_U_k_bij_hifi_Solve_only_Omega_with_Rterm_NEW_ProdLimiter2_Rtermclipped_1_Selectedfoam, renderView1)

# show color bar/color legend
slice1_4Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# set active source
SetActiveSource(frozen_sq_cyl_U_k_bij_hifi_Solve_only_Omega_with_Rterm_NEW_ProdLimiter2_Rtermclipped_1_Selectedfoam)

# create a new 'Slice'
slice1_5 = Slice(Input=frozen_sq_cyl_U_k_bij_hifi_Solve_only_Omega_with_Rterm_NEW_ProdLimiter2_Rtermclipped_1_Selectedfoam)
slice1_5.SliceType = 'Plane'
slice1_5.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice1_5.SliceType.Origin = [4.0, 0.0, 6.0]

# rename source object
RenameSource('x=5', slice1_5)

# Properties modified on slice1_5.SliceType
slice1_5.SliceType.Origin = [5.0, 0.0, 6.0]

# Properties modified on slice1_5.SliceType
slice1_5.SliceType.Origin = [5.0, 0.0, 6.0]

# show data in view
slice1_5Display = Show(slice1_5, renderView1)
# trace defaults for the display properties.
slice1_5Display.Representation = 'Surface'
slice1_5Display.ColorArrayName = ['POINTS', 'k']
slice1_5Display.LookupTable = kLUT
slice1_5Display.OSPRayScaleArray = 'k'
slice1_5Display.OSPRayScaleFunction = 'PiecewiseFunction'
slice1_5Display.SelectOrientationVectors = 'None'
slice1_5Display.ScaleFactor = 1.6
slice1_5Display.SelectScaleArray = 'None'
slice1_5Display.GlyphType = 'Arrow'
slice1_5Display.GlyphTableIndexArray = 'None'
slice1_5Display.DataAxesGrid = 'GridAxesRepresentation'
slice1_5Display.PolarAxes = 'PolarAxesRepresentation'

# hide data in view
Hide(frozen_sq_cyl_U_k_bij_hifi_Solve_only_Omega_with_Rterm_NEW_ProdLimiter2_Rtermclipped_1_Selectedfoam, renderView1)

# show color bar/color legend
slice1_5Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# set active source
SetActiveSource(frozen_sq_cyl_U_k_bij_hifi_Solve_only_Omega_with_Rterm_NEW_ProdLimiter2_Rtermclipped_1_Selectedfoam)

# create a new 'Slice'
slice1_6 = Slice(Input=frozen_sq_cyl_U_k_bij_hifi_Solve_only_Omega_with_Rterm_NEW_ProdLimiter2_Rtermclipped_1_Selectedfoam)
slice1_6.SliceType = 'Plane'
slice1_6.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice1_6.SliceType.Origin = [4.0, 0.0, 6.0]

# rename source object
RenameSource('x=6', slice1_6)

# Properties modified on slice1_6.SliceType
slice1_6.SliceType.Origin = [6.0, 0.0, 6.0]

# Properties modified on slice1_6.SliceType
slice1_6.SliceType.Origin = [6.0, 0.0, 6.0]

# show data in view
slice1_6Display = Show(slice1_6, renderView1)
# trace defaults for the display properties.
slice1_6Display.Representation = 'Surface'
slice1_6Display.ColorArrayName = ['POINTS', 'k']
slice1_6Display.LookupTable = kLUT
slice1_6Display.OSPRayScaleArray = 'k'
slice1_6Display.OSPRayScaleFunction = 'PiecewiseFunction'
slice1_6Display.SelectOrientationVectors = 'None'
slice1_6Display.ScaleFactor = 1.6
slice1_6Display.SelectScaleArray = 'None'
slice1_6Display.GlyphType = 'Arrow'
slice1_6Display.GlyphTableIndexArray = 'None'
slice1_6Display.DataAxesGrid = 'GridAxesRepresentation'
slice1_6Display.PolarAxes = 'PolarAxesRepresentation'

# hide data in view
Hide(frozen_sq_cyl_U_k_bij_hifi_Solve_only_Omega_with_Rterm_NEW_ProdLimiter2_Rtermclipped_1_Selectedfoam, renderView1)

# show color bar/color legend
slice1_6Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# set active source
SetActiveSource(slice1_5)

# set active source
SetActiveSource(slice1_6)

# set active source
SetActiveSource(frozen_sq_cyl_U_k_bij_hifi_Solve_only_Omega_with_Rterm_NEW_ProdLimiter2_Rtermclipped_1_Selectedfoam)

# create a new 'Slice'
slice1_7 = Slice(Input=frozen_sq_cyl_U_k_bij_hifi_Solve_only_Omega_with_Rterm_NEW_ProdLimiter2_Rtermclipped_1_Selectedfoam)
slice1_7.SliceType = 'Plane'
slice1_7.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice1_7.SliceType.Origin = [4.0, 0.0, 6.0]

# set active source
SetActiveSource(slice1)

# set active source
SetActiveSource(slice1_7)

# rename source object
RenameSource('x=7', slice1_7)

# Properties modified on slice1_7.SliceType
slice1_7.SliceType.Origin = [7.0, 0.0, 6.0]

# Properties modified on slice1_7.SliceType
slice1_7.SliceType.Origin = [7.0, 0.0, 6.0]

# show data in view
slice1_7Display = Show(slice1_7, renderView1)
# trace defaults for the display properties.
slice1_7Display.Representation = 'Surface'
slice1_7Display.ColorArrayName = ['POINTS', 'k']
slice1_7Display.LookupTable = kLUT
slice1_7Display.OSPRayScaleArray = 'k'
slice1_7Display.OSPRayScaleFunction = 'PiecewiseFunction'
slice1_7Display.SelectOrientationVectors = 'None'
slice1_7Display.ScaleFactor = 1.6
slice1_7Display.SelectScaleArray = 'None'
slice1_7Display.GlyphType = 'Arrow'
slice1_7Display.GlyphTableIndexArray = 'None'
slice1_7Display.DataAxesGrid = 'GridAxesRepresentation'
slice1_7Display.PolarAxes = 'PolarAxesRepresentation'

# hide data in view
Hide(frozen_sq_cyl_U_k_bij_hifi_Solve_only_Omega_with_Rterm_NEW_ProdLimiter2_Rtermclipped_1_Selectedfoam, renderView1)

# show color bar/color legend
slice1_7Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# set active source
SetActiveSource(frozen_sq_cyl_U_k_bij_hifi_Solve_only_Omega_with_Rterm_NEW_ProdLimiter2_Rtermclipped_1_Selectedfoam)

# create a new 'Slice'
slice1_8 = Slice(Input=frozen_sq_cyl_U_k_bij_hifi_Solve_only_Omega_with_Rterm_NEW_ProdLimiter2_Rtermclipped_1_Selectedfoam)
slice1_8.SliceType = 'Plane'
slice1_8.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice1_8.SliceType.Origin = [4.0, 0.0, 6.0]

# rename source object
RenameSource('x=8', slice1_8)

# Properties modified on slice1_8.SliceType
slice1_8.SliceType.Origin = [8.0, 0.0, 6.0]

# Properties modified on slice1_8.SliceType
slice1_8.SliceType.Origin = [8.0, 0.0, 6.0]

# show data in view
slice1_8Display = Show(slice1_8, renderView1)
# trace defaults for the display properties.
slice1_8Display.Representation = 'Surface'
slice1_8Display.ColorArrayName = ['POINTS', 'k']
slice1_8Display.LookupTable = kLUT
slice1_8Display.OSPRayScaleArray = 'k'
slice1_8Display.OSPRayScaleFunction = 'PiecewiseFunction'
slice1_8Display.SelectOrientationVectors = 'None'
slice1_8Display.ScaleFactor = 1.6
slice1_8Display.SelectScaleArray = 'None'
slice1_8Display.GlyphType = 'Arrow'
slice1_8Display.GlyphTableIndexArray = 'None'
slice1_8Display.DataAxesGrid = 'GridAxesRepresentation'
slice1_8Display.PolarAxes = 'PolarAxesRepresentation'

# hide data in view
Hide(frozen_sq_cyl_U_k_bij_hifi_Solve_only_Omega_with_Rterm_NEW_ProdLimiter2_Rtermclipped_1_Selectedfoam, renderView1)

# show color bar/color legend
slice1_8Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# set active source
SetActiveSource(slice1_7)

# set active source
SetActiveSource(slice1_8)

# set active source
SetActiveSource(frozen_sq_cyl_U_k_bij_hifi_Solve_only_Omega_with_Rterm_NEW_ProdLimiter2_Rtermclipped_1_Selectedfoam)

# create a new 'Slice'
slice1_9 = Slice(Input=frozen_sq_cyl_U_k_bij_hifi_Solve_only_Omega_with_Rterm_NEW_ProdLimiter2_Rtermclipped_1_Selectedfoam)
slice1_9.SliceType = 'Plane'
slice1_9.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice1_9.SliceType.Origin = [4.0, 0.0, 6.0]

# rename source object
RenameSource('x=9', slice1_9)

# Properties modified on slice1_9.SliceType
slice1_9.SliceType.Origin = [9.0, 0.0, 6.0]

# Properties modified on slice1_9.SliceType
slice1_9.SliceType.Origin = [9.0, 0.0, 6.0]

# show data in view
slice1_9Display = Show(slice1_9, renderView1)
# trace defaults for the display properties.
slice1_9Display.Representation = 'Surface'
slice1_9Display.ColorArrayName = ['POINTS', 'k']
slice1_9Display.LookupTable = kLUT
slice1_9Display.OSPRayScaleArray = 'k'
slice1_9Display.OSPRayScaleFunction = 'PiecewiseFunction'
slice1_9Display.SelectOrientationVectors = 'None'
slice1_9Display.ScaleFactor = 1.6
slice1_9Display.SelectScaleArray = 'None'
slice1_9Display.GlyphType = 'Arrow'
slice1_9Display.GlyphTableIndexArray = 'None'
slice1_9Display.DataAxesGrid = 'GridAxesRepresentation'
slice1_9Display.PolarAxes = 'PolarAxesRepresentation'

# hide data in view
Hide(frozen_sq_cyl_U_k_bij_hifi_Solve_only_Omega_with_Rterm_NEW_ProdLimiter2_Rtermclipped_1_Selectedfoam, renderView1)

# show color bar/color legend
slice1_9Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# set active source
SetActiveSource(frozen_sq_cyl_U_k_bij_hifi_Solve_only_Omega_with_Rterm_NEW_ProdLimiter2_Rtermclipped_1_Selectedfoam)

# create a new 'Slice'
slice1_10 = Slice(Input=frozen_sq_cyl_U_k_bij_hifi_Solve_only_Omega_with_Rterm_NEW_ProdLimiter2_Rtermclipped_1_Selectedfoam)
slice1_10.SliceType = 'Plane'
slice1_10.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice1_10.SliceType.Origin = [4.0, 0.0, 6.0]

# rename source object
RenameSource('x=10', slice1_10)

# Properties modified on slice1_10.SliceType
slice1_10.SliceType.Origin = [10.0, 0.0, 6.0]

# Properties modified on slice1_10.SliceType
slice1_10.SliceType.Origin = [10.0, 0.0, 6.0]

# show data in view
slice1_10Display = Show(slice1_10, renderView1)
# trace defaults for the display properties.
slice1_10Display.Representation = 'Surface'
slice1_10Display.ColorArrayName = ['POINTS', 'k']
slice1_10Display.LookupTable = kLUT
slice1_10Display.OSPRayScaleArray = 'k'
slice1_10Display.OSPRayScaleFunction = 'PiecewiseFunction'
slice1_10Display.SelectOrientationVectors = 'None'
slice1_10Display.ScaleFactor = 1.6
slice1_10Display.SelectScaleArray = 'None'
slice1_10Display.GlyphType = 'Arrow'
slice1_10Display.GlyphTableIndexArray = 'None'
slice1_10Display.DataAxesGrid = 'GridAxesRepresentation'
slice1_10Display.PolarAxes = 'PolarAxesRepresentation'

# hide data in view
Hide(frozen_sq_cyl_U_k_bij_hifi_Solve_only_Omega_with_Rterm_NEW_ProdLimiter2_Rtermclipped_1_Selectedfoam, renderView1)

# show color bar/color legend
slice1_10Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# set active source
SetActiveSource(frozen_sq_cyl_U_k_bij_hifi_Solve_only_Omega_with_Rterm_NEW_ProdLimiter2_Rtermclipped_1_Selectedfoam)

# create a new 'Slice'
slice1_11 = Slice(Input=frozen_sq_cyl_U_k_bij_hifi_Solve_only_Omega_with_Rterm_NEW_ProdLimiter2_Rtermclipped_1_Selectedfoam)
slice1_11.SliceType = 'Plane'
slice1_11.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice1_11.SliceType.Origin = [4.0, 0.0, 6.0]

# rename source object
RenameSource('x=11', slice1_11)

# Properties modified on slice1_11.SliceType
slice1_11.SliceType.Origin = [11.0, 0.0, 6.0]

# Properties modified on slice1_11.SliceType
slice1_11.SliceType.Origin = [11.0, 0.0, 6.0]

# show data in view
slice1_11Display = Show(slice1_11, renderView1)
# trace defaults for the display properties.
slice1_11Display.Representation = 'Surface'
slice1_11Display.ColorArrayName = ['POINTS', 'k']
slice1_11Display.LookupTable = kLUT
slice1_11Display.OSPRayScaleArray = 'k'
slice1_11Display.OSPRayScaleFunction = 'PiecewiseFunction'
slice1_11Display.SelectOrientationVectors = 'None'
slice1_11Display.ScaleFactor = 1.6
slice1_11Display.SelectScaleArray = 'None'
slice1_11Display.GlyphType = 'Arrow'
slice1_11Display.GlyphTableIndexArray = 'None'
slice1_11Display.DataAxesGrid = 'GridAxesRepresentation'
slice1_11Display.PolarAxes = 'PolarAxesRepresentation'

# hide data in view
Hide(frozen_sq_cyl_U_k_bij_hifi_Solve_only_Omega_with_Rterm_NEW_ProdLimiter2_Rtermclipped_1_Selectedfoam, renderView1)

# show color bar/color legend
slice1_11Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# set active source
SetActiveSource(frozen_sq_cyl_U_k_bij_hifi_Solve_only_Omega_with_Rterm_NEW_ProdLimiter2_Rtermclipped_1_Selectedfoam)

# create a new 'Slice'
slice1_12 = Slice(Input=frozen_sq_cyl_U_k_bij_hifi_Solve_only_Omega_with_Rterm_NEW_ProdLimiter2_Rtermclipped_1_Selectedfoam)
slice1_12.SliceType = 'Plane'
slice1_12.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice1_12.SliceType.Origin = [4.0, 0.0, 6.0]

# rename source object
RenameSource('x=12', slice1_12)

# Properties modified on slice1_12.SliceType
slice1_12.SliceType.Origin = [12.0, 0.0, 6.0]

# Properties modified on slice1_12.SliceType
slice1_12.SliceType.Origin = [12.0, 0.0, 6.0]

# show data in view
slice1_12Display = Show(slice1_12, renderView1)
# trace defaults for the display properties.
slice1_12Display.Representation = 'Surface'
slice1_12Display.ColorArrayName = ['POINTS', 'k']
slice1_12Display.LookupTable = kLUT
slice1_12Display.OSPRayScaleArray = 'k'
slice1_12Display.OSPRayScaleFunction = 'PiecewiseFunction'
slice1_12Display.SelectOrientationVectors = 'None'
slice1_12Display.ScaleFactor = 1.6
slice1_12Display.SelectScaleArray = 'None'
slice1_12Display.GlyphType = 'Arrow'
slice1_12Display.GlyphTableIndexArray = 'None'
slice1_12Display.DataAxesGrid = 'GridAxesRepresentation'
slice1_12Display.PolarAxes = 'PolarAxesRepresentation'

# hide data in view
Hide(frozen_sq_cyl_U_k_bij_hifi_Solve_only_Omega_with_Rterm_NEW_ProdLimiter2_Rtermclipped_1_Selectedfoam, renderView1)

# show color bar/color legend
slice1_12Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# set active source
SetActiveSource(frozen_sq_cyl_U_k_bij_hifi_Solve_only_Omega_with_Rterm_NEW_ProdLimiter2_Rtermclipped_1_Selectedfoam)

# create a new 'Slice'
slice1_13 = Slice(Input=frozen_sq_cyl_U_k_bij_hifi_Solve_only_Omega_with_Rterm_NEW_ProdLimiter2_Rtermclipped_1_Selectedfoam)
slice1_13.SliceType = 'Plane'
slice1_13.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice1_13.SliceType.Origin = [4.0, 0.0, 6.0]

# rename source object
RenameSource('x=13', slice1_13)

# Properties modified on slice1_13.SliceType
slice1_13.SliceType.Origin = [13.0, 0.0, 6.0]

# Properties modified on slice1_13.SliceType
slice1_13.SliceType.Origin = [13.0, 0.0, 6.0]

# show data in view
slice1_13Display = Show(slice1_13, renderView1)
# trace defaults for the display properties.
slice1_13Display.Representation = 'Surface'
slice1_13Display.ColorArrayName = ['POINTS', 'k']
slice1_13Display.LookupTable = kLUT
slice1_13Display.OSPRayScaleArray = 'k'
slice1_13Display.OSPRayScaleFunction = 'PiecewiseFunction'
slice1_13Display.SelectOrientationVectors = 'None'
slice1_13Display.ScaleFactor = 1.6
slice1_13Display.SelectScaleArray = 'None'
slice1_13Display.GlyphType = 'Arrow'
slice1_13Display.GlyphTableIndexArray = 'None'
slice1_13Display.DataAxesGrid = 'GridAxesRepresentation'
slice1_13Display.PolarAxes = 'PolarAxesRepresentation'

# hide data in view
Hide(frozen_sq_cyl_U_k_bij_hifi_Solve_only_Omega_with_Rterm_NEW_ProdLimiter2_Rtermclipped_1_Selectedfoam, renderView1)

# show color bar/color legend
slice1_13Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# set active source
SetActiveSource(frozen_sq_cyl_U_k_bij_hifi_Solve_only_Omega_with_Rterm_NEW_ProdLimiter2_Rtermclipped_1_Selectedfoam)

# create a new 'Slice'
slice1_14 = Slice(Input=frozen_sq_cyl_U_k_bij_hifi_Solve_only_Omega_with_Rterm_NEW_ProdLimiter2_Rtermclipped_1_Selectedfoam)
slice1_14.SliceType = 'Plane'
slice1_14.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice1_14.SliceType.Origin = [4.0, 0.0, 6.0]

# rename source object
RenameSource('x=14', slice1_14)

# Properties modified on slice1_14.SliceType
slice1_14.SliceType.Origin = [14.0, 0.0, 6.0]

# Properties modified on slice1_14.SliceType
slice1_14.SliceType.Origin = [14.0, 0.0, 6.0]

# show data in view
slice1_14Display = Show(slice1_14, renderView1)
# trace defaults for the display properties.
slice1_14Display.Representation = 'Surface'
slice1_14Display.ColorArrayName = ['POINTS', 'k']
slice1_14Display.LookupTable = kLUT
slice1_14Display.OSPRayScaleArray = 'k'
slice1_14Display.OSPRayScaleFunction = 'PiecewiseFunction'
slice1_14Display.SelectOrientationVectors = 'None'
slice1_14Display.ScaleFactor = 1.6
slice1_14Display.SelectScaleArray = 'None'
slice1_14Display.GlyphType = 'Arrow'
slice1_14Display.GlyphTableIndexArray = 'None'
slice1_14Display.DataAxesGrid = 'GridAxesRepresentation'
slice1_14Display.PolarAxes = 'PolarAxesRepresentation'

# hide data in view
Hide(frozen_sq_cyl_U_k_bij_hifi_Solve_only_Omega_with_Rterm_NEW_ProdLimiter2_Rtermclipped_1_Selectedfoam, renderView1)

# show color bar/color legend
slice1_14Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# set active source
SetActiveSource(frozen_sq_cyl_U_k_bij_hifi_Solve_only_Omega_with_Rterm_NEW_ProdLimiter2_Rtermclipped_1_Selectedfoam)

# create a new 'Slice'
slice1_15 = Slice(Input=frozen_sq_cyl_U_k_bij_hifi_Solve_only_Omega_with_Rterm_NEW_ProdLimiter2_Rtermclipped_1_Selectedfoam)
slice1_15.SliceType = 'Plane'
slice1_15.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice1_15.SliceType.Origin = [4.0, 0.0, 6.0]

# rename source object
RenameSource('x=15', slice1_15)

# Properties modified on slice1_15.SliceType
slice1_15.SliceType.Origin = [15.0, 0.0, 6.0]

# Properties modified on slice1_15.SliceType
slice1_15.SliceType.Origin = [15.0, 0.0, 6.0]

# show data in view
slice1_15Display = Show(slice1_15, renderView1)
# trace defaults for the display properties.
slice1_15Display.Representation = 'Surface'
slice1_15Display.ColorArrayName = ['POINTS', 'k']
slice1_15Display.LookupTable = kLUT
slice1_15Display.OSPRayScaleArray = 'k'
slice1_15Display.OSPRayScaleFunction = 'PiecewiseFunction'
slice1_15Display.SelectOrientationVectors = 'None'
slice1_15Display.ScaleFactor = 1.6
slice1_15Display.SelectScaleArray = 'None'
slice1_15Display.GlyphType = 'Arrow'
slice1_15Display.GlyphTableIndexArray = 'None'
slice1_15Display.DataAxesGrid = 'GridAxesRepresentation'
slice1_15Display.PolarAxes = 'PolarAxesRepresentation'

# hide data in view
Hide(frozen_sq_cyl_U_k_bij_hifi_Solve_only_Omega_with_Rterm_NEW_ProdLimiter2_Rtermclipped_1_Selectedfoam, renderView1)

# show color bar/color legend
slice1_15Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# set active source
SetActiveSource(frozen_sq_cyl_U_k_bij_hifi_Solve_only_Omega_with_Rterm_NEW_ProdLimiter2_Rtermclipped_1_Selectedfoam)

# create a new 'Slice'
slice1_16 = Slice(Input=frozen_sq_cyl_U_k_bij_hifi_Solve_only_Omega_with_Rterm_NEW_ProdLimiter2_Rtermclipped_1_Selectedfoam)
slice1_16.SliceType = 'Plane'
slice1_16.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice1_16.SliceType.Origin = [4.0, 0.0, 6.0]

# rename source object
RenameSource('x=16', slice1_16)

# Properties modified on slice1_16.SliceType
slice1_16.SliceType.Origin = [16.0, 0.0, 6.0]

# Properties modified on slice1_16.SliceType
slice1_16.SliceType.Origin = [16.0, 0.0, 6.0]

# show data in view
slice1_16Display = Show(slice1_16, renderView1)
# trace defaults for the display properties.
slice1_16Display.Representation = 'Surface'
slice1_16Display.ColorArrayName = ['POINTS', 'k']
slice1_16Display.LookupTable = kLUT
slice1_16Display.OSPRayScaleArray = 'k'
slice1_16Display.OSPRayScaleFunction = 'PiecewiseFunction'
slice1_16Display.SelectOrientationVectors = 'None'
slice1_16Display.ScaleFactor = 1.6
slice1_16Display.SelectScaleArray = 'None'
slice1_16Display.GlyphType = 'Arrow'
slice1_16Display.GlyphTableIndexArray = 'None'
slice1_16Display.DataAxesGrid = 'GridAxesRepresentation'
slice1_16Display.PolarAxes = 'PolarAxesRepresentation'

# hide data in view
Hide(frozen_sq_cyl_U_k_bij_hifi_Solve_only_Omega_with_Rterm_NEW_ProdLimiter2_Rtermclipped_1_Selectedfoam, renderView1)

# show color bar/color legend
slice1_16Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

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
slice1_16Display = Show(slice1_16, spreadSheetView1)

# show data in view
slice1Display = Show(slice1, spreadSheetView1)

# export view
ExportView('/home/unimelb.edu.au/lcampoli/UoM/Testcases/frozen_sq_cyl_U_k_bij_hifi_Solve_only_Omega_with_Rterm_NEW_ProdLimiter2_Rtermclipped_1_Selected/frozen_sq_cyl_U_k_bij_hifi_Solve_only_Omega_with_Rterm_NEW_ProdLimiter2_Rtermclipped_1_Selected/k_x_0.csv', view=spreadSheetView1, FilterColumnsByVisibility=1)

# show data in view
slice1_1Display = Show(slice1_1, spreadSheetView1)

# export view
ExportView('/home/unimelb.edu.au/lcampoli/UoM/Testcases/frozen_sq_cyl_U_k_bij_hifi_Solve_only_Omega_with_Rterm_NEW_ProdLimiter2_Rtermclipped_1_Selected/frozen_sq_cyl_U_k_bij_hifi_Solve_only_Omega_with_Rterm_NEW_ProdLimiter2_Rtermclipped_1_Selected/k_x_1.csv', view=spreadSheetView1, FilterColumnsByVisibility=1)

# show data in view
slice1_2Display = Show(slice1_2, spreadSheetView1)

# set active source
SetActiveSource(slice1_2)

# export view
ExportView('/home/unimelb.edu.au/lcampoli/UoM/Testcases/frozen_sq_cyl_U_k_bij_hifi_Solve_only_Omega_with_Rterm_NEW_ProdLimiter2_Rtermclipped_1_Selected/frozen_sq_cyl_U_k_bij_hifi_Solve_only_Omega_with_Rterm_NEW_ProdLimiter2_Rtermclipped_1_Selected/k_x_2.csv', view=spreadSheetView1, FilterColumnsByVisibility=1)

# show data in view
slice1_3Display = Show(slice1_3, spreadSheetView1)

# export view
ExportView('/home/unimelb.edu.au/lcampoli/UoM/Testcases/frozen_sq_cyl_U_k_bij_hifi_Solve_only_Omega_with_Rterm_NEW_ProdLimiter2_Rtermclipped_1_Selected/frozen_sq_cyl_U_k_bij_hifi_Solve_only_Omega_with_Rterm_NEW_ProdLimiter2_Rtermclipped_1_Selected/k_x_3.csv', view=spreadSheetView1, FilterColumnsByVisibility=1)

# show data in view
slice1_4Display = Show(slice1_4, spreadSheetView1)

# export view
ExportView('/home/unimelb.edu.au/lcampoli/UoM/Testcases/frozen_sq_cyl_U_k_bij_hifi_Solve_only_Omega_with_Rterm_NEW_ProdLimiter2_Rtermclipped_1_Selected/frozen_sq_cyl_U_k_bij_hifi_Solve_only_Omega_with_Rterm_NEW_ProdLimiter2_Rtermclipped_1_Selected/k_x_4.csv', view=spreadSheetView1, FilterColumnsByVisibility=1)

# show data in view
slice1_5Display = Show(slice1_5, spreadSheetView1)

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).