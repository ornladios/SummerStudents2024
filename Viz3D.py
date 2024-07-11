import sys

# Get the filename from the command line arguments
if len(sys.argv) < 2:
    print("Usage: pvpython 3DViz.py <filename.bp>")
    sys.exit(1)

filename = sys.argv[1]

# Ensure the rest of your imports and initializations are in place
import paraview
paraview.compatibility.major = 5
paraview.compatibility.minor = 13

# Import the simple module from paraview
from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()

# Setup views used in the visualization
materialLibrary1 = GetMaterialLibrary()

renderView1 = CreateView('RenderView')
renderView1.ViewSize = [1240, 1060]
renderView1.AxesGrid = 'Grid Axes 3D Actor'
renderView1.CenterOfRotation = [32.0, 32.0, 32.0]
renderView1.StereoType = 'Crystal Eyes'
renderView1.CameraPosition = [23.065975008580676, 50.53798701959655, 241.31507652216507]
renderView1.CameraFocalPoint = [32.1948708394694, 32.11913858039407, 32.43437535549241]
renderView1.CameraViewUp = [0.21821012304671297, 0.9752824459795023, 0.034763382519218844]
renderView1.CameraFocalDisk = 1.0
renderView1.CameraParallelScale = 55.42562584220407
renderView1.LegendGrid = 'Legend Grid Actor'
renderView1.PolarGrid = 'Polar Grid Actor'
renderView1.BackEnd = 'OSPRay raycaster'
renderView1.OSPRayMaterialLibrary = materialLibrary1

SetActiveView(None)

layout1 = CreateLayout(name='Layout #1')
layout1.AssignView(0, renderView1)
layout1.SetSize(1240, 1060)

SetActiveView(renderView1)

uncomsummerbp = ADIOS2VTXReader(registrationName='input.bp', FileName=filename)

contour1 = Contour(registrationName='Contour1', Input=uncomsummerbp)
contour1.ContourBy = ['POINTS', 'F']
contour1.Isosurfaces = [5]
contour1.PointMergeMethod = 'Uniform Binning'

# Setup the visualization in view 'renderView1'
uncomsummerbpDisplay = Show(uncomsummerbp, renderView1, 'UniformGridRepresentation')

uncomsummerbpDisplay.Representation = 'Outline'
uncomsummerbpDisplay.ColorArrayName = [None, '']
uncomsummerbpDisplay.SelectNormalArray = 'None'
uncomsummerbpDisplay.SelectTangentArray = 'None'
uncomsummerbpDisplay.SelectTCoordArray = 'None'
uncomsummerbpDisplay.TextureTransform = 'Transform2'
uncomsummerbpDisplay.OSPRayScaleArray = 'F'
uncomsummerbpDisplay.OSPRayScaleFunction = 'Piecewise Function'
uncomsummerbpDisplay.Assembly = 'Hierarchy'
uncomsummerbpDisplay.SelectedBlockSelectors = ['']
uncomsummerbpDisplay.SelectOrientationVectors = 'None'
uncomsummerbpDisplay.ScaleFactor = 6.4
uncomsummerbpDisplay.SelectScaleArray = 'F'
uncomsummerbpDisplay.GlyphType = 'Arrow'
uncomsummerbpDisplay.GlyphTableIndexArray = 'F'
uncomsummerbpDisplay.GaussianRadius = 0.32
uncomsummerbpDisplay.SetScaleArray = ['POINTS', 'F']
uncomsummerbpDisplay.ScaleTransferFunction = 'Piecewise Function'
uncomsummerbpDisplay.OpacityArray = ['POINTS', 'F']
uncomsummerbpDisplay.OpacityTransferFunction = 'Piecewise Function'
uncomsummerbpDisplay.DataAxesGrid = 'Grid Axes Representation'
uncomsummerbpDisplay.PolarAxes = 'Polar Axes Representation'
uncomsummerbpDisplay.ScalarOpacityUnitDistance = 1.7320508075688772
uncomsummerbpDisplay.OpacityArrayName = ['POINTS', 'F']
uncomsummerbpDisplay.ColorArray2Name = ['POINTS', 'F']
uncomsummerbpDisplay.SliceFunction = 'Plane'
uncomsummerbpDisplay.Slice = 32
uncomsummerbpDisplay.SelectInputVectors = [None, '']
uncomsummerbpDisplay.WriteLog = ''

uncomsummerbpDisplay.ScaleTransferFunction.Points = [0.020257519693470085, 0.0, 0.5, 0.0, 382.29231102133963, 1.0, 0.5, 0.0]
uncomsummerbpDisplay.OpacityTransferFunction.Points = [0.020257519693470085, 0.0, 0.5, 0.0, 382.29231102133963, 1.0, 0.5, 0.0]
uncomsummerbpDisplay.SliceFunction.Origin = [32.0, 32.0, 32.0]

contour1Display = Show(contour1, renderView1, 'GeometryRepresentation')

fTF2D = GetTransferFunction2D('F')

fLUT = GetColorTransferFunction('F')
fLUT.TransferFunction2D = fTF2D
fLUT.RGBPoints = [0.22207867914628324, 0.231373, 0.298039, 0.752941, 95.7048045983622, 0.865003, 0.865003, 0.865003, 191.18753051757812, 0.705882, 0.0156863, 0.14902]
fLUT.ScalarRangeInitialized = 1.0

contour1Display.Representation = 'Surface'
contour1Display.ColorArrayName = ['POINTS', 'F']
contour1Display.LookupTable = fLUT
contour1Display.SelectNormalArray = 'Normals'
contour1Display.SelectTangentArray = 'None'
contour1Display.SelectTCoordArray = 'None'
contour1Display.TextureTransform = 'Transform2'
contour1Display.OSPRayScaleArray = 'F'
contour1Display.OSPRayScaleFunction = 'Piecewise Function'
contour1Display.Assembly = 'Hierarchy'
contour1Display.SelectedBlockSelectors = ['']
contour1Display.SelectOrientationVectors = 'None'
contour1Display.ScaleFactor = 5.074971914291382
contour1Display.SelectScaleArray = 'F'
contour1Display.GlyphType = 'Arrow'
contour1Display.GlyphTableIndexArray = 'F'
contour1Display.GaussianRadius = 0.2537485957145691
contour1Display.SetScaleArray = ['POINTS', 'F']
contour1Display.ScaleTransferFunction = 'Piecewise Function'
contour1Display.OpacityArray = ['POINTS', 'F']
contour1Display.OpacityTransferFunction = 'Piecewise Function'
contour1Display.DataAxesGrid = 'Grid Axes Representation'
contour1Display.PolarAxes = 'Polar Axes Representation'
contour1Display.SelectInputVectors = ['POINTS', 'Normals']
contour1Display.WriteLog = ''

contour1Display.ScaleTransferFunction.Points = [191.15628051757812, 0.0, 0.5, 0.0, 191.18753051757812, 1.0, 0.5, 0.0]
contour1Display.OpacityTransferFunction.Points = [191.15628051757812, 0.0, 0.5, 0.0, 191.18753051757812, 1.0, 0.5, 0.0]

fLUTColorBar = GetScalarBar(fLUT, renderView1)
fLUTColorBar.WindowLocation = 'Any Location'
fLUTColorBar.Position = [0.8275201612903226, 0.05094339622641508]
fLUTColorBar.Title = 'F'
fLUTColorBar.ComponentTitle = ''

fLUTColorBar.Visibility = 1
contour1Display.SetScalarBarVisibility(renderView1, True)

fPWF = GetOpacityTransferFunction('F')
fPWF.Points = [0.22207867914628324, 0.0, 0.5, 0.0, 191.18753051757812, 1.0, 0.5, 0.0]
fPWF.ScalarRangeInitialized = 1

timeAnimationCue1 = GetTimeTrack()
timeKeeper1 = GetTimeKeeper()
animationScene1 = GetAnimationScene()
animationScene1.ViewModules = renderView1
animationScene1.Cues = timeAnimationCue1
animationScene1.AnimationTime

pNG1 = CreateExtractor('PNG', renderView1, registrationName='PNG1')
# trace defaults for the extractor.
pNG1.Trigger = 'Time Step'

# init the 'PNG' selected for 'Writer'
pNG1.Writer.FileName = 'RenderView1_{timestep:06d}{camera}.png'
pNG1.Writer.ImageResolution = [1240, 1060]
pNG1.Writer.Format = 'PNG'

# ----------------------------------------------------------------
# restore active source
SetActiveSource(pNG1)
# ----------------------------------------------------------------

# ------------------------------------------------------------------------------
# Catalyst options
from paraview import catalyst
options = catalyst.Options()
options.GlobalTrigger = 'Time Step'
options.CatalystLiveTrigger = 'Time Step'

# ------------------------------------------------------------------------------
if __name__ == '__main__':
    from paraview.simple import SaveExtractsUsingCatalystOptions
    # Code for non in-situ environments; if executing in post-processing
    # i.e. non-Catalyst mode, let's generate extracts using Catalyst options
    SaveExtractsUsingCatalystOptions(options)