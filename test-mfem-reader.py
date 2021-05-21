# This is a simple ParaView macro to test the paraview-mfem-reader
# by loading some example meshes.

from paraview.simple import *
import os
import time

#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

datapath = "/home/ben/projects/mfem/mfem/data"

filenames = [
    "amr-hex.mesh",
    "amr-quad.mesh",
    "ball-nurbs.mesh",
    "beam-hex-nurbs.mesh",
    "beam-hex.mesh",
    "beam-hex.vtk",
    "beam-quad-amr.mesh",
    "beam-quad-nurbs.mesh",
    "beam-quad.mesh",
    "beam-quad.vtk",
    "beam-tet.mesh",
    "beam-tet.vtk",
    "beam-tri.mesh",
    "beam-tri.vtk",
    "beam-wedge.mesh",
    "beam-wedge.vtk",
    "cube-nurbs.mesh",
    "disc-nurbs.mesh",
    "escher-p2.mesh",
    "escher-p2.vtk",
    "escher-p3.mesh",
    "escher.mesh",
    "escher.vtk",
    "fichera-amr.mesh",
    "fichera-mixed-p2.mesh",
    "fichera-mixed-p2.vtk",
    "fichera-mixed.mesh",
    "fichera-q2.mesh",
    "fichera-q2.vtk",
    "fichera-q3.mesh",
    "fichera.mesh",
    "fichera.vtk",
    "inline-hex.mesh",
    "inline-quad.mesh",
    "inline-segment.mesh",
    "inline-tet.mesh",
    "inline-tri.mesh",
    "inline-wedge.mesh",
    "klein-bottle.mesh",
    "klein-donut.mesh",
    "l-shape.mesh",
    "mobius-strip.mesh",
    "periodic-annulus-sector.msh",
    "periodic-cube.mesh",
    "periodic-cube.msh",
    "periodic-hexagon.mesh",
    "periodic-segment.mesh",
    "periodic-square.mesh",
    "periodic-square.msh",
    "periodic-torus-sector.msh",
    "pipe-nurbs-2d.mesh",
    "pipe-nurbs.mesh",
    "rt-2d-p4-tri.mesh",
    "rt-2d-q3.mesh",
    "square-disc-nurbs-patch.mesh",
    "square-disc-nurbs.mesh",
    "square-disc-p2.mesh",
    "square-disc-p2.vtk",
    "square-disc-p3.mesh",
    "square-disc-surf.mesh",
    "square-disc.mesh",
    "square-disc.vtk",
    "square-mixed.mesh",
    "square-nurbs.mesh",
    "star-hilbert.mesh",
    "star-mixed-p2.mesh",
    "star-mixed-p2.vtk",
    "star-mixed.mesh",
    "star-q2.mesh",
    "star-q2.vtk",
    "star-q3.mesh",
    "star-surf.mesh",
    "star.mesh",
    "star.vtk",
    "toroid-hex.mesh",
    "toroid-wedge.mesh",
]

for fname in filenames:
    print(f"Reading mesh: {fname}")
    source = MFEMreader(FileName=os.path.join(datapath, fname))
    RenameSource(fname, source)
    # source = GetActiveSource()
    # SetActiveSource(source)
    renderView1 = GetActiveViewOrCreate('RenderView')
    display = Show(source, renderView1, 'UnstructuredGridRepresentation')
    display.SetRepresentationType('Surface With Edges')
    display.NonlinearSubdivisionLevel = 3
    # time.sleep(0.8)
    Hide(source, renderView1)
    # RenderAllViews()
