import numpy as np
from paraview.util.vtkAlgorithm import (
    VTKPythonAlgorithmBase,
    smdomain,
    smhint,
    smproperty,
    smproxy,
)
from vtkmodules.numpy_interface import dataset_adapter as dsa
from vtkmodules.vtkCommonDataModel import vtkUnstructuredGrid
import vtk

__author__ = "Benjamin F. Zwick"
__email__ = "benzwick@gmail.com"
__copyright__ = "Copyright (c) 2021 {} <{}>".format(__author__, __email__)
__license__ = "License :: OSI Approved :: BSD-3-Clause License"
__version__ = "0.1.0"

extensions = [
    # MFEM collection
    "mfem_root",
    # MFEM mesh format
    "mesh",
    "mesh.gz",
    # Other meshes supported by MFEM
    "msh",                      # Gmsh
    "vtk",                      # VTK
    # TODO: Add more mesh formats
]

file_description = "MFEM files"

# VTK element type based on MFEM geometry and nodes FE space order
# MFEM: - https://github.com/mfem/mfem/blob/master/fem/geom.hpp
#       - https://github.com/mfem/mfem/blob/master/mesh/vtk.cpp
# VTK:  - https://vtk.org/doc/nightly/html/vtkCellType_8h_source.html
# NOTE: Linear elements load much faster than arbitrary order Lagrange elements

VTKGeometry_Map = [                       # MFEM:
    vtk.VTK_VERTEX,                       # POINT
    vtk.VTK_LINE,                         # SEGMENT
    vtk.VTK_TRIANGLE,                     # TRIANGLE
    vtk.VTK_QUAD,                         # SQUARE
    vtk.VTK_TETRA,                        # TETRAHEDRON
    vtk.VTK_HEXAHEDRON,                   # CUBE
    vtk.VTK_WEDGE]                        # PRISM

VTKGeometry_QuadraticMap = [              # MFEM:
    vtk.VTK_VERTEX,                       # POINT
    vtk.VTK_QUADRATIC_EDGE,               # SEGMENT
    vtk.VTK_QUADRATIC_TRIANGLE,           # TRIANGLE
    vtk.VTK_BIQUADRATIC_QUAD,             # SQUARE
    vtk.VTK_QUADRATIC_TETRA,              # TETRAHEDRON
    vtk.VTK_TRIQUADRATIC_HEXAHEDRON,      # CUBE
    vtk.VTK_BIQUADRATIC_QUADRATIC_WEDGE]  # PRISM

VTKGeometry_HighOrderMap = [              # MFEM:
    vtk.VTK_VERTEX,                       # POINT
    vtk.VTK_LAGRANGE_CURVE,               # SEGMENT
    vtk.VTK_LAGRANGE_TRIANGLE,            # TRIANGLE
    vtk.VTK_LAGRANGE_QUADRILATERAL,       # SQUARE
    vtk.VTK_LAGRANGE_TETRAHEDRON,         # TETRAHEDRON
    vtk.VTK_LAGRANGE_HEXAHEDRON,          # CUBE
    vtk.VTK_LAGRANGE_WEDGE]               # PRISM

VTKGeometry_PrismMap = [0, 2, 1, 3, 5, 4]

VTKGeometry_VertexPermutation = [
    None, None, None, None, None, None, VTKGeometry_PrismMap]

# https://github.com/mfem/mfem/blob/master/mesh/mesh_readers.cpp
Mesh_vtk_quadratic_tet = [
    0, 1, 2, 3, 4, 7, 5, 6, 8, 9]
Mesh_vtk_quadratic_wedge = [
    0, 2, 1, 3, 5, 4, 8, 7, 6, 11, 10, 9, 12, 14, 13, 17, 16, 15]
Mesh_vtk_quadratic_hex = [
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
    24, 22, 21, 23, 20, 25, 26]

VTKGeometry_QuadraticVertexPermutation = [
    None, None, None, None,
    Mesh_vtk_quadratic_tet,
    Mesh_vtk_quadratic_hex,
    Mesh_vtk_quadratic_wedge]

@smproxy.reader(
    label="MFEM reader",
    extensions=extensions,
    file_description=file_description)
class MFEMReader(VTKPythonAlgorithmBase):
    def __init__(self):
        VTKPythonAlgorithmBase.__init__(
            self, nInputPorts=0, nOutputPorts=1,
            outputType="vtkUnstructuredGrid")
        self._filename = None

    @smproperty.stringvector(name="FileName")
    @smdomain.filelist()
    @smhint.filechooser(extensions=extensions,
                        file_description=file_description)
    def SetFileName(self, filename):
        if self._filename != filename:
            self._filename = filename
            self.Modified()

    def RequestData(self, request, inInfo, outInfo):
        output = dsa.WrapDataObject(vtkUnstructuredGrid.GetData(outInfo))

        import os
        import json
        import mfem.ser as mfem

        _, ext = os.path.splitext(self._filename)
        cwd = os.path.dirname(self._filename)

        if ext == ".mfem_root":
            with open(self._filename) as f:
                root = json.load(f)
                cycle = int(root['dsets']['main']['cycle'])
                mesh_filename = root['dsets']['main']['mesh']['path']
                mesh_filename = format(mesh_filename % cycle)
                mesh_filename = os.path.join(cwd, mesh_filename)
        else:
            mesh_filename = self._filename

        mesh = mfem.Mesh(mesh_filename)
        dim = mesh.SpaceDimension()
        nelem = mesh.GetNE()

        # Points
        nodes = mesh.GetNodes()
        if nodes is None:
            nnode = mesh.GetNV()
            points = np.array(mesh.GetVertexArray())
        else:
            nnode = mesh.GetNodes().Size() // dim
            points = np.array(mesh.GetNodes().GetDataArray())
            mesh_fes = mesh.GetNodalFESpace()
            if mesh_fes.GetOrdering() == 0:
                points = points.reshape((dim, nnode)).T
            elif mesh_fes.GetOrdering() == 1:
                points = points.reshape((nnode, dim))
            else:
                raise NotImplementedError
        if dim == 1:
            points = np.hstack([points, np.zeros((nnode, 2))])
        elif dim == 2:
            points = np.hstack([points, np.zeros((nnode, 1))])
        output.SetPoints(points)

        # Cells
        cell_attributes = np.empty((nelem), dtype=int)
        cell_types = np.empty((nelem), dtype=np.ubyte)
        cell_offsets = np.empty((nelem), dtype=int)
        cell_conn = []

        # Element order for each mesh element
        elem_order = np.ones((nelem), dtype=int)

        offset = 0
        if nodes is None:
            for i in range(nelem):
                v = mesh.GetElementVertices(i)
                geom = mesh.GetElementBaseGeometry(i)
                perm = VTKGeometry_VertexPermutation[geom]
                if perm: v = [v[i] for i in perm]
                cell_types[i] = VTKGeometry_Map[geom]
                cell_offsets[i] = offset
                offset += len(v) + 1
                cell_conn.append(len(v))
                cell_conn.extend(v)
        else:
            for i in range(nelem):
                v = mesh_fes.GetElementDofs(i)
                geom = mesh.GetElementBaseGeometry(i)
                order = mesh_fes.GetOrder(i)
                if order == 1:
                    perm = VTKGeometry_VertexPermutation[geom]
                    cell_types[i] = VTKGeometry_Map[geom]
                elif order == 2:
                    perm = VTKGeometry_QuadraticVertexPermutation[geom]
                    cell_types[i] = VTKGeometry_QuadraticMap[geom]
                else:
                    print(f"Failed to read mesh file {mesh_filename}")
                    print(f"ERROR: Only elements of order 1 and 2 are supported")
                    return 1
                if perm: v = [v[i] for i in perm]
                cell_offsets[i] = offset
                offset += len(v) + 1
                cell_conn.append(len(v))
                cell_conn.extend(v)

        output.SetCells(cell_types, cell_offsets, cell_conn)

        # Attributes
        for i in range(nelem):
            cell_attributes[i] = mesh.GetAttribute(i)
        output.CellData.append(cell_attributes, "attribute")

        # Read fields
        if ext == ".mfem_root":
            fields = root['dsets']['main']['fields']
            for name, prop in fields.items():
                filename = prop['path']
                filename = format(filename % cycle)
                filename = os.path.join(cwd, filename)
                gf = mfem.GridFunction(mesh, filename)
                gf_fes = gf.FESpace()
                gf_vdim = gf.VectorDim()
                gf_nnode = gf_fes.GetNDofs() // gf_vdim
                gf_nelem = gf_fes.GetNBE()
                data = np.array(gf.GetDataArray())
                if prop['tags']['assoc'] == 'nodes':
                    assert gf_nnode == nnode, "Mesh and grid function have different number of nodes"
                    if gf_fes.GetOrdering() == 0:
                        data = data.reshape((gf_vdim, gf_nnode)).T
                    elif gf_fes.GetOrdering() == 1:
                        data = data.reshape((gf_nnode, gf_vdim))
                    else:
                        raise NotImplementedError
                    output.PointData.append(data, name)
                elif prop['tags']['assoc'] == 'elements':
                    assert gf_nelem == nelem, "Mesh and grid function have different number of elements"
                    if gf_fes.GetOrdering() == 0:
                        data = data.reshape((gf_vdim, gf_nelem)).T
                    elif gf_fes.GetOrdering() == 1:
                        data = data.reshape((gf_nelem, gf_vdim))
                    else:
                        raise NotImplementedError
                    output.CellData.append(data, name)
                else:
                    raise NotImplementedError("assoc: '{}'".format(prop['tags']['assoc']))

        return 1
