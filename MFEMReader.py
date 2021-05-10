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

__author__ = "Benjamin F. Zwick"
__email__ = "benzwick@gmail.com"
__copyright__ = "Copyright (c) 2021 {} <{}>".format(__author__, __email__)
__license__ = "License :: OSI Approved :: BSD-3-Clause License"
__version__ = "0.1.0"

# VTK element type based on MFEM geometry and nodes FE space order
# MFEM: https://github.com/mfem/mfem/blob/master/fem/geom.hpp
# VTK:  https://vtk.org/doc/nightly/html/vtkCellType_8h_source.html
# NOTE: Linear elements load much faster than arbitrary order Lagrange elements
mfem_to_vtk_type = [
    [ 1,  1,  1], # POINT = 0
    [-1,  3, 68], # SEGMENT = 1
    [-1,  5, 69], # TRIANGLE = 2
    [-1,  9, 70], # SQUARE = 3
    [-1, 10, 71], # TETRAHEDRON = 4
    [-1, 12, 72], # CUBE = 5
    [-1, 13, 73], # PRISM = 6
]

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
        dim = mesh.Dimension()
        nelem = mesh.GetNE()

        # FIXME: read field grid functions first to determine the correct order
        #        to use for the mesh curvature.

        # Points
        if mesh.GetNodes() is None:
            mesh_is_curved = False
            nnode = mesh.GetNV()
            points = mesh.GetVertexArray()
        else:
            print("WARNING: Reading curved meshes is untested and probably broken!")
            mesh_is_curved = True
            nnode = mesh.GetNodes().Size() // dim
            points = mesh.GetNodes().GetDataArray()
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
        del points

        # Cells
        cell_attributes = np.empty((nelem), dtype=int)
        cell_types = np.empty((nelem), dtype=np.ubyte)
        cell_offsets = np.empty((nelem), dtype=int)
        cell_conn = []

        # Function used to get element nodal connectivity
        if mesh_is_curved:
            # FIXME: MFEM ordering is different to VTK
            get_conn = mesh_fes.GetElementDofs
        else:
            get_conn = mesh.GetElementVertices

        # Element order for each mesh element
        elem_order = np.ones((nelem), dtype=int)
        if mesh_is_curved:
            for i in range(nelem):
                elem_order[i] = mesh_fes.GetOrder(i)

        offset = 0
        for i in range(nelem):
            a = mesh.GetAttribute(i)
            t = mesh.GetElementType(i)
            v = get_conn(i)
            cell_attributes[i] = a
            # FIXME: add support for mesh order > 1
            cell_types[i] = mfem_to_vtk_type[t][min(elem_order[i], 2)]
            cell_offsets[i] = offset
            offset += len(v) + 1
            cell_conn.append(len(v))
            cell_conn.extend(v)

        output.SetCells(cell_types, cell_offsets, cell_conn)
        del cell_types, cell_offsets, cell_conn, elem_order

        output.CellData.append(cell_attributes, "attribute")
        del cell_attributes

        # Read fields
        if ext == ".mfem_root":
            fields = root['dsets']['main']['fields']
            for name, prop in fields.items():
                filename = prop['path']
                filename = format(filename % cycle)
                filename = os.path.join(cwd, filename)
                gf = mfem.GridFunction(mesh, filename)
                vdim = gf.VectorDim()
                # FIXME: add support for vector fields which will required reshape
                #        of data based on ordering similar to nodes.
                if vdim != 1:
                    raise NotImplementedError("Only scalar fields (VDim == 1) are allowed.")
                data = gf.GetDataArray()
                if prop['tags']['assoc'] == 'nodes':
                    output.PointData.append(data, name)
                    del data
                elif prop['tags']['assoc'] == 'elements':
                    output.CellData.append(data, name)
                    del data
                else:
                    raise NotImplementedError("assoc: '{}'".format(prop['tags']['assoc']))

        return 1
