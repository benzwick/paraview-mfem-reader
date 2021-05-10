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

# MFEM: https://github.com/mfem/mfem/blob/master/fem/geom.hpp
# VTK:  https://vtk.org/doc/nightly/html/vtkCellType_8h_source.html
mfem_to_vtk_type = [
    1,  # POINT = 0
    3,  # SEGMENT = 1
    5,  # TRIANGLE = 2
    9,  # SQUARE = 3
    10, # TETRAHEDRON = 4
    12, # CUBE = 5
    13, # PRISM = 6
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
                mesh_filename = os.path.join(cwd, root['dsets']['main']['mesh']['path'])
                # FIXME: create mesh filename using 'cycle'
                mesh_filename = mesh_filename[:-4] + "000000"
        else:
            mesh_filename = self._filename

        mesh = mfem.Mesh(mesh_filename)
        dim = mesh.Dimension()
        nelem = mesh.GetNE()
        nvert = mesh.GetNV()

        # Points
        points = mesh.GetVertexArray()
        if dim == 1:
            points = np.hstack([points, np.zeros((nvert, 2))])
        elif dim == 2:
            points = np.hstack([points, np.zeros((nvert, 1))])
        output.SetPoints(points)
        del points

        # Cells
        cell_attributes = np.empty((nelem), dtype=int)
        cell_types = np.empty((nelem), dtype=np.ubyte)
        cell_offsets = np.empty((nelem), dtype=int)
        cell_conn = []

        offset = 0
        for i in range(nelem):
            a = mesh.GetAttribute(i)
            t = mesh.GetElementType(i)
            v = mesh.GetElementVertices(i)
            cell_attributes[i] = a
            cell_types[i] = mfem_to_vtk_type[t]
            cell_offsets[i] = offset
            offset += len(v) + 1
            cell_conn.append(len(v))
            cell_conn.extend(v)

        output.SetCells(cell_types, cell_offsets, cell_conn)
        del cell_types, cell_offsets, cell_conn

        output.CellData.append(cell_attributes, "attribute")
        del cell_attributes

        # Read fields
        if ext == ".mfem_root":
            fields = root['dsets']['main']['fields']
            for name, prop in fields.items():
                filename = os.path.join(cwd, prop['path'])
                # FIXME: create field filename using 'cycle'
                filename = filename[:-4] + "000000"
                data = mfem.GridFunction(mesh, filename).GetDataArray()
                if prop['tags']['assoc'] == 'nodes':
                    output.PointData.append(data, name)
                    del data
                if prop['tags']['assoc'] == 'elements':
                    output.CellData.append(data, name)
                    del data
                else:
                    raise NotImplementedError("assoc: '{}'".format(prop['tags']['assoc']))

        return 1
