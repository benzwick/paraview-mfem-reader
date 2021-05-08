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

extensions = ["mesh"]
file_description = "MFEM mesh files"

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

        with open(self._filename, "r") as f:
            lines = f.readlines()

        mesh_format = lines[0].strip()
        assert mesh_format == "MFEM mesh v1.0", \
            "Unsupported mesh format: {}".format(mesh_format)

        for i, l in enumerate(lines):
            if l.startswith('dimension'):
                dim = int(lines[i+1])
            elif l.startswith('elements'):
                nelem = int(lines[i + 1])
                elem_start = i + 2
                elem_end = elem_start + nelem
            elif l.startswith('vertices'):
                nvert = int(lines[i + 1])
                vdim = int(lines[i + 2])
                vert_start = i + 3
                vert_end = vert_start + nvert

        elem = lines[elem_start:elem_end]
        vert = lines[vert_start:vert_end]
        assert len(elem) == nelem
        assert len(vert) == nvert
        del lines

        # Points
        points = np.array([np.array(v.split(), dtype=float) for v in vert])
        if dim == 2:
            points = np.hstack([points, np.zeros((nvert, 1))])
        del vert
        output.SetPoints(points)
        del points

        # Cells
        cell_attributes = np.empty((nelem), dtype=int)
        cell_types = np.empty((nelem), dtype=np.ubyte)
        cell_offsets = np.empty((nelem), dtype=int)
        cell_conn = []

        offset = 0
        for i, e in enumerate(elem):
            s = [int(x) for x in e.split()]
            a, t, v = s[0], s[1], s[2:]
            cell_attributes[i] = a
            cell_types[i] = mfem_to_vtk_type[t]
            cell_offsets[i] = offset
            offset += len(v) + 1
            cell_conn.append(len(v))
            cell_conn.extend(v)
        del elem

        output.SetCells(cell_types, cell_offsets, cell_conn)
        del cell_types, cell_offsets, cell_conn

        output.CellData.append(cell_attributes, "attribute")
        del cell_attributes

        return 1
