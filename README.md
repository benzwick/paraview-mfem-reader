# ParaView MFEM Reader

ParaView plugin for reading [MFEM](https://mfem.org) files.

This plugin can be used to read
[MFEM mesh](https://mfem.org/mesh-formats)
(`.mesh` and `.mesh.gz`) files
as well as other mesh formats supported by MFEM such as Gmsh (`.msh`) files.
It can also be used to read
[MFEM data collection](https://mfem.github.io/doxygen/html/classmfem_1_1VisItDataCollection.html)
(`.mfem_root`) files that contain a mesh and associated fields.

## Installing

This plugin requires [PyMFEM](https://github.com/mfem/PyMFEM)
to be installed within the ParaView Python environment.
Please see the references below for additional details.

1. Set variables:

        export PV_PY_VER=3.7.4
        export PV_NP_VER=1.16.4
        export PV_SIX_VER=1.11.0
        export PV_PY_PREFIX=/opt/paraview/pvpy-${PV_PY_VER}
        export PV_PY_VENV=${PV_PY_PREFIX}-venv

2. Install Python:

        wget https://www.python.org/ftp/python/${PV_PY_VER}/Python-${PV_PY_VER}.tgz
        tar xzvf Python-${PV_PY_VER}.tgz
        cd Python-${PV_PY_VER}
        ./configure --enable-optimizations --prefix=${PV_PY_PREFIX}
        make
        make install
        cd ..

3. Create virtual environment:

        ${PV_PY_PREFIX}/bin/python3 -m venv ${PV_PY_VENV}

4. Install packages into virtual environment:

        source pvpy-${PV_PY_VER}-venv/bin/activate
        pip uninstall -y numpy
        pip install --ignore-installed numpy==$PV_NP_VER
        pip uninstall -y six
        pip install --ignore-installed six==$PV_SIX_VER
        pip install mfem

5. Before starting ParaView (note that python3.7 is used here):

        export PYTHONPATH=${PV_PY_VENV}/lib/python3.7/site-packages

6. Load the `MFEMReader.py` module in ParaView
   via `Tools` &rarr; `Manage Plugins` &rarr; `Load New`.
   You can optionally check the `Auto Load` option.

You should now be able to open MFEM files using `File` &rarr; `Open`.

## References

- <https://github.com/DaanVanVugt/paraview-python-file-reader>
- <https://www.litianyi.me/2019/10/27/paraview-meshio-plugin/>
- <https://github.com/tianyikillua/paraview-dolfinx-reader>
- <https://github.com/nschloe/meshio/blob/main/tools/paraview-meshio-plugin.py>
- <https://blog.kitware.com/easy-customization-of-the-paraview-python-programmable-filter-property-panel>
- <https://www.paraview.org/Wiki/Python_calculator_and_programmable_filter>
- <https://public.kitware.com/pipermail/paraview-developers/2017-January/005051.html>
- <https://discourse.paraview.org/t/how-can-i-install-and-import-other-modules-inside-pvpython/3067/6>
- <https://gitlab.kitware.com/paraview/paraview/-/issues/17891>
- <https://docs.python.org/3/library/venv.html>
