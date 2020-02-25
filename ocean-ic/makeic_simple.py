#!/usr/bin/env python

from __future__ import print_function

import sys, os
import argparse
import subprocess as sp
import netCDF4 as nc

from regridder import regrid

"""
Create ocean model IC based on reanalysis data.
"""

grid_defs_error = \
"""
Grid definitions directory {} not found.
please download it with:
wget http://s3-ap-southeast-2.amazonaws.com/dp-drop/ocean-regrid/grid_defs.tar.gz
and unzip into the same directory as this executable.
"""

def grid_defs_dir():
    """
    Get path to directory where MOM, NEMO, GODAS and ORAS4 grid definitions are
    found.
    """

    if getattr(sys, 'frozen', False):
        basedir = sys._MEIPASS
    else:
        basedir = os.path.dirname(os.path.realpath(__file__))

    return os.path.join(basedir, 'grid_defs')

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('reanalysis_name', help="""
                        Name of src data/grid, must be GODAS or ORAS4""")
    parser.add_argument('temp_reanalysis_file', help='Temperature file from reanalysis.')
    parser.add_argument('salt_reanalysis_file', help='Salt file from reanalysis.')

    parser.add_argument('model_name', help="""
                        Name of model, must be MOM, MOM1 or NEMO""")
    parser.add_argument('output_file', help='Name of the destination/output file.')
    args = parser.parse_args()

    assert args.model_name == 'MOM' or args.model_name == 'MOM1' or \
        args.model_name == 'NEMO'
    assert args.reanalysis_name == 'GODAS' or args.reanalysis_name == 'ORAS4'

    if os.path.exists(args.output_file):
        print("Output file {} already exists, ".format(args.output_file) + \
               "please move or delete.", file=sys.stderr)
        return 1

    # Set up the reanalysis model and grid definitions.
    grid_defs = grid_defs_dir()
    if not os.path.exists(grid_defs):
        print(grid_defs_error.format(grid_defs), file=sys.stderr)
        return 1

    if args.reanalysis_name == 'GODAS':
        reanalysis_hgrids = (os.path.join(grid_defs, 'pottmp.2016.nc'),)
        reanalysis_vgrid = os.path.join(grid_defs, 'pottmp.2016.nc')
    else:
        reanalysis_hgrids = (os.path.join(grid_defs, 'coordinates_grid_T.nc'),
                             os.path.join(grid_defs, 'coordinates_grid_U.nc'),
                             os.path.join(grid_defs, 'coordinates_grid_V.nc'))
        reanalysis_vgrid = os.path.join(grid_defs, 'coordinates_grid_T.nc')

    if args.model_name == 'MOM':
        model_hgrid = os.path.join(grid_defs, 'ocean_hgrid.nc')
        model_vgrid = os.path.join(grid_defs, 'ocean_vgrid.nc')
        model_mask = os.path.join(grid_defs, 'ocean_mask.nc')
        mm_arg = ['--model_mask', model_mask]
    elif args.model_name == 'MOM1':
        model_hgrid = os.path.join(grid_defs, 'grid_spec.nc')
        model_vgrid = os.path.join(grid_defs, 'grid_spec.nc')
        model_mask = os.path.join(grid_defs, 'grid_spec.nc')
        mm_arg = ['--model_mask', model_mask]
    else:
        model_hgrid = os.path.join(grid_defs, 'coordinates.nc')
        model_vgrid = os.path.join(grid_defs, 'data_1m_potential_temperature_nomask.nc')
        model_mask = None
        mm_arg = []

    args = [args.reanalysis_name, reanalysis_hgrids[0], reanalysis_vgrid,
            args.temp_reanalysis_file, args.salt_reanalysis_file, args.model_name,
            model_hgrid, model_vgrid, args.output_file] + mm_arg
    exe = os.path.join(os.path.dirname(__file__), 'makeic.py')
    return sp.call([exe] + args)


if __name__ == '__main__':
    sys.exit(main())
