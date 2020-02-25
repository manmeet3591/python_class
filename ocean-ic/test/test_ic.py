
from __future__ import print_function

import pytest
import os
import subprocess as sp
import sh
import netCDF4 as nc
import numpy as np

data_tarball = 'test_data.tar.gz'
data_tarball_url = 'http://s3-ap-southeast-2.amazonaws.com/dp-drop/ocean-ic/test/test_data.tar.gz'

def check_output_grid(model_name, output):

    with nc.Dataset(output) as f:
        if model_name == 'MOM':
            lats = f.variables['GRID_Y_T'][:]
        else:
            assert model_name == 'NEMO'
            lats = f.variables['nav_lat'][:]

    assert np.min(lats) < -76.0
    assert np.max(lats) > 80.0

def check_output_fields(model_name, output):

    with nc.Dataset(output) as f:
        if model_name == 'MOM':
            assert f.variables['temp'].units == 'C' or \
                'celsius' in f.variables['temp'].units.lower()
            assert f.variables['salt'].units == 'psu'

            temp = f.variables['temp'][:]
            salt = f.variables['salt'][:]
        else:
            assert model_name == 'NEMO'
            temp = f.variables['votemper'][:]
            salt = f.variables['vosaline'][:]

    assert np.max(temp) < 40.0
    assert np.min(temp) > -10.0

    assert np.max(salt) < 50.0
    assert np.min(salt) > 0.0

def mom_godas_simple(which_mom, input_dir, output_dir):

    assert which_mom == 'MOM' or which_mom == 'MOM1'

    output = os.path.join(output_dir, '{}_godas_simple_ic.nc'.format(which_mom.lower()))
    if os.path.exists(output):
        os.remove(output)

    src_name = 'GODAS'
    src_temp_file = os.path.join(input_dir, 'pottmp.2001.nc')
    src_salt_file = os.path.join(input_dir, 'salt.2001.nc')
    dest_name = which_mom
    dest_data_file = output

    args = [src_name, src_temp_file, src_salt_file,
            dest_name, dest_data_file]

    my_dir = os.path.dirname(os.path.realpath(__file__))
    cmd = [os.path.join(my_dir, '../', 'makeic_simple.py')] + args
    ret = sp.call(cmd)
    assert(ret == 0)

    # Check that outputs exist.
    assert(os.path.exists(output))

    check_output_fields('MOM', output)
    check_output_grid('MOM', output)

class TestRegrid():

    @pytest.fixture
    def input_dir(self):
        test_dir = os.path.dirname(os.path.realpath(__file__))
        test_data_dir = os.path.join(test_dir, 'test_data')
        test_data_tarball = os.path.join(test_dir, data_tarball)

        if not os.path.exists(test_data_dir):
            if not os.path.exists(test_data_tarball):
                sh.wget('-P', test_dir, data_tarball_url)
            sh.tar('zxvf', test_data_tarball, '-C', test_dir)

        return os.path.join(test_data_dir, 'input')

    @pytest.fixture
    def output_dir(self):
        test_dir = os.path.dirname(os.path.realpath(__file__))
        test_data_dir = os.path.join(test_dir, 'test_data')

        return os.path.join(test_data_dir, 'output')

    @pytest.fixture
    def grid_def_dir(self):
        test_dir = os.path.dirname(os.path.realpath(__file__))
        grid_data_dir = os.path.join(test_dir, '../', 'grid_defs')

        return grid_data_dir


    @pytest.mark.godas
    def test_mom_godas_simple(self, input_dir, output_dir):
        mom_godas_simple('MOM', input_dir, output_dir)

    @pytest.mark.godas
    @pytest.mark.mom1
    def test_mom1_godas_simple(self, input_dir, output_dir):
        mom_godas_simple('MOM1', input_dir, output_dir)

    @pytest.mark.mom
    def test_mom_godas(self, input_dir, output_dir, grid_def_dir):

        output = os.path.join(output_dir, 'mom_godas_ic.nc')
        if os.path.exists(output):
            os.remove(output)

        src_name = 'GODAS'
        src_hgrid = os.path.join(grid_def_dir, 'pottmp.2016.nc')
        src_vgrid = os.path.join(grid_def_dir, 'pottmp.2016.nc')
        src_temp_file = os.path.join(input_dir, 'pottmp.2001.nc')
        src_salt_file = os.path.join(input_dir, 'salt.2001.nc')
        dest_name = 'MOM'
        dest_hgrid = os.path.join(grid_def_dir, 'ocean_hgrid.nc')
        dest_vgrid = os.path.join(grid_def_dir, 'ocean_vgrid.nc')
        mask_file = os.path.join(grid_def_dir, 'ocean_mask.nc')

        dest_data_file = output

        args = [src_name, src_hgrid, src_vgrid, src_temp_file, src_salt_file,
                dest_name, dest_hgrid, dest_vgrid, dest_data_file,
                '--model_mask', mask_file]

        my_dir = os.path.dirname(os.path.realpath(__file__))
        cmd = [os.path.join(my_dir, '../', 'makeic.py')] + args
        ret = sp.call(cmd)
        assert(ret == 0)

        # Check that outputs exist.
        assert(os.path.exists(output))

        check_output_fields('MOM', output)
        check_output_grid('MOM', output)


    @pytest.mark.grid_def
    def test_godas_change_grid_def(self, input_dir, output_dir, grid_def_dir):
        """
        Use different GODAS grid definitions by the same source files and
        expect the result to be the same.
        """

        outputA = os.path.join(output_dir, 'nemo_godas_2004_A_ic.nc')
        outputB = os.path.join(output_dir, 'nemo_godas_2004_B_ic.nc')
        if os.path.exists(outputA):
            os.remove(outputA)
        if os.path.exists(outputB):
            os.remove(outputB)

        src_name = 'GODAS'
        src_hgrid = os.path.join(grid_def_dir, 'pottmp.2016.nc')
        src_vgrid = os.path.join(grid_def_dir, 'pottmp.2016.nc')
        src_temp_file = os.path.join(input_dir, 'pottmp.2004.nc')
        src_salt_file = os.path.join(input_dir, 'salt.2004.nc')
        dest_name = 'NEMO'
        dest_hgrid = os.path.join(grid_def_dir, 'coordinates.nc')
        dest_vgrid = os.path.join(grid_def_dir,
                                    'data_1m_potential_temperature_nomask.nc')
        dest_data_file = outputA

        args = [src_name, src_hgrid, src_vgrid, src_temp_file, src_salt_file,
                dest_name, dest_hgrid, dest_vgrid, dest_data_file]

        my_dir = os.path.dirname(os.path.realpath(__file__))
        cmd = [os.path.join(my_dir, '../', 'makeic.py')] + args
        ret = sp.call(cmd)
        assert(ret == 0)

        # Now modify the src grid.
        src_hgrid = os.path.join(input_dir, 'pottmp.2004.nc')
        src_vgrid = os.path.join(input_dir, 'pottmp.2004.nc')
        dest_data_file = outputB

        args = [src_name, src_hgrid, src_vgrid, src_temp_file, src_salt_file,
                dest_name, dest_hgrid, dest_vgrid, dest_data_file]

        my_dir = os.path.dirname(os.path.realpath(__file__))
        cmd = [os.path.join(my_dir, '../', 'makeic.py')] + args
        ret = sp.call(cmd)
        assert(ret == 0)

        # Check that outputs exist.
        assert(os.path.exists(outputA))
        assert(os.path.exists(outputB))

        check_output_fields('NEMO', outputA)
        check_output_grid('NEMO', outputA)
        check_output_fields('NEMO', outputB)
        check_output_grid('NEMO', outputB)

        # Check that fields are the same.
        with nc.Dataset(outputA) as f:
            temp_A = f.variables['votemper'][:]
            salt_A = f.variables['vosaline'][:]

        with nc.Dataset(outputB) as f:
            temp_B = f.variables['votemper'][:]
            salt_B = f.variables['vosaline'][:]

            np.array_equal(temp_A, temp_B)
            np.array_equal(salt_A, salt_B)

    @pytest.mark.nemo
    def test_nemo_godas(self, input_dir, output_dir):

        output = os.path.join(output_dir, 'nemo_godas_ic.nc')
        if os.path.exists(output):
            os.remove(output)

        src_name = 'GODAS'
        src_temp_file = os.path.join(input_dir, 'pottmp.2001.nc')
        src_salt_file = os.path.join(input_dir, 'salt.2001.nc')
        dest_name = 'NEMO'
        dest_data_file = output

        args = [src_name, src_temp_file, src_salt_file,
                dest_name, dest_data_file]

        my_dir = os.path.dirname(os.path.realpath(__file__))
        cmd = [os.path.join(my_dir, '../', 'makeic_simple.py')] + args
        ret = sp.call(cmd)
        assert(ret == 0)

        # Check that outputs exist.
        assert(os.path.exists(output))

        check_output_fields('NEMO', output)
        check_output_grid('NEMO', output)

    @pytest.mark.oras4
    def test_mom_oras4(self, input_dir, output_dir):

        output = os.path.join(output_dir, 'mom_oras4_ic.nc')
        if os.path.exists(output):
            os.remove(output)

        src_name = 'ORAS4'
        src_temp_file = os.path.join(input_dir, 'thetao_oras4_1m_2014_grid_T.nc')
        src_salt_file = os.path.join(input_dir, 'so_oras4_1m_2014_grid_T.nc')
        dest_name = 'MOM'
        dest_data_file = output

        args = [src_name, src_temp_file, src_salt_file,
                dest_name, dest_data_file]

        my_dir = os.path.dirname(os.path.realpath(__file__))
        cmd = [os.path.join(my_dir, '../', 'makeic_simple.py')] + args
        ret = sp.call(cmd)
        assert(ret == 0)

        # Check that outputs exist.
        assert(os.path.exists(output))

        check_output_fields('MOM', output)
        check_output_grid('MOM', output)

    @pytest.mark.nemo
    def test_nemo_oras4(self, input_dir, output_dir):

        output = os.path.join(output_dir, 'nemo_oras4_ic.nc')
        if os.path.exists(output):
            os.remove(output)

        src_name = 'ORAS4'
        src_temp_file = os.path.join(input_dir, 'thetao_oras4_1m_2014_grid_T.nc')
        src_salt_file = os.path.join(input_dir, 'so_oras4_1m_2014_grid_T.nc')
        dest_name = 'NEMO'
        dest_data_file = output

        args = [src_name, src_temp_file, src_salt_file,
                dest_name, dest_data_file]

        my_dir = os.path.dirname(os.path.realpath(__file__))
        cmd = [os.path.join(my_dir, '../', 'makeic_simple.py')] + args
        ret = sp.call(cmd)
        assert(ret == 0)

        # Check that outputs exist.
        assert(os.path.exists(output))

        check_output_fields('NEMO', output)
        check_output_grid('NEMO', output)
