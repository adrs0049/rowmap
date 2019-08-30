#!/usr/bin/python
# -*- coding: utf-8 -*-
# Author: Andreas Buttenschoen

import h5py as h5


class MOLFile:
    def __init__(self, shape, eqns, fname=None, name='MOL_data'):
        self.shape = shape
        self.eqns  = eqns

        # filename to write to
        self.fname = fname

        # name of the object
        self.name = name

        # groups
        self.groups = {}

        # file handler
        self.h5f = None

        # field name
        self.field_name = 'data'

        # setup output
        self._setup_h5_output()


    def _setup_h5_output(self):
        self.__open_file(mode='w')
        for i in range(self.eqns):
            grp_name = self.name + '_%d' % i
            grp = self.h5f.create_group(grp_name)
            self.groups[grp.name] = i

        # add a group for timing information
        self.h5f.create_dataset('time', (0, ), maxshape=(None, ), chunks=True)

        # create dataset
        self.new_dset()

        self.__close_file()


    def __open_file(self, mode='a'):
        if self.h5f is None:
            self.h5f = h5.File(self.fname, mode)


    def __close_file(self):
        if self.h5f is not None:
            self.h5f.close()
            self.h5f = None


    def new_dset(self):
        self.__open_file()

        for grp in self.h5f.values():
            if isinstance(grp, h5.Dataset):
                continue

            gid = self.groups[grp.name]
            dset_shape = (0, ) + tuple(self.shape)
            dset_max   = (None, ) + tuple(self.shape)
            grp.create_dataset(self.field_name, dset_shape,
                               maxshape=dset_max, dtype='f', chunks=True)

            print('Found group %s with id %d.\nCreating dataset with shape: %s.' % (grp.name, gid, dset_shape))


    def add(self, t, y):
        # write directly to the file
        self.__open_file()
        try:
            tdset = self.h5f['time']
            tdset.resize(len(tdset)+1, axis=0)
            tdset[len(tdset)-1] = t

            for grp in self.h5f.values():
                if isinstance(grp, h5.Dataset):
                    continue

                gid  = self.groups[grp.name]
                dset = grp[self.field_name]
                shape = dset.shape

                cidx = shape[0]
                dset.resize(cidx+1, axis=0)

                # write data to the file
                dset[cidx, :] = y[gid, :]

        except Exception as e:
            print('MOLFile Error:', e)
            self.__close_file()
        else:
            self.__close_file()


    def finish(self):
        pass


    def __del__(self):
        self.__close_file()

