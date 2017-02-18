#!/usr/bin/env python
#
# LSST Data Management System
#
# Copyright 2008-2016  AURA/LSST.
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <https://www.lsstcorp.org/LegalNotices/>.
#
import unittest

import lsst.utils.tests
import lsst.meas.base.tests
from lsst.meas.base.tests import AlgorithmTestCase
from lsst.meas.extensions.simpleShape import SimpleShapeResultKey


class SimpleShapeMFTestCase(AlgorithmTestCase, lsst.utils.tests.TestCase):

    def setUp(self):
        self.bbox = lsst.afw.geom.Box2I(lsst.afw.geom.Point2I(-20, -30),
                                        lsst.afw.geom.Extent2I(240, 160))
        self.dataset = lsst.meas.base.tests.TestDataset(self.bbox)
        self.dataset.addSource(100000.0, lsst.afw.geom.Point2D(149.9, 50.3),
                               lsst.afw.geom.ellipses.Quadrupole(8, 9, 3))

        self.expectedKeySet = set(['ext_simpleShape_SimpleShape_IxSigma',
                                   'ext_simpleShape_SimpleShape_Ix_Iy_Cov',
                                   'ext_simpleShape_SimpleShape_IxxSigma',
                                   'ext_simpleShape_SimpleShape_Ixx_Ix_Cov',
                                   'ext_simpleShape_SimpleShape_Ixx_Ixy_Cov',
                                   'ext_simpleShape_SimpleShape_Ixx_Iy_Cov',
                                   'ext_simpleShape_SimpleShape_Ixx_Iyy_Cov',
                                   'ext_simpleShape_SimpleShape_IxySigma',
                                   'ext_simpleShape_SimpleShape_Ixy_Ix_Cov',
                                   'ext_simpleShape_SimpleShape_Ixy_Iy_Cov',
                                   'ext_simpleShape_SimpleShape_IySigma',
                                   'ext_simpleShape_SimpleShape_IyySigma',
                                   'ext_simpleShape_SimpleShape_Iyy_Ix_Cov',
                                   'ext_simpleShape_SimpleShape_Iyy_Ixy_Cov',
                                   'ext_simpleShape_SimpleShape_Iyy_Iy_Cov',
                                   'ext_simpleShape_SimpleShape_flag',
                                   'ext_simpleShape_SimpleShape_x',
                                   'ext_simpleShape_SimpleShape_xx',
                                   'ext_simpleShape_SimpleShape_xy',
                                   'ext_simpleShape_SimpleShape_y',
                                   'ext_simpleShape_SimpleShape_yy'])

    def tearDown(self):
        del self.bbox
        del self.dataset

    def testKeys(self):
        """
        Test that Simple shape will run as a plugin, a functor key can be
        created, and result object can be created. Testing of the algorithm
        is left to testSimpleShape.py
        """

        task = self.makeSingleFrameMeasurementTask("ext_simpleShape_SimpleShape")
        exposure, catalog = self.dataset.realize(10.0, task.schema)
        # Check the keys are in the schema
        foundKeySet = set(catalog.schema.getNames()) & self.expectedKeySet
        missingKeySet = self.expectedKeySet - foundKeySet
        self.assertEqual(missingKeySet, set())

        # check that SimpleShapeResultKey can be created and that it returns a results object
        key = SimpleShapeResultKey(catalog.schema["ext_simpleShape_SimpleShape"])
        task.run(catalog, exposure)

        result = catalog[0].get(key)
        for attr in ("center", "ellipse", "covariance", "getFlag"):
            self.assertTrue(hasattr(result, attr), "Result Object missing {}".format(attr))


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()

if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
