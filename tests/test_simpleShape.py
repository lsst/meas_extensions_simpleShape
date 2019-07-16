#
# LSST Data Management System
#
# Copyright 2008-2017  AURA/LSST.
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

import numpy as np

import lsst.utils.tests
import lsst.geom
import lsst.afw.geom.ellipses as el
import lsst.afw.image
from lsst.meas.extensions.simpleShape import SimpleShape


class SimpleShapeTestCase(lsst.utils.tests.TestCase):

    def setUp(self):
        self.ellipseCores = [
            el.Quadrupole(25.0, 25.0, 0.0),
            el.Quadrupole(27.0, 22.0, -5.0),
            el.Quadrupole(23.0, 28.0, 2.0),
        ]
        self.centers = [
            lsst.geom.Point2D(0.0, 0.0),
            lsst.geom.Point2D(2.0, 3.0),
            lsst.geom.Point2D(-1.0, 2.5),
        ]
        for ellipseCore in self.ellipseCores:
            ellipseCore.scale(2)
        self.bbox = lsst.geom.Box2I(lsst.geom.Point2I(-500, -500), lsst.geom.Point2I(50, 50))
        self.xg, self.yg = np.meshgrid(
            np.arange(self.bbox.getBeginX(), self.bbox.getEndX(), dtype=float),
            np.arange(self.bbox.getBeginY(), self.bbox.getEndY(), dtype=float)
        )

    def evaluateGaussian(self, ellipse):
        '''
        Create an elliptical Guassian as a Numpy array. Does not normalize the Gaussian
        '''
        gt = ellipse.getGridTransform()
        xt = gt[gt.XX] * self.xg + gt[gt.XY] * self.yg + gt[gt.X]
        yt = gt[gt.YX] * self.xg + gt[gt.YY] * self.yg + gt[gt.Y]
        return np.exp(-0.5 * (xt**2 + yt**2))

    def buildImageAndMoments(self, dEllipseCore, dCenter, wEllipseCore, wCenter):
        '''
        Generates a elliptical Gaussian image from input ellipse (dEllipse)
        and uses an elliptical Gaussian (wEllipse) to calculate the shape
        of the generated image.
        '''
        dEllipse = el.Ellipse(dEllipseCore, dCenter)
        image = lsst.afw.image.MaskedImageF(self.bbox)
        image.getImage().getArray()[:, :] = self.evaluateGaussian(dEllipse)
        wEllipse = el.Ellipse(wEllipseCore, wCenter)
        result = SimpleShape.computeMoments(wEllipse, image)
        return result

    def testCorrectWeightedMoments(self):
        '''
        Test that the measured moments can be corrected for the fact that the measurement
        contains information on the moments of the weight function used to make the
        measurement. The results should only contain the objects shape and not the moments
        used to make the measurement. This is a subset of the functionality in testNoNoiseGaussians.
        Because the other test tests a broader scope, it is useful to have this test happen under
        more restrictive conditions.
        '''
        for dEllipseCore in self.ellipseCores:
            for dCenter in self.centers:
                dEllipse = el.Ellipse(dEllipseCore, dCenter)
                dArray = self.evaluateGaussian(dEllipse)
                for wEllipseCore in self.ellipseCores:
                    for wCenter in self.centers:
                        wEllipse = el.Ellipse(wEllipseCore, wCenter)
                        wArray = self.evaluateGaussian(wEllipse)
                        product = dArray * wArray
                        i0 = np.sum(product)
                        ix = np.sum(product * self.xg) / i0
                        iy = np.sum(product * self.yg) / i0
                        ixx = np.sum(product * (self.xg - ix)**2) / i0
                        iyy = np.sum(product * (self.yg - iy)**2) / i0
                        ixy = np.sum(product * (self.xg - ix) * (self.yg - iy)) / i0
                        mEllipseCore = el.Quadrupole(ixx, iyy, ixy)
                        mCenter = lsst.geom.Point2D(ix, iy)
                        SimpleShape.correctWeightedMoments(wEllipseCore, mEllipseCore, mCenter)
                        self.assertFloatsAlmostEqual(mEllipseCore.getParameterVector(),
                                                     dEllipseCore.getParameterVector(),
                                                     rtol=1E-8, atol=1E-11)

    def testNoNoiseGaussians(self):
        '''
        Test that shape moments can be actuately determined for Gaussian images with no noise
        '''
        for ellipseCore in self.ellipseCores:
            for center in self.centers:
                result = self.buildImageAndMoments(ellipseCore, center, ellipseCore, center)
                self.assertFloatsAlmostEqual(result.ellipse.getParameterVector(),
                                             ellipseCore.getParameterVector(),
                                             rtol=3E-3, atol=1E-15)
                self.assertFloatsAlmostEqual(np.array(result.center),
                                             np.array(center),
                                             rtol=1E-8, atol=1E-15)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
