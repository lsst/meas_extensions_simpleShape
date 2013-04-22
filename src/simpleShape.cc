// -*- lsst-c++ -*-
/*
 * LSST Data Management System
 * Copyright 2008-2013 LSST Corporation.
 *
 * This product includes software developed by the
 * LSST Project (http://www.lsst.org/).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the LSST License Statement and
 * the GNU General Public License along with this program.  If not,
 * see <http://www.lsstcorp.org/LegalNotices/>.
 */

#include "lsst/meas/extensions/simpleShape.h"

namespace lsst { namespace meas { namespace extensions { namespace simpleShape {

SimpleShape::SimpleShape(SimpleShapeControl const & ctrl, afw::table::Schema & schema) :
    algorithms::ShapeAlgorithm(ctrl, schema, "shape measured with fixed (non-adaptive) Gaussian weight")
{
    if (ctrl.recentroid) {
        _centroidKeys = addCentroidFields( // addCentroidFields is in afw::table, via Koenig lookup on schema
            schema, ctrl.name + ".centroid", "centroid measured with fixed circular Gaussian weight"
        );
    }
}

template<typename PixelT>
void SimpleShape::_apply(
    afw::table::SourceRecord & source,
    afw::image::Exposure<PixelT> const & exposure,
    afw::geom::Point2D const & center
) const {
    SimpleShapeControl const & ctrl = static_cast<SimpleShapeControl const &>(this->getControl());
    // set flags so an exception throw produces a flagged source
    source.set(this->getKeys().flag, true);
    if (ctrl.recentroid) {
        source.set(_centroidKeys.flag, true);
    }
    afw::geom::ellipses::Ellipse weight(afw::geom::ellipses::Axes(ctrl.sigma), center);
    Result result = measure(weight, exposure.getMaskedImage(), ctrl.nSigmaRegion, ctrl.recentroid);
    source.set(this->getKeys().meas, result.ellipse);
    source.set(this->getKeys().err, result.covariance.block<3,3>(0,0));
    source.set(this->getKeys().flag, false);
    if (ctrl.recentroid) {
        source.set(_centroidKeys.meas, result.center);
        source.set(_centroidKeys.err, result.covariance.block<2,2>(3,3));
        source.set(_centroidKeys.flag, false);
    }
}

PTR(algorithms::AlgorithmControl) SimpleShapeControl::_clone() const {
    return boost::make_shared<SimpleShapeControl>(*this);
}

PTR(algorithms::Algorithm) SimpleShapeControl::_makeAlgorithm(
    afw::table::Schema & schema,
    PTR(daf::base::PropertyList) const &
) const {
    return boost::make_shared<SimpleShape>(*this, boost::ref(schema));
}

LSST_MEAS_ALGORITHM_PRIVATE_IMPLEMENTATION(SimpleShape);

}}}} // namespace lsst::meas::extensions::simpleShape
