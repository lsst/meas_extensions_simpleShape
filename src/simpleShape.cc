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

#include "lsst/pex/exceptions.h"
#include "lsst/meas/algorithms/ImageMoments.h"
#include "lsst/meas/extensions/simpleShape.h"

namespace lsst { namespace meas { namespace extensions { namespace simpleShape {

/*
 * All of the actual work is done by meas::algorithms::ImageMoments; this file just involves
 * meeting the requirements of the pluggable algorithm framework:
 *
 *  - Construction via the control object's makeAlgorithm method and the algorithm constructor
 *
 *  - Registration of the Schema fields the algorithm will fill, while saving the Key objects
 *    that will be used to refer to those fields later.  Most of this is done by the ShapeAlgorithm
 *    base class constructor, but we optionally also allocate fields to save a centroid.
 *
 *  - Calling the algorithmic code in measure() from the all-important _apply() member function
 *
 *  - Instantiating templates via the LSST_MEAS_ALGORITHM_PRIVATE_IMPLEMENTATION macro
 */

namespace {

class SimpleShape : public algorithms::ShapeAlgorithm {
public:

    SimpleShape(SimpleShapeControl const & ctrl, afw::table::Schema & schema);

private:

    template <typename PixelT>
    void _apply(
        afw::table::SourceRecord & source,
        afw::image::Exposure<PixelT> const & exposure,
        afw::geom::Point2D const & center
    ) const;

    LSST_MEAS_ALGORITHM_PRIVATE_INTERFACE(SimpleShape);

    afw::table::KeyTuple<afw::table::Centroid> _centroidKeys;
};

SimpleShape::SimpleShape(SimpleShapeControl const & ctrl, afw::table::Schema & schema) :
    algorithms::ShapeAlgorithm(ctrl, schema, "shape measured with fixed (non-adaptive) Gaussian weight"),
    _centroidKeys(
        addCentroidFields( // addCentroidFields is in afw::table, via Koenig lookup on schema
            schema, ctrl.name + ".centroid", "centroid measured with fixed circular Gaussian weight"
        )
    )
{}

template<typename PixelT>
void SimpleShape::_apply(
    afw::table::SourceRecord & source,
    afw::image::Exposure<PixelT> const & exposure,
    afw::geom::Point2D const & center
) const {
    SimpleShapeControl const & ctrl = static_cast<SimpleShapeControl const &>(this->getControl());
    // set flags so an exception throw produces a flagged source
    source.set(this->getKeys().flag, true);
    source.set(_centroidKeys.flag, true);
    algorithms::ImageMoments moments(ctrl.nSigmaRegion, ctrl.useApproximateExp);
    algorithms::ImageMoments::EllipseResult result = moments.measureEllipse(
        exposure.getMaskedImage(),
        afw::geom::ellipses::Axes(ctrl.sigma, ctrl.sigma),
        center
    );
    source.set(this->getKeys().meas, result.getQuadrupole());
    source.set(this->getKeys().err, result.covariance->block<3,3>(0,0));
    source.set(this->getKeys().flag, false);
    source.set(_centroidKeys.meas, result.getCentroid() + afw::geom::Extent2D(center));
    source.set(_centroidKeys.err, result.covariance->block<2,2>(3,3));
    source.set(_centroidKeys.flag, false);
}

LSST_MEAS_ALGORITHM_PRIVATE_IMPLEMENTATION(SimpleShape);

} // anonymous

PTR(algorithms::AlgorithmControl) SimpleShapeControl::_clone() const {
    return boost::make_shared<SimpleShapeControl>(*this);
}

PTR(algorithms::Algorithm) SimpleShapeControl::_makeAlgorithm(
    afw::table::Schema & schema,
    PTR(daf::base::PropertyList) const &
) const {
    return boost::make_shared<SimpleShape>(*this, boost::ref(schema));
}

}}}} // namespace lsst::meas::extensions::simpleShape
