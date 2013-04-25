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
#ifndef LSST_MEAS_EXTENSIONS_simpleShape_h_INCLUDED
#define LSST_MEAS_EXTENSIONS_simpleShape_h_INCLUDED

#include "lsst/pex/config.h"
#include "lsst/afw/geom/ellipses.h"
#include "lsst/meas/algorithms/ShapeControl.h"

namespace lsst { namespace meas { namespace extensions { namespace simpleShape {

class SimpleShapeControl : public algorithms::ShapeControl {
public:

    LSST_CONTROL_FIELD(sigma, double, "Sigma of circular Gaussian used as weight function, in pixels");
    LSST_CONTROL_FIELD(nSigmaRegion, double, "Maximum radius for pixels to include, in units of sigma");
    LSST_CONTROL_FIELD(useApproximateExp, bool, "Whether to use fast approximate exponential function");

    SimpleShapeControl() :
        algorithms::ShapeControl("shape.simple"), sigma(1.5), nSigmaRegion(3), useApproximateExp(true) {}

private:
    virtual PTR(algorithms::AlgorithmControl) _clone() const;
    virtual PTR(algorithms::Algorithm) _makeAlgorithm(
        afw::table::Schema & schema, PTR(daf::base::PropertyList) const & metadata
    ) const;
};

}}}} // namespace lsst::meas::extensions::simpleShape

#endif // !LSST_MEAS_EXTENSIONS_simpleShape_h_INCLUDED
