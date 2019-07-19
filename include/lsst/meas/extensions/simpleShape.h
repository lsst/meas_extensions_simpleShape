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

#include <bitset>

#include "lsst/pex/config.h"
#include "lsst/geom.h"
#include "lsst/afw/geom/ellipses.h"
#include "lsst/afw/table.h"
#include "lsst/afw/image.h"
#include "lsst/meas/base/Algorithm.h"
#include "lsst/meas/base.h"

namespace lsst { namespace meas { namespace extensions { namespace simpleShape {

class SimpleShapeResult;

/**
 *  @brief A C++ control class to handle SdssShapeAlgorithm's configuration
 *
 */
class SimpleShapeControl {
 public:
    LSST_CONTROL_FIELD(sigma, double, "Sigma of circular Gaussian used as weight function, in pixels");
    LSST_CONTROL_FIELD(nSigmaRegion, double, "Maximum radius for pixels to include, in units of sigma");

    SimpleShapeControl() : sigma(1.5), nSigmaRegion(3) {}
};

class SimpleShapeResultKey : public afw::table::FunctorKey<SimpleShapeResult> {
 public:
    static SimpleShapeResultKey addFields(
            afw::table::Schema & schema,
            std::string const & name
    );

    // Default constructor
    SimpleShapeResultKey() {}

     /**
     *  @brief Construct from a subschema, assuming _xx, _yy, etc. subfields
     *
     *  If a schema has "a_xx", "a_yy", etc. fields, this constructor allows you to construct
     *  a SimpleShapeResultKey via:
     *  @code
     *  SimpleShapeResultKey k(schema["a"]);
     *  @endcode
     */
    SimpleShapeResultKey(lsst::afw::table::SubSchema const & s);

    virtual SimpleShapeResult get(afw::table::BaseRecord  const & record) const;
    virtual void set(afw::table::BaseRecord & record, SimpleShapeResult const & value) const;

    //@{
    /// Compare the FunctorKey for equality with another, using the underlying Keys
    bool operator==(SimpleShapeResultKey const & other) const;
    bool operator!=(SimpleShapeResultKey const & other) const { return !(*this == other); }
    //@}

    /// Return True if the key is valid
    bool isValid() const;

    lsst::meas::base::FlagHandler const & getFlagHandler() const { return _flagHandler; }

 private:
    lsst::afw::table::QuadrupoleKey _shapeResult;
    lsst::afw::table::Point2DKey _centroidResult;
    lsst::afw::table::CovarianceMatrixKey<double, 5 > _uncertantyResult;
    lsst::meas::base::FlagHandler _flagHandler;
};


class SimpleShape : public lsst::meas::base::SimpleAlgorithm {
public:

    // Structures and routines to manage flaghandler
    static base::FlagDefinitionList const & getFlagDefinitions();
    static unsigned int const N_FLAGS = 1;
    static base::FlagDefinition const FAILURE;

    typedef SimpleShapeControl Control;

    SimpleShape(Control const & ctrl, std::string const & name, afw::table::Schema & schema);

    /**
     * Compute the Gaussian-weighted moments of an image.
     *
     * @param[in] weight        An ellipse object of Gaussian weights to apply to
     *                          the measurement.
     * @param[in] image         A Masked image instance with int float or double
     *                          pixels.
     * @param[in] nSigmaRegion  Maximum radius for pixels to include, in units
     *                          of sigma
     */
    template <typename T>
    static SimpleShapeResult computeMoments(
        afw::geom::ellipses::Ellipse const & weight,
        afw::image::MaskedImage<T> const & image,
        double nSigmaRegion=3.0
    );

    /**
     *  @brief Convert linear raw moments into an ellipse and centroid, and return the derivative
     *         of the conversion.
     *
     *  @note This function is mainly intended for internal use, and is only exposed publically
     *        so it can be unit-tested in Python.
     *
     *  For weight function @f$w@f$ and data @f$p@f$, the "raw" moments @f$Q@f$ are defined as:
     *  @f{eqnarray*}{
     *    Q_0 &=& \sum_n w(x_n, y_n) p_n \\
     *    Q_{xx} &=& \sum_n w(x_n, y_n) x_n^2 p_n \\
     *    Q_{yy} &=& \sum_n w(x_n, y_n) y_n^2 p_n \\
     *    Q_{xy} &=& \sum_n w(x_n, y_n) x_n y_n p_n \\
     *    Q_x &=& \sum_n w(x_n, y_n) x_n p_n \\
     *    Q_y &=& \sum_n w(x_n, y_n) y_n p_n
     *  @f}
     *  whereas the converted ellipse and centroid moments are:
     *  @f{eqnarray*}{
     *    M_{xx} &=& Q_{xx} / Q_0 - Q_x^2 \\
     *    M_{xx} &=& Q_{yy} / Q_0 - Q_y^2 \\
     *    M_{xx} &=& Q_{xy} / Q_0 - Q_x Q_y \\
     *    M_x &=& Q_x / Q_0 \\
     *    M_y &=& Q_y / Q_0
     *  @f}
     *
     *  Note the slightly unusual ordering; this is for consistency with afw::geom::ellipses::Ellipse.
     */
    static Eigen::Matrix<double,5,6> convertRawMoments(
        Eigen::Matrix<double,6,1> const & q,
        afw::geom::ellipses::Quadrupole & quadrupole,
        geom::Point2D & center
    );

    /**
     *  @brief Correct moments measured with a Gaussian weight function by assuming the data was also
     *         an elliptical Gaussian, and return the derivative of the correction.
     *
     *  @note This function is mainly intended for internal use, and is only exposed publically
     *        so it can be unit-tested in Python.
     *
     *  If we naively measure Gaussian-weighted moments, we'll measure the moments of the product
     *  of the weight function and the data.  What we want is the moments of the data, as if we
     *  had measured them with no weight function (but without sacrificing the S/N benefit that
     *  comes from using a weight function).  To do that, we assume the data is also an elliptical
     *  Gaussian, and "divide" the weight function from the measured moments to compute it.
     *
     *  If @f$W@f$ and @f$M@f$ are the quadruple matrices of the weight function and measurement,
     *  and @f$\eta@f$ is the measured centroid (we work in a coordinate system where the weight
     *  function is centered at the origin), then the corrected quadrupole matrix @f$C@f$ and
     *  centroid are @f$\nu@f$ are:
     *  @f{eqnarray*}{
     *    C &=& \left(M^{-1} - W^{-1}\right)^{-1} \\
     *    \nu &=& C M^{-1} \eta
     *  @f}
     */
    static Eigen::Matrix<double, 5, 5> correctWeightedMoments(
        afw::geom::ellipses::Quadrupole const & weight,
        afw::geom::ellipses::Quadrupole & ellipse,
        geom::Point2D & center
    );

    virtual void measure(
            afw::table::SourceRecord & measRecord,
            afw::image::Exposure<float> const & exposure
            ) const;

    virtual void fail(
        afw::table::SourceRecord & measRecord,
        lsst::meas::base::MeasurementError * error=NULL
    ) const;

 private:
    Control _ctrl;
    SimpleShapeResultKey _resultKey;
    lsst::meas::base::SafeCentroidExtractor _centroidExtractor;
};

/**
 *  Struct to hold the results of SimpleShape when we don't run it as a plugin.
 */
class SimpleShapeResult {
 public:
    afw::geom::ellipses::Quadrupole ellipse; ///< Measured second moments.
    geom::Point2D center; ///< Measured first moments, or the input center if !recentroid
    Eigen::Matrix<double,5,5> covariance; ///< Matrix of uncertainties; ordered Ixx, Iyy, Ixy, Ix, Iy.

#ifndef SWIG
    std::bitset<SimpleShape::N_FLAGS> flags;
#endif

    // Flag getter for Swig which doesn't understand std::bitset
    bool getFlag(int index) const { return flags[index]; }

    SimpleShapeResult(); ///< Constructor; initializes everything to Nan
};


}}}} // namespace lsst::meas::extensions::simpleShape

#endif // !LSST_MEAS_EXTENSIONS_simpleShape_h_INCLUDED
