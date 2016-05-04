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

#include <array>

#include "lsst/pex/exceptions.h"
#include "lsst/afw/geom/ellipses/PixelRegion.h"
#include "lsst/meas/extensions/simpleShape.h"
#include "lsst/afw/math.h"

namespace lsst { namespace meas { namespace extensions { namespace simpleShape {

/* ---- Machinery for iterating over elliptical regions ----------------------------------------------------
 *
 * In contrast to FootprintFunctor (which is what's used by most other measurement algorithms), this
 * machinery doesn't require a virtual function call at each pixel, and hence allows the functor to
 * be inlined.
 * Eventually, something like this will replace FootprintFunctor (see #1836)
 */

namespace {

using afw::geom::Span;

Span clipSpan(Span const & span, afw::geom::Box2I const & box) {
    if (span.getY() < box.getMinY() || span.getY() > box.getMaxY()) return Span();
    return Span(span.getY(),
                std::min(std::max(span.getMinX(), box.getMinX()), box.getMaxX()),
                std::max(std::min(span.getMaxX(), box.getMaxX()), box.getMinX())
        );
}

template <typename Function, typename Iterator>
void iterateSpan(Function function, Iterator pixIter, Span const & span) {
    for (
        Span::Iterator pointIter = span.begin(), pointEnd = span.end();
        pointIter != pointEnd;
        ++pointIter, ++pixIter
    ) {
        boost::unwrap_ref(function)(*pointIter, *pixIter);
    }
}

template <typename Function, typename Image, typename Region>
void iterateRegion(Function function, Image const & image, Region const & region) {
    afw::geom::Box2I bbox = image.getBBox(afw::image::PARENT);
    if (bbox.contains(region.getBBox())) {
        // if the box contains the region, there's no need to check each span to make sure it's entirely
        // within the image
        for (
            typename Region::Iterator spanIter = region.begin(), spanEnd = region.end();
            spanIter != spanEnd;
            ++spanIter
        ) {
            iterateSpan(
                function,
                image.x_at(spanIter->getMinX() - image.getX0(), spanIter->getY() - image.getY0()),
                *spanIter
            );
        }
    } else {
        for (
            typename Region::Iterator spanIter = region.begin(), spanEnd = region.end();
            spanIter != spanEnd;
            ++spanIter
        ) {
            Span span = clipSpan(*spanIter, bbox);
            iterateSpan(
                function,
                image.x_at(span.getMinX() - image.getX0(), span.getY() - image.getY0()),
                span
            );
        }
    }
}

} // anonymous

/* ---- Implementation for SimpleShape algorithm ------------------------------------------------------------
 *
 * All of the pixel-based work is done by an accumulating functor that we invoke using the above machinery.
 * The rest is in the static member function measure(), which provides a way to call the algorithm
 * outside the context of the pluggable source measurement algorithm framework.
 */

namespace {

// Enums used to index arrays for code clarity
enum { Q0=0, QXX, QYY, QXY, QX, QY, N_Q };
enum { MXX=0, MYY, MXY, MX, MY, N_M };

typedef Eigen::Matrix<double,N_Q,1> VectorQ;
typedef Eigen::Matrix<double,N_Q,N_Q> MatrixQ;
typedef Eigen::Matrix<double,N_M,1> VectorM;
typedef Eigen::Matrix<double,N_M,N_M> MatrixM;
typedef Eigen::Matrix<double,N_M,N_Q> MatrixMQ;

struct RawMomentAccumulator {

    template <typename PixelT>
    void operator()(afw::geom::Point2I const & pos, PixelT const & pixel) {
        afw::geom::Extent2D d = afw::geom::Point2D(pos) - _center;
        afw::geom::Extent2D gtd = _gt(d);
        double w = std::exp(-0.5 * (gtd.getX()*gtd.getX() + gtd.getY()*gtd.getY()));
        VectorQ q = VectorQ::Constant(w);
        q[QXX] *= d.getX() * d.getX();
        q[QYY] *= d.getY() * d.getY();
        q[QXY] *= d.getX() * d.getY();
        q[QX] *= d.getX();
        q[QY] *= d.getY();
        moments += q * pixel.image();
        covariance.selfadjointView<Eigen::Lower>().rankUpdate(q, pixel.variance());
    }

    explicit RawMomentAccumulator(afw::geom::ellipses::Ellipse const & weightEllipse) :
        moments(VectorQ::Zero()),
        covariance(MatrixQ::Zero()),
        _center(weightEllipse.getCenter()),
        _gt(weightEllipse.getCore().getGridTransform())
    {}

    VectorQ moments;
    MatrixQ covariance;
private:
    afw::geom::Point2D _center;
    afw::geom::LinearTransform _gt;
};

} // anonymous

template <typename T>
SimpleShapeResult SimpleShape::computeMoments(
    afw::geom::ellipses::Ellipse const & weight,
    afw::image::MaskedImage<T> const & image,
    double nSigmaRegion
) {
    // Define the pixel region we want to accumulate over.
    afw::geom::ellipses::Ellipse regionEllipse(weight);
    regionEllipse.getCore().scale(nSigmaRegion);
    afw::geom::ellipses::PixelRegion region(regionEllipse);
    SimpleShapeResult result;

    // We work in the coordinate system where the origin is the center of the weight function.
    // This will make things easier in the various transformations we'll have to apply to the raw
    // moments.

    // All the pixel operations take place in the next two lines, when we use the above machinery
    // to accumulate the raw moments in a single pass.
    RawMomentAccumulator functor(weight);
    iterateRegion(boost::ref(functor), image, region);

    // Then we convert the raw moments to an ellipse and centroid, and propagate the uncertainty
    // through the covariance matrix.
    MatrixMQ dm_dq = convertRawMoments(functor.moments, result.ellipse, result.center);
    MatrixM mCov = dm_dq * functor.covariance.selfadjointView<Eigen::Lower>() * dm_dq.adjoint();

    // Next, we correct the measured moments for the fact that what we just measured was the
    // moments of the product of the weight function and the data, and we want just the moments
    // of the data.
    MatrixM dc_dm = correctWeightedMoments(
        weight.getCore().as<afw::geom::ellipses::Quadrupole>(),
        result.ellipse,
        result.center
    );
    result.covariance = dc_dm * mCov.selfadjointView<Eigen::Lower>() * dc_dm.adjoint();

    // Finally, we switch back to the native image coordinate system.
    result.center += afw::geom::Extent2D(weight.getCenter());

    return result;
}


MatrixMQ SimpleShape::convertRawMoments(
    VectorQ const & q,
    afw::geom::ellipses::Quadrupole & ellipse,
    afw::geom::Point2D & center
) {
    VectorM m = q.segment<N_M>(1) / q[Q0];
    MatrixMQ dm_dq = MatrixMQ::Zero();;
    dm_dq.block<N_M,N_M>(0,1) = MatrixM::Identity() / q[Q0];
    dm_dq.col(Q0) = -m / q[Q0];
    // now we account for the centroid^2 terms in the second moments
    m[MXX] -= m[MX] * m[MX];
    m[MYY] -= m[MY] * m[MY];
    m[MXY] -= m[MX] * m[MY];
    dm_dq(MXX, QX) = -2.0 * m[MX] / q[Q0];
    dm_dq(MYY, QY) = -2.0 * m[MY] / q[Q0];
    dm_dq(MXY, QX) = m[MY] / q[Q0];
    dm_dq(MXY, QY) = m[MX] / q[Q0];
    double tmp = 2.0 / (q[Q0] * q[Q0] * q[Q0]); // == d(-Q_0^{-2})/d(Q_0)
    dm_dq(MXX, Q0) += tmp * q[QX] * q[QX];
    dm_dq(MYY, Q0) += tmp * q[QY] * q[QY];
    dm_dq(MXY, Q0) += tmp * q[QX] * q[QY];

    ellipse.setIxx(m[MXX]);
    ellipse.setIyy(m[MYY]);
    ellipse.setIxy(m[MXY]);
    center.setX(m[MX]);
    center.setY(m[MY]);

    return dm_dq;
}

MatrixM SimpleShape::correctWeightedMoments(
    afw::geom::ellipses::Quadrupole const & weight,
    afw::geom::ellipses::Quadrupole & ellipse,
    afw::geom::Point2D & center
) {
    Eigen::Matrix2d wMat = weight.getMatrix();
    Eigen::Vector2d mVec = center.asEigen();
    Eigen::Matrix2d mMat = ellipse.getMatrix();
    if (wMat.determinant() <= 0.0) {
        throw LSST_EXCEPT(
            pex::exceptions::RuntimeError,
            "Weight moments matrix is singular"
        );
    }
    if (mMat.determinant() <= 0.0) {
        throw LSST_EXCEPT(
            pex::exceptions::RuntimeError,
            "Measured moments matrix is singular"
        );
    }
    Eigen::Matrix2d mInv = mMat.inverse();
    Eigen::Matrix2d cInv = mInv - wMat.inverse();
    if (cInv.determinant() <= 0.0) {
        throw LSST_EXCEPT(
            pex::exceptions::RuntimeError,
            "Corrected moments matrix is singular"
        );
    }
    Eigen::Matrix2d cMat = cInv.inverse();
    ellipse.setIxx(cMat(0, 0));
    ellipse.setIyy(cMat(1, 1));
    ellipse.setIxy(cMat(0, 1));
    Eigen::Matrix2d cMat_mInv = cMat * mInv;
    Eigen::Vector2d mInv_mVec = mInv * mVec;
    Eigen::Vector2d cVec = cMat_mInv * mVec;
    center.setX(cVec[0]);
    center.setY(cVec[1]);
    Eigen::Matrix2d dcMat_dmxx = cMat_mInv.col(0) * cMat_mInv.col(0).adjoint();
    Eigen::Matrix2d dcMat_dmyy = cMat_mInv.col(1) * cMat_mInv.col(1).adjoint();
    Eigen::Matrix2d dcMat_dmxy = cMat_mInv.col(0) * cMat_mInv.col(1).adjoint()
        + cMat_mInv.col(1) * cMat_mInv.col(0).adjoint();
    Eigen::Vector2d dcVec_dmxx = cMat_mInv.col(0) * mInv_mVec[0];
    Eigen::Vector2d dcVec_dmyy = cMat_mInv.col(1) * mInv_mVec[1];
    Eigen::Vector2d dcVec_dmxy = cMat_mInv.col(0) * mInv_mVec[1] + cMat_mInv.col(1) * mInv_mVec[0];
    Eigen::Matrix2d const & dcVec_dmVec = cMat_mInv;

    // calculations done - now we just need to put these elements back into a 5x5 matrix to return

    MatrixM dc_dm = MatrixM::Zero();
    dc_dm(MXX, MXX) = dcMat_dmxx(0, 0);
    dc_dm(MYY, MXX) = dcMat_dmxx(1, 1);
    dc_dm(MXY, MXX) = dcMat_dmxx(0, 1);
    dc_dm(MXX, MYY) = dcMat_dmyy(0, 0);
    dc_dm(MYY, MYY) = dcMat_dmyy(1, 1);
    dc_dm(MXY, MYY) = dcMat_dmyy(0, 1);
    dc_dm(MXX, MXY) = dcMat_dmxy(0, 0);
    dc_dm(MYY, MXY) = dcMat_dmxy(1, 1);
    dc_dm(MXY, MXY) = dcMat_dmxy(0, 1);
    dc_dm(MX, MXX) = dcVec_dmxx[0];
    dc_dm(MY, MXX) = dcVec_dmxx[1];
    dc_dm(MX, MYY) = dcVec_dmyy[0];
    dc_dm(MY, MYY) = dcVec_dmyy[1];
    dc_dm(MX, MXY) = dcVec_dmxy[0];
    dc_dm(MY, MXY) = dcVec_dmxy[1];
    dc_dm(MX, MX) = dcVec_dmVec(0, 0);
    dc_dm(MX, MY) = dcVec_dmVec(0, 1);
    dc_dm(MY, MX) = dcVec_dmVec(1, 0);
    dc_dm(MY, MY) = dcVec_dmVec(1, 1);

    return dc_dm;
}


SimpleShapeResult::SimpleShapeResult() : ellipse(std::numeric_limits<lsst::meas::base::ErrElement>::quiet_NaN(),
                                     std::numeric_limits<lsst::meas::base::ErrElement>::quiet_NaN(),
                                     std::numeric_limits<lsst::meas::base::ErrElement>::quiet_NaN()),
                             center(std::numeric_limits<lsst::meas::base::ErrElement>::quiet_NaN(),
                                    std::numeric_limits<lsst::meas::base::ErrElement>::quiet_NaN()),
                             covariance(Eigen::Matrix<double,5,5>::Constant(std::numeric_limits<lsst::meas::base::ErrElement>::quiet_NaN()))
{}

static std::array<lsst::meas::base::FlagDefinition, SimpleShape::N_FLAGS> const flagDefs = {{
    {"flag", "general failure flag, set if anything went wrong"}
}};

SimpleShapeResultKey SimpleShapeResultKey::addFields(
        afw::table::Schema & schema,
        std::string const & name
) {
            SimpleShapeResultKey r;
            r._shapeResult = lsst::afw::table::QuadrupoleKey::addFields(schema,
                                                            name, "elliptical Gaussian moments");
            r._centroidResult = lsst::afw::table::Point2DKey::addFields(schema,
                                                            name, "elliptical Gaussian moments", "pixels");
            r._uncertantyResult = lsst::afw::table::CovarianceMatrixKey<double, 5>::addFields(schema,name, 
                                                            std::vector<std::string> ({"Ixx", "Iyy", "Ixy",
                                                                                      "Ix", "Iy"}), "pixels");
            r._flagHandler = lsst::meas::base::FlagHandler::addFields(schema,
                                                            name, flagDefs.begin(), flagDefs.end());
            return r;
}

SimpleShapeResultKey::SimpleShapeResultKey(lsst::afw::table::SubSchema const & s) :
    _shapeResult(s),
    _centroidResult(s),
    _uncertantyResult(s, std::vector<std::string> ({"Ixx", "Iyy", "Ixy", "Ix", "Iy"})),
    _flagHandler(s, flagDefs.begin(), flagDefs.end())
{}

SimpleShapeResult SimpleShapeResultKey::get(lsst::afw::table::BaseRecord const & record) const {
    SimpleShapeResult result;
    result.ellipse = record.get(_shapeResult);
    result.center = record.get(_centroidResult);
    result.covariance = record.get(_uncertantyResult);
    for (int n = 0; n < SimpleShape::N_FLAGS; ++n) {
        result.flags[n] = _flagHandler.getValue(record, n);
    }
    return result;
}

void SimpleShapeResultKey::set(afw::table::BaseRecord & record, SimpleShapeResult const & value) const {
    record.set(_shapeResult, value.ellipse);
    record.set(_centroidResult, value.center);
    record.set(_uncertantyResult, value.covariance);
    for (int n = 0; n < SimpleShape::N_FLAGS; ++n) {
        _flagHandler.setValue(record, n, value.flags[n]);
    }
}

bool SimpleShapeResultKey::operator==(SimpleShapeResultKey const & other) const {
    return _shapeResult == other._shapeResult &&
        _centroidResult == other._centroidResult &&
        _uncertantyResult == other._uncertantyResult;
    //don't bother with flags - if we've gotten this far, it's basically impossible the flags don't match
}

bool SimpleShapeResultKey::isValid() const {
    return _shapeResult.isValid() &&
        _centroidResult.isValid() &&
        _uncertantyResult.isValid();
    //don't bother with flags - if we've gotten this far, it's basically impossible the flags don't match
}


SimpleShape::SimpleShape(
        Control const & ctrl,
        std::string const & name,
        afw::table::Schema & schema
    )
    : _ctrl(ctrl),
      _resultKey(SimpleShapeResultKey::addFields(schema, name)),
      _centroidExtractor(schema, name)
    {}

void SimpleShape::measure(
    afw::table::SourceRecord & source,
    afw::image::Exposure<float> const & exposure
) const {
    afw::geom::Point2D center = _centroidExtractor(source, _resultKey.getFlagHandler());
    afw::geom::ellipses::Ellipse weight(afw::geom::ellipses::Axes(_ctrl.sigma), center);
    // set flags so an exception throw produces a flagged source
    SimpleShapeResult result = computeMoments(weight, exposure.getMaskedImage(), _ctrl.nSigmaRegion);
    source.set(_resultKey, result);
}

void SimpleShape::fail(
        afw::table::SourceRecord & measRecord,
        lsst::meas::base::MeasurementError * error
) const {
    _resultKey.getFlagHandler().handleFailure(measRecord, error);
}

#define INSTANTIATE_IMAGE(IMAGE) \
    template SimpleShapeResult SimpleShape::computeMoments(\
        afw::geom::ellipses::Ellipse const &,              \
        IMAGE const &,                                     \
        double                                             \
    )

INSTANTIATE_IMAGE(lsst::afw::image::MaskedImage<int>);
INSTANTIATE_IMAGE(lsst::afw::image::MaskedImage<float>);
INSTANTIATE_IMAGE(lsst::afw::image::MaskedImage<double>);

}}}} // namespace lsst::meas::extensions::simpleShape
