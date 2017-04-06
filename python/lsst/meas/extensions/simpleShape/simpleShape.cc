/*
 * LSST Data Management System
 *
 * This product includes software developed by the
 * LSST Project (http://www.lsst.org/).
 * See the COPYRIGHT file
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
 * see <https://www.lsstcorp.org/LegalNotices/>.
 */
#include "pybind11/pybind11.h"

#include "numpy/arrayobject.h"
#include "ndarray/pybind11.h"

#include "lsst/meas/extensions/simpleShape.h"
#include "lsst/pex/config/python.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace meas {
namespace extensions {
namespace simpleShape {

namespace {

template <typename T, typename PyClass>
static void declareMoments(PyClass cls) {
    cls.def_static("computeMoments", (SimpleShapeResult(*)(afw::geom::ellipses::Ellipse const &,
                                                           afw::image::MaskedImage<T> const &, double)) &
                                             SimpleShape::computeMoments,
                   "weight"_a, "image"_a, "nSigmaRegion"_a = 3.0);
}

}  // <anonymous>

PYBIND11_PLUGIN(simpleShape) {
    py::module::import("lsst.afw.geom");
    py::module::import("lsst.afw.image");
    py::module::import("lsst.afw.table");
    py::module::import("lsst.meas.base");

    py::module mod("simpleShape");

    if (_import_array() < 0) {
        PyErr_SetString(PyExc_ImportError, "numpy.core.multiarray failed to import");
        return nullptr;
    }

    /* Module level */
    py::class_<SimpleShape, std::shared_ptr<SimpleShape>, base::SimpleAlgorithm> clsSimpleShape(
            mod, "SimpleShape");
    py::class_<SimpleShapeControl> clsSimpleShapeControl(mod, "SimpleShapeControl");
    py::class_<SimpleShapeResult> clsSimpleShapeResult(mod, "SimpleShapeResult");
    py::class_<SimpleShapeResultKey> clsSimpleShapeResultKey(mod, "SimpleShapeResultKey");

    /* Constructors */
    clsSimpleShape.def(py::init<SimpleShape::Control const &, std::string const &, afw::table::Schema &>(),
                       "ctrl"_a, "name"_a, "schema"_a);
    clsSimpleShapeControl.def(py::init<>());

    clsSimpleShapeResultKey.def(py::init<afw::table::SubSchema const &>(), "s"_a);

    /* Members */
    LSST_DECLARE_CONTROL_FIELD(clsSimpleShapeControl, SimpleShapeControl, sigma);
    LSST_DECLARE_CONTROL_FIELD(clsSimpleShapeControl, SimpleShapeControl, nSigmaRegion);

    declareMoments<float>(clsSimpleShape);
    declareMoments<double>(clsSimpleShape);
    clsSimpleShape.def_static("correctWeightedMoments", &SimpleShape::correctWeightedMoments, "weight"_a,
                              "ellipse"_a, "center"_a);

    clsSimpleShapeResult.def_readwrite("ellipse", &SimpleShapeResult::ellipse);
    clsSimpleShapeResult.def_readwrite("center", &SimpleShapeResult::center);
    clsSimpleShapeResult.def_readwrite("covariance", &SimpleShapeResult::covariance);
    clsSimpleShapeResult.def("getFlag", &SimpleShapeResult::getFlag);

    clsSimpleShapeResultKey.def("get", &SimpleShapeResultKey::get, "record"_a);

    return mod.ptr();
}

}  // shapeHSM
}  // extensions
}  // meas
}  // lsst
