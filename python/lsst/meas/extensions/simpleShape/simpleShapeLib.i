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

%define simpleShape_DOCSTRING
"
Access to the classes from the meas_extensions_simpleShape library
"
%enddef

%feature("autodoc", "1");
%module(package="lsst.meas.extensions.simpleShape", docstring=simpleShape_DOCSTRING) simpleShapeLib

%{
#include "lsst/pex/logging.h"
#include "lsst/meas/base.h"
#include "lsst/afw/math.h"
#include "lsst/afw/math/FunctionLibrary.h"
#include "lsst/afw/detection.h"
#include "lsst/meas/extensions/simpleShape.h"

#define PY_ARRAY_UNIQUE_SYMBOL LSST_MEAS_EXT_SIMPLESHAPE_NUMPY_ARRAY_API
#include "numpy/arrayobject.h"
#include "ndarray/swig.h"
#include "ndarray/swig/eigen.h"
%}

%init %{
    import_array();
%}

%include "lsst/p_lsstSwig.i"

%lsst_exceptions();

%include "ndarray.i"

%declareNumPyConverters(Eigen::Matrix<double,5,5>)

%import "lsst/meas/base/baseLib.i"

%shared_ptr(lsst::meas::extensions::simpleShape::SimpleShapeControl);
%shared_ptr(lsst::meas::extensions::simpleShape::SimpleShape);

%include "lsst/meas/extensions/simpleShape.h"

%extend lsst::meas::extensions::simpleShape::SimpleShape {
// add class-scope typedefs as typeobject class attributes, just
// to make C++ and Python interfaces more similar
%pythoncode %{
    Control = SimpleShapeControl
%}
}

%template(computeMoments) lsst::meas::extensions::simpleShape::SimpleShape::computeMoments<float>;
%template(computeMoments) lsst::meas::extensions::simpleShape::SimpleShape::computeMoments<double>;
