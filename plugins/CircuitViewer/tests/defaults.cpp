/* Copyright (c) 2016-2018, EPFL/Blue Brain Project
 * All rights reserved. Do not distribute without permission.
 * Responsible Author: Juan Hernando juan.hernando@epfl.ch
 *
 * This file is part of Brayns <https://github.com/BlueBrain/Brayns>
 *
 * This library is free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License version 3.0 as published
 * by the Free Software Foundation.
 *
 * This library is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
 * details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this library; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */

#include "../CircuitParameters.h"
#include "../MorphologyParameters.h"

#define BOOST_TEST_MODULE brayns
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE(morphologyDefaults)
{
    brayns::MorphologyParameters parameters;
    BOOST_CHECK_EQUAL(parameters.getRadiusCorrection(), 0.f);
    BOOST_CHECK_EQUAL(
        brayns::enumsToBitmask(parameters.getSectionTypes()),
        size_t(brayns::MorphologyParameters::SectionType::all));
    BOOST_CHECK_EQUAL(parameters.getLayout().nbColumns, 0);
    BOOST_CHECK_EQUAL(parameters.useSDF(), false);
}


BOOST_AUTO_TEST_CASE(circuitDefaults)
{
    brayns::CircuitParameters parameters;
    BOOST_CHECK_EQUAL(parameters.getTargets(), "");
    BOOST_CHECK_EQUAL(parameters.getReport(), "");
    BOOST_CHECK_EQUAL(parameters.getStartSimulationTime(), 0.f);
    BOOST_CHECK_EQUAL(parameters.getEndSimulationTime(),
                      std::numeric_limits<float>::max());
    BOOST_CHECK_EQUAL(parameters.getSimulationValuesRange().x(),
                      std::numeric_limits<double>::max());
    BOOST_CHECK_EQUAL(parameters.getSimulationValuesRange().y(),
                      std::numeric_limits<double>::min());
}
