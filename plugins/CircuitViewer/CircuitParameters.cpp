/* Copyright (c) 2015-2018, EPFL/Blue Brain Project
 * All rights reserved. Do not distribute without permission.
 * Responsible Author: Juan Hernando <juan.hernando@epfl.ch>
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

#include "CircuitParameters.h"

#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>

namespace brayns
{
namespace
{
const std::string PARAM_DENSITY = "circuit-density";
const std::string PARAM_BOUNDING_BOX = "circuit-bounding-box";
const std::string PARAM_MESH_FOLDER = "circuit-mesh-folder";
const std::string PARAM_MESH_FILE_PATTERN = "circuit-mesh-filename-pattern";
const std::string PARAM_TRANSFORM_MESHES =
    "circuit-transforma-meshes";
const std::string PARAM_TARGETS = "circuit-targets";
const std::string PARAM_REPORT = "circuit-report";
const std::string PARAM_START_SIMULATION_TIME =
    "circuit-start-simulation-time";
const std::string PARAM_END_SIMULATION_TIME =
    "circuit-end-simulation-time";
const std::string PARAM_SIMULATION_STEP = "circuit-simulation-step";
const std::string PARAM_SIMULATION_RANGE =
    "circuit-simulation-values-range";
const std::string PARAM_RANDOM_SEED = "circuit-random-seed";
}

CircuitParameters::CircuitParameters()
    : AbstractParameters("Circuit loading")
{
    _parameters.add_options() //
        (PARAM_TARGETS.c_str(), po::value<std::string>(),
         "Circuit targets [comma separated strings]")
        //
        (PARAM_DENSITY.c_str(), po::value<float>(),
         "Density of cells in the circuit in percent [float]")
        //
        (PARAM_MESH_FOLDER.c_str(), po::value<std::string>(),
         "Folder containing meshed morphologies [string]")
        //
        (PARAM_REPORT.c_str(), po::value<std::string>(),
         "Circuit report [string]")
        //
        (PARAM_START_SIMULATION_TIME.c_str(), po::value<double>(),
         "Start simulation timestamp [double]")
        //
        (PARAM_END_SIMULATION_TIME.c_str(), po::value<double>(),
         "End simulation timestamp [double]")
        //
        (PARAM_SIMULATION_STEP.c_str(), po::value<double>(),
         "Step between simulation frames [double]")
        //
        (PARAM_SIMULATION_RANGE.c_str(),
         po::value<floats>()->multitoken(),
         "Minimum and maximum values for the simulation [float float]")
        //
        (PARAM_RANDOM_SEED.c_str(), po::value<size_t>(),
         "Random seed for circuit [int]")
        //
        (PARAM_BOUNDING_BOX.c_str(), po::value<floats>()->multitoken(),
         "Does not load circuit geometry outside of the specified "
         "bounding box [float float float float float float]")
        //
        (PARAM_MESH_FILE_PATTERN.c_str(), po::value<std::string>(),
         "Pattern used to determine the name of the file containing a "
         "meshed morphology [string]")
        //
        (PARAM_TRANSFORM_MESHES.c_str(),
         po::bool_switch()->default_value(false),
         "Enable mesh transformation according to circuit information.");
}

void CircuitParameters::parse(const po::variables_map& vm)
{

    if (vm.count(PARAM_TARGETS))
        _targets = vm[PARAM_TARGETS].as<std::string>();
    if (vm.count(PARAM_DENSITY))
        _density = vm[PARAM_DENSITY].as<float>();
    if (vm.count(PARAM_RANDOM_SEED))
        _randomSeed = vm[PARAM_RANDOM_SEED].as<size_t>();
    if (vm.count(PARAM_REPORT))
        _report = vm[PARAM_REPORT].as<std::string>();
    if (vm.count(PARAM_START_SIMULATION_TIME))
        _startSimulationTime = vm[PARAM_START_SIMULATION_TIME].as<double>();
    if (vm.count(PARAM_END_SIMULATION_TIME))
        _endSimulationTime = vm[PARAM_END_SIMULATION_TIME].as<double>();
    if (vm.count(PARAM_SIMULATION_STEP))
        _simulationStep = vm[PARAM_SIMULATION_STEP].as<double>();
    if (vm.count(PARAM_SIMULATION_RANGE))
    {
        floats values = vm[PARAM_SIMULATION_RANGE].as<floats>();
        if (values.size() == 2)
            _simulationValuesRange = Vector2f(values[0], values[1]);
    }
    if (vm.count(PARAM_BOUNDING_BOX))
    {
        const floats values = vm[PARAM_BOUNDING_BOX].as<floats>();
        if (values.size() == 6)
        {
            _boundingBox.reset();
            _boundingBox.merge(Vector3f(values[0], values[1], values[2]));
            _boundingBox.merge(Vector3f(values[3], values[4], values[5]));
        }
        else
            BRAYNS_ERROR << "Invalid number of values for "
                         << PARAM_BOUNDING_BOX << std::endl;
    }
    if (vm.count(PARAM_MESH_FOLDER))
        _meshFolder = vm[PARAM_MESH_FOLDER].as<std::string>();
    if (vm.count(PARAM_MESH_FILE_PATTERN))
        _meshFilePattern = vm[PARAM_MESH_FILE_PATTERN].as<std::string>();
    _transformMeshes = vm[PARAM_TRANSFORM_MESHES].as<bool>();
}

void CircuitParameters::print()
{
    AbstractParameters::print();
    BRAYNS_INFO << "Targets                 : "
                << _targets << std::endl;
    BRAYNS_INFO << "Density                 : "
                << _density << std::endl;
    BRAYNS_INFO << "Report                  : "
                << _report << std::endl;
    BRAYNS_INFO << "Start simulation time   : "
                << _startSimulationTime << std::endl;
    BRAYNS_INFO << "End simulation time     : "
                << _endSimulationTime << std::endl;
    BRAYNS_INFO << "Simulation step         : "
                << _simulationStep << std::endl;
    BRAYNS_INFO << "Simulation values range : "
                << _simulationValuesRange << std::endl;
    BRAYNS_INFO << "Bounding box            : "
                << _boundingBox << std::endl;
    BRAYNS_INFO << "Mesh folder             : "
                << _meshFolder << std::endl;
    BRAYNS_INFO << "Mesh filename pattern   : "
                << _meshFilePattern << std::endl;
    BRAYNS_INFO << "Mesh transformation     : "
                << (_transformMeshes ? "Yes" : "No")
                << std::endl;
}

strings CircuitParameters::getTargetsAsStrings() const
{
    strings targets;
    boost::char_separator<char> separator(",");
    boost::tokenizer<boost::char_separator<char>> tokens(_targets, separator);
    for_each(tokens.begin(), tokens.end(),
             [&targets](const std::string& s) { targets.push_back(s); });
    return targets;
}

float CircuitParameters::getDensity() const
{
    return std::max(0.f, std::min(100.f, _density));
}

}
