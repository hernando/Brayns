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

#pragma once

#include <brayns/common/types.h>
#include <brayns/parameters/AbstractParameters.h>

namespace brayns
{
class CircuitParameters;
}

SERIALIZATION_ACCESS(CircuitParameters)

namespace brayns
{
class CircuitParameters : public AbstractParameters
{
public:
    CircuitParameters();

    /** @copydoc AbstractParameters::parse */
    void parse(const po::variables_map& vm) final;

    /** @copydoc AbstractParameters::print */
    void print() final;

    const std::string& getTargets() const { return _targets; }
    strings getTargetsAsStrings() const;
    /** Density of cells in the circuit in percent (Mainly for testing
     * purposes) */
    float getDensity() const;
    /** Random seed of the circuit */
    size_t getRandomSeed() const { return _randomSeed; }
    /**
     * Defines a bounding box outside of which geometry of a circuit will not be
     * loaded
     */
    const Boxd& getBoundingBox() const { return _boundingBox; }
    void setBoundingBox(const Boxd& value)
    {
        _updateValue(_boundingBox, value);
    }

    const std::string& getReport() const { return _report; }
    /** Defines the range of frames to be loaded for the simulation */
    double getEndSimulationTime() const { return _endSimulationTime; }
    double getStartSimulationTime() const { return _startSimulationTime; }
    double getSimulationStep() const { return _simulationStep; }
    Vector2d getSimulationValuesRange() const { return _simulationValuesRange; }
    bool transformMeshes() const { return _transformMeshes; }
    /** Defines the folder where morphologies meshes are stored. Meshes must
     * have the same name as the h5/SWC morphology file, suffixed with an
     * extension supported by the assimp library.
     */
    const std::string& getMeshFolder() const { return _meshFolder; }
    /**
     * Return the filename pattern use to load meshes
     */
    const std::string& getMeshFilePattern() const { return _meshFilePattern; }
private:
    std::string _targets;
    float _density;
    size_t _randomSeed;

    std::string _report;
    double _startSimulationTime;
    double _endSimulationTime{std::numeric_limits<double>::max()};
    double _simulationStep;
    Vector2d _simulationValuesRange{-std::numeric_limits<double>::max(),
                                    std::numeric_limits<double>::max()};

    Boxd _boundingBox{{0, 0, 0}, {0, 0, 0}};

    std::string _meshFolder;
    std::string _meshFilePattern;
    bool _transformMeshes{true};

    SERIALIZATION_FRIEND(CircuitParameters)
};
}
