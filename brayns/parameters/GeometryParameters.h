/* Copyright (c) 2015-2016, EPFL/Blue Brain Project
 * All rights reserved. Do not distribute without permission.
 * Responsible Author: Cyrille Favreau <cyrille.favreau@epfl.ch>
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

#ifndef GEOMETRYPARAMETERS_H
#define GEOMETRYPARAMETERS_H

#include "AbstractParameters.h"

#include <brayns/common/types.h>

SERIALIZATION_ACCESS(GeometryParameters)

namespace brayns
{
/**
 * Defines how morphologies should be organized in space when the layout mode is
 * selected. The idea is to present the morphology in a grid with a given number
 * of columns, and a spacing in between. The spacing scale is the same as the
 * one from the morphologies.
 */
struct MorphologyLayout
{
    int32_t nbColumns{0};
    int32_t verticalSpacing{0};
    int32_t horizontalSpacing{0};
};

struct CircuitConfiguration
{
    bool useSimulationModel{false};
    Boxd boundingBox{{0, 0, 0}, {0, 0, 0}};
    double density{100};
    std::string meshFilenamePattern;
    std::string meshFolder;
    std::string targets;
    std::string report;
    double startSimulationTime{0};
    double endSimulationTime{std::numeric_limits<float>::max()};
    double simulationStep{0};
    Vector2d simulationValuesRange{std::numeric_limits<double>::max(),
                                   std::numeric_limits<double>::min()};
    size_t randomSeed = 0;
    bool meshTransformation{false};
};

/** Manages geometry parameters
 */
class GeometryParameters : public AbstractParameters
{
public:
    /**
       Parse the command line parameters and populates according class members
     */
    GeometryParameters();

    /** @copydoc AbstractParameters::print */
    void print() final;

    /** Binary representation of a scene to load */
    const std::string& getLoadCacheFile() const { return _loadCacheFile; }
    /** Binary representation of a scene to save */
    const std::string& getSaveCacheFile() const { return _saveCacheFile; }
    /** Circuit targets */
    const std::string& getCircuitTargets() const
    {
        return _circuitConfiguration.targets;
    }
    strings getCircuitTargetsAsStrings() const;
    /** Circuit compartment report */
    const std::string& getCircuitReport() const
    {
        return _circuitConfiguration.report;
    }
    /** Defines the folder where morphologies meshes are stored. Meshes must
     * have the same name as the h5/SWC morphology file, suffixed with an
     * extension supported by the assimp library
     */
    const std::string& getCircuitMeshFolder() const
    {
        return _circuitConfiguration.meshFolder;
    }
    /** ensity of cells in the circuit in percent (Mainly for testing
     * purposes) */
    double getCircuitDensity() const;

    /**
     * Defines a bounding box outside of which geometry of a circuit will not be
     * loaded
     */
    const Boxd& getCircuitBoundingBox() const
    {
        return _circuitConfiguration.boundingBox;
    }
    void setCircuitBoundingBox(const Boxd& value)
    {
        _updateValue(_circuitConfiguration.boundingBox, value);
    }

    /**
     * Defines if a different model is used to handle the simulation geometry.
     * If set to True, the shading of the main geometry model will be done
     * using information stored in a secondary model that contains the
     * simulation information. See OSPRay simulation renderer for more details.
     */
    bool getCircuitUseSimulationModel() const
    {
        return _circuitConfiguration.useSimulationModel;
    }
    void setCircuitUseSimulationModel(const bool value)
    {
        _updateValue(_circuitConfiguration.useSimulationModel, value);
    }
    /**
     * Return the filename pattern use to load meshes
     */
    const std::string& getCircuitMeshFilenamePattern() const
    {
        return _circuitConfiguration.meshFilenamePattern;
    }
    /** Radius multiplier applied to spheres, cones and cylinders.
     * @param value Radius multiplier. Multiplies the radius contained in the
     *        data source by the specified value.
     */
    void setRadiusMultiplier(const float value)
    {
        _updateValue(_radiusMultiplier, value);
    }
    float getRadiusMultiplier() const { return _radiusMultiplier; }
    /** Radius correction applied to spheres and cylinders.
     * @param value Radius value. The radius contained in the data source is
     *        ignored and all geometries use the specified value.
     */
    void setRadiusCorrection(const float value)
    {
        _updateValue(_radiusCorrection, value);
    }
    float getRadiusCorrection() const { return _radiusCorrection; }
    /** Enables a different color for every molecule/morphology/mesh when
     * loading them from a given folder
     */
    ColorScheme getColorScheme() const { return _colorScheme; }
    const std::string& getColorSchemeAsString(const ColorScheme value) const;
    void setColorScheme(const ColorScheme value)
    {
        _updateValue(_colorScheme, value);
    }

    /** Morphology quality */
    GeometryQuality getGeometryQuality() const { return _geometryQuality; }
    const std::string& getGeometryQualityAsString(
        const GeometryQuality value) const;

    /** Morphology section types*/
    const MorphologySectionTypes& getMorphologySectionTypes() const
    {
        return _morphologySectionTypes;
    }
    /** Morphology layout */
    const MorphologyLayout& getMorphologyLayout() const
    {
        return _morphologyLayout;
    }

    /** Defines the range of frames to be loaded for the simulation */
    double getCircuitEndSimulationTime() const
    {
        return _circuitConfiguration.endSimulationTime;
    }
    double getCircuitStartSimulationTime() const
    {
        return _circuitConfiguration.startSimulationTime;
    }
    double getCircuitSimulationStep() const
    {
        return _circuitConfiguration.simulationStep;
    }
    Vector2d getCircuitSimulationValuesRange() const
    {
        return _circuitConfiguration.simulationValuesRange;
    }

    size_t getCircuitMeshTransformation() const
    {
        return _circuitConfiguration.meshTransformation;
    }

    /** Random seed of the circuit */
    size_t getCircuitRandomSeed() const
    {
        return _circuitConfiguration.randomSeed;
    }

    /** Metaballs grid size */
    size_t getMetaballsGridSize() const { return _metaballsGridSize; }
    /** Metaballs threshold */
    float getMetaballsThreshold() const { return _metaballsThreshold; }
    /** Metaballs samples from soma */
    size_t getMetaballsSamplesFromSoma() const
    {
        return _metaballsSamplesFromSoma;
    }

    /** Realistic somas enabled? */
    bool useRealisticSomas() const { return _metaballsGridSize != 0; }
    /**
     * Defines what memory mode should be used between Brayns and the
     * underlying renderer
     */
    MemoryMode getMemoryMode() const { return _memoryMode; };
    bool getMorphologyDampenBranchThicknessChangerate() const
    {
        return _morphologyDampenBranchThicknessChangerate;
    }

    bool getMorphologyUseSDFGeometries() const
    {
        return _morphologyUseSDFGeometries;
    }

    const std::set<BVHFlag>& getDefaultBVHFlags() const
    {
        return _defaultBVHFlags;
    }

protected:
    void parse(const po::variables_map& vm) final;

    // Circuit
    CircuitConfiguration _circuitConfiguration;

    // Scene
    std::string _loadCacheFile;
    std::string _saveCacheFile;

    // Morphology
    float _radiusMultiplier;
    float _radiusCorrection;
    ColorScheme _colorScheme;
    GeometryQuality _geometryQuality;
    std::set<BVHFlag> _defaultBVHFlags;
    MorphologySectionTypes _morphologySectionTypes;
    MorphologyLayout _morphologyLayout;
    size_t _metaballsGridSize;
    float _metaballsThreshold;
    size_t _metaballsSamplesFromSoma;
    bool _morphologyDampenBranchThicknessChangerate;
    bool _morphologyUseSDFGeometries;

    // System parameters
    MemoryMode _memoryMode;

    SERIALIZATION_FRIEND(GeometryParameters)
};
}
#endif // GEOMETRYPARAMETERS_H
