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

    /** Radius multiplier applied to spheres, cones and cylinders.
     * @param value Radius multiplier. Multiplies the radius contained in the
     *        data source by the specified value.
     */
    void setRadiusMultiplier(const float value)
    {
        _updateValue(_radiusMultiplier, value);
    }
    float getRadiusMultiplier() const { return _radiusMultiplier; }

    /** Enables a different color for every scene element whenloading them
     * from a given folder.
     */
    ColorScheme getColorScheme() const { return _colorScheme; }
    const std::string& getColorSchemeAsString(const ColorScheme value) const;
    void setColorScheme(const ColorScheme value)
    {
        _updateValue(_colorScheme, value);
    }

    /** Mesh quality */
    GeometryQuality getGeometryQuality() const { return _geometryQuality; }
    const std::string& getGeometryQualityAsString(
        const GeometryQuality value) const;

    /**
     * Defines what memory mode should be used between Brayns and the
     * underlying renderer
     */
    MemoryMode getMemoryMode() const { return _memoryMode; };

    const std::set<BVHFlag>& getDefaultBVHFlags() const
    {
        return _defaultBVHFlags;
    }

protected:
    void parse(const po::variables_map& vm) final;

    // Scene
    std::string _loadCacheFile;
    std::string _saveCacheFile;
    std::set<BVHFlag> _defaultBVHFlags;

    // Geometry
    float _radiusMultiplier;
    ColorScheme _colorScheme;
    GeometryQuality _geometryQuality;

    // System parameters
    MemoryMode _memoryMode;

    SERIALIZATION_FRIEND(GeometryParameters)
};
}
#endif // GEOMETRYPARAMETERS_H
