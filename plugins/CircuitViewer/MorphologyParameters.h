/* Copyright (c) 2015-2018, EPFL/Blue Brain Project
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

#pragma once

#include <brayns/parameters/AbstractParameters.h>

namespace brayns
{
class MorphologyParameters : public AbstractParameters
{
public:
    /** Morphology element types */
    enum class SectionType
    {
        soma = 0x01,
        axon = 0x02,
        dendrite = 0x04,
        apical_dendrite = 0x08,
        all = 0xff
    };
    using SectionTypes = std::vector<SectionType>;

    enum class ColorScheme
    {
        none = 1,
        by_id = 1,
        by_type = 2,
        by_segment_type = 3,
        by_layer = 4,
        by_mtype = 5,
        by_etype = 6,
        by_target = 7,
    };

    /**
     * Defines how morphologies should be organized in space when the layout
     * mode is selected. The idea is to present the morphology in a grid with a
     * given number of columns, and a spacing in between. The spacing scale is
     * the same as the one from the morphologies.
     */

    struct Layout
    {
        size_t nbColumns{0};
        size_t verticalSpacing{0};
        size_t horizontalSpacing{0};
    };

    MorphologyParameters();

    /** @copydoc AbstractParameters::parse */
    void parse(const po::variables_map& vm) final;

    /** @copydoc AbstractParameters::print */
    void print() final;

    const SectionTypes& getSectionTypes() const { return _sectionTypes; }
    ColorScheme getColorScheme() const { return _colorScheme; }
    const std::string& getColorSchemeAsString(const ColorScheme value) const;
    void setColorScheme(const ColorScheme value)
    {
        _updateValue(_colorScheme, value);
    }

    const Layout& getLayout() const { return _layout; }
    bool useSomaOnly() const
    {
        return (_sectionTypes.size() == 1 &&
                _sectionTypes[0] == SectionType::soma);
    }

    bool useRealisticSomas() const { return _metaballsGridSize != 0; }
    /** Radius correction applied to spheres and cylinders.
     * @param value Radius value. The radius contained in the data source is
     *        ignored and all geometries use the specified value.
     */
    void setRadiusCorrection(const float value)
    {
        _updateValue(_radiusCorrection, value);
    }
    float getRadiusCorrection() const { return _radiusCorrection; }
    bool useSDF() const { return _useSDF; }
    bool getDampenBranchThicknessChangerate() const
    {
        return _dampenBranchThicknessChangerate;
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
    /**
     * Defines if a different model is used to handle the simulation geometry.
     * If set to True, the shading of the main geometry model will be done
     * using information stored in a secondary model that contains the
     * simulation information. See OSPRay simulation renderer for more details.
     */
    bool useSimulationModel() const { return _useSimulationModel; }
    void setUseSimulationModel(const bool value)
    {
        _updateValue(_useSimulationModel, value);
    }

private:
    ColorScheme _colorScheme{ColorScheme::none};
    float _radiusCorrection{0.f};

    SectionTypes _sectionTypes{{SectionType::all}};
    Layout _layout;

    size_t _metaballsGridSize{0};
    float _metaballsThreshold{1.f};
    size_t _metaballsSamplesFromSoma{3};

    bool _useSDF{false};
    bool _dampenBranchThicknessChangerate{false};

    bool _useSimulationModel{true};
};
}
