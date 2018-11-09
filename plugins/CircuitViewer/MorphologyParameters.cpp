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

#include "MorphologyParameters.h"

namespace brayns
{
namespace
{
const std::string PARAM_COLOR_SCHEME = "morphology-color-scheme";
const std::string PARAM_RADIUS_CORRECTION = "radius-correction";
const std::string PARAM_SECTION_TYPES = "morphology-section-types";
const std::string PARAM_LAYOUT = "morphology-layout";
const std::string PARAM_METABALLS_GRIDSIZE = "metaballs-grid-size";
const std::string PARAM_METABALLS_THRESHOLD = "metaballs-threshold";
const std::string PARAM_METABALLS_SAMPLES_FROM_SOMA =
    "metaballs-samples-from-soma";
const std::string PARAM_DAMPEN_BRANCH_THICKNESS_CHANGERATE =
    "morphology-dampen-branch-thickness-changerate";
const std::string PARAM_USE_SDF_GEOMETRIES = "morphology-use-sdf-geometries";
const std::string PARAM_USE_SIMULATION_MODEL =
    "morphology-use-simulation-model";

const std::array<std::string, 12> COLOR_SCHEMES = {
    {"none", "neuron-by-id", "neuron-by-type", "neuron-by-segment-type",
     "neuron-by-layer", "neuron-by-mtype", "neuron-by-etype",
     "neuron-by-target"}};
}

void MorphologyParameters::print()
{
    AbstractParameters::print();
    BRAYNS_INFO << "Color scheme         : "
                << getColorSchemeAsString(_colorScheme) << std::endl;
    BRAYNS_INFO << "Radius correction    : " << _radiusCorrection << std::endl;
    BRAYNS_INFO << "Section types        : " << enumsToBitmask(_sectionTypes)
                << std::endl;
    BRAYNS_INFO << "Use simulation model : "
                << (_useSimulationModel ? "Yes" : "No")
                << std::endl;
    BRAYNS_INFO << "Layout" << std::endl;
    BRAYNS_INFO << " - Columns            : " << _layout.nbColumns << std::endl;
    BRAYNS_INFO << " - Vertical spacing   : " << _layout.verticalSpacing
                << std::endl;
    BRAYNS_INFO << " - Horizontal spacing : " << _layout.horizontalSpacing
                << std::endl;
    BRAYNS_INFO << "Metaballs" << std::endl;
    BRAYNS_INFO << " - Grid size         : " << _metaballsGridSize << std::endl;
    BRAYNS_INFO << " - Threshold         : " << _metaballsThreshold
                << std::endl;
    BRAYNS_INFO << " - Samples from soma : " << _metaballsSamplesFromSoma
                << std::endl;
    BRAYNS_INFO << "SDF geometry" << std::endl;
    BRAYNS_INFO << " - Use SDF                            : " << _useSDF
                << std::endl;
    BRAYNS_INFO << " - Dampen branch thickess change rate : " << _useSDF
                << std::endl;
}

MorphologyParameters::MorphologyParameters()
    : AbstractParameters("Morphology loading")
{
    _parameters.add_options() //
        (PARAM_COLOR_SCHEME.c_str(), po::value<std::string>(),
         "Color scheme to be applied to the geometry "
         "[none|neuron-by-id|neuron-by-type|neuron-by-segment-type|"
         "neuron-by-layer|neuron-by-mtype|neuron-by-etype|neuron-by-"
         "target]")
        //
        (PARAM_RADIUS_CORRECTION.c_str(), po::value<float>(),
         "Forces radius of spheres and cylinders to the specified value "
         "[float]")
        //
        (PARAM_SECTION_TYPES.c_str(), po::value<size_t>(),
         "Morphology section types (1: soma, 2: axon, 4: dendrite, "
         "8: apical dendrite). Values can be added to select more than "
         "one type of section")
        //
        (PARAM_LAYOUT.c_str(), po::value<size_ts>()->multitoken(),
         "Morphology layout defined by number of "
         "columns, vertical spacing, horizontal spacing "
         "[int int int]")
        //
        (PARAM_METABALLS_GRIDSIZE.c_str(), po::value<size_t>(),
         "Metaballs grid size [int]. Activates automated meshing of somas "
         "if different from 0")
        //
        (PARAM_METABALLS_THRESHOLD.c_str(), po::value<float>(),
         "Metaballs threshold [float]")
        //
        (PARAM_METABALLS_SAMPLES_FROM_SOMA.c_str(), po::value<size_t>(),
         "Number of morphology samples (or segments) from soma used by "
         "automated meshing [int]")
        //
        (PARAM_DAMPEN_BRANCH_THICKNESS_CHANGERATE.c_str(),
         po::bool_switch()->default_value(false),
         "Dampen the thickness rate of change for branches in the morphology.")
        //
        (PARAM_USE_SDF_GEOMETRIES.c_str(),
         po::bool_switch()->default_value(false),
         "Use SDF geometries for drawing the morphologies.")
        //
        (PARAM_USE_SIMULATION_MODEL.c_str(),
         po::bool_switch()->default_value(false),
         "Defines if a different model is used to handle the simulation "
         "geometry.");
}

void MorphologyParameters::parse(const po::variables_map& vm)
{
    if (vm.count(PARAM_COLOR_SCHEME))
    {
        _colorScheme = ColorScheme::none;
        const auto& colorScheme = vm[PARAM_COLOR_SCHEME].as<std::string>();
        if (!colorScheme.empty())
        {
            auto it = std::find(COLOR_SCHEMES.begin(), COLOR_SCHEMES.end(),
                                colorScheme);
            if (it == COLOR_SCHEMES.end())
                throw po::error("No match for color scheme '" + colorScheme);

            const auto index = std::distance(COLOR_SCHEMES.begin(), it);
            _colorScheme = static_cast<ColorScheme>(index);
        }
    }
    if (vm.count(PARAM_RADIUS_CORRECTION))
        _radiusCorrection = vm[PARAM_RADIUS_CORRECTION].as<float>();
    if (vm.count(PARAM_SECTION_TYPES))
    {
        _sectionTypes.clear();
        const auto bits = vm[PARAM_SECTION_TYPES].as<size_t>();
        if (bits & size_t(SectionType::soma))
            _sectionTypes.push_back(SectionType::soma);
        if (bits & size_t(SectionType::axon))
            _sectionTypes.push_back(SectionType::axon);
        if (bits & size_t(SectionType::dendrite))
            _sectionTypes.push_back(SectionType::dendrite);
        if (bits & size_t(SectionType::apical_dendrite))
            _sectionTypes.push_back(SectionType::apical_dendrite);
    }

    if (vm.count(PARAM_LAYOUT))
    {
        size_ts values = vm[PARAM_LAYOUT].as<size_ts>();
        if (values.size() == 3)
        {
            _layout.nbColumns = values[0];
            _layout.verticalSpacing = values[1];
            _layout.horizontalSpacing = values[2];
        }
    }
    if (vm.count(PARAM_METABALLS_GRIDSIZE))
        _metaballsGridSize = vm[PARAM_METABALLS_GRIDSIZE].as<size_t>();
    if (vm.count(PARAM_METABALLS_THRESHOLD))
        _metaballsThreshold = vm[PARAM_METABALLS_THRESHOLD].as<float>();
    if (vm.count(PARAM_METABALLS_SAMPLES_FROM_SOMA))
        _metaballsSamplesFromSoma =
            vm[PARAM_METABALLS_SAMPLES_FROM_SOMA].as<size_t>();
    _dampenBranchThicknessChangerate =
        vm[PARAM_DAMPEN_BRANCH_THICKNESS_CHANGERATE].as<bool>();
    _useSDF = vm[PARAM_USE_SDF_GEOMETRIES].as<bool>();
    _useSimulationModel = vm[PARAM_USE_SIMULATION_MODEL].as<bool>();

    markModified();
}

const std::string& MorphologyParameters::getColorSchemeAsString(
    const ColorScheme value) const
{
    return COLOR_SCHEMES[static_cast<size_t>(value)];
}
}
