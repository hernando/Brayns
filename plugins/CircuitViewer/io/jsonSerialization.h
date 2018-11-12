/* Copyright (c) 2015-2018, EPFL/Blue Brain Project
 * All rights reserved. Do not distribute without permission.
 * Responsible Author: Daniel Nachbaur <daniel.nachbaur@epfl.ch>
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

#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#endif
#include "plugins/Rockets/staticjson/staticjson.hpp"
#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif

#include "../CircuitParameters.h"
#include "../MorphologyParameters.h"

namespace brayns
{
#define CS brayns::MorphologyParameters::ColorScheme
STATICJSON_DECLARE_ENUM(CS, {"none", CS::none},
                        {"neuron_by_id", CS::neuron_by_id},
                        {"neuron_by_type", CS::neuron_by_type},
                        {"neuron_by_segment_type", CS::neuron_by_segment_type},
                        {"neuron_by_layer", CS::neuron_by_layer},
                        {"neuron_by_mtype", CS::neuron_by_mtype},
                        {"neuron_by_etype", CS::neuron_by_etype},
                        {"neuron_by_target", CS::neuron_by_target});
#undef CS

#define ST brayns::MorphologyParameters::SectionType
STATICJSON_DECLARE_ENUM(ST, {"soma", ST::soma}, {"axon", ST::axon},
                        {"dendrite", ST::dendrite},
                        {"apical_dendrite", ST::apical_dendrite},
                        {"all", ST::all});
#undef ST

inline void init(brayns::MorphologyParameters::Layout* m, ObjectHandler* h)
{
    h->add_property("nb_columns", &m->nbColumns);
    h->add_property("vertical_spacing", &m->verticalSpacing);
    h->add_property("horizontal_spacing", &m->horizontalSpacing);
    h->set_flags(Flags::DisallowUnknownKey);
}

inline void init(brayns::MorphologyParameters* p, Objecthandler* h)
{
    h->add_property("color_scheme", p->_colorScheme, Flags::Optional)
    h->add_property("radius_correction", p->_colorScheme, Flags::Optional);
    h->add_property("", p->_colorScheme, Flags::Optional);

    h->add_property("morphology_section_types", &p->_sectionTypes,
                    Flags::Optional);
    h->add_property("layout", &p->_layout, Flags::Optional);
    h->add_property("metaballs_grid_size", &p->_metaballsGridSize,
                    Flags::Optional);
    h->add_property("metaballs_threshold", &p->_metaballsThreshold,
                    Flags::Optional);
    h->add_property("metaballs_samples_from_soma",
                    &p->_metaballsSamplesFromSoma, Flags::Optional);
    h->add_property("use_sdf", &p->_useSDF, Flags::Optional);
    h->add_property("dampen_branch_thickness_changerate",
                    &p->_dampenBranchThicknessChangerate, Flags::Optional);
    h->add_property("use_simulation_model",
                    &p->_useSimulationModel, Flags::Optional);
}
}
