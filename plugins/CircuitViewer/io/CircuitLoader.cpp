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

#include "CircuitLoader.h"
#include "ModelData.h"
#include "MorphologyLoader.h"
#include "SimulationHandler.h"

#include <brayns/common/scene/Model.h>
#include <brayns/common/scene/Scene.h>
#include <brayns/common/utils/Utils.h>

#include <brayns/parameters/ApplicationParameters.h>

#include <brain/brain.h>
#include <brion/brion.h>

#if BRAYNS_USE_ASSIMP
#include <brayns/io/MeshLoader.h>
#endif

#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>

using namespace brayns;

namespace
{
constexpr size_t LOAD_BATCH_SIZE = 100;
constexpr size_t DEFAULT_SOMA_RADIUS = 5;

using Property = brayns::PropertyMap::Property;
const Property PROP_DENSITY = {"density", "Density", 100.0};
const Property PROP_RANDOM_SEED = {"randomSeed", "Random seed", 0};
const Property PROP_REPORT = {"report", "Report", std::string("")};
const Property PROP_TARGETS = {"targets", "Targets", std::string("")};
const Property PROP_MESH_FILENAME_PATTERN = {"meshFilenamePattern",
                                             "Mesh filename pattern",
                                             std::string("")};
const Property PROP_MESH_FOLDER = {"meshFolder", "Mesh folder",
                                   std::string("")};
const Property PROP_USE_SIMULATION_MODEL = {"useSimulationModel",
                                            "Use simulation model", false};
const Property PROP_TRANSFORM_MESHES = {"transformMeshes", "Transform meshes",
                                        false};
const Property PROP_COLOR_SCHEME = {"colorScheme", "Color scheme",
                                    brayns::enumToString(ColorScheme::none),
                                    brayns::enumNames<ColorScheme>()};
const Property PROP_START_SIMULATION_TIME = {"startSimulationTime",
                                             "Start simulation time", 0.0};
const Property PROP_END_SIMULATION_TIME = {"endSimulationTime",
                                           "End simulation time", 0.0};
const Property PROP_SIMULATION_STEP = {"simulationStep", "Simulation step",
                                       0.0};
const Property PROP_SYNCHRONOUS_MODE = {"synchronousMode", "Synchronous mode",
                                        false};
const Property PROP_GEOMETRY_QUALITY = {
    "geometryQuality", "Geometry quality", int32_t(GeometryQuality::high),
    brayns::enumerateGeometryQualityEnums()};
const auto LOADER_NAME = "circuit";

using SimulationHandlerPtr = std::shared_ptr<SimulationHandler>;

/**
 * @brief getMeshFilenameFromGID Returns the name of the mesh file according
 * to the --circuit-mesh-folder, --circuit-mesh-filename-pattern command
 * line arguments and a GID
 * @param gid GID of the cell
 * @return A string with the full path of the mesh file
 */
std::string _getMeshFilenameFromGID(const uint64_t gid,
                                    const std::string& pattern,
                                    const std::string& folder)
{
    const std::string gidAsString = std::to_string(gid);
    const std::string GID = "{gid}";

    auto filenamePattern = pattern;
    if (!filenamePattern.empty())
        filenamePattern.replace(filenamePattern.find(GID), GID.length(),
                                gidAsString);
    else
        filenamePattern = gidAsString;
    return folder + "/" + filenamePattern;
}

std::vector<float> _getSomaRadii(const brain::Circuit& circuit,
                                 const float defaultRadius)
{
    std::map<std::string, float> radii;

    // Old nomemclature
    radii["NGC"] = 7.46;
    radii["CRC"] = 7.28;
    radii["BP"] = 6.74;
    radii["L6CTPC"] = 7.00;
    radii["ChC"] = 8.01;
    radii["L4PC"] = 7.84;
    radii["L2PC"] = 8.36;
    radii["LBC"] = 9.28;
    radii["NBC"] = 8.61;
    radii["L5STPC"] = 8.61;
    radii["L6CCPC"] = 7.71;
    radii["L3PC"] = 7.70;
    radii["ADC"] = 7.72;
    radii["L5TTPC"] = 11.17;
    radii["AHC"] = 7.54;
    radii["L5UTPC"] = 9.75;
    radii["MC"] = 9.37;
    radii["DBC"] = 6.58;
    radii["SBC"] = 7.75;
    radii["L4SS"] = 8.87;
    radii["BTC"] = 9.15;
    radii["L4SP"] = 7.76;
    radii["L6CLPC"] = 7.32;

    // New nomemclature
    radii["L1_DAC"] = 7.77242212296;
    radii["L1_NGC-DA"] = 7.71745554606;
    radii["L1_NGC-SA"] = 7.53378756841;
    radii["L1_HAC"] = 7.48774570227;
    radii["L1_DLAC"] = 7.63075244427;
    radii["L1_SLAC"] = 6.54029955183;
    radii["L23_PC"] = 7.7618568695;
    radii["L23_MC"] = 8.25832914484;
    radii["L23_BTC"] = 9.30396906535;
    radii["L23_DBC"] = 6.67491123753;
    radii["L23_BP"] = 6.47212839127;
    radii["L23_NGC"] = 8.40507149696;
    radii["L23_LBC"] = 9.0722948432;
    radii["L23_NBC"] = 8.64244471205;
    radii["L23_SBC"] = 7.77334165573;
    radii["L23_ChC"] = 8.38753128052;
    radii["L4_PC"] = 7.5360350558;
    radii["L4_SP"] = 7.97097187241;
    radii["L4_SS"] = 8.95526583989;
    radii["L4_MC"] = 9.64632482529;
    radii["L4_BTC"] = 8.13213920593;
    radii["L4_DBC"] = 9.06448200771;
    radii["L4_BP"] = 5.80337953568;
    radii["L4_NGC"] = 7.61350679398;
    radii["L4_LBC"] = 8.92902739843;
    radii["L4_NBC"] = 9.02906776877;
    radii["L4_SBC"] = 9.16592723673;
    radii["L4_ChC"] = 8.54107379913;
    radii["L5_TTPC1"] = 12.4378146535;
    radii["L5_TTPC2"] = 12.9319117383;
    radii["L5_UTPC"] = 7.55000627041;
    radii["L5_STPC"] = 8.84197418645;
    radii["L5_MC"] = 9.02451186861;
    radii["L5_BTC"] = 9.24175042372;
    radii["L5_DBC"] = 8.92543775895;
    radii["L5_BP"] = 5.98329114914;
    radii["L5_NGC"] = 6.62666320801;
    radii["L5_LBC"] = 10.0915957915;
    radii["L5_NBC"] = 9.00336164747;
    radii["L5_SBC"] = 10.0064823627;
    radii["L5_ChC"] = 7.85296694438;
    radii["L6_TPC_L1"] = 7.11837192363;
    radii["L6_TPC_L4"] = 7.08293668837;
    radii["L6_UTPC"] = 7.74583534818;
    radii["L6_IPC"] = 7.71743569988;
    radii["L6_BPC"] = 7.14034402714;
    radii["L6_MC"] = 9.55710192858;
    radii["L6_BTC"] = 9.97289562225;
    radii["L6_DBC"] = 9.88600463867;
    radii["L6_BP"] = 5.98329114914;
    radii["L6_NGC"] = 6.62666320801;
    radii["L6_LBC"] = 9.49501685154;
    radii["L6_NBC"] = 9.33885664259;
    radii["L6_SBC"] = 8.06322860718;
    radii["L6_ChC"] = 8.15173133214;

    const auto& mtypes = circuit.getMorphologyTypeNames();
    std::vector<float> output;
    for (auto mtype : mtypes)
    {
        auto iter = radii.find(mtype);
        if (iter == radii.end())
            output.push_back(defaultRadius);
        else
            output.push_back(iter->second);
    }

    return output;
}

/**
 * Return a vector of gids.size() elements with the layer index of each gid
 */
size_ts _getLayerIds(const brain::Circuit& circuit, const brain::GIDSet& gids)
{
    std::vector<brion::GIDSet> layers;
    try
    {
        for (auto layer : {"1", "2", "3", "4", "5", "6"})
            layers.push_back(circuit.getGIDs(std::string("Layer") + layer));
    }
    catch (...)
    {
        BRAYNS_ERROR << "Not all layer targets were found for neuron_by_layer"
                        " color scheme"
                     << std::endl;
    }

    std::vector<brion::GIDSet::const_iterator> iters;
    for (const auto& layer : layers)
        iters.push_back(layer.begin());

    size_ts materialIds;
    materialIds.reserve(gids.size());
    // For each GID find in which layer GID list it appears, assign
    // this index to its material id list.
    for (auto gid : gids)
    {
        size_t i = 0;
        for (i = 0; i != layers.size(); ++i)
        {
            while (iters[i] != layers[i].end() && *iters[i] < gid)
                ++iters[i];
            if (iters[i] != layers[i].end() && *iters[i] == gid)
            {
                materialIds.push_back(i);
                break;
            }
        }
        if (i == layers.size())
            throw std::runtime_error("Layer Id not found for GID " +
                                     std::to_string(gid));
    }
    return materialIds;
}

size_t _getMaterialId(const ColorScheme colorScheme, const uint64_t index,
                      const brain::neuron::SectionType sectionType,
                      const size_ts& perCellMaterialIds)
{
    size_t materialId = 0;
    switch (colorScheme)
    {
    case ColorScheme::neuron_by_id:
        materialId = index;
        break;
    case ColorScheme::neuron_by_segment_type:
        switch (sectionType)
        {
        case brain::neuron::SectionType::soma:
            materialId = 1;
            break;
        case brain::neuron::SectionType::axon:
            materialId = 2;
            break;
        case brain::neuron::SectionType::dendrite:
            materialId = 3;
            break;
        case brain::neuron::SectionType::apicalDendrite:
            materialId = 4;
            break;
        default:
            materialId = 0;
            break;
        }
        break;
    case ColorScheme::neuron_by_target: // no break
    case ColorScheme::neuron_by_etype:  // no break
    case ColorScheme::neuron_by_mtype:  // no break
    case ColorScheme::neuron_by_layer:  // no break
        if (index < perCellMaterialIds.size())
            materialId = perCellMaterialIds[index];
        else
        {
            materialId = NO_MATERIAL;
            BRAYNS_DEBUG << "Failed to get per cell material index"
                         << std::endl;
        }
        break;
    default:
        materialId = NO_MATERIAL;
    }
    return materialId;
}

struct CircuitProperties
{
    CircuitProperties() = default;
    CircuitProperties(const PropertyMap& properties)
    {
        const auto setVariable = [&](auto& variable, const std::string& name,
                                     auto defaultVal) {
            using T = typename std::remove_reference<decltype(variable)>::type;
            variable = properties.getProperty<T>(name, defaultVal);
        };

        const auto setEnumVariable =
            [&](auto& variable, const std::string& name, auto defaultVal) {
                using T = decltype(defaultVal);
                const auto enumStr =
                    properties.getProperty<std::string>(name, enumToString<T>(
                                                                  defaultVal));
                variable = stringToEnum<T>(enumStr);
            };

        setVariable(density, PROP_DENSITY.name, 100.0);
        setVariable(randomSeed, PROP_RANDOM_SEED.name, 0);
        setVariable(report, PROP_REPORT.name, "");
        setVariable(targets, PROP_TARGETS.name, "");
        setVariable(meshFilenamePattern, PROP_MESH_FILENAME_PATTERN.name, "");
        setVariable(meshFolder, PROP_MESH_FOLDER.name, "");
        setVariable(useSimulationModel, PROP_USE_SIMULATION_MODEL.name, false);
        setVariable(transformMeshes, PROP_TRANSFORM_MESHES.name, 0);
        setEnumVariable(colorScheme, PROP_COLOR_SCHEME.name, ColorScheme::none);
        setVariable(startSimulationTime, PROP_START_SIMULATION_TIME.name, 0.0);
        setVariable(endSimulationTime, PROP_END_SIMULATION_TIME.name, 0.0);
        setVariable(simulationStep, PROP_SIMULATION_STEP.name, 0.0);
        setVariable(synchronousMode, PROP_SYNCHRONOUS_MODE.name, false);
        setEnumVariable(geometryQuality, PROP_GEOMETRY_QUALITY.name,
                        GeometryQuality::high);

        boost::char_separator<char> separator(",");
        boost::tokenizer<boost::char_separator<char>> tokens(targets,
                                                             separator);
        for_each(tokens.begin(), tokens.end(),
                 [& targetList = targetList](const std::string& s) {
                     targetList.push_back(s);
                 });
    }

    double density = 0.0;
    int32_t randomSeed = 0;
    std::string report;
    std::vector<std::string> targetList;
    std::string targets;
    std::string meshFilenamePattern;
    std::string meshFolder;
    bool useSimulationModel = false;
    bool transformMeshes = 0;
    ColorScheme colorScheme = ColorScheme::none;

    double startSimulationTime = 0;
    double endSimulationTime = 0;
    double simulationStep = 0;

    bool synchronousMode = false;

    GeometryQuality geometryQuality = GeometryQuality::high;
};

CompartmentReportPtr _openCompartmentReport(const brain::Simulation* simulation,
                                            const std::string& name,
                                            const brion::GIDSet& gids)
{
    try
    {
        auto source = blueConfig.getReportSource(name);
        if (source.getPath().empty())
        {
            BRAYNS_ERROR << "Compartment report not found: " << name
                         << std::endl;
            return {};
        }
        return std::make_shared<brion::CompartmentReport>(source,
                                                          brion::MODE_READ,
                                                          gids);
    }
    catch (const std::exception& e)
    {
        BRAYNS_ERROR << e.what() << std::endl;
        return {};
    }
}

SimulationHandlerPtr _createSimulationHandler(
    const CompartmentReportPtr& report, const CircuitProperties& properties)
{
    return std::make_shared<SimulationHandler>(report,
                                               properties.synchronousMode,
                                               properties.startSimulationTime,
                                               properties.endSimulationTime,
                                               properties.simulationStep);
}

class Impl
{
public:
    Impl(Scene& scene, const PropertyMap& properties)
        : _scene(scene)
        , _properties(properties)
        , _morphologyParams(properties)
    {
        _morphologyParams.useSimulationModel = _properties.useSimulationModel;
    }

    ModelDescriptorPtr importCircuit(const std::string& source,
                                     const LoaderProgress& callback) const
    {
        ModelDescriptorPtr modelDesc;
        // Model (one for the whole circuit)
        ModelMetadata metadata = {{"density",
                                   std::to_string(_properties.density)},
                                  {"report", _properties.report},
                                  {"targets", _properties.targets},
                                  {"mesh-filename-pattern",
                                   _properties.meshFilenamePattern},
                                  {"mesh-folder", _properties.meshFolder}};
        auto model = _scene.createModel();

        // Open Circuit and select GIDs according to specified target
        const brion::BlueConfig bc(source);
        const brain::Circuit circuit(bc);

        brain::GIDSet allGids;
        size_ts targetSizes;

        const strings localTargets = _properties.targetList.empty()
                                         ? strings{{""}}
                                         : _properties.targetList;

        auto resolveTarget = [&](const std::string& name) {
            const float fraction = _properties.density / 100.0f;
            if (simulation)
                return simulation->getGIDs(name, fraction,
                                           _properties.randomSeed);
            return circuit.getRandomGIDs(fraction, name,
                                         _properties.randomSeed);
        };

        for (const auto& target : localTargets)
        {
            const auto gids = resolveTarget(target);
            if (gids.empty())
            {
                BRAYNS_ERROR << "Target " << target
                             << " does not contain any cells" << std::endl;
                continue;
            }

            BRAYNS_INFO << "Target " << target << ": " << gids.size()
                        << " cells" << std::endl;
            allGids.insert(gids.begin(), gids.end());
            targetSizes.push_back(gids.size());
        }

        if (allGids.empty())
            return {};

        auto compartmentReport =
            _openCompartmentReport(bc, _properties.report, allGids);
        if (compartmentReport)
        {
            model->setSimulationHandler(
                _createSimulationHandler(compartmentReport, _properties));
            // Only keep GIDs from the report
            allGids = compartmentReport->getGIDs();
        }

        const auto perCellMaterialIds =
            _createPerCellMaterialIds(circuit, allGids, targetSizes);

        if (_morphologyParams.sectionTypes ==
            std::vector<MorphologySectionType>{MorphologySectionType::soma})

        {
            // Import soma only morphologies as spheres
            _addSomaSpheres(circuit, allGids, callback, *model,
                            reportMapping, perCellMaterialIds);
        }
        else
        {
            // Import meshes
            _importMeshes(callback, *model, allGids,
                          circuit.getTransforms(allGids),
                          perCellMaterialIds);

            // Import morphologies
            const auto useSimulationModel = _properties.useSimulationModel;
            model->useSimulationModel(useSimulationModel);
            if (_properties.meshFolder.empty() || useSimulationModel)
            {
                _importMorphologies(circuit, allGids, callback, *model,
                                    reportMapping, perCellMaterialIds);
            }
        }

        // Create materials
        model->createMissingMaterials();

        // Compute circuit center
        auto positions = circuit.getPositions(allGids);
        Boxf center;
        for (const auto& position : positions)
            center.merge(position);

        Transformation transformation;
        transformation.setRotationCenter(center.getCenter());
        modelDesc =
            std::make_shared<ModelDescriptor>(std::move(model), "Circuit",
                                              source, metadata);
        modelDesc->setTransformation(transformation);

        return modelDesc;
    }

private:
    size_ts _createPerCellMaterialIds(const brain::Circuit& circuit,
                                      const brain::GIDSet& gids,
                                      const size_ts& targetSizes) const
    {
        switch (_properties.colorScheme)
        {
        case ColorScheme::neuron_by_target:
        {
            size_ts ids;
            ids.reserve(gids.size());
            size_t id = 0;
            for (const auto size : targetSizes)
                std::fill_n(std::back_inserter(ids), size, id);
            return ids;
        }
        case ColorScheme::neuron_by_etype:
            return circuit.getElectrophysiologyTypes(gids);
        case ColorScheme::neuron_by_mtype:
            return circuit.getMorphologyTypes(gids);
        case ColorScheme::neuron_by_layer:
            return _getLayerIds(circuit, gids);
        default:
            return size_ts();
        }
    }

    void _addSomaSpheres(
        const brain::Circuit& circuit, const brain::GIDSet& gids,
        const LoaderProgress& callback, Model& model,
        const brain::CompartmentReportMapping* reportMapping,
        const size_ts& perCellMaterialIds) const
    {
        std::stringstream message;
        message << "Adding " << gids.size() << " spherical somas...";
        const auto positions = circuit.getPositions(gids);
        const auto mtypes = circuit.getMorphologyTypes(gids);

        std::vector<float> mtypeRadii = _getSomaRadii(circuit,
                                                      DEFAULT_SOMA_RADIUS);

        for (size_t i = 0; i != gids.size(); ++i)
        {
            const auto materialId =
                _getMaterialId(_properties.colorScheme,
                               i, brain::neuron::SectionType::soma,
                               perCellMaterialIds);
            const auto offset =
                reportMapping ? reportMapping->getOffsets()[i][0] : 0;
            model.addSphere(materialId,
                            {positions[i], mtypeRadii[mtypes[i]], offset});
            callback.updateProgress(message.str(), i / float(gids.size()));
        }

    }

    void _importMeshes(const LoaderProgress& callback BRAYNS_UNUSED,
                       Model& model BRAYNS_UNUSED,
                       const brain::GIDSet& gids BRAYNS_UNUSED,
                       const Matrix4fs& transformations BRAYNS_UNUSED,
                       const size_ts& perCellMaterialIds BRAYNS_UNUSED) const
    {
#if BRAYNS_USE_ASSIMP
        const auto colorScheme = _properties.colorScheme;
        const auto geometryQuality = _properties.geometryQuality;
        MeshLoader meshLoader(_scene);
        size_t loadingFailures = 0;
        const auto meshedMorphologiesFolder = _properties.meshFolder;
        if (meshedMorphologiesFolder.empty())
            return;

        size_t meshIndex = 0;
        // Loading meshes is currently sequential. TODO: Make it parallel!!!
        std::stringstream message;
        message << "Loading " << gids.size() << " meshes...";
        for (const auto& gid : gids)
        {
            size_t materialId = 0;
            if (colorScheme == ColorScheme::neuron_by_id)
                materialId = meshIndex;
            else if (colorScheme != ColorScheme::neuron_by_segment_type)
                materialId = perCellMaterialIds[meshIndex];

            // Load mesh from file
            const auto transformation = _properties.transformMeshes
                                            ? transformations[meshIndex]
                                            : Matrix4f();

            const auto filename =
                _getMeshFilenameFromGID(gid, _properties.meshFilenamePattern,
                                        _properties.meshFolder);
            try
            {
                meshLoader.importMesh(filename, callback, model, meshIndex,
                                      transformation, materialId, colorScheme,
                                      geometryQuality);
            }
            catch (...)
            {
                ++loadingFailures;
            }
            ++meshIndex;
            callback.updateProgress(message.str(),
                                    meshIndex / float(gids.size()));
        }
        if (loadingFailures != 0)
            BRAYNS_WARN << "Failed to import " << loadingFailures << " meshes"
                        << std::endl;
#else
        throw std::runtime_error(
            "assimp dependency is required to load meshes");
#endif
    }

    void _importMorphologies(
        const brain::Circuit& circuit, const brain::GIDSet& gids,
        const LoaderProgress& callback, Model& model,
        const brain::CompartmentReportMapping* reportMapping,
        const size_ts& perCellMaterialIds) const
    {
        MorphologyLoader morphLoader(_scene);

        std::stringstream message;
        message << "Loading " << gids.size() << " morphologies...";

        auto next = gids.begin();
        size_t index = 0;
        std::atomic_size_t current{0};

        std::atomic_size_t loadingFailures{0};
        std::exception_ptr cancelException;

        while (next != gids.end())
        {
            brain::GIDSet batch;
            for (; batch.size() < LOAD_BATCH_SIZE && next != gids.end(); ++next)
                batch.insert(*next);

            const auto morphologies =
                circuit.loadMorphologies(batch,
                                         brain::Circuit::Coordinates::global);

#pragma omp parallel for
            for (uint64_t j = 0; j < morphologies.size(); ++j)
            {
                const auto morphologyIndex = index + j;
                auto materialFunc =
                    [&](const brain::neuron::SectionType type)
                    {
                        if (_properties.useSimulationModel)
                            return size_t{0};
                        else
                            return
                                _getMaterialId(_properties.colorScheme,
                                               morphologyIndex,
                                               type, perCellMaterialIds);
                    };

                const auto& morphology = morphologies[j];

                auto data = morphLoader.processMorphology(
                    *morphology, morphologyIndex, materialFunc,
                    reportMapping, _morphologyParams);
#pragma omp critical
                data.dump(model);

                ++current;
                // Since throwing from inside an parallel for is not a good
                // idea. If the loading gets interrupted the exception will
                // be caught here and rethrown later after exiting the loop.
                try
                {
                    callback.updateProgress(message.str(),
                                            current / float(gids.size()));
                }
                catch (...)
                {
                    cancelException = std::current_exception();
                }
            }
            index += batch.size();
        }

        if (cancelException)
            std::rethrow_exception(cancelException);
    }

private:
    Scene& _scene;
    CircuitProperties _properties;
    MorphologyLoaderParams _morphologyParams;
};
}

namespace brayns
{
CircuitLoader::CircuitLoader(Scene& scene)
    : Loader(scene)
{
}

CircuitLoader::~CircuitLoader()
{
}
bool CircuitLoader::isSupported(const std::string& filename,
                                const std::string& extension
                                    BRAYNS_UNUSED) const
{
    const auto ends_with = [](const std::string& value,
                              const std::string& ending) {
        if (ending.size() > value.size())
            return false;
        return std::equal(ending.rbegin(), ending.rend(), value.rbegin());
    };

    const std::set<std::string> names = {"BlueConfig", "BlueConfig3",
                                         "CircuitConfig", "circuit"};

    for (const auto& name : names)
        if (ends_with(filename, name))
            return true;

    return false;
}

ModelDescriptorPtr CircuitLoader::importFromBlob(
    Blob&& /*blob*/, const LoaderProgress& /*callback*/,
    const PropertyMap& properties BRAYNS_UNUSED, const size_t /*index*/,
    const size_t /*materialID*/) const
{
    throw std::runtime_error("Loading circuit from blob is not supported");
}

ModelDescriptorPtr CircuitLoader::importFromFile(
    const std::string& filename, const LoaderProgress& callback,
    const PropertyMap& propertiesTmp, const size_t /*index*/,
    const size_t /*materialID*/) const
{
    // Fill property map since the actual property types are known now.
    PropertyMap properties = getProperties();
    properties.merge(propertiesTmp);
    auto impl = Impl(_scene, properties);
    try
    {
        return impl.importCircuit(filename, callback);
    }
    catch (const std::exception& error)
    {
        BRAYNS_ERROR << "Failed to open " << filename << ": " << error.what()
                     << std::endl;
        return {};
    }
}

std::string CircuitLoader::getName() const
{
    return LOADER_NAME;
}

std::vector<std::string> CircuitLoader::getSupportedExtensions() const
{
    return {"BlueConfig", "BlueConfig3", "CircuitConfig", "circuit"};
}

PropertyMap CircuitLoader::getProperties() const
{
    PropertyMap pm;
    pm.setProperty(PROP_DENSITY);
    pm.setProperty(PROP_RANDOM_SEED);
    pm.setProperty(PROP_REPORT);
    pm.setProperty(PROP_TARGETS);
    pm.setProperty(PROP_MESH_FILENAME_PATTERN);
    pm.setProperty(PROP_MESH_FOLDER);
    pm.setProperty(PROP_USE_SIMULATION_MODEL);
    pm.setProperty(PROP_TRANSFORM_MESHES);
    pm.setProperty(PROP_COLOR_SCHEME);
    pm.setProperty(PROP_START_SIMULATION_TIME);
    pm.setProperty(PROP_END_SIMULATION_TIME);
    pm.setProperty(PROP_SIMULATION_STEP);
    pm.setProperty(PROP_SYNCHRONOUS_MODE);
    pm.setProperty(PROP_GEOMETRY_QUALITY);

    { // Add all morphology loader properties
        const auto mlpm = MorphologyLoader(_scene).getProperties();
        for (const auto& prop : mlpm.getProperties())
            if (prop && !pm.hasProperty(prop->name))
                pm.setProperty(*prop);
    }

#if BRAYNS_USE_ASSIMP
    { // Add all mesh loader properties
        const auto mlpm = MeshLoader(_scene).getProperties();
        for (const auto& prop : mlpm.getProperties())
            if (prop && !pm.hasProperty(prop->name))
                pm.setProperty(*prop);
    }
#endif

    return pm;
}
}
