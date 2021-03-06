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

#include "Model.h"

#include <brayns/common/Transformation.h>
#include <brayns/common/log.h>
#include <brayns/common/material/Material.h>
#include <brayns/common/material/Texture2D.h>
#include <brayns/common/volume/Volume.h>

#include <boost/filesystem.hpp>

#include <set>

namespace brayns
{
ModelParams::ModelParams(const std::string& path)
    : _name(boost::filesystem::basename(path))
    , _path(path)
{
}

ModelParams::ModelParams(const std::string& name, const std::string& path)
    : _name(name)
    , _path(path)
{
}

ModelDescriptor::ModelDescriptor(ModelPtr model, const std::string& path)
    : ModelParams(path)
    , _model(std::move(model))
{
}

ModelDescriptor::ModelDescriptor(ModelPtr model, const std::string& path,
                                 const ModelMetadata& metadata)
    : ModelParams(path)
    , _metadata(metadata)
    , _model(std::move(model))
{
}

ModelDescriptor::ModelDescriptor(ModelPtr model, const std::string& name,
                                 const std::string& path,
                                 const ModelMetadata& metadata)
    : ModelParams(name, path)
    , _metadata(metadata)
    , _model(std::move(model))
{
}

ModelDescriptor& ModelDescriptor::operator=(const ModelParams& rhs)
{
    if (this == &rhs)
        return *this;
    _updateValue(_boundingBox, rhs.getBoundingBox());
    if (rhs.getName().empty())
        _updateValue(_name, boost::filesystem::basename(rhs.getPath()));
    else
        _updateValue(_name, rhs.getName());
    _updateValue(_path, rhs.getPath());
    _updateValue(_visible, rhs.getVisible());

    // Transformation
    const auto oldRotationCenter = _transformation.getRotationCenter();
    const auto newRotationCenter = rhs.getTransformation().getRotationCenter();
    _updateValue(_transformation, rhs.getTransformation());
    if (newRotationCenter == Vector3f())
        // If no rotation center is specified in the model params, the one set
        // by the model loader is used
        _transformation.setRotationCenter(oldRotationCenter);

    return *this;
}

void ModelDescriptor::addInstance(const ModelInstance& instance)
{
    _instances.push_back(instance);
    _instances.rbegin()->setInstanceID(_nextInstanceID++);
    if (_model)
        _model->markInstancesDirty();
}

void ModelDescriptor::removeInstance(const size_t id)
{
    auto i = std::remove_if(_instances.begin(), _instances.end(),
                            [id](const auto& instance) {
                                return id == instance.getInstanceID();
                            });
    if (i == _instances.end())
        return;

    _instances.erase(i, _instances.end());

    if (_model)
        _model->markInstancesDirty();
}

ModelInstance* ModelDescriptor::getInstance(const size_t id)
{
    auto i = std::find_if(_instances.begin(), _instances.end(),
                          [id](const auto& instance) {
                              return id == instance.getInstanceID();
                          });
    return i == _instances.end() ? nullptr : &(*i);
}

Boxf ModelDescriptor::getInstancesBounds() const
{
    Boxf bounds;
    for (const auto& instance : getInstances())
    {
        if (!instance.getVisible() || !_model)
            continue;

        bounds.merge(
            transformBox(getModel().getBounds(), instance.getTransformation()));
    }
    return bounds;
}

bool Model::empty() const
{
    return _spheres.empty() && _cylinders.empty() && _cones.empty() &&
           _trianglesMeshes.empty() && _sdf.geometries.empty() &&
           _streamlines.empty() && _volumes.empty() &&
           _bounds.isEmpty(); // _bounds only used for volumes, will disappear.
}

uint64_t Model::addSphere(const size_t materialId, const Sphere& sphere)
{
    _spheresDirty = true;
    _spheres[materialId].push_back(sphere);
    _bounds.merge(sphere.center);
    return _spheres[materialId].size() - 1;
}

uint64_t Model::addCylinder(const size_t materialId, const Cylinder& cylinder)
{
    _cylindersDirty = true;
    _cylinders[materialId].push_back(cylinder);
    _bounds.merge(cylinder.center);
    _bounds.merge(cylinder.up);
    return _cylinders[materialId].size() - 1;
}

uint64_t Model::addCone(const size_t materialId, const Cone& cone)
{
    _conesDirty = true;
    _cones[materialId].push_back(cone);
    _bounds.merge(cone.center);
    _bounds.merge(cone.up);
    return _cones[materialId].size() - 1;
}

void Model::addStreamline(const size_t materialId, const Streamline& streamline)
{
    if (streamline.position.size() < 2)
        throw std::runtime_error(
            "Number of vertices is less than two which is minimum needed for a "
            "streamline.");

    if (streamline.position.size() != streamline.color.size())
        throw std::runtime_error("Number of vertices and colors do not match.");

    if (streamline.position.size() != streamline.radius.size())
        throw std::runtime_error("Number of vertices and radii do not match.");

    // Calculate bounds
    for (size_t index = 0; index < streamline.position.size(); ++index)
    {
        const auto& pos = streamline.position[index];
        const float radius = streamline.radius[index];
        const auto radiusVec = Vector3f(radius, radius, radius);
        _bounds.merge(pos + radiusVec);
        _bounds.merge(pos - radiusVec);
    }

    auto& streamlinesData = _streamlines[materialId];

    const size_t startIndex = streamlinesData.vertex.size();
    const size_t endIndex = startIndex + streamline.position.size() - 1;

    for (size_t index = startIndex; index < endIndex; ++index)
        streamlinesData.indices.push_back(index);

    for (size_t i = 0; i < streamline.position.size(); i++)
        streamlinesData.vertex.push_back(
            Vector4f(streamline.position[i], streamline.radius[i]));

    for (const auto& color : streamline.color)
        streamlinesData.vertexColor.push_back(color);

    _streamlinesDirty = true;
}

uint64_t Model::addSDFGeometry(const size_t materialId, const SDFGeometry& geom,
                               const std::vector<size_t>& neighbourIndices)
{
    const uint32_t geomIdx = _sdf.geometries.size();
    _sdf.geometryIndices[materialId].push_back(geomIdx);
    _sdf.neighbours.push_back(neighbourIndices);
    _sdf.geometries.push_back(geom);

    _bounds.merge(getSDFBoundingBox(geom));
    _sdfGeometriesDirty = true;

    return geomIdx;
}

void Model::updateSDFGeometryNeighbours(
    size_t geometryIdx, const std::vector<size_t>& neighbourIndices)
{
    _sdf.neighbours[geometryIdx] = neighbourIndices;
    _sdfGeometriesDirty = true;
}

void Model::addVolume(VolumePtr volume)
{
    _volumes.push_back(volume);
    _bounds.merge(volume->getBounds());
    _volumesDirty = true;
}

void Model::removeVolume(VolumePtr volume)
{
    auto i = std::find(_volumes.begin(), _volumes.end(), volume);
    if (i == _volumes.end())
        return;

    _volumes.erase(i);
    _volumesDirty = true;
}

bool Model::dirty() const
{
    return _spheresDirty || _cylindersDirty || _conesDirty ||
           _trianglesMeshesDirty || _sdfGeometriesDirty || _instancesDirty;
}

void Model::setMaterialsColorMap(const MaterialsColorMap colorMap)
{
    size_t index = 0;
    for (auto material : _materials)
    {
        material.second->setSpecularColor(Vector3f(0.f));
        material.second->setOpacity(1.f);
        material.second->setReflectionIndex(0.f);
        material.second->setEmission(0.f);

        switch (colorMap)
        {
        case MaterialsColorMap::none:
            switch (index)
            {
            case 0: // Default
            case 1: // Soma
                material.second->setDiffuseColor(Vector3f(0.9f, 0.9f, 0.9f));
                break;
            case 2: // Axon
                material.second->setDiffuseColor(Vector3f(0.2f, 0.2f, 0.8f));
                break;
            case 3: // Dendrite
                material.second->setDiffuseColor(Vector3f(0.8f, 0.2f, 0.2f));
                break;
            case 4: // Apical dendrite
                material.second->setDiffuseColor(Vector3f(0.8f, 0.2f, 0.8f));
                break;
            default:
                material.second->setDiffuseColor(
                    Vector3f(float(std::rand() % 255) / 255.f,
                             float(std::rand() % 255) / 255.f,
                             float(std::rand() % 255) / 255.f));
            }
            break;
        case MaterialsColorMap::gradient:
        {
            const float a = float(index) / float(_materials.size());
            material.second->setDiffuseColor(Vector3f(a, 0.f, 1.f - a));
            break;
        }
        case MaterialsColorMap::pastel:
            material.second->setDiffuseColor(
                Vector3f(0.5f + float(std::rand() % 127) / 255.f,
                         0.5f + float(std::rand() % 127) / 255.f,
                         0.5f + float(std::rand() % 127) / 255.f));
            break;
        case MaterialsColorMap::random:
            material.second->setDiffuseColor(
                Vector3f(float(rand() % 255) / 255.f,
                         float(rand() % 255) / 255.f,
                         float(rand() % 255) / 255.f));
            switch (rand() % 10)
            {
            case 0:
                // Transparency only
                material.second->setOpacity(float(std::rand() % 100) / 100.f);
                material.second->setRefractionIndex(1.2f);
                material.second->setSpecularColor(Vector3f(1.f));
                material.second->setSpecularExponent(10.f);
                break;
            case 1:
                // Light emission
                material.second->setEmission(std::rand() % 20);
                break;
            case 2:
                // Reflection only
                material.second->setReflectionIndex(float(std::rand() % 100) /
                                                    100.f);
                material.second->setSpecularColor(Vector3f(1.f));
                material.second->setSpecularExponent(10.f);
                break;
            case 3:
                // Reflection and refraction
                material.second->setReflectionIndex(float(std::rand() % 100) /
                                                    100.f);
                material.second->setOpacity(float(std::rand() % 100) / 100.f);
                material.second->setRefractionIndex(1.2f);
                material.second->setSpecularColor(Vector3f(1.f));
                material.second->setSpecularExponent(10.f);
                break;
            case 4:
                // Reflection and glossiness
                material.second->setReflectionIndex(float(std::rand() % 100) /
                                                    100.f);
                material.second->setSpecularColor(Vector3f(1.f));
                material.second->setSpecularExponent(10.f);
                material.second->setGlossiness(float(std::rand() % 100) /
                                               100.f);
                break;
            case 5:
                // Transparency and glossiness
                material.second->setOpacity(float(std::rand() % 100) / 100.f);
                material.second->setRefractionIndex(1.2f);
                material.second->setSpecularColor(Vector3f(1.f));
                material.second->setSpecularExponent(10.f);
                material.second->setGlossiness(float(std::rand() % 100) /
                                               100.f);
                break;
            }
            break;
        case MaterialsColorMap::shades_of_grey:
            float value = float(std::rand() % 255) / 255.f;
            material.second->setDiffuseColor(Vector3f(value, value, value));
            break;
        }
        material.second->commit();
        ++index;
    }
}

void Model::logInformation()
{
    updateSizeInBytes();

    uint64_t nbSpheres = 0;
    uint64_t nbCylinders = 0;
    uint64_t nbCones = 0;
    uint64_t nbMeshes = _trianglesMeshes.size();
    for (const auto& spheres : _spheres)
        nbSpheres += spheres.second.size();
    for (const auto& cylinders : _cylinders)
        nbCylinders += cylinders.second.size();
    for (const auto& cones : _cones)
        nbCones += cones.second.size();

    BRAYNS_DEBUG << "Spheres: " << nbSpheres << ", Cylinders: " << nbCylinders
                 << ", Cones: " << nbCones << ", Meshes: " << nbMeshes
                 << ", Memory: " << _sizeInBytes << " bytes ("
                 << _sizeInBytes / 1048576 << " MB), Bounds: " << _bounds
                 << std::endl;
}

MaterialPtr Model::getMaterial(const size_t materialId) const
{
    const auto it = _materials.find(materialId);
    if (it == _materials.end())
        throw std::runtime_error("Material " + std::to_string(materialId) +
                                 " is not registered in the model");
    return it->second;
}

void Model::updateSizeInBytes()
{
    _sizeInBytes = 0;
    for (const auto& spheres : _spheres)
        _sizeInBytes += spheres.second.size() * sizeof(Sphere);
    for (const auto& cylinders : _cylinders)
        _sizeInBytes += cylinders.second.size() * sizeof(Cylinder);
    for (const auto& cones : _cones)
        _sizeInBytes += cones.second.size() * sizeof(Cones);
    for (const auto& trianglesMesh : _trianglesMeshes)
    {
        const auto& mesh = trianglesMesh.second;
        _sizeInBytes += mesh.indices.size() * sizeof(Vector3f);
        _sizeInBytes += mesh.normals.size() * sizeof(Vector3f);
        _sizeInBytes += mesh.colors.size() * sizeof(Vector4f);
        _sizeInBytes += mesh.indices.size() * sizeof(Vector3ui);
        _sizeInBytes += mesh.textureCoordinates.size() * sizeof(Vector2f);
    }
    for (const auto& volume : _volumes)
        _sizeInBytes += volume->getSizeInBytes();
}

void Model::createMissingMaterials(const bool castSimulationData)
{
    std::set<size_t> materialIds;
    for (auto& spheres : _spheres)
        materialIds.insert(spheres.first);
    for (auto& cylinders : _cylinders)
        materialIds.insert(cylinders.first);
    for (auto& cones : _cones)
        materialIds.insert(cones.first);
    for (auto& meshes : _trianglesMeshes)
        materialIds.insert(meshes.first);
    for (auto& sdfGeometries : _sdf.geometryIndices)
        materialIds.insert(sdfGeometries.first);

    for (const auto materialId : materialIds)
    {
        const auto it = _materials.find(materialId);
        if (it == _materials.end())
        {
            auto material =
                createMaterial(materialId, std::to_string(materialId));
            material->setCastSimulationData(castSimulationData);
        }
    }
}
}
