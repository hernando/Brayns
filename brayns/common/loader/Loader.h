/* Copyright (c) 2015-2018, EPFL/Blue Brain Project
 *
 * Responsible Author: Daniel.Nachbaur@epfl.ch
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

#include <brayns/common/PropertyMap.h>
#include <brayns/common/types.h>

#include <functional>

#ifdef BRAYNS_USE_OPENMP
#include <omp.h>
#endif

struct LoaderPropertyMap
{
    std::map<std::string, std::vector<brayns::PropertyMap::Property>>
        properties;

    template <typename T>
    inline T getProperty(const std::string& name, const T& valIfNotFound) const
    {
        if (properties.find(name) == properties.end())
            return valIfNotFound;

        for (const auto& property : properties.at(name))
        {
            try
            {
                T result = property.get<T>();
                return result;
            }
            catch (const boost::bad_any_cast& /*e*/)
            {
            }
        }
        throw std::runtime_error("Could not match " + std::string(name));
        return valIfNotFound;
    }

    template <typename T>
    inline T getEnumProperty(const std::string& name,
                             const std::map<std::string, int32_t> enumValues,
                             const T& valIfNotFound) const
    {
        if (properties.find(name) == properties.end())
            return valIfNotFound;

        for (const auto& property : properties.at(name))
        {
            try
            {
                // If property is storing enum as string
                const auto enumName =
                    static_cast<std::string>(property.get<std::string>());

                for (const auto& kv : enumValues)
                    if (enumName == kv.first)
                        return static_cast<T>(kv.second);
            }
            catch (const boost::bad_any_cast& /*e*/)
            {
            }
            try
            {
                // If property is storing enum as int
                T result = static_cast<T>(property.get<int32_t>());
                return result;
            }
            catch (const boost::bad_any_cast& /*e*/)
            {
            }
        }
        throw std::runtime_error("Could not match " + std::string(name));
        return valIfNotFound;
    }

    void addProperty(const brayns::PropertyMap::Property& property)
    {
        properties[property.name].push_back(property);
    }
};

namespace brayns
{
/**
 * A class for providing progress feedback
 */
class LoaderProgress
{
public:
    /**
     * The callback for each progress update with the signature (message,
     * fraction of progress in 0..1 range)
     */
    using CallbackFn = std::function<void(const std::string&, float)>;

    LoaderProgress(CallbackFn callback)
        : _callback(std::move(callback))
    {
    }

    LoaderProgress() = default;
    ~LoaderProgress() = default;

    /**
     * Update the current progress of an operation and call the callback
     */
    void updateProgress(const std::string& message, const float fraction) const
    {
#ifdef BRAYNS_USE_OPENMP
        if (omp_get_thread_num() == 0)
#endif
            if (_callback)
                _callback(message, fraction);
    }

    CallbackFn _callback;
};

struct LoaderSupport
{
    std::string name;
    std::vector<std::string> extensions;
};

/**
 * A base class for data loaders to unify loading data from blobs and files, and
 * provide progress feedback.
 */
class Loader
{
public:
    Loader(Scene& scene)
        : _scene(scene)
    {
    }

    virtual ~Loader() = default;

    virtual LoaderSupport getLoaderSupport() const = 0;
    virtual std::pair<std::string, PropertyMap> getLoaderProperties() const = 0;

    /**
     * Import the data from the blob and return the created model.
     *
     * @param blob the blob containing the data to import
     * @param index Index of the element, mainly used for material
     * assignment
     * @param defaultMaterialId the default material to use
     * @return the model that has been created by the loader
     */
    virtual ModelDescriptorPtr importFromBlob(
        Blob&& blob, const LoaderProgress& callback,
        const LoaderPropertyMap& properties, const size_t index = 0,
        const size_t defaultMaterialId = NO_MATERIAL) const = 0;

    /**
     * Import the data from the given file and return the created model.
     *
     * @param filename the file containing the data to import
     * @param index Index of the element, mainly used for material assignment
     * @param defaultMaterialId the default material to use
     * @return the model that has been created by the loader
     */
    virtual ModelDescriptorPtr importFromFile(
        const std::string& filename, const LoaderProgress& callback,
        const LoaderPropertyMap& properties, const size_t index = 0,
        const size_t defaultMaterialId = NO_MATERIAL) const = 0;

    /**
     * Query the loader if it can load the given file
     */
    virtual bool isSupported(const std::string& filename,
                             const std::string& extension) const = 0;

protected:
    Scene& _scene;
};
}
