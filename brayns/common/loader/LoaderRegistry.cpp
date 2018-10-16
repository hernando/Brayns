/* Copyright (c) 2015-2018, EPFL/Blue Brain Project
 * All rights reserved. Do not distribute without permission.
 * Responsible Author: Daniel.Nachbaur <daniel.nachbaur@epfl.ch>
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

#include "LoaderRegistry.h"

#include <brayns/common/utils/Utils.h>

#include <boost/filesystem.hpp>
#include <boost/range/iterator_range_core.hpp>
namespace fs = boost::filesystem;

namespace brayns
{
void LoaderRegistry::registerLoader(std::unique_ptr<Loader> loader)
{
    _loaders.push_front(std::move(loader));
}

std::vector<LoaderSupport> LoaderRegistry::getLoaderSupport() const
{
    std::vector<LoaderSupport> ret(_loaders.size());

    std::transform(_loaders.begin(), _loaders.end(), ret.begin(),
                   [](const auto& loader) {
                       return loader->getLoaderSupport();
                   });

    ret.erase(std::remove_if(ret.begin(), ret.end(),
                             [](const auto& p) { return p.name == "archive"; }),
              ret.end());
    return ret;
}

std::vector<std::pair<std::string, PropertyMap>>
    LoaderRegistry::generateLoaderPropertyMaps() const
{
    std::vector<std::pair<std::string, PropertyMap>> ret(_loaders.size());

    std::transform(_loaders.begin(), _loaders.end(), ret.begin(),
                   [](const auto& loader) {
                       return loader->getLoaderProperties();
                   });

    ret.erase(std::remove_if(ret.begin(), ret.end(),
                             [](const auto& p) {
                                 return p.first == "archive";
                             }),
              ret.end());
    return ret;
}

bool LoaderRegistry::isSupportedFile(const std::string& filename) const
{
    if (fs::is_directory(filename))
        return false;

    const auto extension = extractExtension(filename);
    for (const auto& loader : _loaders)
        if (loader->isSupported(filename, extension))
            return true;
    return false;
}

bool LoaderRegistry::isSupportedType(const std::string& type) const
{
    for (const auto& loader : _loaders)
        if (loader->isSupported("", type))
            return true;
    return false;
}

const Loader& LoaderRegistry::getSuitableLoader(
    const std::string& filename, const std::string& filetype,
    const std::string& loaderName) const
{
    if (fs::is_directory(filename))
        throw std::runtime_error("'" + filename + "' is a directory");

    const auto extension =
        filetype.empty() ? extractExtension(filename) : filetype;

    // If we have a specific loader we first check if it is an archive
    if (!loaderName.empty())
    {
        Loader* archiveLoader = nullptr;
        for (const auto& loader : _loaders)
            if (loader->getLoaderSupport().name == "archive")
                archiveLoader = loader.get();

        if (archiveLoader && archiveLoader->isSupported(filename, extension))
            return *archiveLoader;

        // Do not use archive loader, find specified one
        for (const auto& loader : _loaders)
            if (loader->getLoaderSupport().name == loaderName)
                return *loader.get();

        throw std::runtime_error("No loader found with name '" + loaderName +
                                 "'");
    }

    for (const auto& loader : _loaders)
        if (loader->isSupported(filename, extension))
            return *loader;

    throw std::runtime_error("No loader found for filename '" + filename +
                             "' and filetype '" + filetype + "'");
}
}
