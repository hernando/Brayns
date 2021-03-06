/* Copyright (c) 2015-2018, EPFL/Blue Brain Project
 * All rights reserved. Do not distribute without permission.
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

#include <brayns/common/tasks/TaskFunctor.h>
#include <brayns/common/types.h>

namespace brayns
{
/**
 * A task functor which loads data from blob or file path and adds the loaded
 * data to the scene.
 */
class LoadModelFunctor : public TaskFunctor
{
public:
    LoadModelFunctor(EnginePtr engine);
    LoadModelFunctor(LoadModelFunctor&&) = default;
    ModelDescriptorPtr operator()(Blob&& blob);
    ModelDescriptorPtr operator()(const std::string& path);

private:
    ModelDescriptorPtr _performLoad(
        const std::function<ModelDescriptorPtr()>& loadData);

    ModelDescriptorPtr _loadData(Blob&& blob);
    ModelDescriptorPtr _loadData(const std::string& path);

    void _updateProgress(const std::string& message, const size_t increment);

    std::function<void(std::string, float)> _getProgressFunc();

    EnginePtr _engine;
    size_t _currentProgress{0};
    size_t _nextTic{0};
};
}
