/* Copyright (c) 2018, EPFL/Blue Brain Project
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

#define BOOST_TEST_MODULE braynsClipPlanes

#include <jsonSerialization.h>

#include <tests/paths.h>

#include "ClientServer.h"
#include "PDiffHelpers.h"

const std::string ADD_CLIP_PLANE("add-clip-plane");
const std::string GET_CLIP_PLANES("get-clip-planes");
const std::string REMOVE_CLIP_PLANES("remove-clip-planes");
const std::string UPDATE_CLIP_PLANE("update-clip-plane");

BOOST_AUTO_TEST_SUITE(clip_plane_api)

BOOST_GLOBAL_FIXTURE(ClientServer);

BOOST_AUTO_TEST_CASE(add_plane)
{
    BOOST_REQUIRE(getScene().getClipPlanes().empty());
    const brayns::Plane equation{{1.0, 2.0, 3.0, 4.0}};
    const auto result =
        makeRequest<brayns::Plane, brayns::ClipPlane>(ADD_CLIP_PLANE, equation);
    BOOST_CHECK_EQUAL(result.getID(), 0);
    BOOST_CHECK(result.getPlane() == equation);
    BOOST_REQUIRE_EQUAL(getScene().getClipPlanes().size(), 1);
    BOOST_CHECK(getScene().getClipPlane(0)->getPlane() == equation);
    BOOST_CHECK(getScene().getClipPlane(1) == brayns::ClipPlanePtr());

    getScene().removeClipPlane(0);
    BOOST_CHECK(getScene().getClipPlanes().empty());
}

BOOST_AUTO_TEST_CASE(get_planes)
{
    const brayns::Plane equation1{{1.0, 1.0, 1.0, 1.0}};
    const brayns::Plane equation2{{2.0, 2.0, 2.0, 2.0}};

    const auto id1 = getScene().addClipPlane(equation1);
    const auto id2 = getScene().addClipPlane(equation2);

    const auto result = makeRequest<brayns::ClipPlanes>(GET_CLIP_PLANES);
    BOOST_CHECK_EQUAL(result.size(), 2);
    BOOST_CHECK(result[0]->getPlane() == equation1);
    BOOST_CHECK(result[1]->getPlane() == equation2);

    getScene().removeClipPlane(id1);
    getScene().removeClipPlane(id2);
}

BOOST_AUTO_TEST_CASE(update_plane)
{
    Client client(ClientServer::instance());

    const brayns::Plane equation1{{1.0, 1.0, 1.0, 1.0}};
    const brayns::Plane equation2{{2.0, 2.0, 2.0, 2.0}};

    bool called = false;
    const auto id1 = getScene().addClipPlane(equation1);
    client.client.connect<brayns::ClipPlane>(
        UPDATE_CLIP_PLANE,
        [id1, equation2, &called](const brayns::ClipPlane& plane) {
            BOOST_CHECK_EQUAL(plane.getID(), id1);
            BOOST_CHECK(plane.getPlane() == equation2);
            called = true;
        });
    process();

    makeRequest<brayns::ClipPlane, bool>(UPDATE_CLIP_PLANE,
                                         brayns::ClipPlane(id1, equation2));
    client.process(); // This is for the scene update
    client.process();

    getScene().removeClipPlane(id1);
}

BOOST_AUTO_TEST_CASE(remove_planes)
{
    const brayns::Plane equation{{1.0, 2.0, 3.0, 4.0}};
    const auto id1 = getScene().addClipPlane(equation);
    const auto id2 = getScene().addClipPlane(equation);
    const auto id3 = getScene().addClipPlane(equation);
    makeRequest<size_ts, bool>(REMOVE_CLIP_PLANES, {id2});
    makeRequest<size_ts, bool>(REMOVE_CLIP_PLANES, {id1, id3});
    BOOST_CHECK(getScene().getClipPlanes().empty());
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(clip_plane_rendering)

void testClipping(bool orthographic = false)
{
    const char* argv[] = {"clipPlanes",    "demo", "--disable-accumulation",
                          "--window-size", "200",  "200"};
    const int argc = sizeof(argv) / sizeof(char*);

    brayns::Brayns brayns(argc, argv);
    const std::string original =
        orthographic ? "demo_ortho.png" : "snapshot.png";

    const std::string clipped = orthographic ? "demo_clipped_ortho.png"
                                             : "demo_clipped_perspective.png";

    auto& engine = brayns.getEngine();
    auto& scene = engine.getScene();
    auto& camera = engine.getCamera();

    camera.setInitialState(scene.getBounds());
    if (orthographic)
        camera.setCurrentType("orthographic");
    brayns.commitAndRender();
    BOOST_CHECK(compareTestImage(original, engine.getFrameBuffer()));

    auto id1 = scene.addClipPlane({1.0, 0.0, 0.0, -0.5});
    auto id2 = scene.addClipPlane({0.0, -1.0, 0.0, 0.5});
    brayns.commitAndRender();
    BOOST_CHECK(compareTestImage(clipped, engine.getFrameBuffer()));

    scene.removeClipPlane(id1);
    scene.removeClipPlane(id2);
    brayns.commitAndRender();
    BOOST_CHECK(
        compareTestImage(original, engine.getFrameBuffer()));

    id1 = scene.addClipPlane({1.0, 0.0, 0.0, -0.5});
    id2 = scene.addClipPlane({0.0, 1.0, 0.0, 0.5});
    scene.getClipPlane(id2)->setPlane({0.0, -1.0, 0.0, 0.5});
    brayns.commitAndRender();
    BOOST_CHECK(compareTestImage(clipped, engine.getFrameBuffer()));
}

BOOST_AUTO_TEST_CASE(perspective)
{
    testClipping();
}

BOOST_AUTO_TEST_CASE(orthographic)
{
    testClipping(true);
}

BOOST_AUTO_TEST_SUITE_END()
