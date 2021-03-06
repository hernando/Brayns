/* Copyright (c) 2015-2016, EPFL/Blue Brain Project
 * All rights reserved. Do not distribute without permission.
 * Author: Jafet Villafranca Diaz <jafet.villafrancadiaz@epfl.ch>
 *
 * Ray-cone intersection:
 * based on Ching-Kuang Shene (Graphics Gems 5, p. 227-230)
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

// ospray
#include "ospray/SDK/common/Model.ih"
#include "ospray/SDK/common/Ray.ih"
#include "ospray/SDK/geometry/Geometry.ih"
#include "ospray/SDK/math/box.ih"
#include "ospray/SDK/math/vec.ih"
// embree
#include "embree2/rtcore.isph"
#include "embree2/rtcore_geometry_user.isph"
#include "embree2/rtcore_scene.isph"

struct ExtendedCones
{
    uniform Geometry geometry;

    uniform uint8 *uniform data;

    float radius;
    float length;
    int materialID;
    int offset_center;
    int offset_up;
    int offset_centerRadius;
    int offset_upRadius;
    int offset_timestamp;
    int offset_texture_coords;
    int offset_materialID;
    int32 numExtendedCones;
    int32 bytesPerCone;
};

void ExtendedCones_bounds(uniform ExtendedCones *uniform geometry,
                          uniform size_t primID, uniform box3fa &bbox)
{
    uniform uint8 *uniform conePtr =
        geometry->data + geometry->bytesPerCone * primID;
    uniform float extent = geometry->radius;
    if (geometry->offset_centerRadius >= 0)
        extent = *((uniform float *)(conePtr + geometry->offset_centerRadius));

    if (geometry->offset_upRadius >= 0)
    {
        uniform float upRadius =
            *((uniform float *)(conePtr + geometry->offset_upRadius));
        if (upRadius > extent)
            extent = upRadius;
    }
    uniform vec3f v0 = *((uniform vec3f *)(conePtr + geometry->offset_center));
    uniform vec3f v1 = *((uniform vec3f *)(conePtr + geometry->offset_up));
    bbox = make_box3fa(min(v0, v1) - make_vec3f(extent),
                       max(v0, v1) + make_vec3f(extent));
}

void ExtendedCones_intersect(uniform ExtendedCones *uniform geometry,
                             varying Ray &ray, uniform size_t primID)
{
    uniform uint8 *uniform conePtr =
        geometry->data + geometry->bytesPerCone * primID;

    uniform float radius0 = geometry->radius;
    if (geometry->offset_centerRadius >= 0)
        radius0 = *((uniform float *)(conePtr + geometry->offset_centerRadius));

    uniform float radius1 = geometry->radius;
    if (geometry->offset_upRadius >= 0)
        radius1 = *((uniform float *)(conePtr + geometry->offset_upRadius));

    uniform float timestamp =
        *((uniform float *)(conePtr + geometry->offset_timestamp));

    if (ray.time > 0 && timestamp > ray.time)
        return;

    uniform vec3f v0 = *((uniform vec3f *)(conePtr + geometry->offset_center));
    uniform vec3f v1 = *((uniform vec3f *)(conePtr + geometry->offset_up));

    if (radius0 < radius1)
    {
        // swap radii and positions, so radius0 and v0 are always at the bottom
        uniform float tmpRadius = radius1;
        radius1 = radius0;
        radius0 = tmpRadius;

        uniform vec3f tmpPos = v1;
        v1 = v0;
        v0 = tmpPos;
    }

    const vec3f upVector = v1 - v0;
    const float upLength = length(upVector);

    // Compute the height of the full cone, in order to obtain its vertex
    const float deltaRadius = radius0 - radius1;
    const float tanA = deltaRadius / upLength;
    const float coneHeight = radius0 / tanA;
    const float squareTanA = tanA * tanA;
    const float div = sqrtf(1.f + squareTanA);
    if (div == 0.f)
        return;
    const float cosA = 1.f / div;

    const vec3f V = v0 + normalize(upVector) * coneHeight;
    const vec3f v = normalize(v0 - V);

    // Normal of the plane P determined by V and ray
    vec3f n = normalize(cross(ray.dir, V - ray.org));
    const float dotNV = dot(n, v);
    if (dotNV > 0.f)
        n = neg(n);

    const float squareCosTheta = 1.f - dotNV * dotNV;
    const float cosTheta = sqrtf(squareCosTheta);
    if (cosTheta < cosA)
        return; // no intersection

    if (squareCosTheta == 0.f)
        return;

    const float squareTanTheta = (1.f - squareCosTheta) / squareCosTheta;
    const float tanTheta = sqrtf(squareTanTheta);

    // Compute u-v-w coordinate system
    const vec3f u = normalize(cross(v, n));
    const vec3f w = normalize(cross(u, v));

    // Circle intersection of cone with plane P
    const vec3f uComponent = sqrtf(squareTanA - squareTanTheta) * u;
    const vec3f vwComponent = v + tanTheta * w;
    const vec3f delta1 = vwComponent + uComponent;
    const vec3f delta2 = vwComponent - uComponent;
    const vec3f rayApex = V - ray.org;

    const vec3f normal1 = cross(ray.dir, delta1);
    const float length1 = length(normal1);

    if (length1 == 0.f)
        return;

    const float r1 = dot(cross(rayApex, delta1), normal1) / (length1 * length1);

    const vec3f normal2 = cross(ray.dir, delta2);
    const float length2 = length(normal2);

    if (length2 == 0.f)
        return;

    const float r2 = dot(cross(rayApex, delta2), normal2) / (length2 * length2);

    float t_in = r1;
    float t_out = r2;
    if (r2 > 0.f)
    {
        if (r1 > 0.f)
        {
            if (r1 > r2)
            {
                t_in = r2;
                t_out = r1;
            }
        }
        else
            t_in = r2;
    }

    if (t_in > ray.t0 && t_in < ray.t)
    {
        const vec3f p1 = ray.org + t_in * ray.dir;
        // consider only the parts within the extents of the truncated cone
        if (dot(p1 - v1, v) > 0.f && dot(p1 - v0, v) < 0.f)
        {
            ray.primID = primID;
            ray.geomID = geometry->geometry.geomID;
            ray.t = t_in;
            const vec3f surfaceVec = normalize(p1 - V);
            ray.Ng = cross(cross(v, surfaceVec), surfaceVec);
            return;
        }
    }
    if (t_out > ray.t0 && t_out < ray.t)
    {
        const vec3f p2 = ray.org + t_out * ray.dir;
        // consider only the parts within the extents of the truncated cone
        if (dot(p2 - v1, v) > 0.f && dot(p2 - v0, v) < 0.f)
        {
            ray.primID = primID;
            ray.geomID = geometry->geometry.geomID;
            ray.t = t_out;
            const vec3f surfaceVec = normalize(p2 - V);
            ray.Ng = cross(cross(v, surfaceVec), surfaceVec);
        }
    }
    return;
}

static void ExtendedCones_postIntersect(uniform Geometry *uniform geometry,
                                        uniform Model *uniform model,
                                        varying DifferentialGeometry &dg,
                                        const varying Ray &ray,
                                        uniform int64 flags)
{
    uniform ExtendedCones *uniform this =
        (uniform ExtendedCones * uniform)geometry;
    dg.geometry = geometry;
#if ((OSPRAY_VERSION_MAJOR == 1) && (OSPRAY_VERSION_MINOR < 5))
    dg.material = geometry->material;
#endif
    vec3f Ng = ray.Ng;
    vec3f Ns = Ng;

    uniform uint8 *conePtr = this->data + this->bytesPerCone * ray.primID;
    // Store value as texture coordinate
    const vec2f tex =
        *((uniform vec2f *)(conePtr + this->offset_texture_coords));

    dg.st.x = tex.x;
    dg.st.y = tex.y;

    if (flags & DG_NORMALIZE)
    {
        Ng = normalize(Ng);
        Ns = normalize(Ns);
    }
    if (flags & DG_FACEFORWARD)
    {
        if (dot(ray.dir, Ng) >= 0.f)
            Ng = neg(Ng);
        if (dot(ray.dir, Ns) >= 0.f)
            Ns = neg(Ns);
    }
    if ((flags & DG_MATERIALID) && (this->offset_materialID >= 0))
        dg.materialID =
            *((uniform uint32 * varying)(conePtr + this->offset_materialID));
    dg.Ng = Ng;
    dg.Ns = Ns;
}

export void *uniform ExtendedCones_create(void *uniform cppEquivalent)
{
    uniform ExtendedCones *uniform geom = uniform new uniform ExtendedCones;
    Geometry_Constructor(&geom->geometry, cppEquivalent,
                         ExtendedCones_postIntersect, 0, 0, 0);
    return geom;
}

export void ExtendedConesGeometry_set(
    void *uniform _geom, void *uniform _model, void *uniform data,
    int uniform numExtendedCones, int uniform bytesPerCone,
    float uniform radius, float uniform length, int uniform materialID,
    int uniform offset_center, int uniform offset_up,
    int uniform offset_centerRadius, int uniform offset_upRadius,
    int uniform offset_timestamp, int uniform offset_texture_coords,
    int uniform offset_materialID)
{
    uniform ExtendedCones *uniform geom =
        (uniform ExtendedCones * uniform)_geom;
    uniform Model *uniform model = (uniform Model * uniform)_model;

    uniform uint32 geomID =
        rtcNewUserGeometry(model->embreeSceneHandle, numExtendedCones);

    geom->geometry.model = model;
    geom->geometry.geomID = geomID;
    geom->numExtendedCones = numExtendedCones;
    geom->radius = radius;
    geom->length = length;
    geom->data = (uniform uint8 * uniform)data;
    geom->materialID = materialID;
    geom->bytesPerCone = bytesPerCone;

    geom->offset_center = offset_center;
    geom->offset_up = offset_up;
    geom->offset_centerRadius = offset_centerRadius;
    geom->offset_upRadius = offset_upRadius;
    geom->offset_timestamp = offset_timestamp;
    geom->offset_texture_coords = offset_texture_coords;
    geom->offset_materialID = offset_materialID;

    rtcSetUserData(model->embreeSceneHandle, geomID, geom);
    rtcSetBoundsFunction(model->embreeSceneHandle, geomID,
                         (uniform RTCBoundsFunc)&ExtendedCones_bounds);
    rtcSetIntersectFunction(
        model->embreeSceneHandle, geomID,
        (uniform RTCIntersectFuncVarying)&ExtendedCones_intersect);
    rtcSetOccludedFunction(
        model->embreeSceneHandle, geomID,
        (uniform RTCOccludedFuncVarying)&ExtendedCones_intersect);
    rtcEnable(model->embreeSceneHandle, geomID);
}
