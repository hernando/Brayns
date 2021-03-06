/* Copyright (c) 2018, EPFL/Blue Brain Project
 * All rights reserved. Do not distribute without permission.
 * Responsible Author: Cyrille Favreau <cyrille.favreau@epfl.ch>
 *                     Grigori Chevtchenko <grigori.chevtchenko@epfl.ch>
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

#include "StereoCamera.ih"
#include "math/LinearSpace.ih"
#include "math/math.ih"
#include "math/sampling.ih"

void StereoCamera_initRay(uniform Camera *uniform _self, varying Ray &ray,
                          const varying CameraSample &sample)
{
    uniform StereoCamera *uniform self = (uniform StereoCamera * uniform)_self;

    vec2f screen = sample.screen;
    screen = Camera_subRegion(_self, screen);
    vec3f pixelPos;
    vec3f org = make_vec3f(0.f);
    if (self->side_by_side)
    {
        if (screen.x < 0.5f)
        {
            pixelPos.x = (2.f * screen.x - 0.5f) * self->imgPlane_size_x;
            org.x = -self->interpupillaryDistance * 0.5f;
        }
        else
        {
            pixelPos.x = (2.f * screen.x - 1.5f) * self->imgPlane_size_x;
            org.x = self->interpupillaryDistance * 0.5f;
        }
    }
    else
    {
        pixelPos.x = (screen.x - 0.5f) * self->imgPlane_size_x;
        org.x = self->interpupillaryDistance * 0.5f;
    }

    pixelPos.y = (screen.y - 0.5f) * self->imgPlane_size_y;
    pixelPos.z = -self->distanceToPlane;

    vec3f dir = normalize(pixelPos - org);

    const LinearSpace3f cameraSpace =
        make_LinearSpace3f(self->dir_du, self->dir_dv, self->dir_cam);
    dir = normalize(cameraSpace * dir);
    org = cameraSpace * org;

    const float time = Camera_shutterTime(_self, sample.time);
    org = self->org + org;
    setRay(ray, org, normalize(dir), self->super.nearClip, infinity, time);
}

/*! create a new ispc-side version of a StereoCamera - with given
    pointer to the c-class version - and return both class pointer and
    pointer to internal data back via the supplied reference
    parameters */
export void *uniform StereoCamera_create(void *uniform cppE)
{
    uniform StereoCamera *uniform self = uniform new uniform StereoCamera;
    self->super.cppEquivalent = cppE;
    self->super.initRay = StereoCamera_initRay;
    self->super.doesDOF = false;
    return self;
}

export void StereoCamera_set(
    void *uniform _self, const uniform vec3f &org, const uniform vec3f &dir_cam,
    const uniform vec3f &dir_du, const uniform vec3f &dir_dv,
    const uniform float distanceToPlane, const uniform float imgPlane_size_y,
    const uniform float imgPlane_size_x,
    const uniform float interpupillaryDistance, const uniform bool side_by_side)
{
    uniform StereoCamera *uniform self = (uniform StereoCamera * uniform)_self;
    self->org = org;
    self->dir_cam = dir_cam;
    self->dir_du = dir_du;
    self->dir_dv = dir_dv;
    self->distanceToPlane = distanceToPlane;
    self->imgPlane_size_y = imgPlane_size_y;
    self->imgPlane_size_x = imgPlane_size_x;
    self->interpupillaryDistance = interpupillaryDistance;
    self->side_by_side = side_by_side;
}
