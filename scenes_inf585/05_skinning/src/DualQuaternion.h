//
// Created by barth on 27/01/24.
//

#ifndef INC_05_SKINNING_DUALQUATERNION_H
#define INC_05_SKINNING_DUALQUATERNION_H

#include "cgp/cgp.hpp"

using namespace cgp;

class DualQuaternion {
public:
    DualQuaternion(quaternion real, quaternion dual) : real(real), dual(dual) {}

    static DualQuaternion FromRotationTranslation(quaternion rotation, vec3 translation) {
        auto real = rotation;
        auto dual = 0.5 * quaternion(translation.x, translation.y, translation.z, 0) * rotation;

        return {real, dual};
    }

    DualQuaternion operator*(float scalar) {
        return {real * scalar, dual * scalar};
    }

    DualQuaternion operator+(DualQuaternion other) {
        return {real + other.real, dual + other.dual};
    }

    quaternion getRotation() {
        auto rotation = real;
        rotation = rotation / cgp::norm(rotation);

        return rotation;
    }

    vec3 getTranslation() {
        auto translation = 2 * dual * cgp::conjugate(getRotation());
        return vec3(translation.x, translation.y, translation.z);
    }

    matrix_stack<float, 4, 4> getMatrix() {
        auto rotation = getRotation();
        auto translation = getTranslation();

        // this was done by ChatGPT
        return {
            1 - 2 * rotation.y * rotation.y - 2 * rotation.z * rotation.z, 2 * rotation.x * rotation.y - 2 * rotation.z * rotation.w, 2 * rotation.x * rotation.z + 2 * rotation.y * rotation.w, translation.x,
            2 * rotation.x * rotation.y + 2 * rotation.z * rotation.w, 1 - 2 * rotation.x * rotation.x - 2 * rotation.z * rotation.z, 2 * rotation.y * rotation.z - 2 * rotation.x * rotation.w, translation.y,
            2 * rotation.x * rotation.z - 2 * rotation.y * rotation.w, 2 * rotation.y * rotation.z + 2 * rotation.x * rotation.w, 1 - 2 * rotation.x * rotation.x - 2 * rotation.y * rotation.y, translation.z,
            0, 0, 0, 1
        };
    }

private:
    quaternion real;
    quaternion dual;
};

#endif //INC_05_SKINNING_DUALQUATERNION_H
