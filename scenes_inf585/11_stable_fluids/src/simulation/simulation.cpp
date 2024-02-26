#include "simulation.hpp"

using namespace cgp;


void divergence_free(grid_2D<vec2> &new_velocity, grid_2D<vec2> const &velocity, grid_2D<float> &divergence,
                     grid_2D<float> &gradient_field) {
    // v = projection of v0 on divergence free vector field
    //
    // v : Final vector field to be filled
    // v0: Initial vector field (non divergence free)
    // divergence: temporary buffer used to compute the divergence of v0
    // gradient_field: temporary buffer used to compute v = v0 - nabla(gradient_field)


    // TO do:
    // 1. Compute divergence of v0
    // 2. Compute gradient_field such that nabla(gradient_field)^2 = div(v0)
    // 3. Compute v = v0 - nabla(gradient_field)

    int const N = int(velocity.dimension.x);

    for (int k = 0; k < 10; k++) {
        // compute divergences
        for(int x = 1; x < N - 1; ++x) {
            for(int y = 1; y < N - 1; ++y) {
                divergence(x, y) = (velocity(x + 1, y).x - velocity(x - 1, y).x +
                                    velocity(x, y + 1).y - velocity(x, y - 1).y) / 2.0f;
            }
        }

        // compute gradients
        for (int x = 1; x < N - 1; ++x) {
            for (int y = 1; y < N - 1; ++y) {
                float divergence_v0 = divergence(x, y);

                gradient_field(x, y) = (gradient_field(x - 1, y) + gradient_field(x + 1, y) +
                                        gradient_field(x, y - 1) + gradient_field(x, y + 1) -
                                        divergence_v0) / 4.0f;
            }
        }
    }

    // compute final velocity
    for (int x = 1; x < N - 1; ++x) {
        for (int y = 1; y < N - 1; ++y) {
            vec2 nabla_gradient_field = vec2(gradient_field(x + 1, y) - gradient_field(x - 1, y),
                                             gradient_field(x, y + 1) - gradient_field(x, y - 1)) / 2.0f;
            new_velocity(x, y) = velocity(x, y) - nabla_gradient_field;
        }
    }

}