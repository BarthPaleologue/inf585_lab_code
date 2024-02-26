#include "simulation.hpp"

using namespace cgp;

// Convert a density value to a pressure
float density_to_pressure(float rho, float rho0, float stiffness) {
    return stiffness * (rho - rho0);
}

float W_laplacian_viscosity(vec3 const &p_i, vec3 const &p_j, float h) {
    return 15.0f * std::pow(h - norm(p_i - p_j), 3.0) / (2.0f * 3.14159f * std::pow(h, 6.0));
}

vec3 W_gradient_pressure(vec3 const &p_i, vec3 const &p_j, float h) {
    vec3 gradient = (p_i - p_j) / norm(p_i - p_j);

    float const r = norm(p_i - p_j);

    gradient *= -3.0f * 315.0f * 2.0f * r * std::pow(h * h - r * r, 2.0f) / (64.0f * 3.14159f * std::pow(h, 9));

    return gradient;
}

float W_density(vec3 const &p_i, const vec3 &p_j, float h) {
    float const r = norm(p_i - p_j);
    //assert_cgp_no_msg(r <= h);
    return 315.0f * std::pow(h * h - r * r, 3.0f) / (64.0f * 3.14159f * std::pow(h, 9));
}


void update_density(numarray<particle_element> &particles, float h, float m) {
    // To do: Compute the density value (particles[i].rho) at each particle position
    //  rho_i = \sum_j m W_density(pi,pj)
    int const N = particles.size();
    for (int i = 0; i < N; ++i) {
        float rho = 0;
        for (int j = 0; j < N; ++j) {
            if (norm(particles[i].p - particles[j].p) > h) continue;
            rho += m * W_density(particles[i].p, particles[j].p, h);
        }
        particles[i].rho = rho;
    }
}

// Convert the particle density to pressure
void update_pressure(numarray<particle_element> &particles, float rho0, float stiffness) {
    const int N = particles.size();
    for (int i = 0; i < N; ++i) {
        particles[i].pressure = density_to_pressure(particles[i].rho, rho0, stiffness);
    }
}

// Compute the forces and update the acceleration of the particles
void update_force(numarray<particle_element> &particles, float h, float m, float nu) {
    const int N = particles.size();
    for (int i = 0; i < N; ++i) {
        // gravity
        particles[i].f = m * vec3{0, -9.81f, 0};

        // pressure
        vec3 pressure_force = vec3();
        for (int j = 0; j < N; ++j) {
            if (i == j) continue;

            vec3 const &p_i = particles[i].p;
            vec3 const &p_j = particles[j].p;
            float const rho_j = particles[j].rho;
            float const p_i_pressure = particles[i].pressure;
            float const p_j_pressure = particles[j].pressure;

            if (norm(p_i - p_j) > h) continue;
            pressure_force += m * (p_i_pressure + p_j_pressure) * W_gradient_pressure(p_i, p_j, h) / (2.0f * rho_j);
        }

        particles[i].f += -m * pressure_force / particles[i].rho;

        // viscosity
        vec3 viscosity_force = vec3();
        for (int j = 0; j < N; ++j) {
            if (i == j) continue;

            vec3 const p_i_viscosity = particles[i].v;
            vec3 const p_j_viscosity = particles[j].v;

            float const rho_j = particles[j].rho;

            if (norm(particles[i].p - particles[j].p) > h) continue;
            viscosity_force +=
                    m * (p_j_viscosity - p_i_viscosity) * W_laplacian_viscosity(particles[i].p, particles[j].p, h) /
                    rho_j;
        }

        particles[i].f += m * nu * viscosity_force;
    }
}

void simulate(float dt, numarray<particle_element> &particles, sph_parameters_structure const &sph_parameters) {

    // Update values
    update_density(particles, sph_parameters.h, sph_parameters.m);                   // First compute updated density
    update_pressure(particles, sph_parameters.rho0, sph_parameters.stiffness);       // Compute associated pressure
    update_force(particles, sph_parameters.h, sph_parameters.m, sph_parameters.nu);  // Update forces

    // Numerical integration
    float const damping = 0.005f;
    int const N = particles.size();
    float const m = sph_parameters.m;
    for (int k = 0; k < N; ++k) {
        vec3 &p = particles[k].p;
        vec3 &v = particles[k].v;
        vec3 &f = particles[k].f;

        v = (1 - damping) * v + dt * f / m;
        p = p + dt * v;
    }


    // Collision
    float const epsilon = 1e-3f;
    for (int k = 0; k < N; ++k) {
        vec3 &p = particles[k].p;
        vec3 &v = particles[k].v;

        // small perturbation to avoid alignment
        if (p.y < -1) {
            p.y = -1 + epsilon * rand_uniform();
            v.y *= -0.5f;
        }
        if (p.x < -1) {
            p.x = -1 + epsilon * rand_uniform();
            v.x *= -0.5f;
        }
        if (p.x > 1) {
            p.x = 1 - epsilon * rand_uniform();
            v.x *= -0.5f;
        }
    }

}