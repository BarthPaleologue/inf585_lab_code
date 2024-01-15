#include "ffd.hpp"



using namespace cgp;

// Helper to compute permutation(n, k) - k among n
//   avoids the recursive formula (too costly)
int binomial_coeff(int n, int k)
{
    int res = 1;
    if(k>n-k)
        k = n - k;
    for(int i=0; i<k; ++i) {
        res *= (n-i);
        res /= (i+1);
    }
    return res;
}



// Precompute the FFD elements that depends only on the initial position of the shape (its barycentric coordinates in the grid)
numarray<grid_3D<float> > precompute_weights(numarray<vec3>& position, int Nx, int Ny, int Nz)
{
	numarray<grid_3D<float>> weights;
	weights.resize(position.size());

	// TO DO: compute the weights such that
	//   weights[k](kx,ky,kz) = b_kx(x_k) * b_ky(y_k) * b_zk(z_k)
	//    with (x_k,y_k,z_k) <= position[k]
	//         (kx,ky,kz) varying on the grid dimension
	//         b_kx, b_ky, b_kz being the Bezier basis functions
	//
	// Note: a grid 3D is a (array-like) structure provided by the library in order to conveniently store and access elements ordered with 3 paramters:
	// - Resizing a grid_3D:
	//   grid_3D G; G.resize(Nx, Ny, Nz);
	// - Accessing/Modifying an element at index (kx,ky,kz):
	//   element=G(kx,ky,kz); or G(kx,ky,kz)=element;

    for (int k = 0; k < position.size(); k++) {
        weights[k].resize(Nx, Ny, Nz);
        for (int kx = 0; kx < Nx; kx++) {
            for (int ky = 0; ky < Ny; ky++) {
                for (int kz = 0; kz < Nz; kz++) {
                    float b_kx = binomial_coeff(Nx-1, kx) * pow(1-position[k].x, Nx-1-kx) * pow(position[k].x, kx);
                    float b_ky = binomial_coeff(Ny-1, ky) * pow(1-position[k].y, Ny-1-ky) * pow(position[k].y, ky);
                    float b_kz = binomial_coeff(Nz-1, kz) * pow(1-position[k].z, Nz-1-kz) * pow(position[k].z, kz);
                    weights[k](kx, ky, kz) = b_kx * b_ky * b_kz;
                }
            }
        }
    }

	return weights;
}


// Computation of the FFD deformation on the position with respect to the grid
void ffd_deform(numarray<vec3>& position, grid_3D<vec3> const& grid, numarray<grid_3D<float> > const& weights)
{
	int const Nx = grid.dimension.x;
	int const Ny = grid.dimension.y;
	int const Nz = grid.dimension.z;

	int const N_vertex = position.size();

	// TO DO: compute the deformed position
	// General formulation:
	// For all position k of the shape to be deformed
	//     position[k] = sum_{x,y,z} weights[k] * grid(x,y,z)

    for (int k = 0; k < N_vertex; k++) {
        vec3 sum = vec3(0,0,0);
        for (int kx = 0; kx < Nx; kx++) {
            for (int ky = 0; ky < Ny; ky++) {
                for (int kz = 0; kz < Nz; kz++) {
                    sum += weights[k](kx, ky, kz) * grid(kx, ky, kz);
                }
            }
        }
        position[k] = sum;
    }


}
