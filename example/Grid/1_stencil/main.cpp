#include "Grid/grid_dist_id.hpp"
#include "data_type/aggregate.hpp"
#include "Decomposition/CartDecomposition.hpp"

/*
 * ### WIKI 1 ###
 *
 * ## Simple example
 *
 * This example show how to move grid_key in order to create a Laplacian stencil,
 * be careful, the function move are convenient, but not the fastest implementation
 *
 * ### WIKI END ###
 *
 */

/*
 *
 * ### WIKI 2 ###
 *
 * Define some convenient constants and types
 *
 */
constexpr size_t x = 0;
constexpr size_t y = 1;
constexpr size_t z = 2;

constexpr size_t A = 0;
constexpr size_t B = 0;

typedef aggregate<float[3],float[3]> grid_point;

int main(int argc, char* argv[])
{
	//
	// ### WIKI 3 ###
	//
	// Initialize the library and several objects
	//
	openfpm_init(&argc,&argv);

	//
	// ### WIKI 4 ###
	//
	// Create several object needed later, in particular
	// * A 3D box that define the domain
	// * an array of 3 unsigned integer that define the size of the grid on each dimension
	// * A Ghost object that will define the extension of the ghost part for each sub-domain in physical units

	Box<3,float> domain({0.0,0.0,0.0},{1.0,1.0,1.0});
	size_t sz[3];
	sz[0] = 100;
	sz[1] = 100;
	sz[2] = 100;

	// Ghost
	Ghost<3,float> g(0.03);

	//
	// ### WIKI 4 ###
	//
	// Create a distributed grid in 3D (1° template parameter) space in with float precision (2° template parameter)
	// each grid point contain a vector of dimension 3 (float[3]),
	// using a CartesianDecomposition strategy (4° parameter) (the parameter 1° and 2° inside CartDecomposition must match 1° and 2°
	// of grid_dist_id)
	//
	// Constructor parameters:
	//
	// * sz: size of the grid on each dimension
	// * domain: where the grid is defined
	// * g: ghost extension
	//
	grid_dist_id<3, float, grid_point> g_dist(sz,domain,g);

	// ### WIKI 5 ###
	//
	// Get an iterator that go throught the point of the domain (No ghost)
	//
	auto dom = g_dist.getDomainIterator();
	
	// ### WIKI END ###

	while (dom.isNext())
	{
		//
		// ### WIKI 6 ###
		//
		// Get the local grid key, the local grid key store internally the sub-domain id (each sub-domain contain a grid)
		// and the local grid point id identified by 2 integers in 2D 3 integer in 3D and so on. These two distinct elements are
		// available with key.getSub() and key.getKey()
		//
		auto key = dom.get();

		//
		// ### WIKI 7 ###
		//
		// Here we convert the local grid position, into global position, key_g contain 3 integers that identify the position
		// of the grid point in global coordinates
		//
		//
		auto key_g = g_dist.getGKey(key);

		//
		// ### WIKI 8 ###
		//
		// we write on the grid point of position (i,j,k) the value i*i + j*j + k*k on the component [0] of the vector
		g_dist.template get<0>(key)[0] = key_g.get(0)*key_g.get(0) + key_g.get(1)*key_g.get(1) + key_g.get(2)*key_g.get(2);

		//
		// ### WIKI 9 ###
		//
		// next point

		++dom;

		// ### WIKI END ###
	}

	//
	// ### WIKI 10 ###
	//
	// Each sub-domain has an extended part, that is materially contained from another processor that in general is not synchronized
	// ghost_get<0> synchronize the property 0 (the vector) in the ghost part
	//
	//
	g_dist.template ghost_get<0>();
	
	//
	// ### WIKI 11 ###
	//
	// Get again another iterator, iterate across all the domain points, calculating a Laplace stencil
	//
	//
	auto dom2 = g_dist.getDomainIterator();
	
	while (dom2.isNext())
	{
		auto key = dom2.get();

		// Laplace stencil
		g_dist.template get<B>(key)[1] = g_dist.template get<A>(key.move(x,1))[0] + g_dist.template get<A>(key.move(x,-1))[0] +
		                                 g_dist.template get<A>(key.move(y,1))[0] + g_dist.template get<A>(key.move(y,-1))[0] +
										 g_dist.template get<A>(key.move(z,1))[0] + g_dist.template get<A>(key.move(z,-1))[0] -
										 6*g_dist.template get<A>(key)[0];
		                    

		++dom2;
	}

	//
	// ### WIKI 12 ###
	//
	// Finally we want a nice output to visualize the information stored by the distributed grid
	//
	g_dist.write("output");

	//
	// ### WIKI 14 ###
	//
	// Deinitialize the library
	//
	openfpm_finalize();
}
