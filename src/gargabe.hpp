/*
 * gargabe.hpp
 *
 *  Created on: Jan 13, 2015
 *      Author: i-bird
 */

#ifndef GARGABE_HPP_
#define GARGABE_HPP_



	template <unsigned int j, unsigned int i, typename Graph> void optimize(size_t start_p, Graph & graph)
	{
		// We assume that Graph is the rapresentation of a cartesian graph
		// this mean that the direction d is at the child d

		// Create an Hyper-cube

		HyperCube<dim> hyp;

		// Get the number of wavefronts

		size_t n_wf = hyp.getNumberOfElements_R(0);

		// Get the number of intersecting wavefront



		// Get the number of sub-dimensional common wavefront
		// basically are a list of all the subdomain common to two or more

		// Create n_wf wavefront queue

		openfpm::vector<wavefront> v_w;
		v.reserve(n_wf);

		// direction of expansion

		size_t domain_id = 0;
		int exp_dir = 0;
		bool can_expand = true;

		// while is possible to expand

		while (can_expand)
		{
			// for each direction of expansion expand the wavefront

			for (int d = 0 ; d < n_wf ; d++)
			{
				// get the wavefront at direction d

				openfpm::vector<size_t> & wf_d = v_w.get<wavefront::domains>(d);

				// flag to indicate if the wavefront can expand

				bool w_can_expand = true;

				// for each subdomain

				for (size_t sub = 0 ; sub < wf_d.size() ; sub++)
				{
					// check if the adjacent domain in direction d exist
					// and is of the same id

					// get the starting subdomain
					size_t sub_w = wf_d.get<0>(sub);

					// we get the processor id of the neighborhood sub-domain on direction d
					size_t exp_p = graph.getChild(sub_w,d).get<j>();

					// we check if it is the same processor id
					if (exp_p != domain_id)
					{
						w_can_expand = false;
					}
				}

				// if we can expand the wavefront expand it
				if (w_can_expand == true)
				{
					// for each subdomain
					for (size_t sub = 0 ; sub < wf_d.size() ; sub++)
					{
						// update the position of the wavefront
						wf_d.get<0>(sub) = wf_d.get<0>(sub) + gh.stride(d);
					}

					// here we add sub-domains to all the other queues
					// get the face of the hyper-cube

					SubHyperCube<dim,dim-1> sub_hyp = hyp.getSubHyperCube(d);

					std::vector<comb<dim>> q_comb = sub_hyp.getCombinations_R(dim-2);
				}
			}
		}

		// For each point in the Hyper-cube check if we can move the wave front


	}

#ifndef PARALLEL_DECOMPOSITION
//		CreateSubspaces();
#endif

#ifndef USE_METIS_GP

		// Here we do not use METIS
		// Distribute the divided domains

		// Get the number of processing units
		size_t Np = v_cl.getProcessingUnits();

		// Get the ID of this processing unit
		// and push the subspace is taking this
		// processing unit

		for (size_t p_id = v_cl.getProcessUnitID(); p_id < Np ; p_id += Np)
			id_sub.push_back(p_id);
#else


#endif




		/*
		 * CartDecomposition.cpp
		 *
		 *  Created on: Aug 15, 2014
		 *      Author: Pietro Incardona
		 */

		#include "CartDecomposition.hpp"



		/*! \brief The the bulk part of the data set, or the data that does not depend
		 *  from the ghosts layers
		 *
		 * The the bulk part of the data set, or the data that does not depend from the
		 *  ghosts layers
		 *
		 */

		/*template<typename T> T CartDecomposition<T>::getBulk(T data)
		{
			// for each element in data

			for (size_t i = 0; i < data.size() ; i++)
			{
				if (localSpace.isInside())
			}

		}

		template<typename T> T CartDecomposition<T>::getInternal()
		{

		}*/

		/*! \brief Check if is border or bulk
		 *
		 * \param neighboorhood define the neighboorhood of all the points
		 * \return true if border, false if bulk
		 *
		 */

		bool borderOrBulk(neighborhood & nb)
		{
			device::grid<1,size_t> nbr = nb.next();

			// check the neighborhood

			// get neighborhood iterator

			grid_key_dx_iterator<dim> iterator_nbr = nbr.getIterator();

			while (iterator_nbr.hasNext())
			{
				grid_key_dx key_nbr = iterator_nbr.next();

				// check if the neighboorhood is internal

				if(subspace.isBound(data.template get<Point::x>(key_nbr)) == false)
				{
					// it is border

					return true;

					ret.bord.push_back(key);
					break;
				}
			}

			return false;
		}

		/*! \brief This function divide the data set into bulk, border, external and internal part
		 *
		 * \tparam dim dimensionality of the structure storing your data
		 *         (example if they are in 3D grid, has to be 3)
		 * \tparam T type of object we are dividing
		 * \tparam device type of layout selected
		 * \param data 1-dimensional grid of point
		 * \param nb define the neighborhood of all the points
		 * \return a structure with the set of objects divided
		 *
		 */

		template<unsigned int dim, typename T, template<typename> class layout, typename Memory, template<unsigned int, typename> class Domain, template<typename, typename, typename> class data_s>
		dataDiv<T> CartDecomposition<dim,T,layout>::divide(device::grid<1,Point<dim,T>> & data, neighborhood & nb)
		{
			//! allocate the 3 subset

			dataDiv<T> ret;

			ret.bord = new boost::shared_ptr<T>(new T());
			ret.inte = new boost::shared_ptr<T>(new T());
			ret.ext = new boost::shared_ptr<T>(new T());

			//! get grid iterator

			grid_key_dx_iterator<dim> iterator = data.getIterator();

			//! we iterate trough all the set of objects

			while (iterator.hasNext())
			{
				grid_key_dx<dim> key = iterator.next();

				//! Check if the object is inside the subspace

				if (subspace.isBound(data.template get<Point<3,T>::x>(key)))
				{
					//! Check if the neighborhood is inside the subspace

					if (borderOrBulk(nb) == true)
					{
						// It is border

						ret.bord.push_back(key);
					}
					else
					{
						// It is bulk

						ret.bulk.push_back(key);
					}
				}
				else
				{
					//! it is external

					ret.ext.push_back(key);
				}
			}
		}


#endif /* GARGABE_HPP_ */
