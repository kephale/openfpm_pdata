#ifndef DEC_OPTIMIZER_HPP
#define DEC_OPTIMIZER_HPP

#include "Grid/iterators/grid_key_dx_iterator_sub.hpp"

/*! \brief this class represent a wavefront of dimension dim
 *
 * \dim Dimensionality of the wavefront (dimensionality of the space
 *                                       where it live so the wavefront
 *                                       is dim-1)
 *
 */

template <unsigned int dim>
class wavefront
{
public:

	typedef boost::fusion::vector<size_t[dim],size_t[dim]> type;

	type data;

	static const int start = 0;
	static const int stop = 1;
	static const int max_prop = 2;

	/* \brief Get the key to the point 1
	 *
	 * \return the key to the point 1
	 *
	 */

	grid_key_dx<dim> getKP1()
	{
		// grid key to return
		grid_key_dx<dim> ret(boost::fusion::at_c<start>(data));

		return ret;
	}

	/* \brief Get the key to point 2
	 *
	 * \return the key to the point 2
	 *
	 */

	grid_key_dx<dim> getKP2()
	{
		// grid key to return
		grid_key_dx<dim> ret(boost::fusion::at_c<stop>(data));

		return ret;
	}

	/* \brief produce a box from an encapsulated object
	 *
	 * \param encap encapsulated object
	 *
	 */

	template<typename encap> static Box<dim,size_t> getBox(const encap && enc)
	{
		Box<dim,size_t> bx;

		// Create the object from the encapsulation

		getBox(enc,bx);

		return bx;
	}

	/* \brief produce a box from an encapsulated object
	 *
	 * \param encap encapsulated object
	 *
	 */

	template<typename encap> static void getBox(const encap & enc, Box<dim,size_t> & bx)
	{
		// Create the object from the encapsulation

		for (int i = 0 ; i < dim ; i++)
		{
			bx.setLow(i,enc.template get<wavefront::start>()[i]);
			bx.setHigh(i,enc.template get<wavefront::stop>()[i]);
		}
	}
};

/*! \brief This class take a graph representing the space decomposition and produce a
 *         simplified version
 *
 * Given a Graph_CSR and a seed point produce an alternative decomposition in boxes with
 * less sub-domain. In the following we referee with sub-domain the boxes produced by this
 * algorithm and sub-sub-domain the sub-domain before reduction
 *
 */

template <unsigned int dim, typename Graph>
class dec_optimizer
{
	// create a grid header for helping

	grid_sm<dim,void> gh;

private:

	/*! \brief Expand one wavefront
	 *
	 * \param v_w wavefronts
	 * \param w_comb wavefront expansion combinations
	 * \param d direction of expansion
	 * \param bc boundary condition
	 *
	 */
	void expand_one_wf(openfpm::vector<wavefront<dim>> & v_w, std::vector<comb<dim>> & w_comb , size_t d)
	{
		for (size_t j = 0 ; j < dim ; j++)
		{
			v_w.template get<wavefront<dim>::stop>(d)[j] = v_w.template get<wavefront<dim>::stop>(d)[j] + w_comb[d].c[j];
			v_w.template get<wavefront<dim>::start>(d)[j] = v_w.template get<wavefront<dim>::start>(d)[j] + w_comb[d].c[j];
		}
	}


	/*! \brief Adjust the other wavefronts
	 *
	 * \param d direction
	 *
	 */
	void adjust_others_wf(openfpm::vector<wavefront<dim>> & v_w,  HyperCube<dim> & hyp, std::vector<comb<dim>> & w_comb, size_t d)
	{
		// expand the intersection of the wavefronts

		std::vector<comb<dim>> q_comb = SubHyperCube<dim,dim-1>::getCombinations_R(w_comb[d],dim-2);

		// Eliminate the w_comb[d] direction

		for (size_t k = 0 ; k < q_comb.size() ; k++)
		{
			for (size_t j = 0 ; j < dim ; j++)
			{
				if (w_comb[d].c[j] != 0)
				{
					q_comb[k].c[j] = 0;
				}
			}
		}

		// for all the combinations
		for (size_t j = 0 ; j < q_comb.size() ; j++)
		{
			size_t id = hyp.LinId(q_comb[j]);

			// get the combination of the direction d

			bool is_pos = hyp.isPositive(d);

			// is positive, modify the stop point or the starting point

			for (size_t s = 0 ; s < dim ; s++)
			{
				if (is_pos == true)
				{v_w.template get<wavefront<dim>::stop>(id)[s] = v_w.template get<wavefront<dim>::stop>(id)[s] + w_comb[d].c[s];}
				else
				{v_w.template get<wavefront<dim>::start>(id)[s] = v_w.template get<wavefront<dim>::start>(id)[s] + w_comb[d].c[s];}
			}
		}
	}

	/* \brief Fill the wavefront position
	 *
	 * \tparam prp property to set
	 *
	 * \param graph we are processing
	 * \param Box to fill
	 * \param id value to fill with
	 *
	 */

	template<unsigned int prp> void write_wavefront(Graph & graph,openfpm::vector<wavefront<dim>> & v_w)
	{
		// fill the wall domain with 0

		fill_domain<prp>(graph,gh.getBox(),0);

		// fill the wavefront

		for (int i = 0 ; i < v_w.size() ; i++)
		{
			Box<dim,size_t> box = wavefront<dim>::getBox(v_w.get(i));

			fill_domain<prp>(graph,box,1);
		}
	}

	/* \brief Fill the domain
	 *
	 * \tparam p_sub property to set with the sub-domain id
	 *
	 * \param graph we are processing
	 * \param Box to fill
	 * \param ids value to fill with
	 *
	 */

	template<unsigned int p_sub> void fill_domain(Graph & graph,const Box<dim,size_t> & box, long int ids)
	{
		// Create a subgrid iterator
		grid_key_dx_iterator_sub<dim,do_not_print_warning_on_adjustment<dim>> g_sub(gh,box.getKP1(),box.getKP2());

		// iterate through all grid points

		while (g_sub.isNext())
		{
			// get the actual key
			const grid_key_dx<dim> & gk = g_sub.get();

			// get the vertex and set the sub id

			graph.vertex(gh.LinId(gk)).template get<p_sub>() = ids;

			// next subdomain
			++g_sub;
		}
	}

	/* \brief Add the boundary domain of id p_id to the queue
	 *
	 * \tparam i-property where is stored the decomposition
	 *
	 * \param domains vector with domains to process
	 * \param graph we are processing
	 * \param w_comb hyper-cube combinations
	 * \param p_id processor id
	 * \param box_nn_processor list of neighborhood processors for the box
	 * \param bc Boundary conditions
	 *
	 */
	template<unsigned int p_sub, unsigned int p_id> void add_to_queue(openfpm::vector<size_t> & domains, openfpm::vector<wavefront<dim>> & v_w, Graph & graph,  std::vector<comb<dim>> & w_comb, long int pr_id, openfpm::vector< openfpm::vector<size_t> > & box_nn_processor, const size_t(& bc)[dim])
	{
		// create a new queue
		openfpm::vector<size_t> domains_new;

		// push in the new queue, the old domains of the queue that are not assigned element

		for (size_t j = 0 ; j < domains.size() ; j++)
		{
			long int gs = graph.vertex(domains.get(j)).template get<p_sub>();
			if (gs < 0)
			{
				// not assigned push it

				domains_new.add(domains.get(j));
			}
		}

		// Create an Hyper-cube
		HyperCube<dim> hyp;

		for (size_t d = 0 ; d < v_w.size() ; d++)
		{
			expand_one_wf(v_w,w_comb,d);
			adjust_others_wf(v_w,hyp,w_comb,d);
		}

		// for each expanded wavefront create a sub-grid iterator and add the sub-domain

		for (size_t d = 0 ; d < v_w.size() ; d++)
		{
			// Create a sub-grid iterator
			grid_key_dx_iterator_sub_bc<dim,do_not_print_warning_on_adjustment<dim>> g_sub(gh,v_w.template get<wavefront<dim>::start>(d),v_w.template get<wavefront<dim>::stop>(d),bc);

			// iterate through all grid points

			while (g_sub.isNext())
			{
				// get the actual key
				const grid_key_dx<dim> & gk = g_sub.get();

				// get the vertex and if does not have a sub-id and is assigned ...
				long int pid = graph.vertex(gh.LinId(gk)).template get<p_sub>();

				// Get the processor id of the sub-sub-domain
				long int pp_id = graph.vertex(gh.LinId(gk)).template get<p_id>();

				if (pp_id != pr_id)
				{
					// this box is contiguous to other processors

					box_nn_processor.get(box_nn_processor.size()-1).add(pp_id);
				}

				// if the sub-sub-domain is not assigned
				if (pid < 0)
				{
					// ... and we are not processing the full graph
					if (pr_id != -1)
					{
						// ... and the processor id of the sub-sub-domain match the part we are processing, add to the queue

						if ( pr_id == pp_id)
							domains_new.add(gh.LinId(gk));
					}
					else
						domains_new.add(gh.LinId(gk));
				}

				++g_sub;
			}
		}

		// sort and make it unique the numbers in processor list
	    std::sort(box_nn_processor.get(box_nn_processor.size()-1).begin(), box_nn_processor.get(box_nn_processor.size()-1).end());
	    auto last = std::unique(box_nn_processor.get(box_nn_processor.size()-1).begin(), box_nn_processor.get(box_nn_processor.size()-1).end());
	    box_nn_processor.get(box_nn_processor.size()-1).erase(last, box_nn_processor.get(box_nn_processor.size()-1).end());

		// copy the new queue to the old one (it not copied, C++11 move semantic)
		domains.swap(domains_new);
	}

	/* \brief Find the biggest hyper-cube
	 *
	 * starting from one initial sub-domain find the biggest hyper-cube
	 * output the box, and fill a list of neighborhood processor
	 *
	 * \tparam p_sub id of the property storing the sub-decomposition
	 * \tparam p_id id of the property containing the decomposition
	 *
	 * \param start_p initial domain
	 * \param graph representing the grid of sub-sub-domain
	 * \param box produced box
	 * \param v_w Wavefronts
	 * \param w_comb wavefronts directions (0,0,1) (0,0,-1) (0,1,0) (0,-1,0) ...
	 *
	 */
	template <unsigned int p_sub, unsigned int p_id> void expand_from_point(size_t start_p, Graph & graph, Box<dim,size_t> & box, openfpm::vector<wavefront<dim>> & v_w , std::vector<comb<dim>> & w_comb)
	{
		// We assume that Graph is the rapresentation of a cartesian graph
		// this mean that the direction d is at the child d

		// Get the number of wavefronts
		size_t n_wf = w_comb.size();

		// Create an Hyper-cube
		HyperCube<dim> hyp;

		// direction of expansion

		size_t domain_id = graph.vertex(start_p).template get<p_id>();
		bool can_expand = true;

		// while is possible to expand

		while (can_expand)
		{
			// reset can expand
			can_expand = false;

			// for each direction of expansion expand the wavefront

			for (size_t d = 0 ; d < n_wf ; d++)
			{
				// number of processed sub-domain
				size_t n_proc_sub = 0;

				// flag to indicate if the wavefront can expand
				bool w_can_expand = true;

				// Create an iterator of the expanded wavefront
				grid_key_dx<dim> start = grid_key_dx<dim>(v_w.template get<wavefront<dim>::start>(d)) + w_comb[d];
				grid_key_dx<dim> stop = grid_key_dx<dim>(v_w.template get<wavefront<dim>::stop>(d)) + w_comb[d];
				grid_key_dx_iterator_sub<dim,do_not_print_warning_on_adjustment<dim>> it(gh,start,stop);

				// for each sub-domain in the expanded wavefront
				while (it.isNext())
				{
					// get the wavefront sub-domain id
					size_t sub_w_e = gh.LinId(it.get());

					// we get the processor id of the neighborhood sub-domain on direction d
					// (expanded wavefront)
					size_t exp_p = graph.vertex(sub_w_e).template get<p_id>();

					// Check if already assigned
					long int ass = graph.vertex(sub_w_e).template get<p_sub>();

					// we check if it is the same processor id and is not assigned
					w_can_expand &= ((exp_p == domain_id) & (ass < 0));

					// next domain
					++it;

					// increase the number of processed sub-domain
					n_proc_sub++;
				}

				// if we did not processed sub-domain, we cannot expand
				w_can_expand &= (n_proc_sub != 0);

				// if you can expand one wavefront we did not reach the end
				can_expand |= w_can_expand;

				// if we can expand the wavefront expand it
				if (w_can_expand == true)
				{
					// expand the wavefront
					for (size_t j = 0 ; j < dim ; j++)
					{
						v_w.template get<wavefront<dim>::stop>(d)[j] = v_w.template get<wavefront<dim>::stop>(d)[j] + w_comb[d].c[j];
						v_w.template get<wavefront<dim>::start>(d)[j] = v_w.template get<wavefront<dim>::start>(d)[j] + w_comb[d].c[j];
					}

					// expand the intersection of the wavefronts

					std::vector<comb<dim>> q_comb = SubHyperCube<dim,dim-1>::getCombinations_R(w_comb[d],dim-2);

					// Eliminate the w_comb[d] direction

					for (size_t k = 0 ; k < q_comb.size() ; k++)
					{
						for (size_t j = 0 ; j < dim ; j++)
						{
							if (w_comb[d].c[j] != 0)
							{
								q_comb[k].c[j] = 0;
							}
						}
					}

					// for all the combinations
					for (size_t j = 0 ; j < q_comb.size() ; j++)
					{
						size_t id = hyp.LinId(q_comb[j]);

						// get the combination of the direction d

						bool is_pos = hyp.isPositive(d);

						// is positive, modify the stop point or the starting point

						for (size_t s = 0 ; s < dim ; s++)
						{
							if (is_pos == true)
							{v_w.template get<wavefront<dim>::stop>(id)[s] = v_w.template get<wavefront<dim>::stop>(id)[s] + w_comb[d].c[s];}
							else
							{v_w.template get<wavefront<dim>::start>(id)[s] = v_w.template get<wavefront<dim>::start>(id)[s] + w_comb[d].c[s];}
						}
					}
				}
			}
		}

		// get back the hyper-cube produced

		for (size_t i = 0 ; i < dim ; i++)
		{
			// get the index of the wavefront direction
			size_t p_f = hyp.positiveFace(i);
			size_t n_f = hyp.negativeFace(i);

			// set the box
			box.setHigh(i,v_w.template get<wavefront<dim>::stop>(p_f)[i]);
			box.setLow(i,v_w.template get<wavefront<dim>::start>(n_f)[i]);
		}
	}

	/*! \brief Initialize the wavefronts
	 *
	 * \param starting point of the wavefront set
	 * \param v_w Wavefront to initialize
	 *
	 */
	void InitializeWavefront(grid_key_dx<dim> & start_p, openfpm::vector<wavefront<dim>> & v_w)
	{
		// Wavefront to initialize

		for (size_t i = 0 ; i < v_w.size() ; i++)
		{
			for (size_t j = 0 ; j < dim ; j++)
			{
				v_w.template get<wavefront<dim>::start>(i)[j] = start_p.get(j);
				v_w.template get<wavefront<dim>::stop>(i)[j] = start_p.get(j);
			}
		}
	}

	/*! \brief Get the first seed
	 *
	 * search in the graph for one sub-domain labelled with processor id
	 * to use as seed
	 *
	 * \tparam p_id property in the graph storing the sub-domain id
	 *
	 * \param Graph graph
	 * \param id processor id
	 *
	 */

	template<unsigned int p_id, unsigned int p_sub> grid_key_dx<dim> search_seed(Graph & graph, long int id)
	{
		// if no processor is selected return the first point
		if (id < -1)
		{
			grid_key_dx<dim> key;
			key.zero();

			return key;
		}

		// Create a grid iterator
		grid_key_dx_iterator<dim> g_sub(gh);

		// iterate through all grid points

		while (g_sub.isNext())
		{
			// get the actual key
			const grid_key_dx<dim> & gk = g_sub.get();

			// if the subdomain has the id we are searching stop
			if ((long int)graph.vertex(gh.LinId(gk)).template get<p_id>() == id && graph.vertex(gh.LinId(gk)).template get<p_sub>() == -1)
			{
				return gk;
			}

			++g_sub;
		}

		// If not found return an invalid key
		grid_key_dx<dim> key;
		key.invalid();

		return key;
	}


	/*! \brief optimize the graph
	 *
	 * Starting from a domain (hyper-cubic), it create wavefront at the boundary and expand
	 * the boundary until the wavefronts cannot expand any more.
	 * To the domains inside the hyper-cube one sub-id is assigned. This procedure continue until
	 * all the domain of one p_id has a sub-id
	 *
	 * \tparam j property containing the decomposition
	 * \tparam i property to fill with the sub-decomposition
	 *
	 * \param start_p seed point
	 * \param graph we are processing
	 * \param p_id Processor id (if p_id == -1 the optimization is done for all the processors)
	 * \param list of sub-domain boxes produced by the algorithm
	 * \param box_nn_processor for each box it list all the neighborhood processor
	 * \param bc Boundary condition
	 * \param init_sub_id when true p_sub property is initial set to -1 [default true]
	 * \param sub_id starting sub_id enumeration [default 0]
	 *
	 * \return last assigned sub-id
	 *
	 */
	template <unsigned int p_sub, unsigned int p_id> size_t optimize(grid_key_dx<dim> & start_p, Graph & graph, long int pr_id, openfpm::vector<Box<dim,size_t>> & lb, openfpm::vector< openfpm::vector<size_t> > & box_nn_processor, const size_t (& bc)[dim], bool init_sub_id = true, size_t sub_id = 0)
	{
		// queue
		openfpm::vector<size_t> v_q;

		// Create an hyper-cube
		HyperCube<dim> hyp;

		// Get the wavefront combinations
		std::vector<comb<dim>> w_comb = hyp.getCombinations_R(dim-1);

		// wavefronts
		openfpm::vector<wavefront<dim>> v_w(w_comb.size());

		// fill the sub decomposition with negative number

		if (init_sub_id == true)
			fill_domain<p_sub>(graph,gh.getBox(),-1);

		// push the first domain
		v_q.add(gh.LinId(start_p));

		while (v_q.size() != 0)
		{
			// Box
			Box<dim,size_t> box;

			// Get the grid_key position from the linearized id
			start_p = gh.InvLinId(v_q.get(0));

			// Initialize the wavefronts from the domain start_p
			InitializeWavefront(start_p,v_w);

			// Create a list for the box
			box_nn_processor.add();

			// Create the biggest box containing the domain
			expand_from_point<p_sub,p_id>(v_q.get(0),graph,box,v_w,w_comb);

			// Add the created box to the list of boxes
			lb.add(box);

			// fill the domain
			fill_domain<p_sub>(graph,box,sub_id);

			// add the surrounding sub-domain to the queue
			add_to_queue<p_sub,p_id>(v_q,v_w,graph,w_comb,pr_id,box_nn_processor,bc);

			// increment the sub_id
			sub_id++;
		}

		return sub_id;
	}

public:

	/*! \brief Constructor
	 *
	 * \param g Graph to simplify
	 * \param sz size of the grid on each dimension
	 *
	 */

	dec_optimizer(Graph & g, const size_t (& sz)[dim])
	:gh(sz)
	{
		// The graph g is suppose to represent a cartesian grid
		// No check is performed on g
	}

	/*! \brief optimize the graph
	 *
	 * Starting from a sub-sub-domain, it create wavefronts at the boundary and expand
	 * the boundary until the wavefronts cannot expand any more, creating a sub-domain covering more sub-sub-domain.
	 * This procedure continue until all the domain is covered by a sub-domains
	 *
	 * \tparam j property containing the processor decomposition
	 * \tparam i property to fill with the sub-domain-decomposition id
	 *
	 * \param start_p seed point
	 * \param graph we are processing
	 *
	 */
	template <unsigned int p_sub, unsigned int p_id> void optimize(grid_key_dx<dim> & start_p, Graph & graph, const size_t (& bc)[dim])
	{
		// temporal vector
		openfpm::vector<Box<dim,size_t>> tmp;

		// temporal vector
		openfpm::vector< openfpm::vector<size_t> > box_nn_processor;

		// optimize
		optimize<p_sub,p_id>(start_p,graph,-1,tmp, box_nn_processor,bc);
	}

	/*! \brief optimize the graph
	 *
	 * Starting from a sub-sub-domain, it create wavefronts at the boundary and expand
	 * the boundary until the wavefronts cannot expand any more, creating a sub-domain covering more sub-sub-domain.
	 * This procedure continue until all the sub-domain of the processor p_id are covered by a sub-domains
	 *
	 * \tparam j property containing the decomposition
	 * \tparam i property to fill with the sub-domain-decomposition id
	 *
	 * \param graph we are processing
	 * \param p_id Processor id (if p_id == -1 the optimization is done for all the processors)
	 * \param list of sub-domain boxes
	 *
	 */
	template <unsigned int p_sub, unsigned int p_id> void optimize(Graph & graph, long int pr_id, openfpm::vector<Box<dim,size_t>> & lb, openfpm::vector< openfpm::vector<size_t> > & box_nn_processor, const size_t (& bc)[dim])
	{
		grid_key_dx<dim> key_seed;
		key_seed.zero();

		// if processor is -1 call optimize with -1 to do on all processors and exit
		if (pr_id == -1)
		{
			optimize<p_sub,p_id>(key_seed,graph,pr_id,lb,box_nn_processor,bc);
			return;
		}

		size_t sub_id = 0;

		// fill the sub decomposition with negative number
		fill_domain<p_sub>(graph,gh.getBox(),-1);

		key_seed = search_seed<p_id,p_sub>(graph,pr_id);

		while (key_seed.isValid())
		{
			// optimize
			sub_id = optimize<p_sub,p_id>(key_seed,graph,pr_id,lb,box_nn_processor,bc,false,sub_id);

			// new seed
			key_seed = search_seed<p_id,p_sub>(graph,pr_id);
		}
	}
};

#endif
