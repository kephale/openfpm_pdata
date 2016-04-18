/*
 * ParMetisDistribution.hpp
 *
 *  Created on: Nov 19, 2015
 *      Author: Antonio Leo
 */


#ifndef SRC_DECOMPOSITION_DISTRIBUTION_PARMETISDISTRIBUTION_HPP_
#define SRC_DECOMPOSITION_DISTRIBUTION_PARMETISDISTRIBUTION_HPP_


#include "SubdomainGraphNodes.hpp"
#include "parmetis_util.hpp"
#include "Graph/ids.hpp"

#define PARMETIS_DISTRIBUTION_ERROR 100002

/*! \brief Class that distribute sub-sub-domains across processors using ParMetis Library
 *
 * Given a graph and setting Computational cost, Communication cost (on the edge) and
 * Migration cost or total Communication costs, it produce the optimal balanced distribution
 *
 * In addition to Metis it provide the functionality to refine the previously computed
 * decomposition
 *
 * ### Initialize a Cartesian graph and decompose
 * \snippet Distribution_unit_tests.hpp Initialize a ParMetis Cartesian graph and decompose
 *
 * ### Refine the decomposition with parmetis
 * \snippet Distribution_unit_tests.hpp refine the decomposition with parmetis
 *
 */
template<unsigned int dim, typename T>
class ParMetisDistribution
{
	//! Vcluster
	Vcluster & v_cl;

	//! Structure that store the cartesian grid information
	grid_sm<dim, void> gr;

	//! rectangular domain to decompose
	Box<dim, T> domain;

	//! Global sub-sub-domain graph
	Graph_CSR<nm_v, nm_e> gp;

	//! Convert the graph to parmetis format
	Parmetis<Graph_CSR<nm_v, nm_e>> parmetis_graph;

	//! Init vtxdist needed for Parmetis
	//
	// vtxdist is a common array across processor, it indicate how
	// vertex are distributed across processors
	//
	// Example we have 3 processors
	//
	// processor 0 has 3 vertices
	// processor 1 has 5 vertices
	// processor 2 has 4 vertices
	//
	// vtxdist contain, 0,3,8,12
	//
	// vtx dist is the unique global-id of the vertices
	//
	openfpm::vector<rid> vtxdist;

	//! partitions
	openfpm::vector<openfpm::vector<idx_t>> partitions;

	//! Init data structure to keep trace of new vertices distribution in processors (needed to update main graph)
	openfpm::vector<openfpm::vector<gid>> v_per_proc;

	//! Hashmap to access to the global position given the re-mapped one (needed for access the map)
	std::unordered_map<rid, gid> m2g;

	//! Flag to check if weights are used on vertices
	bool verticesGotWeights = false;

	/*! \brief Update main graph ad subgraph with the received data of the partitions from the other processors
	 *
	 */
	void updateGraphs()
	{
		size_t Np = v_cl.getProcessingUnits();

		// Init n_vtxdist to gather informations about the new decomposition
		openfpm::vector<rid> n_vtxdist(Np + 1);
		for (size_t i = 0; i <= Np; i++)
			n_vtxdist.get(i).id = 0;

		// Update the main graph with received data from processor i
		for (size_t i = 0; i < Np; i++)
		{
			size_t ndata = partitions.get(i).size();
			size_t k = 0;

			// Update the main graph with the received informations
			for (rid l = vtxdist.get(i); k < ndata && l < vtxdist.get(i + 1); k++, ++l)
			{
				// Create new n_vtxdist (just count processors vertices)
				++n_vtxdist.get(partitions.get(i).get(k) + 1);

				// Update proc id in the vertex (using the old map)
				vertexByMapId(l).template get<nm_v::proc_id>() = partitions.get(i).get(k);

				// Add vertex to temporary structure of distribution (needed to update main graph)
				v_per_proc.get(partitions.get(i).get(k)).add(getVertexGlobalId(l));
			}
		}

		// Create new n_vtxdist (accumulate the counters)
		for (size_t i = 2; i <= Np; i++)
			n_vtxdist.get(i) += n_vtxdist.get(i - 1);

		// Copy the new decomposition in the main vtxdist
		for (size_t i = 0; i <= Np; i++)
			vtxdist.get(i) = n_vtxdist.get(i);

		// Renumber the main graph and re-create the map
		for (size_t p = 0; p < (size_t)Np; p++)
		{
			size_t i = 0;
			for (rid j = vtxdist.get(p); j < vtxdist.get(p + 1); ++j, i++)
			{
				setMapId(j, v_per_proc.get(p).get(i));
				gp.vertex(v_per_proc.get(p).get(i).id).template get<nm_v::id>() = j.id;
			}
		}
	}

	/*! \brief operator to access the vertex by mapped position
	 *
	 * operator to access the vertex
	 *
	 * \param id re-mapped id of the vertex to access
	 *
	 */
	inline auto vertexByMapId(rid id) -> decltype( gp.vertex(m2g.find(id)->second.id) )
	{
		return gp.vertex(m2g.find(id)->second.id);
	}

	/*! \brief operator to remap vertex to a new position
	 *
	 * \param n re-mapped position
	 * \param g global position
	 *
	 */
	inline void setMapId(rid n, gid g)
	{
		m2g[n] = g;
	}

	/*! \brief Get the global id of the vertex given the re-mapped one
	 *
	 * \param remapped id
	 * \return global id
	 *
	 */
	gid getVertexGlobalId(rid n)
	{
		return m2g.find(n)->second;
	}

	/*! \brief operator to init ids vector
	 *
	 * operator to init ids vector
	 *
	 */
	void initLocalToGlobalMap()
	{
		gid g;
		rid i;
		i.id = 0;

		m2g.clear();
		for ( ; (size_t)i.id < gp.getNVertex(); ++i)
		{
			g.id = i.id;

			m2g.insert( { i, g });
		}
	}

	/*! \brief Callback of the sendrecv to set the size of the array received
	 *
	 * \param msg_i Index of the message
	 * \param total_msg Total numeber of messages
	 * \param total_p Total number of processors to comunicate with
	 * \param i Processor id
	 * \param ri Request id
	 * \param ptr Void pointer parameter for additional data to pass to the call-back
	 */
	static void * message_receive(size_t msg_i, size_t total_msg, size_t total_p, size_t i, size_t ri, void * ptr)
	{
		openfpm::vector < openfpm::vector < idx_t >> *v = static_cast<openfpm::vector<openfpm::vector<idx_t>> *>(ptr);

		v->get(i).resize(msg_i / sizeof(idx_t));

		return &(v->get(i).get(0));
	}

public:

	/*! Constructor for the ParMetis class
	 *
	 * \param v_cl Vcluster to use as communication object in this class
	 */
	ParMetisDistribution(Vcluster & v_cl)
	:v_cl(v_cl), parmetis_graph(v_cl, v_cl.getProcessingUnits()), vtxdist(v_cl.getProcessingUnits() + 1), partitions(v_cl.getProcessingUnits()), v_per_proc(v_cl.getProcessingUnits())
	{
	}

	/*! Copy constructor
	 *
	 * \param pm Distribution to copy
	 *
	 */
	ParMetisDistribution(const ParMetisDistribution<dim,T> & pm)
	:v_cl(pm.v_cl),parmetis_graph(v_cl, v_cl.getProcessingUnits())
	{
		this->operator=(pm);
	}

	/*! Copy constructor
	 *
	 * \param pm Distribution to copy
	 *
	 */
	ParMetisDistribution(ParMetisDistribution<dim,T> && pm)
	{
		this->operator=(pm);
	}

	/*! \brief Create the Cartesian graph
	 *
	 * \param grid info
	 * \param dom domain
	 */
	void createCartGraph(grid_sm<dim, void> & grid, Box<dim, T> dom)
	{
		size_t bc[dim];

		for (size_t i = 0 ; i < dim ; i++)
			bc[i] = NON_PERIODIC;

		// Set grid and domain
		gr = grid;
		domain = dom;

		// Create a cartesian grid graph
		CartesianGraphFactory<dim, Graph_CSR<nm_v, nm_e>> g_factory_part;
		gp = g_factory_part.template construct<NO_EDGE, nm_v::id, T, dim - 1, 0>(gr.getSize(), domain, bc);
		initLocalToGlobalMap();

		//! Get the number of processing units
		size_t Np = v_cl.getProcessingUnits();

		//! Division of vertices in Np graphs
		//! Put (div+1) vertices in mod graphs
		//! Put div vertices in the rest of the graphs
		size_t mod_v = gr.size() % Np;
		size_t div_v = gr.size() / Np;

		for (size_t i = 0; i <= Np; i++)
		{
			if (i < mod_v)
				vtxdist.get(i).id = (div_v + 1) * i;
			else
				vtxdist.get(i).id = (div_v) * i + mod_v;
		}

		// Init to 0.0 axis z (to fix in graphFactory)
		if (dim < 3)
		{
			for (size_t i = 0; i < gp.getNVertex(); i++)
			{
				gp.vertex(i).template get<nm_v::x>()[2] = 0.0;
			}
		}
		for (size_t i = 0; i < gp.getNVertex(); i++)
		{
			gp.vertex(i).template get<nm_v::global_id>() = i;
		}

	}

	/*! \brief Get the graph of the distribution
	 *
	 */
	Graph_CSR<nm_v, nm_e> & getGraph()
	{
		return gp;
	}

	/*! \brief Create the decomposition
	 *
	 */
	void decompose()
	{

		//! Get the processor id
		size_t p_id = v_cl.getProcessUnitID();

		//! Get the number of processing units
		size_t Np = v_cl.getProcessingUnits();

		// Number of local vertex
		size_t nl_vertex = vtxdist.get(p_id+1).id - vtxdist.get(p_id).id;

		parmetis_graph.initSubGraph(gp, vtxdist, m2g, verticesGotWeights);

		//! Decompose
		parmetis_graph.decompose<nm_v::proc_id>(vtxdist);

		//! Get result partition for this processors
		idx_t *partition = parmetis_graph.getPartition();

		//! Prepare vector of arrays to contain all partitions
		partitions.get(p_id).resize(nl_vertex);
		std::copy(partition, partition + nl_vertex, &partitions.get(p_id).get(0));

		// Communicate the local distribution to the other processors
		// to reconstruct individually the global graph
		openfpm::vector<size_t> prc;
		openfpm::vector<size_t> sz;
		openfpm::vector<void *> ptr;

		for (size_t i = 0; i < Np; i++)
		{
			if (i != v_cl.getProcessUnitID())
			{
				prc.add(i);
				sz.add(nl_vertex * sizeof(idx_t));
				ptr.add(partitions.get(p_id).getPointer());
			}
		}

		v_cl.sendrecvMultipleMessagesNBX(prc.size(), &sz.get(0), &prc.get(0), &ptr.get(0), message_receive, &partitions,
		NONE);

		// Update graphs with the received data
		updateGraphs();
	}

	/*! \brief Refine current decomposition
	 *
	 * It makes a refinement of the current decomposition using Parmetis function RefineKWay
	 * After that it also does the re-mapping of the graph
	 *
	 */
	void refine()
	{
		size_t Np = v_cl.getProcessingUnits();
		size_t p_id = v_cl.getProcessUnitID();

		// Number of local vertex
		rid nl_vertex = vtxdist.get(p_id+1) - vtxdist.get(p_id);

		// Reset parmetis graph and reconstruct it
		parmetis_graph.reset(gp, vtxdist, m2g, verticesGotWeights);

		// Refine
		parmetis_graph.refine<nm_v::proc_id>(vtxdist);

		// Get result partition for this processor
		idx_t * partition = parmetis_graph.getPartition();

		partitions.get(p_id).resize(nl_vertex.id);
		std::copy(partition, partition + nl_vertex.id, &partitions.get(p_id).get(0));

		// Reset data structure to keep trace of new vertices distribution in processors (needed to update main graph)
		for (size_t i = 0; i < Np; ++i)
		{
			v_per_proc.get(i).clear();
		}

		openfpm::vector<size_t> prc;
		openfpm::vector<size_t> sz;
		openfpm::vector<void *> ptr;

		for (size_t i = 0; i < Np; i++)
		{
			if (i != v_cl.getProcessUnitID())
			{
				partitions.get(i).clear();
				prc.add(i);
				sz.add(nl_vertex.id * sizeof(idx_t));
				ptr.add(partitions.get(p_id).getPointer());
			}
		}

		// Exchange informations through processors
		v_cl.sendrecvMultipleMessagesNBX(prc.size(), &sz.get(0), &prc.get(0), &ptr.get(0), message_receive, &partitions,
		NONE);

		// Update graphs with the new distributions
		updateGraphs();
	}

	/*! \brief Compute the unbalance of the processor compared to the optimal balance
	 *
	 * \return the unbalance from the optimal one 0.01 mean 1%
	 */
	float getUnbalance()
	{
		long t_cost = 0;

		long min, max, sum;
		float unbalance;

		t_cost = getProcessorLoad();

		min = t_cost;
		max = t_cost;
		sum = t_cost;

		v_cl.min(min);
		v_cl.max(max);
		v_cl.sum(sum);
		v_cl.execute();

		unbalance = ((float) (max - min)) / (float) (sum / v_cl.getProcessingUnits());

		return unbalance * 100;
	}

	/*! \brief function that return the position of the vertex in the space
	 *
	 * \param id vertex id
	 * \param pos vector that will contain x, y, z
	 *
	 */
	void getSubSubDomainPosition(size_t id, T (&pos)[dim])
	{
		if (id >= gp.getNVertex())
			std::cerr << "Such vertex doesn't exist (id = " << id << ", " << "total size = " << gp.getNVertex() << ")\n";

		// Copy the geometrical informations inside the pos vector
		pos[0] = gp.vertex(id).template get<nm_v::x>()[0];
		pos[1] = gp.vertex(id).template get<nm_v::x>()[1];
		if (dim == 3)
			pos[2] = gp.vertex(id).template get<nm_v::x>()[2];
	}

	/*! \brief Function that set the weight of the vertex
	 *
	 * \param id vertex id
	 * \param weight to give to the vertex
	 *
	 */
	inline void setComputationCost(size_t id, size_t weight)
	{
		if (!verticesGotWeights)
			verticesGotWeights = true;

		if (id >= gp.getNVertex())
			std::cerr << "Such vertex doesn't exist (id = " << id << ", " << "total size = " << gp.getNVertex() << ")\n";

		// Update vertex in main graph
		gp.vertex(id).template get<nm_v::computation>() = weight;
	}

	/*! \brief Checks if weights are used on the vertices
	 *
	 * \return true if weights are used in the decomposition
	 */
	bool weightsAreUsed()
	{
		return verticesGotWeights;
	}

	/*! \brief function that get the weight of the vertex
	 *
	 * \param id vertex id
	 *
	 */
	size_t getSubSubDomainComputationCost(size_t id)
	{
		if (id >= gp.getNVertex())
			std::cerr << "Such vertex doesn't exist (id = " << id << ", " << "total size = " << gp.getNVertex() << ")\n";

		return gp.vertex(id).template get<nm_v::computation>();
	}

	/*! \brief Compute the processor load counting the total weights of its vertices
	 *
	 * \return the computational load of the processor graph
	 */
	size_t getProcessorLoad()
	{
		size_t load = 0;

		// Processor id
		size_t p_id = v_cl.getProcessUnitID();

		for (rid i = vtxdist.get(p_id); i < vtxdist.get(p_id+1) ; ++i)
			load += gp.vertex(m2g.find(i)->second.id).template get<nm_v::computation>();

		return load;
	}

	/*! \brief Set migration cost of the vertex id
	 *
	 * \param id of the vertex to update
	 * \param migration cost of the migration
	 */
	void setMigrationCost(size_t id, size_t migration)
	{
		if (id >= gp.getNVertex())
			std::cerr << "Such vertex doesn't exist (id = " << id << ", " << "total size = " << gp.getNVertex() << ")\n";

		gp.vertex(id).template get<nm_v::migration>() = migration;
	}

	/*! \brief Set communication cost of the edge id
	 *
	 * \param v_id Id of the source vertex of the edge
	 * \param e i child of the vertex
	 * \param communication Communication value
	 */
	void setCommunicationCost(size_t v_id, size_t e, size_t communication)
	{
		size_t e_id = v_id + e;

		if (e_id >= gp.getNEdge())
			std::cerr << "Such edge doesn't exist (id = " << e_id << ", " << "total size = " << gp.getNEdge() << ")\n";

		gp.getChildEdge(v_id, e).template get<nm_e::communication>() = communication;
	}

	/*! \brief Returns total number of sub-sub-domains in the distribution graph
	 *
	 */
	size_t getNSubSubDomains()
	{
		return gp.getNVertex();
	}

	/*! \brief Returns total number of neighbors of the sub-sub-domain id
	 *
	 * \param i id of the sub-sub-domain
	 */
	size_t getNSubSubDomainNeighbors(size_t id)
	{
		if (id >= gp.getNVertex())
			std::cerr << "Such vertex doesn't exist (id = " << id << ", " << "total size = " << gp.getNVertex() << ")\n";

		return gp.getNChilds(id);
	}

	/*! \brief Print the current distribution and save it to VTK file
	 *
	 * \param file filename
	 *
	 */
	void write(const std::string & file)
	{
		VTKWriter<Graph_CSR<nm_v, nm_e>, VTK_GRAPH> gv2(gp);
		gv2.write(std::to_string(v_cl.getProcessUnitID()) + "_" + file + ".vtk");
	}

	/*! \brief Operator definition for DistParMetisDistribution class
	 *
	 * \param dist
	 */
	const ParMetisDistribution<dim,T> & operator=(const ParMetisDistribution<dim,T> & dist)
	{
		gr = dist.gr;
		domain = dist.domain;
		gp = dist.gp;
		vtxdist = dist.vtxdist;
		partitions = dist.partitions;
		v_per_proc = dist.v_per_proc;
		verticesGotWeights = dist.verticesGotWeights;

		return *this;
	}

	/*! \brief Operator definition for DistParMetisDistribution class
	 *
	 * \param dist
	 */
	const ParMetisDistribution<dim,T> & operator=(ParMetisDistribution<dim,T> && dist)
	{
		v_cl = dist.v_cl;
		gr = dist.gr;
		domain = dist.domain;
		gp.swap(dist.gp);
		vtxdist.swap(dist.vtxdist);
		partitions.swap(dist.partitions);
		v_per_proc.swap(dist.v_per_proc);
		verticesGotWeights = dist.verticesGotWeights;

		return *this;
	}
};

#endif /* SRC_DECOMPOSITION_DISTRIBUTION_PARMETISDISTRIBUTION_HPP_ */
