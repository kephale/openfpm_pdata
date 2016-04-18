/*
 * DistParMetisDistribution.hpp
 *
 *  Created on: Nov 19, 2015
 *      Author: Antonio Leo
 */

#ifndef SRC_DECOMPOSITION_DISTRIBUTION_DISTPARMETISDISTRIBUTION_HPP_
#define SRC_DECOMPOSITION_DISTRIBUTION_DISTPARMETISDISTRIBUTION_HPP_

#include "SubdomainGraphNodes.hpp"
#include "parmetis_dist_util.hpp"
#include "Graph/dist_map_graph.hpp"
#include "Graph/DistGraphFactory.hpp"

/*! \brief Class that distribute sub-sub-domains across processors using ParMetis Library (totally distributed version)
 *
 * Given a graph and setting Computational cost, Communication cost (on the edge) and
 * Migration cost or total Communication costs, it produce the optimal balanced distribution
 *
 * In addition to Metis it provide the functionality to refine the previously computed
 * decomposition
 *
 * It uses a distributed graph (DistGraph_CSR) to store the distribution
 *
 * ### Initialize a DistParMetis Cartesian graph and decompose
 * \snippet Distribution_unit_tests.hpp Initialize a DistParMetis Cartesian graph and decompose
 *
 * ### refine the decomposition with dist_parmetis
 * \snippet Distribution_unit_tests.hpp refine with parmetis the decomposition
 *
 */
template<unsigned int dim, typename T>
class DistParMetisDistribution
{
	// Vcluster
	Vcluster & v_cl;

	// Structure storing the Cartesian grid information
	grid_sm<dim, void> gr;

	// Rectangular domain to decompose
	Box<dim, T> domain;

	// Processor sub-sub-domain graph
	DistGraph_CSR<nm_v, nm_e> g;

	// Graph in the Parmetis format
	DistParmetis<DistGraph_CSR<nm_v, nm_e>> parmetis_graph;

	// Distribution vector needed for Parmetis
	openfpm::vector<idx_t> vtxdist;

	// Vector containing the partitions computed by Parmetis
	openfpm::vector<openfpm::vector<idx_t>> partitions;

	// Flag to check if weights are set on vertices
	bool verticesGotWeights = false;

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

	/*! Constructor for the DistParMetis class
	 *
	 * \param v_cl Vcluster to use as communication object in this class
	 */
	DistParMetisDistribution(Vcluster & v_cl) :
			v_cl(v_cl), parmetis_graph(v_cl, v_cl.getProcessingUnits()), vtxdist(v_cl.getProcessingUnits() + 1), partitions(v_cl.getProcessingUnits())
	{
	}

	/*! \brief Initialize the distribution graph
	 *
	 * \param grid Grid of the decomposition
	 * \param dom Domain of the simulation
	 */
	void createCartGraph(grid_sm<dim, void> & grid, Box<dim, T> dom)
	{
		// Set grid and domain
		gr = grid;
		domain = dom;

		// Create sub graph
		DistGraphFactory<dim, DistGraph_CSR<nm_v, nm_e>> dist_g_factory;
		g = dist_g_factory.template construct<NO_EDGE, T, dim - 1, 0>(gr.getSize(), domain);
		g.getDecompositionVector(vtxdist);

		if (dim == 2)
			for (size_t i = 0; i < g.getNVertex(); i++)
				g.vertex(i).template get<nm_v::x>()[2] = 0;

	}

	/*! \brief Get the decomposition graph of the processor
	 *
	 */
	DistGraph_CSR<nm_v, nm_e> & getGraph()
	{
		return g;
	}

	/*! \brief Decompose the graph using Parmetis
	 *
	 *	Decompose the graph using the PartKway algorithm
	 *
	 */
	void decompose()
	{
		// Initialize sub graph in Parmetis format
		parmetis_graph.initSubGraph(g);

		// Decompose
		parmetis_graph.decompose<nm_v::proc_id>(g);

		// Get result partition for this processors
		idx_t *partition = parmetis_graph.getPartition();

		// Send the vertices in the processors accordingly to the partition vector
		for (size_t i = 0, j = g.firstId(); i < g.getNVertex() && j <= g.lastId(); i++, j++)
		{
			if ((size_t)partition[i] != v_cl.getProcessUnitID())
				g.q_move(g.nodeById(j), partition[i]);
		}
		g.redistribute();
	}

	/*! \brief Refine current decomposition using Parmetis
	 *
	 * It makes a refinement of the current decomposition using Parmetis function RefineKWay
	 * After that it also does the remapping of the graph
	 *
	 */
	void refine()
	{
		// Reset parmetis graph and reconstruct it
		parmetis_graph.reset(g);

		// Refine
		parmetis_graph.refine<nm_v::proc_id>(g);

		// Get result partition for this processors
		idx_t *partition = parmetis_graph.getPartition();

		// Send the vertices in the processors accordingly to the partition vector
		for (size_t i = 0, j = g.firstId(); i < g.getNVertex() && j <= g.lastId(); i++, j++)
		{
			if ((size_t)partition[i] != v_cl.getProcessUnitID())
				g.q_move(g.nodeById(j), partition[i]);
		}
		g.redistribute();
	}

	/*! \brief Compute the un-balance value
	 *
	 * \return the un-balance value
	 */
	float getUnbalance()
	{
		long t_cost = getProcessorLoad();

		long min, max, sum;
		float unbalance;

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

	/*! \brief Function that return the position of the vertex in the space
	 *
	 * \param id vertex id
	 * \param pos vector that will contain x, y, z
	 *
	 */
	void getSubSubDomainPosition(size_t id, T (&pos)[dim])
	{
		if (id >= g.getNVertex())
			std::cerr << "Position - Such vertex doesn't exist (id = " << id << ", " << "total size = " << g.getNVertex() << ")\n";

		pos[0] = g.vertex(id).template get<nm_v::x>()[0];
		pos[1] = g.vertex(id).template get<nm_v::x>()[1];
		if (dim == 3)
			pos[2] = g.vertex(id).template get<nm_v::x>()[2];
	}

	/*! \brief Function that set the weight of the vertex
	 *
	 * \param id vertex id
	 * \param weight to give to the vertex
	 *
	 */
	inline void setComputationCost(size_t id, size_t weight)
	{
		verticesGotWeights = true;

		if (id >= g.getNVertex())
			std::cerr << "Weight - Such vertex doesn't exist (id = " << id << ", " << "total size = " << g.getNVertex() << ")\n";

		// If the vertex is inside this processor update the value
		g.vertex(id).template get<nm_v::computation>() = weight;

	}

	/*! \brief Checks if weights are used on the vertices
	 *
	 * \return true if weights are used in the decomposition
	 */
	bool weightsAreUsed()
	{
		return verticesGotWeights;
	}

	/*! \brief Function that get the weight of the vertex
	 *
	 * \param id vertex id
	 *
	 */
	size_t getVertexWeight(size_t id)
	{
		if (id >= g.getNVertex())
			std::cerr << "Such vertex doesn't exist (id = " << id << ", " << "total size = " << g.getTotNVertex() << ")\n";

		return g.vertex(id).template get<nm_v::computation>();
	}

	/*! \brief Compute the processor load counting the total weights of its vertices
	 *
	 * \return The computational load of the processor graph
	 */
	size_t getProcessorLoad()
	{
		size_t load = 0;

		for (size_t i = 0; i < g.getNVertex(); i++)
		{
			load += g.vertex(i).template get<nm_v::computation>();
		}
		return load;
	}

	/*! \brief Set migration cost of the vertex
	 *
	 * \param id of the vertex
	 * \param migration cost of the migration
	 */
	void setMigrationCost(size_t id, size_t migration)
	{
		if (id >= g.getNVertex())
			std::cerr << "Migration - Such vertex doesn't exist (id = " << id << ", " << "total size = " << g.getNVertex() << ")\n";

		g.vertex(id).template get<nm_v::migration>() = migration;
	}

	/*! \brief Set communication cost of the edge id
	 *
	 * \param v_id Id of the source vertex of the edge
	 * \param e i child of the vertex
	 * \param communication Communication value
	 */
	void setCommunicationCost(size_t v_id, size_t e, size_t communication)
	{
		g.getChildEdge(v_id, e).template get<nm_e::communication>() = communication;
	}

	/*! \brief Returns total number of sub-sub-domains in the distribution graph (number of vertices in the graph)
	 *
	 */
	size_t getNSubSubDomains()
	{
		return g.getNVertex();
	}

	/*! \brief Returns total number of neighbors of the sub-sub-domain
	 *
	 * \param i id of the sub-sub-domain
	 */
	size_t getNSubSubDomainNeighbors(size_t id)
	{
		if (id >= g.getNVertex())
			std::cerr << "Neighbors - Such vertex doesn't exist (id = " << id << ", " << "total size = " << g.getNVertex() << ")\n";

		return g.getNChilds(id);
	}

	/*! \brief Print current graph and save it to file
	 *
	 */
	void write(const std::string & file)
	{
		VTKWriter<DistGraph_CSR<nm_v, nm_e>, DIST_GRAPH> gv2(g);
		gv2.write(std::to_string(file + ".vtk"));
	}

	/*! \brief Operator definition for DistParMetisDistribution class
	 *
	 * \param dist
	 */
	const DistParMetisDistribution<dim, T> & operator=(const DistParMetisDistribution<dim, T> & dist)
	{
		v_cl = dist.v_cl;
		gr = dist.gr;
		domain = dist.domain;
		g = dist.g;
		vtxdist = dist.vtxdist;
		partitions = dist.partitions;
		verticesGotWeights = dist.verticesGotWeights;

		return *this;
	}

	/*! \brief Operator definition for DistParMetisDistribution class
	 *
	 *	\param dist
	 */
	const DistParMetisDistribution<dim, T> & operator=(const DistParMetisDistribution<dim, T> && dist)
	{
		v_cl = dist.v_cl;
		gr = dist.gr;
		domain = dist.domain;
		g.swap(dist.g);
		vtxdist.swap(dist.vtxdist);
		partitions.swap(dist.partitions);
		verticesGotWeights = dist.verticesGotWeights;

		return *this;
	}
};

#endif /* SRC_DECOMPOSITION_DISTRIBUTION_DISTPARMETISDISTRIBUTION_HPP_ */
