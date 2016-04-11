/*
 * DLB.hpp
 *
 *  Created on: Nov 20, 2015
 *      Author: Antonio Leo
 */

#ifndef SRC_DECOMPOSITION_DLB_HPP_
#define SRC_DECOMPOSITION_DLB_HPP_

#include "Plot/GoogleChart.hpp"
#include "Plot/util.hpp"

//! Time structure for statistical purposes
typedef struct
{
	size_t simulationStartTime = 0;
	size_t simulationEndTime;
	double timeStep = 0.1;
	size_t iterationStartTime;
	size_t iterationEndTime;
} Times;

/*! Class that implements the two heuristics to determine when a re-balance of the distribution is needed.
 *
 *  Used heuristics are: SAR and Un-balance Threshold (Default)\n
 *
 *  To chose the heuristic use the method setHeuristic(Heuristic)
 *
 *  In the SAR heuristic the following formula is applied:\n
 *  \f$W_{n} = \frac{\sum_{j=1}^{n} (T_{max}(j) - T_{avg}(j)) + C} {n}\f$
 *
 *  \f$T_{max}(j)\f$ – wall-clock time of bottleneck process in time step j\n
 *  \f$T_{avg}(j)\f$ – average wall-clock time for time step j over all processes\n
 *  \f$C\f$ – cost of re-decomposing the problem\n
 *  \f$n\f$ – number of time steps since last re-decomposition\n
 *  \n
 *  For small n, load balance is good and W decreases since C is amortized over an increasing number of time steps.
 *  As the accumulated idle time starts to dominate, W starts to rise. At this point, C has been fully amortized.
 *  Re-decompose when \f$W_{n} > W_{n-1}\f$\n
 *
 *  In the Un-balance Threshold heuristic the re-balance is triggered when the un-balance level exceeds a certain level.
 *  Levels can be chosen in the ThresholdLevel type.
 */
class DLB
{
public:

	//! Type of DLB heuristics
	enum Heuristic
	{
		SAR_HEURISTIC, UNBALANCE_THRLD
	};

	//! Level of un-balance needed to trigger the re-balance
	enum ThresholdLevel
	{
		THRLD_LOW = 5, THRLD_MEDIUM = 7, THRLD_HIGH = 10
	};

private:

	//! Runtime virtual cluster machine
	Vcluster & v_cl;

	//! Structure that will contain all the timings
	Times timeInfo;

	//! Wn for SAR heuristic
	float w_n = -1;

	//! Computation cost for SAR heuristic
	float c_c = 5;

	//! Number of time-steps since the previous DLB
	size_t n_ts = 1;

	//! Idle time accumulated so far, needed for SAR heuristic
	float i_time = 0;

	//! Vector to collect all timings
	openfpm::vector<long> times;

	//! Type of the heuristic to use
	Heuristic heuristic = UNBALANCE_THRLD;

	//! Un-balance value
	float unbalance = -1;

	//! Threshold value
	ThresholdLevel thl = THRLD_MEDIUM;

	// Vector to store all the W values (need for output)
	openfpm::vector<float> ws;

	/*! \brief Function that gather times informations and decides if a rebalance is needed it uses the SAR heuristic
	 *
	 * \param t
	 *
	 */
	inline bool SAR()
	{
		long t = timeInfo.iterationEndTime - timeInfo.iterationStartTime;
		long t_max = t;
		float t_avg = t;

		// Exchange time informations through processors
		v_cl.max(t_max);
		v_cl.sum(t_avg);
		v_cl.execute();

		t_avg /= v_cl.getProcessingUnits();

		// add idle time to vector
		i_time += t_max - t_avg;

		// Compute Wn
		float nw_n = (i_time + c_c) / n_ts;

		if (w_n == -1)
			w_n = nw_n;

		ws.add(nw_n);

		if (nw_n > w_n)
		{
			i_time = 0;
			n_ts = 1;
			w_n = nw_n;
			return true;
		}
		else
		{
			++n_ts;
			w_n = nw_n;
			return false;
		}
	}

	/*! \brief Check if the un-balance has exceeded the threshold
	 *
	 * \return true if re-balance is needed, false otherwise
	 */
	bool unbalanceThreshold()
	{
		if (unbalance == -1)
		{
			std::cerr << "Error: Un-balance value must be set before checking DLB.";
			return false;
		}

		if (unbalance > thl)
		{
			return true;
		}

		return false;
	}

public:

	/*! \brief Constructor for DLB class
	 *
	 * \param v_cl virtual cluster object
	 */
	DLB(Vcluster & v_cl) :
			v_cl(v_cl)
	{
	}

	/*! \brief Set the heuristic to use (default: un-balance threshold)
	 *
	 * \param h
	 */
	void setHeurisitc(Heuristic h)
	{
		heuristic = h;
	}

	/*! \brief Get the heuristic to use
	 *
	 */
	Heuristic getHeurisitc()
	{
		return heuristic;
	}

	/*! \brief check if a re-balance is needed using the SAR heuristic
	 *
	 */
	bool rebalanceNeeded()
	{
		if (heuristic == SAR_HEURISTIC)
		{
			return SAR();
		}
		else
		{
			return unbalanceThreshold();
		}
	}

	/*! \brief Set start time for the simulation
	 *
	 * \param simulationStartTime time when the whole simulation starts
	 */
	void setSimulationStartTime(size_t t)
	{
		timeInfo.simulationStartTime = t;
	}

	/*! \brief Get start time for the simulation
	 *
	 */
	size_t getSimulationStartTime()
	{
		return timeInfo.simulationStartTime;
	}

	/*! \brief Set end time for the simulation
	 *
	 * \param simulationEndTime time when the whole simulation ends
	 */
	void setSimulationEndTime(size_t t)
	{
		timeInfo.simulationEndTime = t;
	}

	/*! \brief Get end time for the simulation
	 *
	 */
	size_t getSimulationEndTime()
	{
		return timeInfo.simulationEndTime;
	}

	/*! \brief Set start time for the single iteration
	 *
	 */
	void startIteration()
	{
		timeInfo.iterationStartTime = clock();
	}

	/*! \brief Set start time for the single iteration
	 *
	 * \param iterationStartTime time when the single iteration starts
	 */
	void startIteration(size_t t)
	{
		timeInfo.iterationStartTime = t;
	}

	/*! \brief Set end time for the single iteration
	 *
	 */
	void endIteration()
	{
		timeInfo.iterationEndTime = clock();
	}

	/*! \brief Set end time for the single iteration
	 *
	 * \param iterationEndTime time when the single iteration ends
	 */
	void endIteration(size_t t)
	{
		timeInfo.iterationEndTime = t;
	}

	/*! \brief Set time step for the single iteration
	 *
	 * \param timestep value, can be also a 0.1 s
	 */
	void setTimeStep(double t)
	{
		timeInfo.timeStep = t;
	}

	/*! \brief Set the cost of the re-balancing
	 *
	 * \param computation value of the computation cost (default: 5)
	 */
	void setComputationCost(size_t computation)
	{
		c_c = computation;
	}

	/*! \brief Get the cost of the re-balancing
	 *
	 */
	float getComputationCost()
	{
		return c_c;
	}

	/*! \brief Get how many time-steps have passed since the last re-balancing
	 *
	 */
	size_t getNTimeStepSinceDLB()
	{
		return n_ts;
	}

	/*! \brief Set un-balance value
	 *
	 * \param computation value of the computation cost (default: 5)
	 */
	void setUnbalance(float u)
	{
		unbalance = u;
	}

	/*! \brief Set un-balance value
	 *
	 * \param computation value of the computation cost (default: 5)
	 */
	void setThresholdLevel(ThresholdLevel t)
	{
		thl = t;
	}

	/*! \brief print google chart of SAR
	 *
	 */
	void write()
	{
		openfpm::vector<std::string> x;
		openfpm::vector<openfpm::vector<float>> y;

		for (size_t i = 0; i < ws.size(); i++)
			x.add(std::to_string(i));

		for (size_t i = 0; i < ws.size(); i++)
			y.add( { ws.get(i) });

		// Google charts options
		GCoptions options;

		options.title = std::string("SAR heuristic");
		options.yAxis = std::string("W");
		options.xAxis = std::string("Timestep");
		options.lineWidth = 1.0;

		GoogleChart cg;
		cg.AddPointsGraph(x, y, options);
		cg.write("SAR.html");
	}

};

#endif /* SRC_DECOMPOSITION_DLB_HPP_ */
