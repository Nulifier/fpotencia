/*
 * File:   Circuit.h
 * Author: Santiago Peñate Vera
 *
 * Created on 6 de agosto de 2014, 10:06
 * Copyright (C) 2014 Santiago Peñate Vera
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */
#pragma once

#include <cmath>
#include <optional>
#include <vector>

#include "fpotencia_libs.h"

#include "Bus.h"
#include "ExternalGrid.h"
#include "Generator.h"
#include "Line.h"
#include "Load.h"
#include "Shunt.h"
#include "Solution.h"
#include "Transformer.h"

namespace fPotencia {
	class Circuit final {
	public:
		using Buses = std::vector<Bus>;
		using Loads = std::vector<Load>;
		using Generators = std::vector<Generator>;
		using Lines = std::vector<Line>;
		using Transformers = std::vector<Transformer>;
		using Shunts = std::vector<Shunt>;
		using ExternalGrids = std::vector<ExternalGrid>;

		/// Addmitance matrix.
		sp_cx_mat Y;

		/// Full impedance matrix
		const cx_mat& Z() {
			if (!m_Z) {
				compose_Z();
			}
			return *m_Z;
		}

		/// Reduced circuit impedance matrix, excluding the slack bus.
		const cx_mat& Zred() {
			if (!m_Zred) {
				compose_Zred();
			}
			return *m_Zred;
		};

		/// Buses contained in the grid
		Buses buses;

		/// All loads in the grid
		Loads loads;

		/// All generators connected to the grid
		Generators generators;

		/// All lines in the grid
		Lines lines;

		/// All transformers in the grid
		Transformers transformers;

		/// All shunts present in the grid
		Shunts shunts;

		/// References external grids connected to this one
		ExternalGrids externalGrids;

		/// Indices of PQ (load) buses
		std::vector<unsigned int> loadBusIndices;

		/// Indices of PV (generator) buses
		std::vector<unsigned int> generatorBusIndices;

		/// Indicies of VD (slack) bus(es)
		std::vector<unsigned int> slackBusIndices;

		/**
		 * The voltage applyed by default to the initial solution.
		 * This default voltage is a complex number and initialized to
		 * $1 + 0j$ by the constructor.
		 */
		cx_double default_voltage = cx_double(1.0, 0.0);

		/**
		 * Prepares the circuit for the solver.
		 * Usually it only needs to be done once, unless the circuit topology
		 * changes. i.e. addition of a line
		 */
		void compile(bool guess_angles);

		/**
		 * Adds a bus to the busses list.
		 * Busses must be added like this, becouse this function asigns them a
		 * number.
		 */
		void add_Bus(Bus& bus);

		/**
		 * Remove bus by object: Removes a bus from the list, but all the buses
		 * indices remain the same.
		 */
		void remove_Bus(Bus bus);

		solution get_initial_solution();

		cx_solution get_initial_cx_solution();

		void set_solution(cx_solution sol);

		/// Gets the real part of a Y element.
		[[nodiscard]] double G(int i, int j) const { return Y.coeff(i, j).real(); }

		/// Gets the imaginary part of a Y element.
		[[nodiscard]] double B(int i, int j) const { return Y.coeff(i, j).imag(); }

		/**
		 * Checks the correctness of the powr flow of the circuit current
		 * solution.
		 * Checks the equation S= Vx(YxV)*
		 */
		void check_solution();

		/**
		 * Performs a DC power flow analysis to calculate an initial guess of the circuit
		 * voltage angles that remain saved at the vector 'dc_angles'
		 */
		void correct_initial_solution();

		/// Print circuit power flows and voltages.
		void print() const;

		void print_buses_state();

		/**
		 * Print complex matrix
		 */
		void printCXMat(const cx_mat& m, const std::string& header) const;

		/**
		 * Print sparse matrix
		 */
		void printMat(const sp_mat& m, const std::string& header) const;

		/* This function sets the load and generation power values (in actual
		 * values not in p.u) and updates the solution objects, keeping the
		 * lattest voltages. This is usefull for continuation power flow, where
		 * the previous soution is recycled.
		 */
		// void setPowerValues(double loadP[], double loadQ[], double genP[],
		// double genQ[]);

	private:
		/*
		 * Removes a column and a row at the index k
		 */
		static void removeCross(cx_mat& matrix, unsigned int k);

		/*
		 * Adds a column and a row at the index k
		 */
		static void expandOnPoint(cx_mat& matrix, unsigned int k);

		/*
		 * Returns the index of a bus given the bus object
		 */
		int find_bus(Bus bus);

		/*
		 * Returns the index of a bus given the bus name
		 */
		int find_bus(std::string bus_name);

		/*
		 * Generates the circuit initial solutions
		 */
		void generate_initial_solution(bool keep_last_solution = false);

		/// Returns the magnitude of the maximum power of any bus bar.
		[[nodiscard]] double get_max_power() const;

		/*
		 * Calculates the current and power flows on the branch elements given
		 * the solution (the solution is calculates by a solver module)
		 */
		void calculate_flows(const cx_solution& sol);

		/*
		 * Composes the circuit admittance matrix
		 */
		void compose_Y();

		/// Composes the reduced circuit impedance matrix by inverting the
		/// admittance matrix.
		void compose_Zred();

		/// Composes the circuit impedance matrix by inverting the admittance
		/// matrix.
		void compose_Z();

		/*
		 * Vector of voltage angles obtained in the function
		 * 'correct_initial_solution()'
		 */
		vec dc_angles;

		/*
		 * Solution for polar voltage sovers like the Newton-Raphson solver
		 */
		solution sol;

		/*
		 * Solution for cartesian voltage solvers like the Gauss (Z-Matrix)
		 * solver
		 */
		cx_solution cx_sol;

		/*
		 * Default base power, but it is updated to a suitable
		 * one when compiling the circuit
		 */
		double Sbase = 100;

		double Zbase;

		double Ybase;

		double Vbase;

		// Laziliy evaluated values
		std::optional<cx_mat> m_Z;
		std::optional<cx_mat> m_Zred;
	};
}
