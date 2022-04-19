/*
 * File:   main.cpp
 * Author: Santiago Peñate Vera
 *
 * Created on 8 de agosto de 2014, 15:17
 * Copyright (C) 2014 Santiago Peñate Vera
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

//Use mpic++ instead of g++ to compile with parallel capabilities

#include <iostream>
#include <vector>

#include "../fpotencia.h"

#include "gridExamples.h"

using namespace fPotencia;

/*
 * Power flow test script
 */
void power_flow_test() {
	//Eigen::initParallel();

	//Circuit model = circuit_Gutierrez_Castillejos();
	Circuit model = circuit_Lynn_Powell(true);
	//Circuit model = ieee_30_bus_noPV();
	//Circuit model = ieee_30_bus();
	bool estimate_angles = false;
	model.compile(estimate_angles);
	model.printCXMat(model.Y, "Ybus:");
	model.print_buses_state();
	//cout << "Ybus:\n" << model.Y << endl;
	//model.printCXMat(model.Zred, "Zbus:");
	Solver_State state;
	solution prev_sol;

	/*
		std::cout << "\n\n\n NEWTON-RAPHSON POLAR\n\n " << std::endl;
		//New solution
		Solver_NRpolar nrs(model);
		Solver_State state = nrs.solve();
		//nrs.Model.print();
		if (state == Solver_State::Converged) {
			std::cout << "Converged in " << nrs.Iterations << " iterations." << std::endl;
			//nrs.Model.check_solution();
			prev_sol = nrs.Model.get_initial_solution();
			prev_sol.print("Solution:");
		} else {
			std::cout << "NOT Converged " << std::endl;
		}
	
	
		std::cout << "\n\n\n NEWTON-RAPHSON RECTANGULAR \n\n " << std::endl;
		//New solution
		Solver_NRrect nrr(model);
		state = nrr.solve();
		//nrs.Model.print();
		if (state == Solver_State::Converged) {
			std::cout << "Converged in " << nrr.Iterations << " iterations." << std::endl;
			//nrs.Model.check_solution();
			prev_sol = nrr.Model.get_initial_solution();
			prev_sol.print("Solution:");
		} else {
			std::cout << "NOT Converged " << std::endl;
		}


		std::cout << "\n\n\n JACOBI \n\n " << std::endl;
		Solver_Jacobi js(model);
		js.Max_Iter = 100;
		js.EPS = 1e-9;
		state = js.solve();
		if (state == Solver_State::Converged) {
			std::cout << "Converged in " << js.Iterations << " iterations." << std::endl;
			prev_sol = js.Model.get_initial_solution();
			prev_sol.print("Solution:");
		} else {
			std::cout << "NOT Converged " << std::endl;
		}


	std::cout << "\n\n\n IWAMOTO \n\n " << std::endl;
	Solver_Iwamoto Iws(model);
	Iws.Max_Iter = 6;
	Iws.EPS = 1e-9;
	state = Iws.solve();
	if (state == Solver_State::Converged) {
		std::cout << "Converged in " << Iws.Iterations << " iterations." << std::endl;
		prev_sol = Iws.Model.get_initial_solution();
		prev_sol.print("Solution:");
	} else {
		std::cout << "NOT Converged " << std::endl;
	}*/

	std::cout << "\n\n\n NR CURRENT \n\n " << std::endl;
	NRcurrentSolver NRcs(model);
	NRcs.maxIterations = 6;
	NRcs.tolerance = 1e-9;
	state = NRcs.solve();
	if (state == Solver_State::Converged) {
		std::cout << "Converged in " << NRcs.Iterations << " iterations." << std::endl;
		prev_sol = NRcs.Model.get_initial_solution();
		prev_sol.print("Solution:");
	} else {
		std::cout << "NOT Converged " << std::endl;
	}
}

void Lines_constructor_test() {
	double freq = 60; //Hz
	double earth_resistivity = 100; //ohm-m
	double air_permitivity = 1.4240e-2; //uF/mile

	/***************************************************************************
	 * Example 1, one line (3 phases + neutral)
	 **************************************************************************/
	cout << "Example 1, one line (3 phases + neutral)" << endl;

	Conductor conductor1("336,400 26/7 ACSR", 0.0244, 0.306, 0.721, 530, Units_mode::US);
	Conductor neutro("4/0 6/1 ACSR", 0.00814, 0.5920, 0.563, 340, Units_mode::US);

	OverheadLine line;
	line.addConductor(conductor1, ConnectionPhase::A, 0.0, 29.0);
	line.addConductor(conductor1, ConnectionPhase::B, 2.5, 29.0);
	line.addConductor(conductor1, ConnectionPhase::C, 7.0, 29.0);

	Tower tower1(freq, earth_resistivity, air_permitivity, Units_mode::US);
	tower1.addLine(line);
	tower1.addNeutral(neutro, 4.0, 25.0);
	cx_mat Z1 = tower1.Zabc();


	cout << Z1.rows() << "," << Z1.cols() << endl;

	cout << Z1 << endl;
	cout << "Sequence:\n" << abc2seq(Z1) << endl;

	cout << "Shunt:\n" << tower1.Yabc_shunt() << endl;

	/***************************************************************************
	 * Example 2, one tower with two lines (3 + 3 phases + 1 neutral)
	 **************************************************************************/

	cout << "\n\nExample 2, one tower with two lines (3 + 3 phases + 1 neutral)" << endl;
	//Conductor conductor1("336,400 26/7 ACSR", 0.0244, 0.306, 0.721, 530)

	Conductor conductor2("250,000 AA", 0.0171, 0.41, 0.567, 329, Units_mode::US);
	//Conductor neutro("4/0 6/1 ACSR", 0.00814, 0.5920, 0.563, 340);

	OverheadLine line1;
	line1.addConductor(conductor1, ConnectionPhase::A, 0.0, 35.0);
	line1.addConductor(conductor1, ConnectionPhase::B, 2.5, 35.);
	line1.addConductor(conductor1, ConnectionPhase::C, 7.0, 35.0);

	OverheadLine line2;
	line2.addConductor(conductor2, ConnectionPhase::A, 2.5, 33.0);
	line2.addConductor(conductor2, ConnectionPhase::B, 7.0, 33.);
	line2.addConductor(conductor2, ConnectionPhase::C, 0.0, 33.0);


	Tower tower2(freq, earth_resistivity, air_permitivity, Units_mode::US);
	tower2.addLine(line1);
	tower2.addLine(line2);
	tower2.addNeutral(neutro, 4.0, 29.0);

	cx_mat Z2 = tower2.Zabc();


	cout << Z2.rows() << "," << Z2.cols() << endl;
	cout << Z2 << endl;

	cout << "Shunt:\n" << tower2.Yabc_shunt() << endl;

	/***************************************************************************
	 * Example 3: three cables in a trench"
	 **************************************************************************/

	cout << "\n\nExample 3: three cables in a trench" << endl;

	Conductor phase3("250,000 AA", 0.0171, 0.41, 0.567, 329, Units_mode::US);
	Conductor neutro3("14 AWG SLD copper", 0.00208, 14.8722, 0.0641, 20, Units_mode::US);

	UndergroundCable cable3(phase3, neutro3, 13, 1.29, 2.3, Units_mode::US); //XLPE cable

	UndergroundLine line3;
	line3.addCable(cable3, ConnectionPhase::A, 0, 0);
	line3.addCable(cable3, ConnectionPhase::B, 0.5, 0);
	line3.addCable(cable3, ConnectionPhase::C, 1, 0);

	Ditch ditch3(freq, earth_resistivity, air_permitivity, Units_mode::Metric);
	ditch3.addLine(line3);

	cx_mat Z3 = ditch3.Zabc();

	cout << Z3.rows() << "," << Z3.cols() << endl;
	cout << Z3 << endl;
	cout << "Shunt:\n" << ditch3.Yabc_shunt() << endl;



	/***************************************************************************
	 * Example 4: single tape shield cable with neutral
	 **************************************************************************/
	cout << "\n\nExample 4: single tape shield cable with neutral" << endl;

	Conductor phase4("1/0 AA", 0.0111, 0.97, 0.368, 202, Units_mode::US);
	Conductor neutro4("1/0 copper 7 strand", 0.01113, 0.607, 0.368, 310, Units_mode::US);

	UndergroundCable cable4(phase4, 0.88, 5, earth_resistivity, 2.3, Units_mode::US);

	UndergroundLine line4;
	line4.addCable(cable4, ConnectionPhase::C, 0, 0);

	Ditch ditch4(freq, earth_resistivity, air_permitivity, Units_mode::US);
	ditch4.addLine(line4);
	ditch4.addNeutral(neutro4, 0.25, 0.0);

	cx_mat Z4 = ditch4.Zabc();

	cout << Z4.rows() << "," << Z4.cols() << endl;
	cout << Z4 << endl;
	cout << "Shunt:\n" << ditch4.Yabc_shunt() << endl;


	/**********/

	cx_mat3 M;
	cout << "\n\nM:\n" << M << endl;
}

/*
 */
void graph_test() {
	/** Creates a graph with 12 vertices */
	Graph g(7);

	g.addEdge(0, 1); //0
	g.addEdge(0, 2); //1
	g.addEdge(1, 3); //2
	g.addEdge(2, 3); //3
	g.addEdge(2, 5); //4
	g.addEdge(3, 6); //5
	g.addEdge(4, 5); //6
	g.addEdge(4, 6); //7

	g.addEdge(4, 3); //8
	clock_t t1;
	t1 = clock();

	g.BFS(0);

	cout << endl;
	for (int i = 0; i < g.nodes_size(); i++)
		cout << i << "->" << g.Visited[i] << endl;

	cout << "\nNodes" << endl;
	for (int i = 0; i < g.nodes_size(); i++)
		cout << g.Node_labels[i] << ", ";

	cout << "\nEdges" << endl;
	for (int i = 0; i < g.edges_size(); i++)
		cout << g.Edge_labels[i] << ", ";

	float diff = (double) (clock() - t1) / CLOCKS_PER_SEC;
	cout << endl << "The time taken for Breadth first search: " << diff << endl;
}

typedef std::complex<double> cx_double;

void example() {
	Eigen::MatrixXcd mat(4, 4);

	int k = 0;
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++) {
			mat(i, j) = cx_double(k, k * k);
			k++;
		}

	cout << mat << endl;
	cout << mat.topLeftCorner(2, 3) << endl;
}

/*
 */
int main(int argc, char** argv) {

	power_flow_test();
	//example();

	//Lines_constructor_test();

	//graph_test();

	//int k;
	//cin >> k;
	return 0;
}


