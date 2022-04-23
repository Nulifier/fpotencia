/* 
 * File:   LineConstructors.h
 * Author: Santiago Peñate Vera
 *
 * Created on 15 de March de 2015, 10:05
 * Copyright (C) 2015 Santiago Peñate Vera
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "LineConstructors.h"

namespace fPotencia {

	/***************************************************************************
	 * Carson equations
	 ***************************************************************************/

	CarsonEquations::CarsonEquations(double frequency, double earth_resistivity, double air_permitivity) {
		freq = frequency;
		ro = earth_resistivity;
		epsylon = air_permitivity;

		recalculateConstants();
	}

	/*
	 */
	void CarsonEquations::recalculateConstants() {
		C1 = PI * PI * freq * G; //Ohm/mile
		C2 = 4 * PI * freq * G; //Ohm/mile
		C3 = 7.6786 + 0.5 * log(ro / freq); //?
		A = 1.0 / (2.0 * PI * epsylon);
	}

	/*
	 * Auto impedance
	 * r_i: resistance of the conductor i
	 * GMR_i: Geometric Mean Radius of the conductor i
	 */
	cx_double CarsonEquations::zii(double r_i, double GMR_i) {
		return cx_double(r_i + C1, C2 * (log(1.0 / GMR_i) + C3));
	}

	/*
	 * Mutual impedance
	 * Dij: Distance from conductor i to conductor j
	 */
	cx_double CarsonEquations::zij(double Dij) {
		return cx_double(C1, C2 * (log(1.0 / Dij) + C3));
	}

	/*
	 * Self potential coefficient
	 * Sii: Distance from the conductor to its image
	 * RDi: Radius of the conductor
	 */
	double CarsonEquations::Pii(double Sii, double RDi) {
		return A * log(Sii / RDi);
	}

	/*
	 * mutual potential coefficient
	 * Sij: Distance from the conductor to the image of the conductor l
	 * Dij: Distance from conductor i to conductor j
	 */
	double CarsonEquations::Pij(double Sij, double Dij) {
		return A * log(Sij / Dij);
	}

	/***************************************************************************
	 * Conductor
	 ***************************************************************************/

	Conductor::Conductor() {
	}

	Conductor::Conductor(std::string name, double geometric_mean_radius, double resistance, double diameter, double capacity, Units_mode unit_mode) {

		Name = name;

		if (unit_mode == Units_mode::Metric) {
			GMR = geometric_mean_radius * m2ft; //m -> ft
			r = resistance / km2mile; //ohm/m -> ohm/mile
			d = diameter * m2in; //m -> inches
		} else {
			GMR = geometric_mean_radius; //ft
			r = resistance; //ohm/mile
			d = diameter; //inches
		}

		Capacity = capacity; //Ampheres
	}

	/***************************************************************************
	 * Underground Cable
	 ***************************************************************************/

	UndergroundCable::UndergroundCable(Conductor phase_conductor, Conductor neutral_conductor,
			uint number_of_neutrals, double cable_diameter,
			double relative_permitivity, Units_mode unit_mode) {

		std::cout << "Concentric neutral cable" << std::endl;
		double k, d_od;
		units = unit_mode;

		//This is the constructor for concentric neutral cables
		type = Cable_Type::ConcentricNeutral;

		phase = phase_conductor;
		neutral = neutral_conductor;
		Eps_r = relative_permitivity;
		k = number_of_neutrals;
		Eps_0 = 0.01420; //uF/mile

		if (units == Units_mode::Metric) {
			d_od = cable_diameter * m2in;
		} else {
			d_od = cable_diameter;
		}

		//calculated parameters

		//Radius that passes by the center of the concentric neutrals (in)
		R = (d_od - neutral.d) / 24.0; //ft ( 24 is 2*12, where 12 is the inch-feet conversion)

		//the neutral cable is recalculated as an equivalent neutral 
		//given the number of neutrals (k)        
		//See Kersting pag 102
		neutral.GMR = pow(neutral.GMR * k * pow(R, k - 1), 1.0 / k); //ft
		neutral.r = neutral.r / k; //ohm/mile

		double RDc = phase_conductor.d / 2.0;
		double RDs = neutral_conductor.d / 2.0;

		Cshunt = cx_double(0, (2.0 * PI * Eps_r * Eps_0) / (log(R * 12 / RDc) - (1.0 / k) * log(k * RDs / (R * 12))));
	}

	UndergroundCable::UndergroundCable(Conductor phase_conductor, double tape_diameter,
			double tape_thickness, double earth_resistivity,
			double relative_permitivity, Units_mode unit_mode) {

		std::cout << "Tape shield cable" << std::endl;
		double d_s, T, ro;
		units = unit_mode;

		//This is the constructor for tape shield cables
		type = Cable_Type::TapeShield;

		phase = phase_conductor;
		Eps_r = relative_permitivity;
		T = tape_thickness; //mm
		ro = earth_resistivity; //Ohm-m

		if (units == Units_mode::Metric) {
			d_s = tape_diameter * m2in; //in            
		} else {
			d_s = tape_diameter; //in
		}
		//calculated parameters        
		GMR_shield = ((d_s / 2.0) - (T / 2000.0)) / 12.0; //ft
		//double r_shield = 7.9385e8 * ro / (T * d_s); //ohm/mile
		double r_shield = 0.18826 * ro / (T * d_s); //ohm/mile

		//generate the neutral equivalent to the tape
		//string name, double geometric_mean_radius, double resistance,  double diameter, double capacity
		neutral = Conductor("tape", GMR_shield, r_shield, 0.0, 0.0, units);
		neutral.type = Conductor_Type::Tape_Shield;

		double RDc = phase_conductor.d / 2.0;

		Cshunt = cx_double(0, (2.0 * PI * Eps_r * Eps_0) / log(GMR_shield * 12 / RDc));
	}

	/***************************************************************************
	 * Underground line
	 ***************************************************************************/

	UndergroundLine::UndergroundLine() {
		Conductor cond;
		cond.type = Conductor_Type::Void;

		/*Intialize all 3 phases with a "blank" cable*/

		for (uint i = 0; i < 3; i++) {
			conductors[i] = cond;
			conductors_pos[i] = cx_double(0.0, 0.0);
			neutrals[i] = cond;
			neutrals_pos[i] = cx_double(0.0, 0.0);
			C_shunt[i] = cx_double(0.0, 0.0);
		}

	}

	void UndergroundLine::addCable(UndergroundCable cbl, ConnectionPhase phase, double x, double y) {
		uint idx = 0;

		if (phase == ConnectionPhase::A)
			idx = 0;
		else if (phase == ConnectionPhase::B)
			idx = 1;
		else if (phase == ConnectionPhase::C)
			idx = 2;


		conductors[idx] = cbl.phase;
		conductors_pos[idx] = cx_double(x, y);

		neutrals[idx] = cbl.neutral;

		if (cbl.type == Cable_Type::ConcentricNeutral) {
			neutrals_pos[idx] = cx_double(x, y + cbl.R);
		} else if (cbl.type == Cable_Type::TapeShield) {
			neutrals_pos[idx] = cx_double(x, y);
		}

		C_shunt[idx] = cbl.Cshunt;
	}

	/***************************************************************************
	 * Overhead line
	 ***************************************************************************/

	OverheadLine::OverheadLine() {

		/*Intialize all 3 phases with a "blank" cable*/
		Conductor c;
		c.type = Conductor_Type::Void;

		for (uint i = 0; i < 3; i++) {
			conductors[i] = c;
			conductors_pos[i] = cx_double(0.0, 0.0);
		}
	}

	OverheadLine::~OverheadLine() {
	}

	void OverheadLine::addConductor(Conductor conductor, ConnectionPhase phase, double x, double y) {
		uint idx = 0;

		if (phase == ConnectionPhase::A)
			idx = 0;
		else if (phase == ConnectionPhase::B)
			idx = 1;
		else if (phase == ConnectionPhase::C)
			idx = 2;


		conductors[idx] = conductor;
		conductors_pos[idx] = cx_double(x, y);
	}

	/***************************************************************************
	 * Underground line grouping: Ditch
	 ***************************************************************************/

	Ditch::Ditch(Units_mode unit_mode) {
		units = unit_mode;
	}

	Ditch::Ditch(double frequency, double earth_resistivity, double air_permitivity, Units_mode unit_mode) {
		freq = frequency;
		ro = earth_resistivity;
		epsylon = air_permitivity;
		units = unit_mode;
	}

	Ditch::~Ditch() {
	}

	void Ditch::addLine(UndergroundLine line) {
		lines.push_back(line);
	}

	void Ditch::addNeutral(Conductor conductor, double x, double y) {
		common_neutrals.push_back(conductor);
		common_neutrals_pos.push_back(cx_double(x, y));
	}

	void Ditch::Compile() {
		neutrals.clear();
		neutrals_pos.clear();
		phases.clear();
		phases_pos.clear();
		conductors.clear();
		conductors_pos.clear();
		C_shunt.clear();

		//Add the common neutral conductors
		for (uint i = 0; i < common_neutrals.size(); i++) {
			neutrals.push_back(common_neutrals[i]);
			neutrals_pos.push_back(common_neutrals_pos[i]);
		}

		//Add the lines conductors
		for (uint i = 0; i < lines.size(); i++) {
			//add the lines neutral conductors
			for (uint j = 0; j < 3; j++) {
				neutrals.push_back(lines[i].neutrals[j]);
				neutrals_pos.push_back(lines[i].neutrals_pos[j]);
			}

			//add the lines phase conductors
			for (uint j = 0; j < 3; j++) {
				phases.push_back(lines[i].conductors[j]);
				phases_pos.push_back(lines[i].conductors_pos[j]);

				C_shunt.push_back(lines[i].C_shunt[j]);
			}
		}

		//Add the conductors in order; first phases, then neutrals
		//this is like this to be able to do the kron reduction
		for (uint i = 0; i < phases.size(); i++) {
			conductors.push_back(phases[i]);
			conductors_pos.push_back(phases_pos[i]);
		}

		for (uint i = 0; i < neutrals.size(); i++) {
			conductors.push_back(neutrals[i]);
			conductors_pos.push_back(neutrals_pos[i]);
		}

	}

	/*
	 * Conductors distance
	 */
	double Ditch::D(uint i, uint j) {
		double d;
		d = abs(conductors_pos[i] - conductors_pos[j]);
		if (d == 0.0) {
			/* Both conductors are in the same place
			 * This might be because the tape shields virtual neutral
			 * are modelled as a neutral cable located concentric with 
			 * the phase conductor. In this case d is the tape shield GMR
			 */
			if (conductors[i].type == Conductor_Type::Tape_Shield) {
				d = conductors[i].GMR;

			} else {
				d = conductors[j].GMR;
			}
			//cout<<"\nTape shielded D:"<<d<<endl;
		}
		return d;
	}

	/*
	 * Calculates the primitive impedance matrix of the Cable group
	 */
	cx_mat Ditch::Zprimitive() {
		Compile();
		CarsonEquations eq(freq, ro, epsylon);
		uint n = conductors.size();
		uint m = neutrals.size();
		uint p = phases.size();

		cx_mat z(n, n);
		z.setZero(n, n);

		for (uint i = 0; i < n; i++) {
			for (uint j = 0; j < n; j++) {
				if (conductors[i].type != Conductor_Type::Void && conductors[j].type != Conductor_Type::Void) {
					if (i == j) {
						z(i, j) = eq.zii(conductors[i].r, conductors[i].GMR);
					}
					else {
						z(i, j) = eq.zij(D(i, j));
					}
				}
			}
		}
		
		return z;
	}

	/*
	 * Calculates the ABC impedance matrix (aplying the Kron reduction when needed)
	 */
	cx_mat Ditch::Zabc() {
		
		cx_mat z = Zprimitive();
		uint n = conductors.size();
		uint m = neutrals.size();
		uint p = phases.size();

		if (units == Units_mode::Metric) {
			if (m > 0)
				return Kron_cx(z, p, m) / mile2km;
			else
				return z / mile2km;
		} else {
			if (m > 0)
				return Kron_cx(z, p, m);
			else
				return z;
		}
	}

	/*
	 * Computes the neutral transformation matrix
	 * Used in computing the neutral cables currents
	 */
	cx_mat Ditch::Neutral_transformation_Matrix() {        
		cx_mat z = Zprimitive();
		uint m = neutrals.size();
		uint n = conductors.size();
		uint p = n - m;
		
		cx_mat znj = z.bottomLeftCorner(m, p);
		cx_mat znn = z.bottomRightCorner(m, m);

		Eigen::FullPivLU<cx_mat> lu(znn);
		return -1*lu.inverse()*znj;
	}

	/*
	 * Return the shunt admittance matrix of the ditch in S/mile
	 */
	cx_mat Ditch::Yabc_shunt() {

		uint p = C_shunt.size();
		cx_mat z(p, p);
		double w = 2.0 * PI * freq;
		z.setZero(p, p);
		for (uint i = 0; i < p; i++)
			z(i, i) = w * C_shunt[i];

		//Up to here z is in uS, we'll return Siemens
		if (units == Units_mode::Metric)
			return z * 1e-6 / mile2km;
		else
			return z * 1e-6;
	}

	/***************************************************************************
	 * Overhead line grouping: Tower
	 ***************************************************************************/

	Tower::Tower(Units_mode unit_mode) {
		units = unit_mode;
	}

	Tower::Tower(double frequency, double earth_resistivity, double air_permitivity, Units_mode unit_mode) {
		freq = frequency;
		ro = earth_resistivity;
		epsylon = air_permitivity;
		units = unit_mode;
	}

	Tower::~Tower() {
	}

	void Tower::addLine(OverheadLine line) {
		lines.push_back(line);
	}

	void Tower::addNeutral(Conductor conductor, double x, double y) {
		neutrals.push_back(conductor);
		neutrals_pos.push_back(cx_double(x, y));
	}

	void Tower::Compile() {
		conductors.clear();
		conductors_pos.clear();

		//Add the lines conductors
		for (uint i = 0; i < lines.size(); i++) {
			//add the lines phase conductors
			for (uint j = 0; j < 3; j++) {
				conductors.push_back(lines[i].conductors[j]);
				conductors_pos.push_back(lines[i].conductors_pos[j]);
			}
		}

		for (uint i = 0; i < neutrals.size(); i++) {
			conductors.push_back(neutrals[i]);
			conductors_pos.push_back(neutrals_pos[i]);
		}

	}

	/*
	 * Conductors distance
	 */
	double Tower::D(uint i, uint j) {
		double d;
		d = abs(conductors_pos[i] - conductors_pos[j]);
		if (d == 0.0) {
			/* Both conductors are in the same place
			 * This might be because the tape shields virtual neutral
			 * are modelled as a neutral cable located concentric with 
			 * the phase conductor. In this case d is the tape shield GMR
			 */
			if (conductors[i].type == Conductor_Type::Tape_Shield)
				d = conductors[i].GMR;
			else
				d = conductors[j].GMR;
		}
		return d;
	}

	/*
	 * Distance from the conductor i to the image of conductor j
	 */
	double Tower::S(uint i, uint j) {
		double d;
		d = abs(conductors_pos[i] - conj(conductors_pos[j]));
		if (d == 0.0) {
			/* Both conductors are in the same place
			 * This might be because the tape shields virtual neutral
			 * are modelled as a neutral cable located concentric with 
			 * the phase conductor. In this case d is the tape shield GMR
			 */
			if (conductors[i].type == Conductor_Type::Tape_Shield)
				d = conductors[i].GMR;
			else
				d = conductors[j].GMR;
		}
		return d;
	}

	/*
	 * Computes the primitive impedance matrix of a line
	 */
	cx_mat Tower::Zprimitive() {
		Compile();
		CarsonEquations eq(freq, ro, epsylon);
		uint n = conductors.size();

		cx_mat z(n, n);
		z.setZero(n, n);

		for (uint i = 0; i < n; i++) {
			for (uint j = 0; j < n; j++) {
				if (conductors[i].type != Conductor_Type::Void && conductors[j].type != Conductor_Type::Void) {
					if (i == j)
						z(i, j) = eq.zii(conductors[i].r, conductors[i].GMR);
					else
						z(i, j) = eq.zij(D(i, j));
				}
			}
		}

		return z;
	}

	/*
	 * Compose the ABC phase impedance matrix
	 */
	cx_mat Tower::Zabc() {
		
		cx_mat z = Zprimitive();
		uint m = neutrals.size();
		uint n = conductors.size();
		
		std::cout<< z <<std::endl;
		
		if (units == Units_mode::Metric) {
			if (m > 0)
				return Kron_cx(z, n - m, m) / mile2km;
			else
				return z / mile2km;
		} else {
			if (m > 0)
				return Kron_cx(z, n - m, m);
			else
				return z;
		}
	}

	/*
	 * Computes the neutral transformation matrix
	 * Used in computing the neutral cables currents
	 */
	cx_mat Tower::Neutral_transformation_Matrix() {
		
		cx_mat z = Zprimitive();
		uint m = neutrals.size();
		uint n = conductors.size();
		uint p = n - m;
		
		cx_mat znj = z.bottomLeftCorner(m, p);
		cx_mat znn = z.bottomRightCorner(m, m);

		Eigen::FullPivLU<cx_mat> lu(znn);
		return -1*lu.inverse() * znj;
	}

	/*
	 * Return the shunt admittance matrix of the tower in S/mile
	 */
	cx_mat Tower::Yabc_shunt() {
		Compile();
		CarsonEquations eq(freq, ro, epsylon);
		uint n = conductors.size();
		uint m = neutrals.size();
		mat z(n, n);

		mat s(n, n);

		for (uint i = 0; i < n; i++) {
			for (uint j = 0; j < n; j++) {
				s(i, j) = S(i, j);
				if (i == j)
					z(i, j) = eq.Pii(S(i, i), conductors[i].GMR);
				else
					z(i, j) = eq.Pij(S(i, j), D(i, j));
			}
		}

		if (m > 0)
			z = Kron(z, n - 1, 1);

		Eigen::FullPivLU<mat> lu(z);
		z = lu.inverse();

		//Up to here z is in uS, we'll return Siemens
		if (units == Units_mode::Metric)
			return cx_double(0, freq * 2.0 * PI) * z * 1e-6 / mile2km;
		else
			return cx_double(0, freq * 2.0 * PI) * z * 1e-6;

	}
}
