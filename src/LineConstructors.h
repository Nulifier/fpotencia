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
#pragma once

#include "Solution.h"
#include "CommonFunctions.h"
#include "fpotencia_libs.h"

/*
 * This file contains the neccesary classes to model 3-phase overhead lines
 * and cables for distribution.
 * 
 * the refference book is:
 * Distribution systems modelling and analisys (3rd Edition)
 * by William H Kersting, CRC press.
 * 
 * The magnitudes are as used in the book, that means a quite heretic mix of 
 * units like they like in the US.
 * 
 * Units converters from and to the international metric system are provided
 */

namespace fPotencia {
	/*
	 * Carsosn's equations to compute a line impedance and capacitance
	 */
	class CarsonEquations {
	public:
		CarsonEquations(double frequency, double earth_resistivity, double air_permitivity);

		virtual ~CarsonEquations();

		/*
		 */
		cx_double zii(double r_i, double GMR_i);

		/*
		 */
		cx_double zij(double Dij);

		/*
		 */
		double Pii(double Sij, double RDi);

		/*
		 */
		double Pij(double Sij, double Dij);



		/*Equations constants that are calculated withinthe class*/
		double A;
		double C1 = .0;
		double C2 = .0;
		double C3 = .0;
	private:

		/*
		 */
		void recalculateConstants();

		double G = 0.1609347e-3; //Ohm/mile
		double freq = 50.0; //Hz
		double ro = 100.0; //Ohm-m
		double epsylon = 1.4240e-2; //uF/mile
	};

	/*
	 * Simple metalic conductor  
	 *
	 * When a conductor is displayes in duplex, triplex or cuadruplex
	 * configurations, the GMR changes and of course the resistance
	 */
	class Conductor final {
	public:

		Conductor();

		/*
		 * name: name of the conductor
		 * geometric_mean_radius: GMR (in ft)
		 * resistance: conductor resistance (in Ohm/mile)
		 * diameter: Conductor diameter (in inch)
		 * capacity: conductor maximum amperes (in A)
		 */
		Conductor(std::string name, double geometric_mean_radius, double resistance, double diameter, double capacity, Units_mode unit_mode);

		std::string Name = "";
		double GMR = 0.0; //ft
		double r = 0.0; //ohm/mile
		double d = 0.0; //inches
		double Capacity = 0.0; //Amperes

		Conductor_Type type = Conductor_Type::Cable;
	};

	/*
	 * Cable containing distributed embeeded neutral   
	 */
	class UndergroundCable {
	public:

		UndergroundCable();
		virtual ~UndergroundCable();


		/*
		 * Concentric neutral cable
		 * phase_conductor,
		 * neutral_conductor,
		 * number_of_neutrals,
		 * cable_diameter: cable diameter (in inches)
		 */
		UndergroundCable(Conductor phase_conductor,
				Conductor neutral_conductor,
				uint number_of_neutrals,
				double cable_diameter,
				double relative_permitivity,
				Units_mode unit_mode);

		/*
		 * Tape shield cable
		 * phase_conductor
		 * tape_diameter: tape diameter (in inches)
		 * tape_thickness: tape thickness (in mm)
		 * earth_resistivity: earth resistivity (in ohm-m, a ussual value is 100)
		 */
		UndergroundCable(Conductor phase_conductor,
				double tape_diameter,
				double tape_thickness,
				double earth_resistivity,
				double relative_permitivity,
				Units_mode unit_mode);


		Cable_Type type;
		Conductor phase;
		Conductor neutral;

		//Especific of COncentric neutral cables
		double R = 0.0; //ft

		//Especific of tape shield cables
		double GMR_shield = 0.0; //ft

		cx_double Cshunt; //cable cappacitance in uF/mile

	private:

		Units_mode units;

		double Eps_0; //permitivity of free space in uF/mile
		double Eps_r; //relative permitivity of the medium between the phase condutor and the neutrals

		/*typical range of values for Eps_r
		 * PVC -> 3.4 ~ 8.0
		 * EPR -> 2.5 ~ 3.5
		 * PE -> 2.5 ~ 2.6
		 * XLPE -> 2.3 ~ 6.0
		 */
	};

	/*
	 * Underground line definition
	 * An underground lineis composed of one or more underground cables
	 */
	class UndergroundLine {
	public:

		UndergroundLine();

		virtual ~UndergroundLine();

		/*
		 * cbl: cable
		 * x: x position (in ft)
		 * y: y position (in ft)
		 * these are compared to the GMR therefore the GMR must be in the same units
		 */
		void addCable(UndergroundCable cbl, ConnectionPhase phase, double x, double y);

		/*All this vectors have 3 items: A, B, C*/
		Conductor conductors[3];
		Conductor neutrals[3];
		cx_double conductors_pos[3];
		cx_double neutrals_pos[3];
		cx_double C_shunt[3];
	};

	/*
	 * Overhead line definition
	 * It is composed by unisolated conductors
	 */
	class OverheadLine {
	public:
		OverheadLine();

		virtual ~OverheadLine();

		void addConductor(Conductor conductor, ConnectionPhase phase, double x, double y);

		Conductor conductors[3];
		cx_double conductors_pos[3];
	};

	/*
	 * A ditch can host one or more underground lines
	 */
	class Ditch {
	public:

		Ditch(Units_mode unit_mode);

		Ditch(double frequency, double earth_resistivity, double air_permitivity, Units_mode unit_mode);

		virtual~Ditch();

		void addLine(UndergroundLine line);

		void addNeutral(Conductor conductor, double x, double y);

		cx_mat Zabc();

		cx_mat Yabc_shunt();
		
		cx_mat Neutral_transformation_Matrix();

		double freq;
		double ro;
		double epsylon = 8.85e-12; // F/m

		std::vector<UndergroundLine> lines;

		//phases of the cables 
		std::vector<Conductor> phases;
		std::vector<cx_double> phases_pos;

		//neutrals of the cables
		std::vector<Conductor> neutrals;
		std::vector<cx_double> neutrals_pos;

		//extra neutrals
		std::vector<Conductor> common_neutrals;
		std::vector<cx_double> common_neutrals_pos;

		//all the conductors and their positions ordered
		std::vector<Conductor> conductors;
		std::vector<cx_double> conductors_pos;

		std::vector<cx_double> C_shunt;

	private:

		double D(uint i, uint j);

		void Compile();
		
		cx_mat Zprimitive();
		
		Units_mode units;
	};

	/*
	 * A Tower can host one or more underground lines
	 */
	class Tower {
	public:

		Tower(Units_mode unit_mode);

		Tower(double frequency, double earth_resistivity, double air_permitivity, Units_mode unit_mode);

		virtual~Tower();

		void addLine(OverheadLine line);

		void addNeutral(Conductor conductor, double x, double y);

		cx_mat Zabc();

		cx_mat Yabc_shunt();
		
		cx_mat Neutral_transformation_Matrix();

		double freq;
		double ro;
		double epsylon = 8.85e-12; // F/m

		std::vector<OverheadLine> lines;

		//neutrals of the cables
		std::vector<Conductor> neutrals;
		std::vector<cx_double> neutrals_pos;

		//all the conductors and their positions ordered
		std::vector<Conductor> conductors;
		std::vector<cx_double> conductors_pos;

	private:

		double D(uint i, uint j);

		double S(uint i, uint j);

		void Compile();
		
		cx_mat Zprimitive();
		
		Units_mode units;
	};
}
