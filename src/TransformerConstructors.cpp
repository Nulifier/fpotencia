#include "TransformerConstructors.h"


namespace fPotencia {

	/***********************************************************************
	 * Transformer_3phase_Yabc
	 ***********************************************************************/
	
	/*
	 */
	Transformer_3phase_Yabc::Transformer_3phase_Yabc(){}
	
	/*
	 */
	Transformer_3phase_Yabc::~Transformer_3phase_Yabc(){}
	
	/*
	 */
	cx_mat Transformer_3phase_Yabc::Yabc(cx_double leakage_impedance, TransformerConnexionType connexionType, double tap_hv, double tap_lv) {
		//3Phase

		mat Y1, Y2, Y3;
		Y1 << 1, 0, 0,
				0, 1, 0,
				0, 0, 1;

		Y2 << 2, -1, -1,
				-1, 2, -1,
				-1, -1, 2;

		Y3 << -1, 1, 0,
				0, -1, 1,
				1, 0, -1;

		cx_mat3 YI, YII, YIII;
		YI = Y1 * leakage_impedance;
		YII = Y2 * leakage_impedance / 3.0;
		YIII = Y3 * leakage_impedance / sqrt(3.0);
		cx_mat3 Ypp;
		cx_mat3 Yps;
		cx_mat3 Ysp;
		cx_mat3 Yss;

		if (connexionType == TransformerConnexionType::Yg_Yg) {
			Ypp = YI / tap_hv / tap_hv;
			Yps = -YI / tap_hv / tap_lv;
			Ysp = -YI / tap_lv / tap_hv;
			Yss = YI / tap_lv / tap_lv;
		} else if (connexionType == TransformerConnexionType::Yg_Y) {
			Ypp = YII / tap_hv / tap_hv;
			Yps = -YII / tap_hv / tap_lv;
			Ysp = -YII / tap_lv / tap_hv;
			Yss = YII / tap_lv / tap_lv;
		} else if (connexionType == TransformerConnexionType::Yg_D) {
			Ypp = YI / tap_hv / tap_hv;
			Yps = YIII / tap_hv / tap_lv;
			Ysp = YIII.transpose() / tap_lv / tap_hv;
			Yss = YII / tap_lv / tap_lv;
		} else if (connexionType == TransformerConnexionType::Y_Yg) {
			Ypp = YII / tap_hv / tap_hv;
			Yps = -YII / tap_hv / tap_lv;
			Ysp = -YII / tap_lv / tap_hv;
			Yss = YII / tap_lv / tap_lv;
		} else if (connexionType == TransformerConnexionType::Y_Y) {
			Ypp = YII / tap_hv / tap_hv;
			Yps = -YII / tap_hv / tap_lv;
			Ysp = -YII / tap_lv / tap_hv;
			Yss = YII / tap_lv / tap_lv;
		} else if (connexionType == TransformerConnexionType::Y_D) {
			Ypp = YII / tap_hv / tap_hv;
			Yps = YIII / tap_hv / tap_lv;
			Ysp = YIII.transpose() / tap_lv / tap_hv;
			Yss = YII / tap_lv / tap_lv;
		} else if (connexionType == TransformerConnexionType::D_Yg) {
			Ypp = YII / tap_hv / tap_hv;
			Yps = YIII.transpose() / tap_hv / tap_lv;
			Ysp = YIII / tap_lv / tap_hv;
			Yss = YII / tap_lv / tap_lv;
		} else if (connexionType == TransformerConnexionType::D_Y) {
			Ypp = YII / tap_hv / tap_hv;
			Yps = YIII.transpose() / tap_hv / tap_lv;
			Ysp = YIII / tap_lv / tap_hv;
			Yss = YII / tap_lv / tap_lv;
		} else if (connexionType == TransformerConnexionType::D_D) {
			Ypp = YII / tap_hv / tap_hv;
			Yps = -YII / tap_hv / tap_lv;
			Ysp = -YII / tap_lv / tap_hv;
			Yss = YII / tap_lv / tap_lv;
		}

		cx_mat Y = cx_mat(6, 6);
		Y << Ypp, Yps,
				Ysp, Yss;

		return Y;
	}

	/***************************************************************************
	 * Transformer impedances calculators based on the short circuit test results
	 ***************************************************************************/

	/*
	 */
	Transformer_ShortCircuit_Constructor::Transformer_ShortCircuit_Constructor(std::string name,
			TransformerConnexionType connexionType,
			double HV_nominal_voltage,
			double LV_nominal_voltage,
			double Nominal_power,
			double Copper_losses,
			double Iron_losses,
			double No_load_current,
			double Short_circuit_voltage,
			double GX_HV1, double GR_HV1) {

		Name = name;

		Uhv = HV_nominal_voltage;

		Ulv = LV_nominal_voltage;

		Sn = Nominal_power;

		Pcu = Copper_losses;

		Pfe = Iron_losses;

		I0 = No_load_current;

		Usc = Short_circuit_voltage;

		connType = connexionType;

		calculate_model();
	}
	
	/*
	 */
	Transformer_ShortCircuit_Constructor::~Transformer_ShortCircuit_Constructor(){}

	/*
	 * This function calculates the transformer magnetizing impedance and
	 * leaking impedance
	 */
	void Transformer_ShortCircuit_Constructor::calculate_model() {

		//Nominal impedance HV (Ohm)
		double Zn_hv = Uhv * Uhv / Sn;

		//Nominal impedance LV (Ohm)
		double Zn_lv = Ulv * Ulv / Sn;

		// Short circuit impedance (p.u.)
		double zsc = Usc / 100;

		// Short circuit resistance (p.u.)
		double rsc = (Pcu / 1000) / Sn;

		// Short circuit reactance (p.u.)
		double xsc = sqrt(zsc * zsc - rsc * rsc);

		// HV resistance (p.u.)
		double rcu_hv = rsc * GR_hv1;

		// LV resistance (p.u.)
		double rcu_lv = rsc * (1 - GR_hv1);

		// HV shunt reactance (p.u.)
		double xs_hv = xsc * GX_hv1;

		// LV shunt reactance (p.u.)
		double xs_lv = xsc * (1 - GX_hv1);

		//Shunt resistance (p.u.)
		double rfe = Sn / (Pfe / 1000);

		// Magnetization impedance (p.u.)
		double zm = 1 / (I0 / 100);

		// Magnetization reactance (p.u.)
		double xm;
		if (rfe > zm)
			xm = 1 / sqrt(1 / (zm * zm) - 1 / (rfe * rfe));
		else
			xm = 0; //the sqare root cannot be computed

		//Calculated parameters in per unit
		leakage_impedance = cx_double(rsc, xsc);
		magnetizing_impedance = cx_double(rfe, xm);

		Transformer_3phase_Yabc threePhaseConstructor;
		Y_abc = threePhaseConstructor.Yabc(leakage_impedance, connType, 1.0, 1.0);
	}



}