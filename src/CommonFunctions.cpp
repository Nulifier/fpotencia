#include "CommonFunctions.h"

namespace fPotencia {

	/*
	 * Performs the Kron reduction of an impedance matrix
	 * Zabc: Primitive matrix
	 * p: number of poles (that are not neutral)
	 * m: number of neutrals
	 * 
	 * p+m=number of conductors
	 */
	cx_mat Kron_cx(cx_mat Z, uint p, uint m) {
		cx_mat zij = Z.topLeftCorner(p, p);
		cx_mat zin = Z.topRightCorner(p, m);
		cx_mat znj = Z.bottomLeftCorner(m, p);
		cx_mat znn = Z.bottomRightCorner(m, m);

		Eigen::FullPivLU<cx_mat> lu(znn);
		return zij - zin * lu.inverse() * znj;
	}

	/*
	 * Performs the Kron reduction of an impedance matrix
	 * Zabc: Primitive matrix
	 * p: number of poles (that are not neutral)
	 * m: number of neutrals
	 * 
	 * p+m=number of conductors
	 */
	mat Kron(mat Z, uint p, uint m) {
		mat zij = Z.topLeftCorner(p, p);
		mat zin = Z.topRightCorner(p, m);
		mat znj = Z.bottomLeftCorner(m, p);
		mat znn = Z.bottomRightCorner(m, m);

		Eigen::FullPivLU<mat> lu(znn);
		return zij - zin * lu.inverse() * znj;
	}

	/*
	 * Converts the ABC impedances to sequence (012) Impedances
	 */
	cx_mat abc2seq(cx_mat Z) {

		assert(Z.cols() == 3 && Z.rows() == 3);

		double ang = 2.0 * PI / 3.0; //120 deg in radians
		cx_double a(cos(ang), sin(ang));
		cx_double a2 = a*a;
		cx_double one(1., 0.);
		cx_mat mat_(3, 3);
		mat_ << one, one, one,
				one, a2, a,
				one, a, a2;
		cx_mat mat_inv(3, 3);
		mat_inv << one, one, one,
				one, a, a2,
				one, a2, a;
		mat_inv /= 3.0;
		return mat_inv * Z * mat_;
	}

	/*
	 * This function adds a 3x3 matrix M to the three phase admittance matrix Y
	 * at the buses location i, j
	 */
	void addMatToY(cx_mat &Y, cx_mat3 M, int i, int j) {
		for (int a = 0; a < 3; a++)
			for (int b = 0; b < 3; b++) {
				Y.coeffRef(i * 3 + a, j * 3 + b) += M.coeff(a, b);
			}
	}
}
