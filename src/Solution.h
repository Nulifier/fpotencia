/*
 * File:   Solution.h
 * Author: Santiago Peñate Vera
 *
 * Created on 8 de agosto de 2014, 10:43
 * Copyright (C) 2014 Santiago Peñate Vera
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */
#pragma once

#include "fpotencia_libs.h"

namespace fPotencia {
	/**
	 * Complex voltage and power solution.
	 * This implementation of the solution object contains both voltage and
	 * power in complex mode. This is the solution "shape" that suits best
	 * algorithms like gauss-seidel or jacobi.
	 */
	class cx_solution final {
	public:
		/// Complex power at each bus.
		cx_vec S;

		/// Complex voltage at each bus.
		cx_vec V;

		/// Get the active power.
		[[nodiscard]] vec P() const { return S.real(); }

		/// Get the reactive power.
		[[nodiscard]] vec Q() const { return S.imag(); }

		bool initialized = false;

		[[nodiscard]] unsigned int length() const { return S.size(); }

		/// Resizes the solution object and initializes it to zero.
		void resize(unsigned int n);

		/// Print the solution values.
		void print(const std::string& title) const;
		
		/// Active power of the bus at index k.
		[[nodiscard]] double Pr(uint k) const { return S.coeff(k).real(); }

		/// Reactive power of the bus at index k.
		[[nodiscard]] double Pi(uint k) const { return S.coeff(k).imag(); }

		/// Real voltage of the bus at index k.
		[[nodiscard]] double Vr(uint k) const { return V.coeff(k).real(); }

		/// Imaginary voltage of the bus at index k.
		[[nodiscard]] double Vi(uint k) const { return V.coeff(k).imag(); }
	};

	/**
	 * Complex power and polar voltage solution.
	 * This is the implementation of the solution object with voltage in polar
	 * form. This way of representation of the voltage is usefull for the Newton
	 * -Raphson algorithm, in order to better represent the PV buses constraint.
	 */
	class solution final {
	public:
		vec P;

		vec Q;

		vec Vmag;

		vec Varg;

		/// Complex power at each bus.
		[[nodiscard]] cx_vec S() const { cx_vec T(length()); T.real() = P; T.imag() = Q; return T; }
		
		/// Complex power of the bus at index k.
		[[nodiscard]] cx_double S(uint k) const { return {P.coeff(k), Q.coeff(k)}; }

		/// Complex voltage at each bus.
		[[nodiscard]] cx_vec V() const;

		/// Complex voltage of the bus at index k.
		[[nodiscard]] cx_double V(uint k) const { return std::polar(Vmag.coeff(k), Varg.coeff(k)); }

		bool initialized = false;

		[[nodiscard]] unsigned int length() const { return P.size(); }

		/// Resizes the solution object and initializes it to zero.
		void resize(unsigned int n);

		/// Print the solution values.
		void print(const std::string& title) const;

		/// Real voltage of the bus at index k.
		[[nodiscard]] double Vr(uint k) const { return Vmag.coeff(k) * std::cos(Varg.coeff(k)); }

		/// Imaginary voltage of the bus at index k.
		[[nodiscard]] double Vi(uint k) const { return Vmag.coeff(k) * std::sin(Varg.coeff(k)); }

		/// Convert to rectangular solution.
		[[nodiscard]] cx_solution toComplex() const;
	};
}
