#include "Types.h"
#include "LineType.h"
#include "Line.h"


namespace fPotencia {
    void Line::calculateLineParameters()
    {
        valueType_ = lineType_.valueType();
        impedance_ = lineType_.impedance() * length_;
        shuntAdmittance_ = lineType_.shuntAdmittance() * length_;
    }


    void Line::calculateAdmittance()
    {
        cx_double y = cx_double(1, 0) / impedance_;
        cx_double ys = shuntAdmittance_ / cx_double(2, 0);

        //Fill the internal matrix
        admittance_ = cx_mat(2, 2);
        admittance_(0, 0) = y + ys;
        admittance_(0, 1) = -y;
        admittance_(1, 1) = y + ys;
        admittance_(1, 0) = -y;

#ifdef DEBUG
        // Check for illegal values:
        for (cx_mat::Index r = 0; r < admittance_.rows(); ++r) {
            for (cx_mat::Index c = 0; c < admittance_.cols(); ++c) {
                assert(! std::isinf(admittance_(r, c).real()));
                assert(! std::isinf(admittance_(r, c).imag()));
            }
        }
#endif
    }


    void Line::placeAdmittance(sp_cx_mat& systemAdmittance) const
    {
        if (valueType_ == fPotencia::pu) {
            systemAdmittance.coeffRef(bus1Id_, bus1Id_) +=
                    admittance_.coeff(0, 0);
            systemAdmittance.coeffRef(bus1Id_, bus2Id_) +=
                    admittance_.coeff(0, 1);
            systemAdmittance.coeffRef(bus2Id_, bus2Id_) +=
                    admittance_.coeff(1, 1);
            systemAdmittance.coeffRef(bus2Id_, bus1Id_) +=
                    admittance_.coeff(1, 0);
        } else {
            systemAdmittance.coeffRef(bus1Id_, bus1Id_) +=
                    admittance_.coeff(0, 0) * impedanceBase_;
            systemAdmittance.coeffRef(bus1Id_, bus2Id_) +=
                    admittance_.coeff(0, 1) * impedanceBase_;
            systemAdmittance.coeffRef(bus2Id_, bus2Id_) +=
                    admittance_.coeff(1, 1) * impedanceBase_;
            systemAdmittance.coeffRef(bus2Id_, bus1Id_) +=
                    admittance_.coeff(1, 0) * impedanceBase_;
        }
    }


    void Line::setPower(cx_solution const& sol)
    {
        cx_mat voltage(2, 1);
        cx_mat current(2, 1);
        cx_mat power(2, 1);

        voltage(0, 0) = sol.V[bus1Id_];
        voltage(1, 0) = sol.V[bus2Id_];

        current = admittance_ * voltage;

        power(0, 0) = voltage(0, 0) * conj(current(0, 0));
        power(1, 0) = voltage(1, 0) * conj(current(1, 0));

        currentFlowBus1ToBus2_ = current(0, 0);
        currentFlowBus2ToBus1_ = current(1, 0);
        powerFlowBus1ToBus2_ = power(0, 0);
        powerFlowBus2ToBus1_ = power(1, 0);

        powerLoss_ = powerFlowBus1ToBus2_ + powerFlowBus2ToBus1_;
    }


    /*******************************************************************************
     *LineType3 class implementation
     ******************************************************************************/

    /*
     */
    LineType3::LineType3(string name, cx_mat3 Z_abc, cx_mat3 Y_abc) {
        Name = name;
        Zabc = Z_abc;
        Yabc = Y_abc;
    }

    LineType3::~LineType3() {
    }

    /*******************************************************************************
     *Line3 class implementation
     ******************************************************************************/

    /*
     */
    Line3::Line3(string name, int connection_bus1, int connection_bus2, LineType3 line_type, double line_lenght) {
        Name = name;
        bus1 = connection_bus1;
        bus2 = connection_bus2;
        SetType(line_type);
    }

    Line3::~Line3() {
    }

    /*
     */
    void Line3::SetType(LineType3 &line_type) {
        cx_mat3 U;
        U(0.0) = cx_double(1.0, 0.0);
        U(1.1) = cx_double(1.0, 0.0);
        U(2.2) = cx_double(1.0, 0.0);
        cx_mat3 Zabc = line_type.Zabc * lenght;
        cx_mat3 Yabc = line_type.Yabc * lenght;

        a = U + 0.5 * Zabc * Yabc;
        b = Zabc;
        c = Yabc + 0.25 * Yabc * Zabc * Yabc;
        d = U + 0.5 * Zabc * Yabc;

        Eigen::FullPivLU<cx_mat> lu(a);
        A = lu.inverse();
        B = A * b;
    }

}
