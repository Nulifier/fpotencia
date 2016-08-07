#include "LineType.h"
#include "Line.h"

namespace fPotencia {
    Line::Line(string name, int connection_bus1, int connection_bus2, LineType line_type, double line_lenght) {
        Name = name;
        bus1 = connection_bus1;
        bus2 = connection_bus2;
        lenght = line_lenght;

        SetType(line_type);

        //cout << Name << "[" << bus1 << ", " << bus2 << "]" << " -> [z:" << impedance << ", y:" << shunt_admittance << "]" << endl;
    }

    /*
     * Line object destructor
     */
    Line::~Line() {
        //type_used = NULL;
        //delete type_used;
    }

    /*
     * This function calculates the line impedance and admittance from a line type
     * model
     */
    void Line::SetType(LineType line_type) {
        shunt_admittance = line_type.shuntAdmittance() * lenght;
        impedance = line_type.impedance() * lenght;
        values_in_per_unit = (line_type.valueType() == fPotencia::pu);

        //create the element admittance matrix
        cx_double y = cx_double(1, 0) / impedance;
        cx_double ys = shunt_admittance / cx_double(2, 0);

        //Fill the internal matrix
        Y_element = cx_mat(2, 2);
        Y_element(0, 0) = y + ys;
        Y_element(0, 1) = -y;
        Y_element(1, 1) = y + ys;
        Y_element(1, 0) = -y;

        //check for inf or nan
        for (int i = 0; i < 2; i++)
            for (int j = 0; j < 2; j++)
                if (std::isinf(Y_element.coeff(i, j).real()) || std::isinf(Y_element.coeff(i, j).imag())) {
                    std::cout << Y_element << std::endl;
                    std::stringstream ss;
                    ss << "Line>>" << Name << ": infinite or nan values in the element Y at: " << i << "," << j;
                    throw std::invalid_argument(ss.str());
                    //std::cout << "Line>>" << Name << ": infinite or nan values in the element Y at: " << i << "," << j << endl;
                }
    }

    /*
     * Returns the component admittance matrix in sparse format, this way
     * the composition of the circuit admittance matrix is straight forward
     */
    void Line::get_element_Y(int n, sp_cx_mat &Yret) {

        //dimension check
        if (bus1 > (n - 1) || bus2 > (n - 1)) {
            std::stringstream ss;
            ss << "Line>>" << Name << ": Wrong Y dimension: " << n;
            throw std::invalid_argument(ss.str());
            //std::cout << "Line>>" << Name << ": Wrong Y dimension: " << n << endl;
            return;
        }



        //set the circuit matrix values
        if (values_in_per_unit) {
            Yret.coeffRef(bus1, bus1) += Y_element.coeff(0, 0);
            Yret.coeffRef(bus1, bus2) += Y_element.coeff(0, 1);
            Yret.coeffRef(bus2, bus2) += Y_element.coeff(1, 1);
            Yret.coeffRef(bus2, bus1) += Y_element.coeff(1, 0);
        } else {
            Yret.coeffRef(bus1, bus1) += Y_element.coeff(0, 0) * Zbase;
            Yret.coeffRef(bus1, bus2) += Y_element.coeff(0, 1) * Zbase;
            Yret.coeffRef(bus2, bus2) += Y_element.coeff(1, 1) * Zbase;
            Yret.coeffRef(bus2, bus1) += Y_element.coeff(1, 0) * Zbase;
        }
    }

    /*
     * This function calculates the amount of current going through the line
     * given a circuit solution
     */
    void Line::calculate_current(cx_solution sol) {
        cx_mat voltage(2, 1);
        cx_mat current(2, 1);
        cx_mat power(2, 1);

        voltage(0, 0) = sol.V[bus1];
        voltage(1, 0) = sol.V[bus2];

        current = Y_element * voltage;

        power(0, 0) = voltage(0, 0) * conj(current(0, 0));
        power(1, 0) = voltage(1, 0) * conj(current(1, 0));

        current_bus1_to_bus2 = current(0, 0);
        current_bus2_to_bus1 = current(1, 0);

        power_bus1_to_bus2 = power(0, 0);
        power_bus2_to_bus1 = power(1, 0);

        power_losses = power_bus1_to_bus2 + power_bus2_to_bus1;
        /*if (power_bus1_to_bus2.real() > power_bus2_to_bus1.real())
                power_losses = power_bus1_to_bus2 + power_bus2_to_bus1;
                else
                power_losses = power_bus2_to_bus1 - power_bus1_to_bus2;*/
    }

    /*
     * This function prints all the line calculated parameters
     */
    void Line::print() {
        cout << Name << endl;
        cout << "\t r:" << impedance.real() << ", x:" << impedance.imag() << ", c: " << shunt_admittance.imag() << endl;
        cout << "\tPower" << endl;
        cout << "\t bus 1 to 2: " << power_bus1_to_bus2 << endl;
        cout << "\t bus 2 to 1: " << power_bus2_to_bus1 << endl;

        cout << "\t Losses: " << power_losses << endl;

        cout << "\tCurrent" << endl;
        cout << "\t bus 1 to 2: " << current_bus1_to_bus2 << endl;
        cout << "\t bus 2 to 1: " << current_bus2_to_bus1 << endl;
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
