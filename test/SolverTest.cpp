#include "fpotencia.h"
#include "SolverTest.h"


using namespace fPotencia;


SolverTest::SolverTest()
{
}


SolverTest::~SolverTest() noexcept
{
}


fPotencia::Circuit SolverTest::generateIeee14Bus() const
{
    Circuit model;


    // Buses:

    Bus b1("Bus1", undefined_bus_type, 10.0);
    Bus b2("Bus2", undefined_bus_type, 10.0);
    Bus b3("Bus3", undefined_bus_type, 10.0);
    Bus b4("Bus4", undefined_bus_type, 10.0);
    Bus b5("Bus5", undefined_bus_type, 10.0);
    Bus b6("Bus6", undefined_bus_type, 10.0);
    Bus b7("Bus7", undefined_bus_type, 10.0);
    Bus b8("Bus8", undefined_bus_type, 10.0);
    Bus b9("Bus9", undefined_bus_type, 10.0);
    Bus b10("Bus10", undefined_bus_type, 10.0);
    Bus b11("Bus11", undefined_bus_type, 10.0);
    Bus b12("Bus12", undefined_bus_type, 10.0);
    Bus b13("Bus13", undefined_bus_type, 10.0);
    Bus b14("Bus14", undefined_bus_type, 10.0);

    model.add_Bus(b1);
    model.add_Bus(b2);
    model.add_Bus(b3);
    model.add_Bus(b4);
    model.add_Bus(b5);
    model.add_Bus(b6);
    model.add_Bus(b7);
    model.add_Bus(b8);
    model.add_Bus(b9);
    model.add_Bus(b10);
    model.add_Bus(b11);
    model.add_Bus(b12);
    model.add_Bus(b13);
    model.add_Bus(b14);


    // External Grid:

    ExternalGrid eg("External1", b1.index);
    model.externalGrids.push_back(eg);


    // Line Types and Lines:

    LineType ltype1, ltype2, ltype3, ltype4;
    ltype1
            .resistance(0.05)
            .reactance(0.11)
            .susceptance(0.02)
            .valueType(fPotencia::pu);
    ltype2
            .resistance(0.03)
            .reactance(0.08)
            .susceptance(0.02)
            .valueType(fPotencia::pu);
    ltype3
            .resistance(0.04)
            .reactance(0.09)
            .susceptance(0.02)
            .valueType(fPotencia::pu);
    ltype4
            .resistance(0.06)
            .reactance(0.13)
            .susceptance(0.03)
            .valueType(fPotencia::pu);

    Line l1, l2, l3, l4, l5, l6, l7;
    l1.buses(b1.index, b2.index)
            .lineType(ltype1)
            .length(1.0);
    l2.buses(b1.index, b3.index)
            .lineType(ltype2)
            .length(1.0);
    l3.buses(b1.index, b5.index)
            .lineType(ltype3)
            .length(1.0);
    l4.buses(b2.index, b3.index)
            .lineType(ltype3)
            .length(1.0);
    l5.buses(b2.index, b5.index)
            .lineType(ltype3)
            .length(1.0);
    l6.buses(b3.index, b4.index)
            .lineType(ltype4)
            .length(1.0);
    l7.buses(b4.index, b5.index)
            .lineType(ltype3)
            .length(1.0);

    model.lines.push_back(l1);
    model.lines.push_back(l2);
    model.lines.push_back(l3);
    model.lines.push_back(l4);
    model.lines.push_back(l5);
    model.lines.push_back(l6);
    model.lines.push_back(l7);


    // Generators:

    Generator g2, g3, g4, g5, g6, g7, g8, g9, g10, g11, g12, g13, g14;
    g2
            .bus(b2)
            .realPower(21.7)
            .voltage(12.7, Generator::volts)
            .reactivePowerLimits(-40, 50);
    g3
            .bus(b3)
            .realPower(94.2)
            .voltage(19.0, Generator::volts)
            .reactivePowerLimits(0, -40);
    g4
            .bus(b4)
            .realPower(47.8)
            .voltage(-3.9, Generator::volts)
            .reactivePowerLimits(0, 0);
    g5
            .bus(b5)
            .realPower(7.6)
            .voltage(1.6, Generator::volts)
            .reactivePowerLimits(0, 0);
    g6
            .bus(b6)
            .realPower(11.2)
            .voltage(7.5, Generator::volts)
            .reactivePowerLimits(-6.0, 24.0);
    g7
            .bus(b7)
            .realPower(0.0)
            .voltage(0.0, Generator::volts)
            .reactivePowerLimits(0.0, 0.0);
    g8
            .bus(b8)
            .realPower(0.0)
            .voltage(0.0, Generator::volts)
            .reactivePowerLimits(-6.0, 24.0);
    g9
            .bus(b9)
            .realPower(29.5)
            .voltage(16.6, Generator::volts)
            .reactivePowerLimits(0.0, 0.0);
    g10
            .bus(b10)
            .realPower(9.0)
            .voltage(5.8, Generator::volts)
            .reactivePowerLimits(0.0, 0.0);
    g11
            .bus(b11)
            .realPower(3.5)
            .voltage(1.8, Generator::volts)
            .reactivePowerLimits(0.0, 0.0);
    g12
            .bus(b12)
            .realPower(6.1)
            .voltage(1.6, Generator::volts)
            .reactivePowerLimits(0.0, 0.0);
    g13
            .bus(b13)
            .realPower(13.5)
            .voltage(5.8, Generator::volts)
            .reactivePowerLimits(0.0, 0.0);
    g14
            .bus(b14)
            .realPower(14.9)
            .voltage(5.0, Generator::volts)
            .reactivePowerLimits(0.0, 0.0);

    model.generators.push_back(g2);
    model.generators.push_back(g3);
    model.generators.push_back(g4);
    model.generators.push_back(g5);
    model.generators.push_back(g6);
    model.generators.push_back(g7);
    model.generators.push_back(g8);
    model.generators.push_back(g9);
    model.generators.push_back(g10);
    model.generators.push_back(g11);
    model.generators.push_back(g12);
    model.generators.push_back(g13);
    model.generators.push_back(g14);


    // Loads:

    Load ld1, ld2, ld3, ld4, ld5, ld6, ld7, ld8, ld9, ld10, ld11, ld12, ld13,
            ld14;
    ld1
            .bus(b1)
            .name("Load 1")
            .realPower(0)
            .reactivePower(0);
    ld2
            .bus(b2)
            .name("Load 2")
            .realPower(21.7)
            .reactivePower(12.7);
    ld3
            .bus(b3)
            .name("Load 3")
            .realPower(94.2)
            .reactivePower(19.0);
    ld4
            .bus(b4)
            .name("Load 4")
            .realPower(47.8)
            .reactivePower(-3.9);
    ld5
            .bus(b5)
            .name("Load 5")
            .realPower(7.6)
            .reactivePower(1.6);
    ld6
            .bus(b6)
            .name("Load 6")
            .realPower(11.2)
            .reactivePower(7.5);
    ld7
            .bus(b7)
            .name("Load 7")
            .realPower(0.0)
            .reactivePower(0.0);
    ld8
            .bus(b8)
            .name("Load 8")
            .realPower(0.0)
            .reactivePower(0.0);
    ld9
            .bus(b9)
            .name("Load 9")
            .realPower(29.5)
            .reactivePower(16.6);
    ld10
            .bus(b10)
            .name("Load 10")
            .realPower(9.0)
            .reactivePower(5.8);
    ld11
            .bus(b11)
            .name("Load 11")
            .realPower(3.5)
            .reactivePower(1.8);
    ld12
            .bus(b12)
            .name("Load 12")
            .realPower(6.1)
            .reactivePower(1.6);
    ld13
            .bus(b13)
            .name("Load 13")
            .realPower(13.5)
            .reactivePower(5.8);
    ld14
            .bus(b14)
            .name("Load 14")
            .realPower(14.9)
            .reactivePower(5.0);

    model.loads.push_back(ld1);
    model.loads.push_back(ld2);
    model.loads.push_back(ld3);
    model.loads.push_back(ld4);
    model.loads.push_back(ld5);
    model.loads.push_back(ld6);
    model.loads.push_back(ld7);
    model.loads.push_back(ld8);
    model.loads.push_back(ld9);
    model.loads.push_back(ld10);
    model.loads.push_back(ld11);
    model.loads.push_back(ld12);
    model.loads.push_back(ld13);
    model.loads.push_back(ld14);


    return model;
}


Circuit SolverTest::generateLynnPowellWithGenerator() const
{
    Circuit model;


    // Buses:

    Bus b1("bus1", undefined_bus_type, 132.0);
    Bus b2("bus2", undefined_bus_type, 132.0);
    Bus b3("bus3", undefined_bus_type, 132.0);
    Bus b4("bus4", undefined_bus_type, 132.0);
    Bus b5("bus5", undefined_bus_type, 132.0);

    model.add_Bus(b1);
    model.add_Bus(b2);
    model.add_Bus(b3);
    model.add_Bus(b4);
    model.add_Bus(b5);


    // External Grids:

    ExternalGrid eg("External1", b1.index);
    model.externalGrids.push_back(eg);


    // Line types and lines:

    LineType ltype1, ltype2, ltype3, ltype4;
    ltype1
            .resistance(0.05)
            .reactance(0.11)
            .susceptance(0.02)
            .valueType(fPotencia::pu);
    ltype2
            .resistance(0.03)
            .reactance(0.08)
            .susceptance(0.02)
            .valueType(fPotencia::pu);
    ltype3
            .resistance(0.04)
            .reactance(0.09)
            .susceptance(0.02)
            .valueType(fPotencia::pu);
    ltype4
            .resistance(0.06)
            .reactance(0.13)
            .susceptance(0.03)
            .valueType(fPotencia::pu);

    Line l1, l2, l3, l4, l5, l6, l7;
    l1
            .buses(b1.index, b2.index)
            .lineType(ltype1)
            .length(1.0);
    l2
            .buses(b1.index, b3.index)
            .lineType(ltype1)
            .length(1.0);
    l3
            .buses(b1.index, b5.index)
            .lineType(ltype2)
            .length(1.0);
    l4
            .buses(b2.index, b3.index)
            .lineType(ltype3)
            .length(1.0);
    l5
            .buses(b2.index, b5.index)
            .lineType(ltype3)
            .length(1.0);
    l6
            .buses(b3.index, b4.index)
            .lineType(ltype4)
            .length(1.0);
    l7
            .buses(b4.index, b5.index)
            .lineType(ltype3)
            .length(1.0);
    model.lines.push_back(l1);
    model.lines.push_back(l2);
    model.lines.push_back(l3);
    model.lines.push_back(l4);
    model.lines.push_back(l5);
    model.lines.push_back(l6);
    model.lines.push_back(l7);


    // Loads:

    Load ld2, ld3, ld5;
    ld2.bus(b2).name("Load 1").realPower(40).reactivePower(20);
    ld3.bus(b3).name("Load 2").realPower(25).reactivePower(15);
    ld5.bus(b5).name("Load 4").realPower(50).reactivePower(20);
    model.loads.push_back(ld2);
    model.loads.push_back(ld3);
    model.loads.push_back(ld5);


    // Generators:

    Generator g1;
    g1
            .name("Generator 1")
            .bus(b4)
            .realPower(30)
            .voltage(1.0, Generator::pu)
            .reactivePowerLimits(-15, 20);
    model.generators.push_back(g1);


    return model;
}



fPotencia::Circuit SolverTest::generateLynnPowellWithoutGenerator() const
{
    Circuit model;


    // Buses:

    Bus b1("bus1", undefined_bus_type, 132.0);
    Bus b2("bus2", undefined_bus_type, 132.0);
    Bus b3("bus3", undefined_bus_type, 132.0);
    Bus b4("bus4", undefined_bus_type, 132.0);
    Bus b5("bus5", undefined_bus_type, 132.0);

    model.add_Bus(b1);
    model.add_Bus(b2);
    model.add_Bus(b3);
    model.add_Bus(b4);
    model.add_Bus(b5);


    // External Grids:

    ExternalGrid eg("External1", b1.index);
    model.externalGrids.push_back(eg);


    // Line types and lines:

    LineType ltype1, ltype2, ltype3, ltype4;
    ltype1
            .resistance(0.05)
            .reactance(0.11)
            .susceptance(0.02)
            .valueType(fPotencia::pu);
    ltype2
            .resistance(0.03)
            .reactance(0.08)
            .susceptance(0.02)
            .valueType(fPotencia::pu);
    ltype3
            .resistance(0.04)
            .reactance(0.09)
            .susceptance(0.02)
            .valueType(fPotencia::pu);
    ltype4
            .resistance(0.06)
            .reactance(0.13)
            .susceptance(0.03)
            .valueType(fPotencia::pu);

    Line l1, l2, l3, l4, l5, l6, l7;
    l1
            .buses(b1.index, b2.index)
            .lineType(ltype1)
            .length(1.0);
    l2
            .buses(b1.index, b3.index)
            .lineType(ltype1)
            .length(1.0);
    l3
            .buses(b1.index, b5.index)
            .lineType(ltype2)
            .length(1.0);
    l4
            .buses(b2.index, b3.index)
            .lineType(ltype3)
            .length(1.0);
    l5
            .buses(b2.index, b5.index)
            .lineType(ltype3)
            .length(1.0);
    l6
            .buses(b3.index, b4.index)
            .lineType(ltype4)
            .length(1.0);
    l7
            .buses(b4.index, b5.index)
            .lineType(ltype3)
            .length(1.0);    model.lines.push_back(l1);
    model.lines.push_back(l2);
    model.lines.push_back(l3);
    model.lines.push_back(l4);
    model.lines.push_back(l5);
    model.lines.push_back(l6);
    model.lines.push_back(l7);

    // Loads:

    Load ld2, ld3, ld4, ld5;
    ld2.bus(b2).name("Load 1").realPower(40).reactivePower(20);
    ld3.bus(b3).name("Load 2").realPower(25).reactivePower(15);
    ld4.bus(b4).name("Load 3").realPower(40).reactivePower(20);
    ld5.bus(b5).name("Load 4").realPower(50).reactivePower(20);
    model.loads.push_back(ld2);
    model.loads.push_back(ld3);
    model.loads.push_back(ld4);
    model.loads.push_back(ld5);


    return model;
}
