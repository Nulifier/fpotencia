#include "test-models.h"

fPotencia::Circuit TestModels::ieee14Model() {
	using namespace fPotencia;

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
	LineType ltype1("line type 1", 0.05, 0.11, 0.02, true);
	LineType ltype2("line type 2", 0.03, 0.08, 0.02, true);
	LineType ltype3("line type 3", 0.04, 0.09, 0.02, true);
	LineType ltype4("line type 4", 0.06, 0.13, 0.03, true);

	Line l1("Line 1-2", b1.index, b2.index, ltype1, 1.0);
	Line l2("Line 1-3", b1.index, b3.index, ltype1, 1.0);
	Line l3("Line 1-5", b1.index, b5.index, ltype2, 1.0);
	Line l4("Line 2-3", b2.index, b3.index, ltype3, 1.0);
	Line l5("Line 2-5", b2.index, b5.index, ltype3, 1.0);
	Line l6("Line 3-4", b3.index, b4.index, ltype4, 1.0);
	Line l7("Line 4-5", b4.index, b5.index, ltype3, 1.0);

	model.lines.push_back(l1);
	model.lines.push_back(l2);
	model.lines.push_back(l3);
	model.lines.push_back(l4);
	model.lines.push_back(l5);
	model.lines.push_back(l6);
	model.lines.push_back(l7);

	// Generators:
	Generator g2("Gen2", b2.index, 21.7, 12.7, -40.0, 50.0, false);
	Generator g3("Gen3", b3.index, 94.2, 19.0, 0.0, 40.0, false);
	Generator g4("Gen4", b4.index, 47.8, -3.9, 0.0, 0.0, false);
	Generator g5("Gen5", b5.index, 7.6, 1.6, 0.0, 0.0, false);
	Generator g6("Gen6", b6.index, 11.2, 7.5, -6.0, 24.0, false);
	Generator g7("Gen7", b7.index, 0.0, 0.0, 0.0, 0.0, false);
	Generator g8("Gen8", b8.index, 0.0, 0.0, -6.0, 24.0, false);
	Generator g9("Gen9", b9.index, 29.5, 16.6, 0.0, 0.0, false);
	Generator g10("Gen10", b10.index, 9.0, 5.8, 0.0, 0.0, false);
	Generator g11("Gen11", b11.index, 3.5, 1.8, 0.0, 0.0, false);
	Generator g12("Gen12", b12.index, 6.1, 1.6, 0.0, 0.0, false);
	Generator g13("Gen13", b13.index, 13.5, 5.8, 0.0, 0.0, false);
	Generator g14("Gen14", b14.index, 14.9, 5.0, 0.0, 0.0, false);

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
	Load ld1("Load1", b1.index, 0.0, 0.0);
	Load ld2("Load2", b2.index, 21.7, 12.7);
	Load ld3("Load3", b3.index, 94.2, 19.0);
	Load ld4("Load4", b4.index, 47.8, -3.9);
	Load ld5("Load5", b5.index, 7.6, 1.6);
	Load ld6("Load6", b6.index, 11.2, 7.5);
	Load ld7("Load7", b7.index, 0.0, 0.0);
	Load ld8("Load8", b8.index, 0.0, 0.0);
	Load ld9("Load9", b9.index, 29.5, 16.6);
	Load ld10("Load10", b10.index, 9.0, 5.8);
	Load ld11("Load11", b11.index, 3.5, 1.8);
	Load ld12("Load12", b12.index, 6.1, 1.6);
	Load ld13("Load13", b13.index, 13.5, 5.8);
	Load ld14("Load14", b14.index, 14.9, 5.0);

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

fPotencia::Circuit TestModels::lynnPowerllWithGenerator() {
	using namespace fPotencia;

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
	LineType ltype1("line type 1", 0.05, 0.11, 0.02, true);
	LineType ltype2("line type 2", 0.03, 0.08, 0.02, true);
	LineType ltype3("line type 3", 0.04, 0.09, 0.02, true);
	LineType ltype4("line type 4", 0.06, 0.13, 0.03, true);

	Line l1("Line 1-2", b1.index, b2.index, ltype1, 1.0);
	Line l2("Line 1-3", b1.index, b3.index, ltype1, 1.0);
	Line l3("Line 1-5", b1.index, b5.index, ltype2, 1.0);
	Line l4("Line 2-3", b2.index, b3.index, ltype3, 1.0);
	Line l5("Line 2-5", b2.index, b5.index, ltype3, 1.0);
	Line l6("Line 3-4", b3.index, b4.index, ltype4, 1.0);
	Line l7("Line 4-5", b4.index, b5.index, ltype3, 1.0);
	model.lines.push_back(l1);
	model.lines.push_back(l2);
	model.lines.push_back(l3);
	model.lines.push_back(l4);
	model.lines.push_back(l5);
	model.lines.push_back(l6);
	model.lines.push_back(l7);

	// Loads:
	Load ld2("Load1", b2.index, 40, 20);
	Load ld3("Load2", b3.index, 25, 15);
	Load ld5("Load4", b5.index, 50, 20);
	model.loads.push_back(ld2);
	model.loads.push_back(ld3);
	model.loads.push_back(ld5);

	// Generators:
	Generator g1("Generator", b4.index, 30, 1.0, -15, 20, true);
	model.generators.push_back(g1);

	return model;
}

fPotencia::Circuit TestModels::lynnPowerllWithoutGenerator() {
	using namespace fPotencia;

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
	LineType ltype1("line type 1", 0.05, 0.11, 0.02, true);
	LineType ltype2("line type 2", 0.03, 0.08, 0.02, true);
	LineType ltype3("line type 3", 0.04, 0.09, 0.02, true);
	LineType ltype4("line type 4", 0.06, 0.13, 0.03, true);

	Line l1("Line 1-2", b1.index, b2.index, ltype1, 1.0);
	Line l2("Line 1-3", b1.index, b3.index, ltype1, 1.0);
	Line l3("Line 1-5", b1.index, b5.index, ltype2, 1.0);
	Line l4("Line 2-3", b2.index, b3.index, ltype3, 1.0);
	Line l5("Line 2-5", b2.index, b5.index, ltype3, 1.0);
	Line l6("Line 3-4", b3.index, b4.index, ltype4, 1.0);
	Line l7("Line 4-5", b4.index, b5.index, ltype3, 1.0);
	model.lines.push_back(l1);
	model.lines.push_back(l2);
	model.lines.push_back(l3);
	model.lines.push_back(l4);
	model.lines.push_back(l5);
	model.lines.push_back(l6);
	model.lines.push_back(l7);

	// Loads:
	Load ld2("Load1", b2.index, 40, 20);
	Load ld3("Load2", b3.index, 25, 15);
	Load ld4("Load3", b4.index, 40, 20);
	Load ld5("Load4", b5.index, 50, 20);
	model.loads.push_back(ld2);
	model.loads.push_back(ld3);
	model.loads.push_back(ld4);
	model.loads.push_back(ld5);

	return model;
}
