#include "../fpotencia.h"

using namespace fPotencia;

Circuit ieee_30_bus() {

	Circuit model("circuit 1");

	//buses creation
	Bus b0("bus0", VD, 135.0);
	model.add_Bus(b0);
	Bus b1("bus1", PV, 135.0);
	model.add_Bus(b1);
	Bus b2("bus2", PQ, 135.0);
	model.add_Bus(b2);
	Bus b3("bus3", PQ, 135.0);
	model.add_Bus(b3);
	Bus b4("bus4", PQ, 135.0);
	model.add_Bus(b4);
	Bus b5("bus5", PQ, 135.0);
	model.add_Bus(b5);
	Bus b6("bus6", PQ, 135.0);
	model.add_Bus(b6);
	Bus b7("bus7", PQ, 135.0);
	model.add_Bus(b7);
	Bus b8("bus8", PQ, 135.0);
	model.add_Bus(b8);
	Bus b9("bus9", PQ, 135.0);
	model.add_Bus(b9);
	Bus b10("bus10", PQ, 135.0);
	model.add_Bus(b10);
	Bus b11("bus11", PQ, 135.0);
	model.add_Bus(b11);
	Bus b12("bus12", PV, 135.0);
	model.add_Bus(b12);
	Bus b13("bus13", PQ, 135.0);
	model.add_Bus(b13);
	Bus b14("bus14", PQ, 135.0);
	model.add_Bus(b14);
	Bus b15("bus15", PQ, 135.0);
	model.add_Bus(b15);
	Bus b16("bus16", PQ, 135.0);
	model.add_Bus(b16);
	Bus b17("bus17", PQ, 135.0);
	model.add_Bus(b17);
	Bus b18("bus18", PQ, 135.0);
	model.add_Bus(b18);
	Bus b19("bus19", PQ, 135.0);
	model.add_Bus(b19);
	Bus b20("bus20", PQ, 135.0);
	model.add_Bus(b20);
	Bus b21("bus21", PV, 135.0);
	model.add_Bus(b21);
	Bus b22("bus22", PV, 135.0);
	model.add_Bus(b22);
	Bus b23("bus23", PQ, 135.0);
	model.add_Bus(b23);
	Bus b24("bus24", PQ, 135.0);
	model.add_Bus(b24);
	Bus b25("bus25", PQ, 135.0);
	model.add_Bus(b25);
	Bus b26("bus26", PV, 135.0);
	model.add_Bus(b26);
	Bus b27("bus27", PQ, 135.0);
	model.add_Bus(b27);
	Bus b28("bus28", PQ, 135.0);
	model.add_Bus(b28);
	Bus b29("bus29", PQ, 135.0);
	model.add_Bus(b29);

	//load creation
	Shunt sh4("Shunt4", 4, 0.0, -5.26315789474);
	model.shunts.push_back(sh4);
	Shunt sh23("Shunt23", 23, 0.0, -25.0);
	model.shunts.push_back(sh23);

	//load creation
	Load ld1("Load1", 1, 21.7, 12.7);
	model.loads.push_back(ld1);
	Load ld2("Load2", 2, 2.4, 1.2);
	model.loads.push_back(ld2);
	Load ld3("Load3", 3, 7.6, 1.6);
	model.loads.push_back(ld3);
	Load ld6("Load6", 6, 22.8, 10.9);
	model.loads.push_back(ld6);
	Load ld7("Load7", 7, 30.0, 30.0);
	model.loads.push_back(ld7);
	Load ld9("Load9", 9, 5.8, 2.0);
	model.loads.push_back(ld9);
	Load ld11("Load11", 11, 11.2, 7.5);
	model.loads.push_back(ld11);
	Load ld13("Load13", 13, 6.2, 1.6);
	model.loads.push_back(ld13);
	Load ld14("Load14", 14, 8.2, 2.5);
	model.loads.push_back(ld14);
	Load ld15("Load15", 15, 3.5, 1.8);
	model.loads.push_back(ld15);
	Load ld16("Load16", 16, 9.0, 5.8);
	model.loads.push_back(ld16);
	Load ld17("Load17", 17, 3.2, 0.9);
	model.loads.push_back(ld17);
	Load ld18("Load18", 18, 9.5, 3.4);
	model.loads.push_back(ld18);
	Load ld19("Load19", 19, 2.2, 0.7);
	model.loads.push_back(ld19);
	Load ld20("Load20", 20, 17.5, 11.2);
	model.loads.push_back(ld20);
	Load ld22("Load22", 22, 3.2, 1.6);
	model.loads.push_back(ld22);
	Load ld23("Load23", 23, 8.7, 6.7);
	model.loads.push_back(ld23);
	Load ld25("Load25", 25, 3.5, 2.3);
	model.loads.push_back(ld25);
	Load ld28("Load28", 28, 2.4, 0.9);
	model.loads.push_back(ld28);
	Load ld29("Load29", 29, 10.6, 1.9);
	model.loads.push_back(ld29);


	//Generation creation
	Generator g0("Generator0", 0, 23.54, 1.0, -20.0, 150.0, true);
	model.generators.push_back(g0);
	Generator g1("Generator1", 1, 60.97, 1.0, -20.0, 60.0, true);
	model.generators.push_back(g1);
	Generator g2("Generator2", 21, 21.59, 1.0, -15.0, 62.5, true);
	model.generators.push_back(g2);
	Generator g3("Generator3", 26, 26.91, 1.0, -15.0, 48.7, true);
	model.generators.push_back(g3);
	Generator g4("Generator4", 22, 19.2, 1.0, -10.0, 40.0, true);
	model.generators.push_back(g4);
	Generator g5("Generator5", 12, 37.0, 1.0, -15.0, 44.7, true);
	model.generators.push_back(g5);


	//Branch creation
	TransformerType TType0("Transformer_type0", cx_double(0.02, 0.06), cx_double(0.0, 33.3333333333));
	Transformer T0("Transformer0", 0, 1, TType0);
	model.transformers.push_back(T0);
	TransformerType TType1("Transformer_type1", cx_double(0.05, 0.19), cx_double(0.0, 50.0));
	Transformer T1("Transformer1", 0, 2, TType1);
	model.transformers.push_back(T1);
	TransformerType TType2("Transformer_type2", cx_double(0.06, 0.17), cx_double(0.0, 50.0));
	Transformer T2("Transformer2", 1, 3, TType2);
	model.transformers.push_back(T2);
	LineType LType0("Line_type0", 0.01, 0.04, 0.0, true);
	Line L0("Line0", 2, 3, LType0, 1.0);
	model.lines.push_back(L0);
	TransformerType TType3("Transformer_type3", cx_double(0.05, 0.2), cx_double(0.0, 50.0));
	Transformer T3("Transformer3", 1, 4, TType3);
	model.transformers.push_back(T3);
	TransformerType TType4("Transformer_type4", cx_double(0.06, 0.18), cx_double(0.0, 50.0));
	Transformer T4("Transformer4", 1, 5, TType4);
	model.transformers.push_back(T4);
	LineType LType1("Line_type1", 0.01, 0.04, 0.0, true);
	Line L1("Line1", 3, 5, LType1, 1.0);
	model.lines.push_back(L1);
	TransformerType TType5("Transformer_type5", cx_double(0.05, 0.12), cx_double(0.0, 100.0));
	Transformer T5("Transformer5", 4, 6, TType5);
	model.transformers.push_back(T5);
	TransformerType TType6("Transformer_type6", cx_double(0.03, 0.08), cx_double(0.0, 100.0));
	Transformer T6("Transformer6", 5, 6, TType6);
	model.transformers.push_back(T6);
	LineType LType2("Line_type2", 0.01, 0.04, 0.0, true);
	Line L2("Line2", 5, 7, LType2, 1.0);
	model.lines.push_back(L2);
	LineType LType3("Line_type3", 0.0, 0.21, 0.0, true);
	Line L3("Line3", 5, 8, LType3, 1.0);
	model.lines.push_back(L3);
	LineType LType4("Line_type4", 0.0, 0.56, 0.0, true);
	Line L4("Line4", 5, 9, LType4, 1.0);
	model.lines.push_back(L4);
	LineType LType5("Line_type5", 0.0, 0.21, 0.0, true);
	Line L5("Line5", 8, 10, LType5, 1.0);
	model.lines.push_back(L5);
	LineType LType6("Line_type6", 0.0, 0.11, 0.0, true);
	Line L6("Line6", 8, 9, LType6, 1.0);
	model.lines.push_back(L6);
	LineType LType7("Line_type7", 0.0, 0.26, 0.0, true);
	Line L7("Line7", 3, 11, LType7, 1.0);
	model.lines.push_back(L7);
	LineType LType8("Line_type8", 0.0, 0.14, 0.0, true);
	Line L8("Line8", 11, 12, LType8, 1.0);
	model.lines.push_back(L8);
	LineType LType9("Line_type9", 0.12, 0.26, 0.0, true);
	Line L9("Line9", 11, 13, LType9, 1.0);
	model.lines.push_back(L9);
	LineType LType10("Line_type10", 0.07, 0.13, 0.0, true);
	Line L10("Line10", 11, 14, LType10, 1.0);
	model.lines.push_back(L10);
	LineType LType11("Line_type11", 0.09, 0.2, 0.0, true);
	Line L11("Line11", 11, 15, LType11, 1.0);
	model.lines.push_back(L11);
	LineType LType12("Line_type12", 0.22, 0.2, 0.0, true);
	Line L12("Line12", 13, 14, LType12, 1.0);
	model.lines.push_back(L12);
	LineType LType13("Line_type13", 0.08, 0.19, 0.0, true);
	Line L13("Line13", 15, 16, LType13, 1.0);
	model.lines.push_back(L13);
	LineType LType14("Line_type14", 0.11, 0.22, 0.0, true);
	Line L14("Line14", 14, 17, LType14, 1.0);
	model.lines.push_back(L14);
	LineType LType15("Line_type15", 0.06, 0.13, 0.0, true);
	Line L15("Line15", 17, 18, LType15, 1.0);
	model.lines.push_back(L15);
	LineType LType16("Line_type16", 0.03, 0.07, 0.0, true);
	Line L16("Line16", 18, 19, LType16, 1.0);
	model.lines.push_back(L16);
	LineType LType17("Line_type17", 0.09, 0.21, 0.0, true);
	Line L17("Line17", 9, 19, LType17, 1.0);
	model.lines.push_back(L17);
	LineType LType18("Line_type18", 0.03, 0.08, 0.0, true);
	Line L18("Line18", 9, 16, LType18, 1.0);
	model.lines.push_back(L18);
	LineType LType19("Line_type19", 0.03, 0.07, 0.0, true);
	Line L19("Line19", 9, 20, LType19, 1.0);
	model.lines.push_back(L19);
	LineType LType20("Line_type20", 0.07, 0.15, 0.0, true);
	Line L20("Line20", 9, 21, LType20, 1.0);
	model.lines.push_back(L20);
	LineType LType21("Line_type21", 0.01, 0.02, 0.0, true);
	Line L21("Line21", 20, 21, LType21, 1.0);
	model.lines.push_back(L21);
	LineType LType22("Line_type22", 0.1, 0.2, 0.0, true);
	Line L22("Line22", 14, 22, LType22, 1.0);
	model.lines.push_back(L22);
	LineType LType23("Line_type23", 0.12, 0.18, 0.0, true);
	Line L23("Line23", 21, 23, LType23, 1.0);
	model.lines.push_back(L23);
	LineType LType24("Line_type24", 0.13, 0.27, 0.0, true);
	Line L24("Line24", 22, 23, LType24, 1.0);
	model.lines.push_back(L24);
	LineType LType25("Line_type25", 0.19, 0.33, 0.0, true);
	Line L25("Line25", 23, 24, LType25, 1.0);
	model.lines.push_back(L25);
	LineType LType26("Line_type26", 0.25, 0.38, 0.0, true);
	Line L26("Line26", 24, 25, LType26, 1.0);
	model.lines.push_back(L26);
	LineType LType27("Line_type27", 0.11, 0.21, 0.0, true);
	Line L27("Line27", 24, 26, LType27, 1.0);
	model.lines.push_back(L27);
	LineType LType28("Line_type28", 0.0, 0.4, 0.0, true);
	Line L28("Line28", 27, 26, LType28, 1.0);
	model.lines.push_back(L28);
	LineType LType29("Line_type29", 0.22, 0.42, 0.0, true);
	Line L29("Line29", 26, 28, LType29, 1.0);
	model.lines.push_back(L29);
	LineType LType30("Line_type30", 0.32, 0.6, 0.0, true);
	Line L30("Line30", 26, 29, LType30, 1.0);
	model.lines.push_back(L30);
	LineType LType31("Line_type31", 0.24, 0.45, 0.0, true);
	Line L31("Line31", 28, 29, LType31, 1.0);
	model.lines.push_back(L31);
	TransformerType TType7("Transformer_type7", cx_double(0.06, 0.2), cx_double(0.0, 50.0));
	Transformer T7("Transformer7", 7, 27, TType7);
	model.transformers.push_back(T7);
	TransformerType TType8("Transformer_type8", cx_double(0.02, 0.06), cx_double(0.0, 100.0));
	Transformer T8("Transformer8", 5, 27, TType8);
	model.transformers.push_back(T8);

	return model;

}

Circuit ieee_30_bus_noPV() {

	Circuit model("circuit 1");

	//buses creation
	Bus b0("bus0", VD, 135.0);
	model.add_Bus(b0);
	Bus b1("bus1", PQ, 135.0);
	model.add_Bus(b1);
	Bus b2("bus2", PQ, 135.0);
	model.add_Bus(b2);
	Bus b3("bus3", PQ, 135.0);
	model.add_Bus(b3);
	Bus b4("bus4", PQ, 135.0);
	model.add_Bus(b4);
	Bus b5("bus5", PQ, 135.0);
	model.add_Bus(b5);
	Bus b6("bus6", PQ, 135.0);
	model.add_Bus(b6);
	Bus b7("bus7", PQ, 135.0);
	model.add_Bus(b7);
	Bus b8("bus8", PQ, 135.0);
	model.add_Bus(b8);
	Bus b9("bus9", PQ, 135.0);
	model.add_Bus(b9);
	Bus b10("bus10", PQ, 135.0);
	model.add_Bus(b10);
	Bus b11("bus11", PQ, 135.0);
	model.add_Bus(b11);
	Bus b12("bus12", PQ, 135.0);
	model.add_Bus(b12);
	Bus b13("bus13", PQ, 135.0);
	model.add_Bus(b13);
	Bus b14("bus14", PQ, 135.0);
	model.add_Bus(b14);
	Bus b15("bus15", PQ, 135.0);
	model.add_Bus(b15);
	Bus b16("bus16", PQ, 135.0);
	model.add_Bus(b16);
	Bus b17("bus17", PQ, 135.0);
	model.add_Bus(b17);
	Bus b18("bus18", PQ, 135.0);
	model.add_Bus(b18);
	Bus b19("bus19", PQ, 135.0);
	model.add_Bus(b19);
	Bus b20("bus20", PQ, 135.0);
	model.add_Bus(b20);
	Bus b21("bus21", PQ, 135.0);
	model.add_Bus(b21);
	Bus b22("bus22", PQ, 135.0);
	model.add_Bus(b22);
	Bus b23("bus23", PQ, 135.0);
	model.add_Bus(b23);
	Bus b24("bus24", PQ, 135.0);
	model.add_Bus(b24);
	Bus b25("bus25", PQ, 135.0);
	model.add_Bus(b25);
	Bus b26("bus26", PQ, 135.0);
	model.add_Bus(b26);
	Bus b27("bus27", PQ, 135.0);
	model.add_Bus(b27);
	Bus b28("bus28", PQ, 135.0);
	model.add_Bus(b28);
	Bus b29("bus29", PQ, 135.0);
	model.add_Bus(b29);

	//load creation
	Shunt sh4("Shunt4", 4, 0.0, -5.26315789474);
	model.shunts.push_back(sh4);
	Shunt sh23("Shunt23", 23, 0.0, -25.0);
	model.shunts.push_back(sh23);

	//load creation
	Load ld1("Load1", 1, 21.7, 12.7);
	model.loads.push_back(ld1);
	Load ld2("Load2", 2, 2.4, 1.2);
	model.loads.push_back(ld2);
	Load ld3("Load3", 3, 7.6, 1.6);
	model.loads.push_back(ld3);
	Load ld6("Load6", 6, 22.8, 10.9);
	model.loads.push_back(ld6);
	Load ld7("Load7", 7, 30.0, 30.0);
	model.loads.push_back(ld7);
	Load ld9("Load9", 9, 5.8, 2.0);
	model.loads.push_back(ld9);
	Load ld11("Load11", 11, 11.2, 7.5);
	model.loads.push_back(ld11);
	Load ld13("Load13", 13, 6.2, 1.6);
	model.loads.push_back(ld13);
	Load ld14("Load14", 14, 8.2, 2.5);
	model.loads.push_back(ld14);
	Load ld15("Load15", 15, 3.5, 1.8);
	model.loads.push_back(ld15);
	Load ld16("Load16", 16, 9.0, 5.8);
	model.loads.push_back(ld16);
	Load ld17("Load17", 17, 3.2, 0.9);
	model.loads.push_back(ld17);
	Load ld18("Load18", 18, 9.5, 3.4);
	model.loads.push_back(ld18);
	Load ld19("Load19", 19, 2.2, 0.7);
	model.loads.push_back(ld19);
	Load ld20("Load20", 20, 17.5, 11.2);
	model.loads.push_back(ld20);
	Load ld22("Load22", 22, 3.2, 1.6);
	model.loads.push_back(ld22);
	Load ld23("Load23", 23, 8.7, 6.7);
	model.loads.push_back(ld23);
	Load ld25("Load25", 25, 3.5, 2.3);
	model.loads.push_back(ld25);
	Load ld28("Load28", 28, 2.4, 0.9);
	model.loads.push_back(ld28);
	Load ld29("Load29", 29, 10.6, 1.9);
	model.loads.push_back(ld29);


	//Generation creation
	/*
	Generator g0("Generator0", 0, 23.54, 1.0, -20.0, 150.0, true);
	model.generators.push_back(g0);
	Generator g1("Generator1", 1, 60.97, 1.0, -20.0, 60.0, true);
	model.generators.push_back(g1);
	Generator g2("Generator2", 21, 21.59, 1.0, -15.0, 62.5, true);
	model.generators.push_back(g2);
	Generator g3("Generator3", 26, 26.91, 1.0, -15.0, 48.7, true);
	model.generators.push_back(g3);
	Generator g4("Generator4", 22, 19.2, 1.0, -10.0, 40.0, true);
	model.generators.push_back(g4);
	Generator g5("Generator5", 12, 37.0, 1.0, -15.0, 44.7, true);
	model.generators.push_back(g5);
	 */

	//Branch creation
	TransformerType TType0("Transformer_type0", cx_double(0.02, 0.06), cx_double(0.0, 33.3333333333));
	Transformer T0("Transformer0", 0, 1, TType0);
	model.transformers.push_back(T0);
	TransformerType TType1("Transformer_type1", cx_double(0.05, 0.19), cx_double(0.0, 50.0));
	Transformer T1("Transformer1", 0, 2, TType1);
	model.transformers.push_back(T1);
	TransformerType TType2("Transformer_type2", cx_double(0.06, 0.17), cx_double(0.0, 50.0));
	Transformer T2("Transformer2", 1, 3, TType2);
	model.transformers.push_back(T2);
	LineType LType0("Line_type0", 0.01, 0.04, 0.0, true);
	Line L0("Line0", 2, 3, LType0, 1.0);
	model.lines.push_back(L0);
	TransformerType TType3("Transformer_type3", cx_double(0.05, 0.2), cx_double(0.0, 50.0));
	Transformer T3("Transformer3", 1, 4, TType3);
	model.transformers.push_back(T3);
	TransformerType TType4("Transformer_type4", cx_double(0.06, 0.18), cx_double(0.0, 50.0));
	Transformer T4("Transformer4", 1, 5, TType4);
	model.transformers.push_back(T4);
	LineType LType1("Line_type1", 0.01, 0.04, 0.0, true);
	Line L1("Line1", 3, 5, LType1, 1.0);
	model.lines.push_back(L1);
	TransformerType TType5("Transformer_type5", cx_double(0.05, 0.12), cx_double(0.0, 100.0));
	Transformer T5("Transformer5", 4, 6, TType5);
	model.transformers.push_back(T5);
	TransformerType TType6("Transformer_type6", cx_double(0.03, 0.08), cx_double(0.0, 100.0));
	Transformer T6("Transformer6", 5, 6, TType6);
	model.transformers.push_back(T6);
	LineType LType2("Line_type2", 0.01, 0.04, 0.0, true);
	Line L2("Line2", 5, 7, LType2, 1.0);
	model.lines.push_back(L2);
	LineType LType3("Line_type3", 0.0, 0.21, 0.0, true);
	Line L3("Line3", 5, 8, LType3, 1.0);
	model.lines.push_back(L3);
	LineType LType4("Line_type4", 0.0, 0.56, 0.0, true);
	Line L4("Line4", 5, 9, LType4, 1.0);
	model.lines.push_back(L4);
	LineType LType5("Line_type5", 0.0, 0.21, 0.0, true);
	Line L5("Line5", 8, 10, LType5, 1.0);
	model.lines.push_back(L5);
	LineType LType6("Line_type6", 0.0, 0.11, 0.0, true);
	Line L6("Line6", 8, 9, LType6, 1.0);
	model.lines.push_back(L6);
	LineType LType7("Line_type7", 0.0, 0.26, 0.0, true);
	Line L7("Line7", 3, 11, LType7, 1.0);
	model.lines.push_back(L7);
	LineType LType8("Line_type8", 0.0, 0.14, 0.0, true);
	Line L8("Line8", 11, 12, LType8, 1.0);
	model.lines.push_back(L8);
	LineType LType9("Line_type9", 0.12, 0.26, 0.0, true);
	Line L9("Line9", 11, 13, LType9, 1.0);
	model.lines.push_back(L9);
	LineType LType10("Line_type10", 0.07, 0.13, 0.0, true);
	Line L10("Line10", 11, 14, LType10, 1.0);
	model.lines.push_back(L10);
	LineType LType11("Line_type11", 0.09, 0.2, 0.0, true);
	Line L11("Line11", 11, 15, LType11, 1.0);
	model.lines.push_back(L11);
	LineType LType12("Line_type12", 0.22, 0.2, 0.0, true);
	Line L12("Line12", 13, 14, LType12, 1.0);
	model.lines.push_back(L12);
	LineType LType13("Line_type13", 0.08, 0.19, 0.0, true);
	Line L13("Line13", 15, 16, LType13, 1.0);
	model.lines.push_back(L13);
	LineType LType14("Line_type14", 0.11, 0.22, 0.0, true);
	Line L14("Line14", 14, 17, LType14, 1.0);
	model.lines.push_back(L14);
	LineType LType15("Line_type15", 0.06, 0.13, 0.0, true);
	Line L15("Line15", 17, 18, LType15, 1.0);
	model.lines.push_back(L15);
	LineType LType16("Line_type16", 0.03, 0.07, 0.0, true);
	Line L16("Line16", 18, 19, LType16, 1.0);
	model.lines.push_back(L16);
	LineType LType17("Line_type17", 0.09, 0.21, 0.0, true);
	Line L17("Line17", 9, 19, LType17, 1.0);
	model.lines.push_back(L17);
	LineType LType18("Line_type18", 0.03, 0.08, 0.0, true);
	Line L18("Line18", 9, 16, LType18, 1.0);
	model.lines.push_back(L18);
	LineType LType19("Line_type19", 0.03, 0.07, 0.0, true);
	Line L19("Line19", 9, 20, LType19, 1.0);
	model.lines.push_back(L19);
	LineType LType20("Line_type20", 0.07, 0.15, 0.0, true);
	Line L20("Line20", 9, 21, LType20, 1.0);
	model.lines.push_back(L20);
	LineType LType21("Line_type21", 0.01, 0.02, 0.0, true);
	Line L21("Line21", 20, 21, LType21, 1.0);
	model.lines.push_back(L21);
	LineType LType22("Line_type22", 0.1, 0.2, 0.0, true);
	Line L22("Line22", 14, 22, LType22, 1.0);
	model.lines.push_back(L22);
	LineType LType23("Line_type23", 0.12, 0.18, 0.0, true);
	Line L23("Line23", 21, 23, LType23, 1.0);
	model.lines.push_back(L23);
	LineType LType24("Line_type24", 0.13, 0.27, 0.0, true);
	Line L24("Line24", 22, 23, LType24, 1.0);
	model.lines.push_back(L24);
	LineType LType25("Line_type25", 0.19, 0.33, 0.0, true);
	Line L25("Line25", 23, 24, LType25, 1.0);
	model.lines.push_back(L25);
	LineType LType26("Line_type26", 0.25, 0.38, 0.0, true);
	Line L26("Line26", 24, 25, LType26, 1.0);
	model.lines.push_back(L26);
	LineType LType27("Line_type27", 0.11, 0.21, 0.0, true);
	Line L27("Line27", 24, 26, LType27, 1.0);
	model.lines.push_back(L27);
	LineType LType28("Line_type28", 0.0, 0.4, 0.0, true);
	Line L28("Line28", 27, 26, LType28, 1.0);
	model.lines.push_back(L28);
	LineType LType29("Line_type29", 0.22, 0.42, 0.0, true);
	Line L29("Line29", 26, 28, LType29, 1.0);
	model.lines.push_back(L29);
	LineType LType30("Line_type30", 0.32, 0.6, 0.0, true);
	Line L30("Line30", 26, 29, LType30, 1.0);
	model.lines.push_back(L30);
	LineType LType31("Line_type31", 0.24, 0.45, 0.0, true);
	Line L31("Line31", 28, 29, LType31, 1.0);
	model.lines.push_back(L31);
	TransformerType TType7("Transformer_type7", cx_double(0.06, 0.2), cx_double(0.0, 50.0));
	Transformer T7("Transformer7", 7, 27, TType7);
	model.transformers.push_back(T7);
	TransformerType TType8("Transformer_type8", cx_double(0.02, 0.06), cx_double(0.0, 100.0));
	Transformer T8("Transformer8", 5, 27, TType8);
	model.transformers.push_back(T8);

	return model;

}

Circuit IEEE_14_bus(bool withPV) {

	Circuit model("circuit 1");

	/*Buses definition*/
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
	//once each bus is added to the list, its Index property is modified

	/*External grids*/
	ExternalGrid eg("External1", b1.index);
	model.external_grids.push_back(eg);

	/*Lines definition*/
	LineType ltype1("line type 1", 0.05, 0.11, 0.02, true);
	LineType ltype2("line type 2", 0.03, 0.08, 0.02, true);
	LineType ltype3("line type 3", 0.04, 0.09, 0.02, true);
	LineType ltype4("line type 4", 0.06, 0.13, 0.03, true);
	//cout << "[" << ltype1.Name << "]" << " -> [" << ltype1.impedance << ", " << ltype1.shunt_admittance << "]" << endl;

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

	/*Generators*/
	if (withPV) {
		Generator g2("Gen2", b2.index, 21.7, 12.7, -40.0, 50.0, false);
		Generator g3("Gen3", b3.index, 94.2, 19.0, 0.0, 40.0, false);
		Generator g4("Gen4", b4.index, 47.8, -3.9, 0.0, 0.0, false);
		Generator g5("Gen5", b5.index, 7.6, 1.6, 0.0, 0.0, false);
		Generator g6("Gen6", b6.index, 11.2, 7.5, -6.0, 24.0, false);
		//Generator g7("Gen7", b7.index, 0.0, 0.0, 0.0, 0.0, false);
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
		model.generators.push_back(g8);
		model.generators.push_back(g9);
		model.generators.push_back(g10);
		model.generators.push_back(g11);
		model.generators.push_back(g12);
		model.generators.push_back(g13);
		model.generators.push_back(g14);
	}

	/*Loads*/
	//Load ld1("Load1", b1.index, 0.0, 0.0);
	Load ld2("Load2", b2.index, 21.7, 12.7);
	Load ld3("Load3", b3.index, 94.2, 19.0);
	Load ld4("Load4", b4.index, 47.8, -3.9);
	Load ld5("Load5", b5.index, 7.6, 1.6);
	Load ld6("Load6", b6.index, 11.2, 7.5);
	//Load ld7("Load7", b7.index, 0.0, 0.0);
	//Load ld8("Load8", b8.index, 0.0, 0.0);
	Load ld9("Load9", b9.index, 29.5, 16.6);
	Load ld10("Load10", b10.index, 9.0, 5.8);
	Load ld11("Load11", b11.index, 3.5, 1.8);
	Load ld12("Load12", b12.index, 6.1, 1.6);
	Load ld13("Load13", b13.index, 13.5, 5.8);
	Load ld14("Load14", b14.index, 14.9, 5.0);

	model.loads.push_back(ld2);
	model.loads.push_back(ld3);
	model.loads.push_back(ld4);
	model.loads.push_back(ld5);
	model.loads.push_back(ld6);
	model.loads.push_back(ld9);
	model.loads.push_back(ld10);
	model.loads.push_back(ld11);
	model.loads.push_back(ld12);
	model.loads.push_back(ld13);
	model.loads.push_back(ld14);
	return model;
}

/*
 * This circuit is the test circuit proposed in the book:
 * Power system load flow analysis by Lynn Powell
 */
Circuit circuit_Lynn_Powell(bool withPV) {

	Circuit model("circuit 1");

	/*Buses definition*/
	Bus b1("bus1", undefined_bus_type, 10.0);
	Bus b2("bus2", undefined_bus_type, 10.0);
	Bus b3("bus3", undefined_bus_type, 10.0);
	Bus b4("bus4", undefined_bus_type, 10.0);
	Bus b5("bus5", undefined_bus_type, 10.0);
	model.add_Bus(b1);
	model.add_Bus(b2);
	model.add_Bus(b3);
	model.add_Bus(b4);
	model.add_Bus(b5);
	//once each bus is added to the list, its Index property is modified

	/*External grids*/
	ExternalGrid eg("External1", b1.index);
	model.external_grids.push_back(eg);

	/*Lines definition*/
	LineType ltype1("line type 1", 0.05, 0.11, 0.02, true);
	LineType ltype2("line type 2", 0.03, 0.08, 0.02, true);
	LineType ltype3("line type 3", 0.04, 0.09, 0.02, true);
	LineType ltype4("line type 4", 0.06, 0.13, 0.03, true);
	//cout << "[" << ltype1.Name << "]" << " -> [" << ltype1.impedance << ", " << ltype1.shunt_admittance << "]" << endl;

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

	/*Generators*/
	if (withPV) {
		Generator g1("Generator", b4.index, 70, 10, -1, 2.3, false);
		model.generators.push_back(g1);
	}

	/*Loads*/
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

/*
 *
 */
Circuit circuit_Gutierrez_Castillejos() {

	Circuit model("circuit 1");

	/*Buses definition*/
	Bus bgen_1("GEN_1", undefined_bus_type, 66.0);
	Bus bgen_2("GEN_2", undefined_bus_type, 66.0);
	Bus bcarga_3("CARGA_3", undefined_bus_type, 66.0);
	Bus bcarga_4("CARGA_4", undefined_bus_type, 66.0);
	Bus bcarga_1("CARGA_1", undefined_bus_type, 66.0);
	model.add_Bus(bgen_1);
	model.add_Bus(bgen_2);


	model.add_Bus(bcarga_3);
	model.add_Bus(bcarga_4);
	model.add_Bus(bcarga_1);
	//once each bus is added to the list, its Index property is modified

	/*External grids*/
	ExternalGrid eg("External1", bgen_1.index);
	model.external_grids.push_back(eg);

	/*Lines definition*/
	//cout << "[" << ltype1.Name << "]" << " -> [" << ltype1.impedance << ", " << ltype1.shunt_admittance << "]" << endl;

	Line l1("Gen_1-Carga_1", bgen_1.index, bcarga_3.index, LineType("", 0.05935, 0.3079, 0.00331, true), 1.0);
	Line l2("GEN_1-CARGA_3", bgen_1.index, bcarga_1.index, LineType("", 0.1511, 0.57677, 0.00594, true), 1.0);
	Line l3("GEN_2-CARGA_3", bgen_2.index, bcarga_3.index, LineType("", 0.03958, 0.20528, 0.00221, true), 1.0);
	Line l4("GEN_2-CARGA_4", bgen_2.index, bcarga_4.index, LineType("", 0.03958, 0.20528, 0.00221, true), 1.0);
	Line l5("GEN_2-CARGA_1", bgen_2.index, bcarga_1.index, LineType("", 0.10995, 0.41956, 0.00432, true), 1.0);
	Line l6("CARGA_3-CARGA_1", bcarga_3.index, bcarga_4.index, LineType("", 0.04947, 0.25659, 0.00276, true), 1.0);
	Line l7("CARGA_4-CARGA_1", bcarga_4.index, bcarga_1.index, LineType("", 0.02968, 0.15397, 0.00166, true), 1.0);
	model.lines.push_back(l1);
	model.lines.push_back(l2);
	model.lines.push_back(l3);
	model.lines.push_back(l4);
	model.lines.push_back(l5);
	model.lines.push_back(l6);
	model.lines.push_back(l7);

	/*Generators*/
	Generator g1("gen1", bgen_2.index, 23, 1, -10, 10, true);
	Generator g2("GEN_2", bgen_2.index, 23, 1, -10, 10, true);
	model.generators.push_back(g1);
	model.generators.push_back(g2);

	/*Loads*/
	Load ld1("GEN_2", bgen_2.index, 4, 1);
	Load ld2("CARGA_3", bcarga_3.index, 25, 5);
	Load ld3("CARGA_4", bcarga_4.index, 22, 6);
	Load ld4("CARGA_1", bcarga_1.index, 29, 8);
	model.loads.push_back(ld1);
	model.loads.push_back(ld2);
	model.loads.push_back(ld3);
	model.loads.push_back(ld4);

	return model;
}