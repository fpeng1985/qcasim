//
// Created by Administrator on 2016/10/14.
//

#define CATCH_CONFIG_MAIN
#include "catch.hpp"

using namespace std;

#include "sim_engine.h"
using namespace hfut;

#include <boost/assign/std/vector.hpp>
using namespace boost::assign;

SCENARIO("majority gate 1", "[majority_gate_1]") {

    shared_ptr<QCACircuit> circuit(new QCACircuit);

    vector<vector<int>> circuit_structure_matrix;
    GIVEN("A circuit majority gate 1") {

        for (int i=0; i<5; ++i) {
            circuit_structure_matrix.push_back(vector<int>());
        }

        circuit_structure_matrix[0] +=  0,  0, -1,  0,  0;
        circuit_structure_matrix[1] +=  0,  1,  1,  1,  0;
        circuit_structure_matrix[2] += -1,  1,  1,  1, -2;
        circuit_structure_matrix[3] +=  0,  1,  1,  1,  0;
        circuit_structure_matrix[4] +=  0,  0, -1,  0,  0;

        circuit->populate_cells(circuit_structure_matrix);

        WHEN("we feed the input to the SimEngine") {
            Polarization input_p;
            input_p.insert(make_pair(make_pair(0, 2), 1));
            input_p.insert(make_pair(make_pair(2, 0), 1));
            input_p.insert(make_pair(make_pair(4, 2), 1));

            SimEngine engine;
            engine.set_circuit(circuit);
            engine.run_simulation(input_p);

            THEN("we get the circuit output") {
                REQUIRE(circuit->get_cell(2, 4)->cell_type == CellType::Output);
                REQUIRE(circuit->get_cell(2, 4)->polarization > 0.5);
            }
        }
    }
}
