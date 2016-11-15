//
// Created by Administrator on 2016/10/13.
//
#define CATCH_CONFIG_MAIN
#include "catch.hpp"

using namespace std;

#include "qca_circuit.h"
using namespace hfut;

#include <boost/assign/std/vector.hpp>
using namespace boost::assign;


SCENARIO("majority gate 1", "[majority_gate_1]") {

    QCACircuit circuit;

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

        WHEN("we populate the circuit using description matrix") {
            circuit.populate_cells(circuit_structure_matrix);

            cout << circuit << endl;

            THEN("we get the circuit circuit_structure") {
                REQUIRE(circuit.get_cell(0, 2)->cell_type == CellType::Input);
                REQUIRE(circuit.get_cell(2, 4)->cell_type == CellType::Output);
            }
        }
    }
}

