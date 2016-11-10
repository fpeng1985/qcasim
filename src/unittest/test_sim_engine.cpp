//
// Created by fpeng on 2016/11/9.
//

#define CATCH_CONFIG_MAIN
#include "catch.hpp"

using namespace std;

#include "sim_engine.h"
using namespace hfut;

#include <boost/assign/std/vector.hpp>
using namespace boost::assign;

SCENARIO("2 cells case 1", "[2 cells case 1]") {
    QCACircuit::CircuitStructure circuit_structure_matrix;
    circuit_structure_matrix.push_back(vector<int>());
    circuit_structure_matrix[0] += -1, 1;

    shared_ptr<QCACircuit> circuit(new QCACircuit);
    circuit->populate_cells(circuit_structure_matrix);

    SimEngine engine;
    engine.set_circuit(circuit);

    GIVEN("a input polarization") {
        Polarization input_p;
        input_p.push_back(make_tuple(0,0,1));

        WHEN("we run simulation") {
            engine.run_simulation(input_p);

            THEN("we get the circuit output") {
                REQUIRE(circuit->get_cell(0,1) != nullptr);
                REQUIRE(abs(circuit->get_cell(0,1)->polarization - 0.952) < 0.01);
            }
        }
    }
}

SCENARIO("2 cells case 2", "[2 cells case 2]") {
    QCACircuit::CircuitStructure circuit_structure_matrix;
    circuit_structure_matrix.push_back(vector<int>());
    circuit_structure_matrix[0] += -1, 0, 1;

    shared_ptr<QCACircuit> circuit(new QCACircuit);
    circuit->populate_cells(circuit_structure_matrix);

    SimEngine engine;
    engine.set_circuit(circuit);

    GIVEN("a input polarization") {
        Polarization input_p;
        input_p.push_back(make_tuple(0,0,1));

        WHEN("we run simulation") {
            engine.run_simulation(input_p);

            THEN("we get the circuit output") {
                REQUIRE(circuit->get_cell(0,2) != nullptr);
                REQUIRE(abs(circuit->get_cell(0,2)->polarization - 0.0961) < 0.01);
            }
        }
    }
}

SCENARIO("2 cells case 3", "[2 cells case 3]") {
    QCACircuit::CircuitStructure circuit_structure_matrix;
    circuit_structure_matrix.push_back(vector<int>());
    circuit_structure_matrix.push_back(vector<int>());
    circuit_structure_matrix[0] += 0, 1;
    circuit_structure_matrix[1] += -1, 0;

    shared_ptr<QCACircuit> circuit(new QCACircuit);
    circuit->populate_cells(circuit_structure_matrix);

    SimEngine engine;
    engine.set_circuit(circuit);

    GIVEN("a input polarization") {
        Polarization input_p;
        input_p.push_back(make_tuple(1,0,1));

        WHEN("we run simulation") {
            engine.run_simulation(input_p);

            THEN("we get the circuit output") {
                REQUIRE(circuit->get_cell(0,1) != nullptr);
                REQUIRE(abs(circuit->get_cell(0,1)->polarization - (-0.562)) < 0.01);
            }
        }
    }
}

SCENARIO("majority gate 1", "[majority_gate_1]") {

    shared_ptr<QCACircuit> circuit(new QCACircuit);

    QCACircuit::CircuitStructure circuit_structure_matrix;
    GIVEN("A circuit majority gate 1") {

        for (int i = 0; i < 5; ++i) {
            circuit_structure_matrix.push_back(vector<int>());
        }

        circuit_structure_matrix[0] +=  0, 0, -1, 0,  0;
        circuit_structure_matrix[1] +=  0, 1,  1, 1,  0;
        circuit_structure_matrix[2] += -1, 1,  1, 1, -2;
        circuit_structure_matrix[3] +=  0, 1,  1, 1,  0;
        circuit_structure_matrix[4] +=  0, 0, -1, 0,  0;

        circuit->populate_cells(circuit_structure_matrix);

        SimEngine engine;
        engine.set_circuit(circuit);

        WHEN("we feed [0,0,0] to the IterativeSimEngine") {
            Polarization input_p;
            input_p.push_back(make_tuple(0, 2, -1));
            input_p.push_back(make_tuple(2, 0, -1));
            input_p.push_back(make_tuple(4, 2, -1));

            engine.run_simulation(input_p);

            THEN("we get the circuit output") {
                REQUIRE(circuit->get_cell(2, 4)->polarization < -0.5);
            }
        }

        WHEN("we feed [0,0,1] to the IterativeSimEngine") {
            Polarization input_p;
            input_p.push_back(make_tuple(0, 2, -1));
            input_p.push_back(make_tuple(2, 0, -1));
            input_p.push_back(make_tuple(4, 2, 1));

            engine.run_simulation(input_p);

            THEN("we get the circuit output") {
                REQUIRE(circuit->get_cell(2, 4)->polarization < -0.5);
            }
        }

        WHEN("we feed [0,1,0] to the IterativeSimEngine") {
            Polarization input_p;
            input_p.push_back(make_tuple(0, 2, -1));
            input_p.push_back(make_tuple(2, 0, 1));
            input_p.push_back(make_tuple(4, 2, -1));

            engine.run_simulation(input_p);

            THEN("we get the circuit output") {
                REQUIRE(circuit->get_cell(2, 4)->polarization < -0.5);
            }
        }

        WHEN("we feed [0,1,1] to the IterativeSimEngine") {
            Polarization input_p;
            input_p.push_back(make_tuple(0, 2, -1));
            input_p.push_back(make_tuple(2, 0, 1));
            input_p.push_back(make_tuple(4, 2, 1));

            engine.run_simulation(input_p);

            THEN("we get the circuit output") {
                REQUIRE(circuit->get_cell(2, 4)->polarization > 0.5);
            }
        }

        WHEN("we feed [1,0,0] to the IterativeSimEngine") {
            Polarization input_p;
            input_p.push_back(make_tuple(0, 2, 1));
            input_p.push_back(make_tuple(2, 0, -1));
            input_p.push_back(make_tuple(4, 2, -1));

            engine.run_simulation(input_p);

            THEN("we get the circuit output") {
                REQUIRE(circuit->get_cell(2, 4)->polarization < -0.5);
            }
        }

        WHEN("we feed [1,0,1] to the IterativeSimEngine") {
            Polarization input_p;
            input_p.push_back(make_tuple(0, 2, 1));
            input_p.push_back(make_tuple(2, 0, -1));
            input_p.push_back(make_tuple(4, 2, 1));

            engine.run_simulation(input_p);

            THEN("we get the circuit output") {
                REQUIRE(circuit->get_cell(2, 4)->polarization > 0.5);
            }
        }

        WHEN("we feed [1,1,0] to the IterativeSimEngine") {
            Polarization input_p;
            input_p.push_back(make_tuple(0, 2, 1));
            input_p.push_back(make_tuple(2, 0, 1));
            input_p.push_back(make_tuple(4, 2, -1));

            engine.run_simulation(input_p);

            THEN("we get the circuit output") {
                REQUIRE(circuit->get_cell(2, 4)->polarization > 0.5);
            }
        }

        WHEN("we feed [1,1,1] to the IterativeSimEngine") {
            Polarization input_p;
            input_p.push_back(make_tuple(0, 2, 1));
            input_p.push_back(make_tuple(2, 0, 1));
            input_p.push_back(make_tuple(4, 2, 1));

            engine.run_simulation(input_p);

            THEN("we get the circuit output") {
                REQUIRE(circuit->get_cell(2, 4)->polarization > 0.5);
            }
        }
    }
}