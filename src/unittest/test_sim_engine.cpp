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

    shared_ptr<QCACircuit> circuit(new QCACircuit);

    vector<vector<int>> circuit_structure_matrix;
    GIVEN("A circuit containing two cells") {

        circuit_structure_matrix.push_back(vector<int>());
        circuit_structure_matrix[0] += -1, 1;
        circuit->populate_cells(circuit_structure_matrix);

        SimEngine engine;
        engine.set_circuit(circuit);

        WHEN("we feed [1] to the SimEngine") {
            Polarization input_p;
            input_p.insert(make_pair(make_pair(0, 0), 1));

            engine.set_input_polarization(input_p);
            engine.set_non_input_polarization_randomly();

//            cout << *circuit << endl;
            auto pola = engine.compute_polarization_from_neighbour_cells(0,1);
            circuit->get_cell(0, 1)->polarization = pola;
//            cout << *circuit << endl;

            THEN("we get the circuit output") {
                cout << pola << endl;
                REQUIRE(abs(pola - 0.952) < 0.01);
            }
        }
    }
}

SCENARIO("2 cells case 2", "[2 cells case 2") {

    shared_ptr<QCACircuit> circuit(new QCACircuit);

    vector<vector<int>> circuit_structure_matrix;
    GIVEN("A circuit containing two cells") {

        circuit_structure_matrix.push_back(vector<int>());
        circuit_structure_matrix[0] += -1, 0, 1;
        circuit->populate_cells(circuit_structure_matrix);

        SimEngine engine;
        engine.set_circuit(circuit);

        WHEN("we feed [1] to the SimEngine") {
            Polarization input_p;
            input_p.insert(make_pair(make_pair(0, 0), 1));

            engine.set_input_polarization(input_p);
            engine.set_non_input_polarization_randomly();

//            cout << *circuit << endl;
            auto pola = engine.compute_polarization_from_neighbour_cells(0,2);
            circuit->get_cell(0, 2)->polarization = pola;
//            cout << *circuit << endl;

            THEN("we get the circuit output") {
                cout << pola << endl;
                REQUIRE(abs(pola - 0.0961) < 0.01);
            }
        }
    }
}

SCENARIO("2 cells case 3", "[2 cells case 3]") {

    shared_ptr<QCACircuit> circuit(new QCACircuit);

    vector<vector<int>> circuit_structure_matrix;
    GIVEN("A circuit containing two cells") {

        circuit_structure_matrix.push_back(vector<int>());
        circuit_structure_matrix.push_back(vector<int>());

        circuit_structure_matrix[0] += 0, 1;
        circuit_structure_matrix[1] += -1, 0;

        circuit->populate_cells(circuit_structure_matrix);


        SimEngine engine;
        engine.set_circuit(circuit);

        WHEN("we feed [1] to the SimEngine") {
            Polarization input_p;
            input_p.insert(make_pair(make_pair(1, 0), 1));

            engine.set_input_polarization(input_p);
            engine.set_non_input_polarization_randomly();

//            cout << *circuit << endl;
            auto pola = engine.compute_polarization_from_neighbour_cells(0,1);
            circuit->get_cell(0, 1)->polarization = pola;
//            cout << *circuit << endl;

            THEN("we get the circuit output") {
                cout << pola << endl;
                REQUIRE(abs(pola - (-0.562)) < 0.01);
            }
        }
    }
}