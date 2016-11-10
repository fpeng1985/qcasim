//
// Created by fpeng on 2016/11/9.
//

#include <iostream>
#include <fstream>
#include <memory>
#include <cassert>
using namespace std;

#include <boost/assign/std/vector.hpp>
using namespace boost::assign;

#include "sim_engine.h"
using namespace hfut;

int main(int argc, char *argv[]) {

    //define circuit structure using vector of vector
    vector<vector<int>> circuit_structure_matrix;

    for (int i = 0; i < 5; ++i) {
        circuit_structure_matrix.push_back(vector<int>());
    }

    circuit_structure_matrix[0] +=  0, 0, -1, 0,  0;
    circuit_structure_matrix[1] +=  0, 1,  1, 1,  0;
    circuit_structure_matrix[2] += -1, 1,  1, 1, -2;
    circuit_structure_matrix[3] +=  0, 1,  1, 1,  0;
    circuit_structure_matrix[4] +=  0, 0, -1, 0,  0;

    //populate circuit
    shared_ptr<QCACircuit> circuit(new QCACircuit);
    circuit->populate_cells(circuit_structure_matrix);

    //setup simulation engine
    IterativeSimEngine engine;
    engine.set_circuit(circuit);

    for (int i=0; i<8; ++i) {
        //initialize input
        int bit1 = i/4;
        int bit2 = (i%4)/2;
        int bit3 = i%2;

        int p1 = 2*bit1-1;
        int p2 = 2*bit2-1;
        int p3 = 2*bit3-1;

        PolarizationList input_p;
        input_p.push_back(make_tuple(0,2,p1));
        input_p.push_back(make_tuple(2,0,p2));
        input_p.push_back(make_tuple(4,2,p3));

        //run simulation
        engine.run_simulation(input_p);

        cout << "(" << p1 << ", " << p2 << ", " << p3 << ") : "  << circuit->get_cell(2,4)->polarization << endl;
        cout << *circuit << endl;
    }

    return 0;
}

