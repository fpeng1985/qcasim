//
// Created by Administrator on 2016/10/14.
//

#include "sim_manager.h"

#include <fstream>
#include <sstream>
#include <list>
#include <iterator>
#include <algorithm>

namespace hfut {

    using namespace std;


    SimManager::SimManager() {
        circuit = make_shared<QCACircuit>();

        engine = make_shared<SimEngine>();
        engine->set_circuit(circuit);
    }

    void SimManager::load_benchmark(const string &path) {
        QCACircuit::CircuitStructure structure;
        //read file into strings
        ifstream ifs(path);
        list<string> lines;
        for (string line; getline(ifs, line);) {
            lines.push_front(line);//insert into the front
        }
        ifs.close();

        //convert string to CircuitStructure
        for (auto &line : lines) {
            structure.emplace_back();
            istringstream iss(line);
            copy(istream_iterator<int>(iss), istream_iterator<int>(), back_inserter(*structure.rbegin()));
        }

        //populate circuit
        circuit->populate_cells(structure);

        //generate all the circuit structures
        structures.clear();
        structures.push_back(structure);
        //TODO
    }

}