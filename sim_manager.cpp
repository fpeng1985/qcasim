//
// Created by Administrator on 2016/10/14.
//

#include "SimManager.h"

#include <fstream>
#include <sstream>
#include <list>
#include <itrator>
#include <algorithm>

namespace hfut {

    using namespace std;

    void SimManager::load_benchmark(const string &path) {
        structures.clear();
        structures.push_back(QCACircuit::CircuitStructure());

        //read file into strings
        ifstream in(path);
        list<string> lines;
        string line;
        while (!in.eof()) {
            get_line(in, line);
            lines.push_front(line);//insert into the front
        }
        in.close();

        //convert string to CircuitStructure and feed into the 0th position
        istringstream iss;
        int i = 0;
        for (auto &line : lines) {
            iss.str(line);

            structures[0].push_back(vector<int>());

            copy(istream_iterator<int>(iss), istream_iterator<int>(), back_inserter(structures[0][i]));

            ++i;
        }

        //populate circuit
        circuit->populate_cells(structures[0]);

        //generate other circuit structures
        //TODO
    }

}