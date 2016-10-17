//
// Created by Administrator on 2016/10/14.
//

#include "sim_manager.h"

#include <fstream>
#include <sstream>
#include <list>
#include <iterator>
#include <algorithm>
#include <cassert>

namespace hfut {

    using namespace std;

    SimManager::SimManager() {
        circuit = make_shared<QCACircuit>();

        engine = make_shared<SimEngine>();
        engine->set_circuit(circuit);
    }

    void SimManager::load_benchmark(const string &path) {
        //read benchmark file into list of string in reverse order
        ifstream ifs(path);
        list<string> lines;
        for (string line; getline(ifs, line);) {
            lines.push_front(line);//insert into the front
        }
        ifs.close();

        //convert the list of string to CircuitStructure, changing row and column order
        QCACircuit::CircuitStructure structure;
        bool first_line = true;
        for (auto &line : lines) {
            istringstream iss(line);
            if (first_line) {
                for (string s; getline(iss, s, '\t');) {
                    structure.emplace_back();
                    structure.rbegin()->push_back(stoi(s));
                }
                first_line = false;
            } else {
                int i=0;
                for (string s; getline(iss, s, '\t'); ++i) {
                    structure[i].push_back(stoi(s));
                }
            }
        }

        //populate circuit
        circuit->populate_cells(structure);

        //generate all the circuit structures
        structures.push_back(structure);

        //define helper types
        typedef std::pair<int, int> Index;
        typedef std::vector<Index> IndexContainer;

        //find all the normal cell's indices
        IndexContainer indices;

        shared_ptr<QCACell> cell;
        int i=0;
        for (auto rit=circuit->row_begin(); rit!=circuit->row_end(); ++rit, ++i) {
            int j=0;
            for (auto cit=circuit->col_begin(rit); cit!=circuit->col_end(rit); ++cit, ++j) {
                cell = *cit;

                if (cell != nullptr && cell->cell_type==CellType::Normal) {
                    indices.push_back(make_pair(i, j));
                }
            }
        }

        //generate all the combinations, in integer representation
        CombinationGenerator combgen;
        vector<vector<int>> combinations;
        combgen.generate_combination(indices.size(), combinations);

        //change the integer representation to index format
        vector<IndexContainer> containers;
        for (auto &comb : combinations) {
            containers.emplace_back();
            for (auto idx : comb) {
                containers.rbegin()->push_back(indices[idx]);
            }
        }

        //fill the structures
        for (auto &container : containers) {
            structures.push_back(structures[0]);

            QCACircuit::CircuitStructure &s = *structures.rbegin();

            for (auto &index : container) {
                int &xidx = index.first;
                int &yidx = index.second;

                s[xidx][yidx] = 0;
            }
        }

        assert(structures.size() == combinations.size()+1);
    }

    void SimManager::CombinationGenerator::generate_combination(int n, vector<vector<int>> &combinations) {
        vector<int> a;
        for (int i=0; i<n; ++i) {
            a.push_back(i);
        }

        vector<int> b;
        b.resize(a.size());

        for (int m=1; m<=n; ++m) {
            combine(a, n, b, m, m, combinations);
        }
    }

    void SimManager::CombinationGenerator::combine(const vector<int> &a, int n, vector<int> &b, int m, const int M, vector<vector<int>> &combinations) {
        for(int i=n; i>=m; i--) {
            b[m-1] = i - 1;
            if (m > 1)
                combine(a, i-1, b, m-1, M, combinations);
            else {
                vector<int> tmp;
                for(int j=M-1; j>=0; j--)
                    tmp.push_back(a[b[j]]);
                combinations.push_back(tmp);
            }
        }
    }

}