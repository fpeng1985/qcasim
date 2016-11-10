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
        //[1]reset input_size
        input_size = 0;

        //[2]read benchmark file into structure, row by row, column by column
        QCACircuit::CircuitStructure structure;

        ifstream ifs(path);
        istringstream iss;
        for (string line; getline(ifs, line);) {
            iss.str(line);
            iss.seekg(0);

            structure.emplace_back();
            copy(istream_iterator<int>(iss), istream_iterator<int>(), back_inserter(*structure.rbegin()));
        }
        ifs.close();
        structures.push_back(structure);

        //[3]find all the non-input cell's indices, and set the input size
        typedef std::pair<int, int> Index;
        typedef std::vector<Index> IndexContainer;

        IndexContainer indices;
        for (size_t i=0; i<structure.size(); ++i) {
            for (size_t j=0; j<structure[i].size(); ++j) {
                if (structure[i][j] != 0) {
                    if (structure[i][j] == -1) {//input cell
                        ++input_size;
                    } else {//non-input cell
                        indices.push_back(make_pair(i, j));
                    }
                }
            }
        }

        //[4]generate all the combinations, in integer representation
        CombinationGenerator combgen;
        vector<vector<int>> combinations;
        combgen.generate_combination(indices.size(), combinations);

        //[5]change the integer representation to index format
        vector<IndexContainer> containers;
        for (auto &comb : combinations) {
            containers.emplace_back();
            for (auto idx : comb) {
                containers.rbegin()->push_back(indices[idx]);
            }
        }

        //[6]generate all the structures
        for (auto &container : containers) {
            structures.push_back(structures[0]);

            QCACircuit::CircuitStructure &s = *structures.rbegin();

            for (auto &index : container) {
                int &ridx = index.first;
                int &cidx = index.second;

                s[ridx][cidx] = 0;
            }
        }

        assert(structures.size() == combinations.size()+1);
    }

    void SimManager::test_benchmark() {
        int input_comb_size = 1<<input_size;

        for (auto &structure : structures) {
            circuit->populate_cells(structure);
            QCATruthTable truth_table;

            //for each input combination compute its output
            for (int i=0; i<input_comb_size; ++i) {
                //initialize input bits
                vector<int> bits;
                int j=i;
                do {
                    bits.push_back(j^1);
                    j = j>>1;
                } while(j!=0);

                //initialize input in both truth value and polarization value form
                QCATruthValueList input_truth_vals;
                PolarizationList input_p = circuit->get_input_polarizations();

                assert(bits.size() == input_size);
                assert(input_p.size() == input_size);

                for (size_t k=0; k<bits.size(); ++k) {
                    input_truth_vals.push_back(make_tuple(get<0>(input_p[k]), get<1>(input_p[k]), bits[k]));
                    get<2>(input_p[k]) = convert_logic_to_polarization(bits[k]);
                }

                //run the simulation
                engine->run_simulation(input_p);

                PolarizationList output_p = circuit->get_output_polarizations();

                QCATruthValueList output_truth_vals;
                for (size_t k=0; k<output_p.size(); ++k) {
                    output_truth_vals.push_back(make_tuple(get<0>(output_p[k]), get<1>(output_p[k]),
                                                           convert_polarization_to_logic(get<2>(output_p[k]))));
                }

                truth_table.insert(make_pair(input_truth_vals, output_truth_vals));
            }

            tables.push_back(truth_table);
        }
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