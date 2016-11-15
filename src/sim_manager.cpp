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

    QCATruthTable::QCATruthTable(const std::vector<std::pair<int, int>> &input_idx, const std::vector<std::pair<int, int>> &output_idx) {
#ifndef NDBUG
        input_size  = input_idx.size();
        output_size = output_idx.size();
#else
        size_t input_size  = input_idx.size();
        size_t output_size = output_idx.size();
#endif

        int input_comb_size  = 1<<input_size;

        QCATruthValueSet inputs;
        QCATruthValueSet outputs;

        for (int i=0; i<input_comb_size; ++i) {
            //get input values
            inputs.clear();

            int j = i;
            for (size_t k=0; k < input_size; ++k) {
                inputs.insert( make_tuple(input_idx[k].first, input_idx[k].second, j^1) );
                j /= 2;
            }

            assert(j == 0);
            assert(inputs.size()  == input_size);

            //get output values
            outputs.clear();

            for (size_t k=0; k<output_size; ++k) {
                outputs.insert( make_tuple(output_idx[k].first, output_idx[k].second, 0) );
            }

            assert(outputs.size() == output_size);

            //insert the mapping from inputs to outputs
            table.insert(make_pair(inputs, outputs));
        }
    }

    std::ostream &operator<<(std::ostream &os, const QCATruthTable &truth_table) {
        for (auto &mapping : truth_table.table) {
            auto &inputs  = mapping.first;
            auto &outputs = mapping.second;

            for (auto &input:inputs) {
                os << "(" << get<0>(input) << "," << get<1>(input) << "," << get<2>(input) << ") ";
            }

            os << " ---> ";

            for (auto &output:outputs) {
                os << "(" << get<0>(output) << "," << get<1>(output) << "," << get<2>(output) << ") ";
            }

            os << endl;
        }

        return os;
    }

    bool operator==(const QCATruthTable &lhs, const QCATruthTable &rhs) {
#ifndef NDBUG
        if (lhs.input_size  != rhs.input_size)  return false;
        if (lhs.output_size != rhs.output_size) return false;
#endif

        if (lhs.table != rhs.table) return false;

        return true;
    }

    bool operator!=(const QCATruthTable &lhs, const QCATruthTable &rhs) {
#ifndef NDBUG
        if (lhs.input_size  != rhs.input_size)  return true;
        if (lhs.output_size != rhs.output_size) return true;
#endif

        if (lhs.table != rhs.table) return true;

        return false;
    }

    SimResultGroup::SimResultGroup(const QCACircuit::CircuitStructure &structure, const QCATruthTable &table,
                                   const std::vector<std::vector<std::pair<int, int>>> &index_combinations) {
        QCACircuit::CircuitStructure cs;

        for (auto &index_combination : index_combinations) {
            cs = structure;

            for (auto &index:index_combination) {
                auto &ridx = index.first;
                auto &cidx = index.second;

                cs[ridx][cidx] = 0;
            }

            results.push_back(SimResult(cs, table));
        }

        _correction_ratial = 0;
    }

    SimManager::SimManager() {
        circuit = make_shared<QCACircuit>();

        engine = make_shared<SimEngine>();
        engine->set_circuit(circuit);
    }

    void SimManager::load_benchmark(const string &path) {
        //[1]read benchmark file into circuit_structure, row by row, column by column
        benchmark_circuit_structure.clear();

        ifstream ifs(path);
        istringstream iss;
        for (string line; getline(ifs, line);) {
            iss.str(line);
            iss.seekg(0);

            benchmark_circuit_structure.emplace_back();
            copy(istream_iterator<int>(iss), istream_iterator<int>(), back_inserter(*benchmark_circuit_structure.rbegin()));
        }
        ifs.close();

        //[2]set cell sizes, fill cell indices
        input_cell_size  = 0;
        output_cell_size = 0;
        normal_cell_size = 0;

        input_idx.clear();
        normal_idx.clear();
        output_idx.clear();

        for (size_t i=0; i<benchmark_circuit_structure.size(); ++i) {
            for (size_t j=0; j<benchmark_circuit_structure[i].size(); ++j) {
                switch (benchmark_circuit_structure[i][j]) {
                    case -1://input cell
                        ++input_cell_size;
                        input_idx.push_back(make_pair(i, j));
                        break;
                    case 1://normal cell
                        ++normal_cell_size;
                        normal_idx.push_back(make_pair(i, j));
                        break;
                    case -2://output cell
                        ++output_cell_size;
                        output_idx.push_back(make_pair(i, j));
                        break;
                    default://no cell
                        continue;
                }
            }
        }

        //[3]initialize and update benchmark truth table
        benchmark_truth_table = QCATruthTable(input_idx, output_idx);
        compute_truth_table(benchmark_circuit_structure, benchmark_truth_table);
    }

    void SimManager::test_benchmark() {
        vector<vector<CellIndex>> combinations;
        CombinationGenerator<CellIndex> combgen;

        for (size_t m=1; m<=normal_cell_size; ++m) {
            combinations.clear();
            combgen.generate_combination(normal_idx, m, combinations);

            SimResultGroup result_grp(benchmark_circuit_structure, benchmark_truth_table, combinations);

            int correct_cnt = 0;
            for (auto &result : result_grp) {
                const QCACircuit::CircuitStructure &structure = result.circuit_structure;
                QCATruthTable &table = result.truth_table;

                compute_truth_table(structure, table);

                //check success
                result.success = true;
                for (auto &mapping : table) {
                    auto &outputs = mapping.second;

                    for (auto &output : outputs) {
                        if (output == -1) {
                            result.success = false;
                            break;
                        }
                    }

                    if (!result.success) break;
                }

                //check correct
                result.correct = (table==benchmark_truth_table);

                if (result.correct) {
                    ++correct_cnt;
                }
            }

            result_grp.set_correction_ratial(correct_cnt * 1.0 / result_grp.size());

            results.insert(make_pair(m, result_grp));
        }
    }

    void SimManager::compute_truth_table(const QCACircuit::CircuitStructure &structure, QCATruthTable &table) {
        circuit->populate_cells(structure);

        PolarizationList input_p;
        PolarizationList output_p;
        for (auto &mapping : table) {
            auto &inputs = mapping.first;

            input_p.clear();
            output_p.clear();

            for (auto &input : inputs) {
                input_p.push_back(make_tuple(get<0>(input), get<1>(input), convert_logic_to_polarization(get<2>(input))));
            }

            engine->run_simulation(input_p);

            output_p = circuit->get_output_polarizations();
            QCATruthValueSet simulated_outputs;
            for (auto &pola : output_p) {
                simulated_outputs.insert(make_tuple(get<0>(pola), get<1>(pola), convert_polarization_to_logic(get<2>(pola))));
            }

            table.set_value(inputs, simulated_outputs);
        }
    }

}