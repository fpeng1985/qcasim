//
// Created by Administrator on 2016/10/13.
//

#include "qca_circuit.h"

#include <iostream>
#include <random>
#include <cstdlib>

namespace hfut {
    using namespace std;

    QCACircuit::QCACircuit() {

    }

    void QCACircuit::populate_cells(const vector<vector<int>> &cell_structure_matrix) {
        for (size_t i=0; i<cell_structure_matrix.size(); ++i) {
            cells.push_back(vector<QCACell*>());

            for (size_t j=0; j<cell_structure_matrix[i].size(); ++j) {
                if (cell_structure_matrix[i][j] == 0) {
                    cells[i].push_back(nullptr);
                } else if (cell_structure_matrix[i][j] == 1) {
                    cells[i].push_back(new QCACell(i, j, 0, CellType::Normal));
                } else if (cell_structure_matrix[i][j] == -1) {
                    cells[i].push_back(new QCACell(i, j, 0, CellType::Input));
                } else if (cell_structure_matrix[i][j] == -2) {
                    cells[i].push_back(new QCACell(i, j, 0, CellType::Output));
                } else {
                        cerr << "Input format error!" << endl;
                        exit(-1);
                }
            }
        }

    }

    QCACell *QCACircuit::get_cell(size_t i, size_t j) {
        if (i>=0 && i<cells.size() && j>=0 && j<cells[i].size()) {
            return cells[i][j];
        } else {
            cout << i << " " << j << endl;
            cerr << "Overflow circuit index boundary!" << endl;
            exit(-2);
        }
    }

    void QCACircuit::clear() {
        for (size_t i=0; i<cells.size(); ++i) {
            for (size_t j=0; j<cells[i].size(); ++j) {
                if (cells[i][j] != nullptr) {
                    delete cells[i][j];
                    cells[i][j] = nullptr;
                }
            }
        }

        cells.clear();
    }

    void QCACircuit::initialize_polarization() {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(-1, 1);

        for (auto &line : cells) {
            for (auto &cell : line) {
                if (cell != nullptr && cell->cell_type != CellType::Input) {
                    cell->polarization = dis(gen);
                }
            }
        }
    }
}