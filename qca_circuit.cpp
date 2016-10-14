//
// Created by Administrator on 2016/10/13.
//

#include "qca_circuit.h"

#include <iostream>

namespace hfut {

    using namespace std;

    void QCACircuit::populate_cells(const vector<vector<int>> &cell_structure_matrix) {
        for (size_t i=0; i<cell_structure_matrix.size(); ++i) {
            cells.push_back(vector<shared_ptr<QCACell>>());

            for (size_t j=0; j<cell_structure_matrix[i].size(); ++j) {
                switch (cell_structure_matrix[i][j]) {
                    case 0:
                        cells[i].push_back(nullptr);
                        break;
                    case 1:
                        cells[i].push_back(make_shared<QCACell>(i, j, 0, CellType::Normal));
                        break;
                    case -1:
                        cells[i].push_back(make_shared<QCACell>(i, j, 0, CellType::Input));
                        break;
                    case -2:
                        cells[i].push_back(make_shared<QCACell>(i, j, 0, CellType::Output));
                        break;
                    default:
                        cerr << "Input format error!" << endl;
                        exit(-1);
                }
            }
        }
    }

    shared_ptr<QCACell> QCACircuit::get_cell(int i, int j) {
        if (i>=0 && i<cells.size() && j>=0 && j<cells[i].size()) {
            return cells[i][j];
        } else {
            return nullptr;
        }
    }

    void QCACircuit::clear() {
        cells.clear();
    }

}