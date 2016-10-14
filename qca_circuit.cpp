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
        if (i>=0 && i<int(cells.size()) && j>=0 && j<int(cells[i].size())) {
            return cells[i][j];
        } else {
            return nullptr;
        }
    }

    void QCACircuit::clear() {
        cells.clear();
    }

    ostream & operator<<(ostream &out, const QCACircuit &circuit) {
        out << "===============================" << endl;
        for (auto rit=circuit.row_begin(); rit!=circuit.row_end(); ++rit) {
            for (auto cit=circuit.col_begin(rit); cit!=circuit.col_end(rit); ++cit) {
                shared_ptr<QCACell> cell = *cit;
                if (cell == nullptr) {
                    out << "X\t";
                } else {
                    out << cell->polarization << "\t";
                }
            }
            out << endl;
        }
        out << "===============================" << endl;
        return out;
    }

}