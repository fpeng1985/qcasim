//
// Created by Administrator on 2016/10/13.
//

#include "qca_circuit.h"

namespace hfut {

    using namespace std;

    void QCACircuit::populate_cells(const CircuitStructure &cell_structure_matrix) {
        clear();

        for (size_t r=0; r<cell_structure_matrix.size(); ++r) {
            cells.push_back(vector<shared_ptr<QCACell>>());

            for (size_t c=0; c<cell_structure_matrix[r].size(); ++c) {
                switch (cell_structure_matrix[r][c]) {
                    case 0:
                        cells[r].push_back(nullptr);
                        break;
                    case 1:
                        cells[r].push_back(make_shared<QCACell>(r, c, 0, CellType::Normal));
                        break;
                    case -1:
                        cells[r].push_back(make_shared<QCACell>(r, c, 0, CellType::Input));
                        break;
                    case -2:
                        cells[r].push_back(make_shared<QCACell>(r, c, 0, CellType::Output));
                        break;
                    default:
                        cerr << "Input format error!" << endl;
                        exit(-1);
                }
            }
        }
    }

    shared_ptr<QCACell> QCACircuit::get_cell(int r, int c) {
        if (r>=0 && r<(cells.size()) && c>=0 && c<(cells[r].size())) {
            return cells[r][c];
        } else {
            return nullptr;
        }
    }

    void QCACircuit::clear() {
        cells.clear();
    }

    ostream & operator<<(ostream &out, const QCACircuit &circuit) {
        out << "=========================================================================================" << endl;
        for (auto rit=circuit.row_begin(); rit!=circuit.row_end(); ++rit) {
            for (auto cit=circuit.col_begin(rit); cit!=circuit.col_end(rit); ++cit) {
                shared_ptr<QCACell> cell = *cit;
                if (cell == nullptr) {
                    out << "X\t\t";
                } else {
                    out << cell->polarization << "\t\t";
                }
            }
            out << endl;
        }
        out << "=========================================================================================" << endl;
        return out;
    }

}