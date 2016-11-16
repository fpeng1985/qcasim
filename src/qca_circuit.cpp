//
// Created by Administrator on 2016/10/13.
//

#include "qca_circuit.h"

#include <queue>

namespace hfut {

    using namespace std;

    void QCACircuit::populate_cells(const CircuitStructure &cell_structure_matrix) {
        clear();

        queue<shared_ptr<QCACell>> queue;
        vector<vector<int>> vflag;
        shared_ptr<QCACell> cell = nullptr;

        for (size_t r=0; r<cell_structure_matrix.size(); ++r) {
            cells.push_back(vector<shared_ptr<QCACell>>());
            vflag.push_back(vector<int>());

            for (size_t c=0; c<cell_structure_matrix[r].size(); ++c) {
                switch (cell_structure_matrix[r][c]) {
                    case 0:
                        cells[r].push_back(nullptr);
                        vflag[r].push_back(0);
                        break;
                    case 1:
                        cells[r].push_back(make_shared<QCACell>(r, c, 0, CellType::Normal));
                        vflag[r].push_back(0);
                        break;
                    case -1:
                        cell = make_shared<QCACell>(r, c, 0, CellType::Input);
                        cells[r].push_back(cell);
                        queue.push(cell);
                        vflag[r].push_back(1);
                        break;
                    case -2:
                        cells[r].push_back(make_shared<QCACell>(r, c, 0, CellType::Output));
                        vflag[r].push_back(0);
                        break;
                    default:
                        cerr << "Input format error!" << endl;
                        exit(-1);
                }

            }
        }

        while (!queue.empty()) {
            cell = queue.front();
            cells_in_bfs.push_back(cell);
            queue.pop();

            auto &ridx = cell->r_index;
            auto &cidx = cell->c_index;

            if ( (get_cell(ridx, cidx-1) != nullptr) && (vflag[ridx][cidx-1] == 0) ) {
                queue.push( cells[ridx][cidx-1] );
                vflag[ridx][cidx-1] = 1;
            }

            if ( (get_cell(ridx, cidx+1) != nullptr) && (vflag[ridx][cidx+1] == 0) ) {
                queue.push( cells[ridx][cidx+1] );
                vflag[ridx][cidx+1] = 1;
            }

            if ( (get_cell(ridx-1, cidx) != nullptr) && (vflag[ridx-1][cidx] == 0) ) {
                queue.push( cells[ridx-1][cidx] );
                vflag[ridx-1][cidx] = 1;
            }

            if ( (get_cell(ridx+1, cidx) != nullptr) && (vflag[ridx+1][cidx] == 0) ) {
                queue.push( cells[ridx+1][cidx] );
                vflag[ridx+1][cidx] = 1;
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

    PolarizationList &&QCACircuit::get_input_polarizations() const {
        PolarizationList input_p;
        shared_ptr<QCACell> cell = nullptr;

        for (auto rit=cells.begin(); rit!=cells.end(); ++rit) {
            for (auto cit=rit->begin(); cit!=rit->end(); ++cit) {
                cell = *cit;
                if (cell != nullptr && cell->cell_type == CellType::Input) {
                    input_p.push_back(make_tuple(cell->r_index, cell->c_index, cell->polarization));
                }
            }
        }

        return move(input_p);
    }

    PolarizationList &&QCACircuit::get_output_polarizations() const {
        PolarizationList input_p;
        shared_ptr<QCACell> cell = nullptr;

        for (auto rit=cells.begin(); rit!=cells.end(); ++rit) {
            for (auto cit=rit->begin(); cit!=rit->end(); ++cit) {
                cell = *cit;
                if (cell != nullptr && cell->cell_type == CellType::Output) {
                    input_p.push_back(make_tuple(cell->r_index, cell->c_index, cell->polarization));
                }
            }
        }

        return move(input_p);
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