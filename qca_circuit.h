//
// Created by Administrator on 2016/10/13.
//

#ifndef QCASIM_QCACIRCUIT_H
#define QCASIM_QCACIRCUIT_H

#include <vector>
#include <string>

#include "qca_cell.h"

namespace hfut {

    using std::vector;
    using std::string;

    class QCACircuit {
    public:
        QCACircuit();

        void populate_cells(const vector<vector<int>> &cell_structure_matrix);
        QCACell *get_cell(size_t i, size_t j);

    private:
        vector<vector<QCACell*>>  cells;
    };

}

#endif //QCASIM_QCACIRCUIT_H
