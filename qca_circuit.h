//
// Created by Administrator on 2016/10/13.
//

#ifndef QCASIM_QCACIRCUIT_H
#define QCASIM_QCACIRCUIT_H

#include <vector>
#include <memory>

#include "qca_cell.h"

namespace hfut {

    class QCACircuit {
    public:
        void populate_cells(const std::vector<std::vector<int>> &cell_structure_matrix);
        std::shared_ptr<QCACell> get_cell(int i, int j);

        void clear();

    private:
        std::vector<std::vector<std::shared_ptr<QCACell>>>  cells;

        friend class SimEngine;
    };

}

#endif //QCASIM_QCACIRCUIT_H
