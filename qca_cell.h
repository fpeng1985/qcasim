//
// Created by Administrator on 2016/10/13.
//

#ifndef QCASIM_QCACELL_H
#define QCASIM_QCACELL_H

#include <cstdlib>

namespace hfut {

    enum class CellType {Input, Output, Normal};

    struct QCACell {
        QCACell(std::size_t i=0, std::size_t j=0, long double p=1, CellType t=CellType::Normal);

        std::size_t x_index;
        std::size_t y_index;
        long double polarization;
        CellType cell_type;
    };

}

#endif //QCASIM_QCACELL_H
