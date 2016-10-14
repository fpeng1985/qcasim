//
// Created by Administrator on 2016/10/13.
//

#include "qca_cell.h"

namespace hfut {

    using namespace std;

    QCACell::QCACell(size_t i, size_t j, long double p, CellType t) : x_index(i), y_index(j), polarization(p), cell_type(t) {}

}
