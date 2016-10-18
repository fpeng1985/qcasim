//
// Created by Administrator on 2016/10/13.
//

#include "qca_cell.h"

namespace hfut {

    using namespace std;

    QCACell::QCACell(size_t r, size_t c, long double p, CellType t) : r_index(r), c_index(c), polarization(p), cell_type(t) {}

}
