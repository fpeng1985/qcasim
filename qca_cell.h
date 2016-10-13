//
// Created by Administrator on 2016/10/13.
//

#ifndef QCASIM_QCACELL_H
#define QCASIM_QCACELL_H


namespace hfut {
    enum class CellType {Input, Output, Normal};

    struct QCACell {
        QCACell(int i=0, int j=0, double p=1, CellType t=CellType::Normal);

        int x_index;
        int y_index;
        double polarization;
        CellType cell_type;
    };

}

#endif //QCASIM_QCACELL_H
