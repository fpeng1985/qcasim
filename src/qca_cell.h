/*!
 * \file qca_cell.h
 * \author Peng Fei
 * \date 2016/10/13
 * \brief QCACell definition
 */

#ifndef QCASIM_QCACELL_H
#define QCASIM_QCACELL_H

#include <cstdlib>

namespace hfut {

    /*!
     * \typedef enum class CellType {Input, Output, Normal};
     * \brief QCA cell type enum
     */
    enum class CellType {Input, Output, Normal};

    //! struct storing the QCA cell infomation
    struct QCACell {
        //! QCACell constructor
        QCACell(std::size_t r=0, std::size_t c=0, long double p=1, CellType t=CellType::Normal);

        std::size_t r_index;//!< QCACell's row index in QCACircuit
        std::size_t c_index;//!< QCACell's column index in QCACircuit
        long double polarization;//!< Polarization value for the cell
        CellType cell_type;//!< Type of the cell, i.e. Input, Output or Normal
    };

}

#endif //QCASIM_QCACELL_H
