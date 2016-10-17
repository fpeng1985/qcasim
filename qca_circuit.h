/*!
 * \file qca_circuit.h
 * \author Peng Fei
 * \date 2016/10/13
 * \brief QCACircuit definition
 */

#ifndef QCASIM_QCACIRCUIT_H
#define QCASIM_QCACIRCUIT_H

#include <iostream>
#include <vector>
#include <memory>

#include "qca_cell.h"

namespace hfut {

    class QCACircuit {
    public:
        /*!
         * \typedef std::vector<std::vector<int>> CircuitStructure;
         * \brief A temporary type for representing the QCACircuit's layout
         */
        typedef std::vector<std::vector<int>> CircuitStructure;

        /*!
         * \typedef std::vector<std::vector<std::shared_ptr<QCACell>>>::iterator RowIterator;
         * \brief row iterator type for the 2d array of QCACircuit containing QCACells
         */
        typedef std::vector<std::vector<std::shared_ptr<QCACell>>>::iterator RowIterator;

        /*!
         * \typedef std::vector<std::vector<std::shared_ptr<QCACell>>>::const_iterator RowConstIterator;
         * \brief const row iterator type for the 2d array of QCACircuit containing QCACells
         */
        typedef std::vector<std::vector<std::shared_ptr<QCACell>>>::const_iterator RowConstIterator;

        /*!
         * \typedef std::vector<std::shared_ptr<QCACell>>::iterator ColIterator;
         * \brief column iterator type for the 2d array of QCACircuit containing QCACells
         */
        typedef std::vector<std::shared_ptr<QCACell>>::iterator ColIterator;

        /*!
         * \typedef std::vector<std::shared_ptr<QCACell>>::const_iterator ColConstIterator;
         * \brief const column iterator type for the 2d array of QCACircuit containing QCACells
         */
        typedef std::vector<std::shared_ptr<QCACell>>::const_iterator ColConstIterator;

        void populate_cells(const CircuitStructure &cell_structure_matrix);
        std::shared_ptr<QCACell> get_cell(int r, int c);

        void clear();

        inline RowIterator row_begin() {
            return cells.begin();
        }
        inline RowIterator row_end() {
            return cells.end();
        }
        inline RowConstIterator row_begin()const {
            return cells.cbegin();
        }
        inline RowConstIterator row_end()const {
            return cells.cend();
        }

        inline ColIterator col_begin(const RowIterator &it) {
            return it->begin();
        }
        inline ColIterator col_end(const RowIterator &it) {
            return it->end();
        }
        inline ColConstIterator col_begin(const RowConstIterator &it)const {
            return it->cbegin();
        }
        inline ColConstIterator col_end(const RowConstIterator &it)const {
            return it->cend();
        }

    private:
        std::vector<std::vector<std::shared_ptr<QCACell>>>  cells;

        friend std::ostream & operator<<(std::ostream &out, const QCACircuit &circuit);
    };

    /*!
     * \fn std::ostream & operator<<(std::ostream &out, const QCACircuit &circuit);
     * \brief A help function used for printing the QCACircuit
     * \param out the output stream object
     * \param circuit the circuit to be printed
     */
    std::ostream & operator<<(std::ostream &out, const QCACircuit &circuit);

}

#endif //QCASIM_QCACIRCUIT_H
