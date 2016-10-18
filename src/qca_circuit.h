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

    //! class representing the QCA circuit
    class QCACircuit {
    public:
        /////////////////////////////////////////////////////////////////////////////////////////////
        /*!
         * \typedef std::vector<std::vector<std::shared_ptr<QCACell>>>::iterator RowIterator
         * \brief row iterator type for the 2d array of QCACircuit containing QCACells
         */
        typedef std::vector<std::vector<std::shared_ptr<QCACell>>>::iterator RowIterator;

        /*!
         * \typedef std::vector<std::vector<std::shared_ptr<QCACell>>>::const_iterator RowConstIterator
         * \brief const row iterator type for the 2d array of QCACircuit containing QCACells
         */
        typedef std::vector<std::vector<std::shared_ptr<QCACell>>>::const_iterator RowConstIterator;

        /*!
         * \typedef std::vector<std::shared_ptr<QCACell>>::iterator ColIterator
         * \brief column iterator type for the 2d array of QCACircuit containing QCACells
         */
        typedef std::vector<std::shared_ptr<QCACell>>::iterator ColIterator;

        /*!
         * \typedef std::vector<std::shared_ptr<QCACell>>::const_iterator ColConstIterator
         * \brief const column iterator type for the 2d array of QCACircuit containing QCACells
         */
        typedef std::vector<std::shared_ptr<QCACell>>::const_iterator ColConstIterator;

        /*!
         * \fn inline RowIterator row_begin()
         * \brief row begin iterator
         */
        inline RowIterator row_begin() {
            return cells.begin();
        }

        /*!
         * \fn inline RowIterator row_end()
         * \brief row end iterator
         */
        inline RowIterator row_end() {
            return cells.end();
        }

        /*!
         * \fn inline RowConstIterator row_begin()const
         * \brief row begin iterator, const version
         */
        inline RowConstIterator row_begin()const {
            return cells.cbegin();
        }

        /*!
         * \fn inline RowConstIterator row_end()const
         * \brief row end iterator, const version
         */
        inline RowConstIterator row_end()const {
            return cells.cend();
        }

        /*!
         * \fn inline Colterator col_begin(const RowIterator &it)
         * \brief column begin iterator
         * \param it row iterator representing the current queried row
         */
        inline ColIterator col_begin(const RowIterator &it) {
            return it->begin();
        }

        /*!
         * \fn inline Colterator col_end(const RowIterator &it)
         * \brief column end iterator
         * \param it row iterator representing the current queried row
         */
        inline ColIterator col_end(const RowIterator &it) {
            return it->end();
        }

        /*!
         * \fn inline ColConstIterator col_begin(const RowConstIterator &it) const
         * \brief column begin iterator, const version
         * \param it row iterator representing the current queried row
         */
        inline ColConstIterator col_begin(const RowConstIterator &it)const {
            return it->cbegin();
        }

        /*!
         * \fn inline ColConstIterator col_end(const RowConstIterator &it) const
         * \brief column end iterator, const version
         * \param it row iterator representing the current queried row
         */
        inline ColConstIterator col_end(const RowConstIterator &it)const {
            return it->cend();
        }
        /////////////////////////////////////////////////////////////////////////////////////////////

        /////////////////////////////////////////////////////////////////////////////////////////////
        /*!
         * \typedef std::vector<std::vector<int>> CircuitStructure
         * \brief A temporary type for representing the QCACircuit's layout
         */
        typedef std::vector<std::vector<int>> CircuitStructure;

        /*!
         * \fn void populate_cells(const CircuitStructure &cell_structure_matrix)
         * \brief setup the 2d array internal data structure for circuit layout
         * \param cell_structure_matrix vector of vector of int representing circuit layout
         */
        void populate_cells(const CircuitStructure &cell_structure_matrix);

        /*!
         * \fn std::shared_ptr<QCACell> get_cell(int r, int c)
         * \brief access method for QCACell
         * \param r the row index for QCACell
         * \param c the column index for QCACell
         */
        std::shared_ptr<QCACell> get_cell(int r, int c);

        /*!
         * \fn void clear()
         * \brief clear the internal data structures
         */
        void clear();
        /////////////////////////////////////////////////////////////////////////////////////////////

    private:
        std::vector<std::vector<std::shared_ptr<QCACell>>>  cells;//!< the internal data structure for circuit layout

        friend std::ostream & operator<<(std::ostream &out, const QCACircuit &circuit);
    };

    /*!
     * \fn std::ostream & operator<<(std::ostream &out, const QCACircuit &circuit)
     * \brief A help function used for printing the QCACircuit
     * \param out the output stream object
     * \param circuit the circuit to be printed
     */
    std::ostream & operator<<(std::ostream &out, const QCACircuit &circuit);

}

#endif //QCASIM_QCACIRCUIT_H
