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

    /*!
     * \typedef std::tuple<int, int, long double> PolarizationValue;
     * \brief tuple data structure containing (row, col) index and polarization information
     */
    typedef std::tuple<int, int, long double> PolarizationValue;

    /*!
     * \typedef std::vector<std::tuple<int, int, long double>> PolarizationList;
     * \brief vector of PolarizationValue representing a group of polarization states
     */
    typedef std::vector<PolarizationValue> PolarizationList;

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
         * \typedef std::vector<std::vector<std::shared_ptr<QCACell>>>::reverse_iterator ReverseRowIterator;
         * \brief reverse row iterator type for the 2d array of QCACircuit containing QCACells
         */
        typedef std::vector<std::vector<std::shared_ptr<QCACell>>>::reverse_iterator ReverseRowIterator;

        /*!
         * \typedef std::vector<std::vector<std::shared_ptr<QCACell>>>::const_reverse_iterator ReverseRowConstIterator;
         * \brief const reverse row iterator type for the 2d array of QCACircuit containing QCACells
         */
        typedef std::vector<std::vector<std::shared_ptr<QCACell>>>::const_reverse_iterator ReverseRowConstIterator;

        /*!
         * \typedef std::vector<std::shared_ptr<QCACell>>::reverse_iterator ReverseColIterator;
         * \brief reverse column iterator type for the 2d array of QCACircuit containing QCACells
         */
        typedef std::vector<std::shared_ptr<QCACell>>::reverse_iterator ReverseColIterator;

        /*!
         * \typedef std::vector<std::shared_ptr<QCACell>>::const_reverse_iterator ReverseColConstIterator;
         * \brief const reverse column iterator type for the 2d array of QCACircuit containing QCACells
         */
        typedef std::vector<std::shared_ptr<QCACell>>::const_reverse_iterator ReverseColConstIterator;

        /*!
         * \typedef std::vector<std::shared_ptr<QCACell>>::iterator BfsIterator;
         * \brief breath first search iterator for qca circuit
         */
        typedef std::vector<std::shared_ptr<QCACell>>::iterator BfsIterator;

        /*!
         * \typedef std::vector<std::shared_ptr<QCACell>>::const_iterator BfsConstIterator;
         * \brief breath first search const iterator for qca circuit
         */
        typedef std::vector<std::shared_ptr<QCACell>>::const_iterator BfsConstIterator;

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

        /*!
         * \fn inline ReverseRowIterator row_rbegin()
         * \brief reverse row begin iterator
         */
        inline ReverseRowIterator row_rbegin() {
            return cells.rbegin();
        }

        /*!
         * \fn inline ReverseRowIterator row_rend()
         * \brief reverse row end iterator
         */
        inline ReverseRowIterator row_rend() {
            return cells.rend();
        }

        /*!
         * \fn inline ReverseRowConstIterator row_rbegin()const
         * \brief reverse row begin iterator, const version
         */
        inline ReverseRowConstIterator row_rbegin()const {
            return cells.crbegin();
        }

        /*!
         * \fn inline ReverseRowConstIterator row_rend()const
         * \brief reverse row end iterator, const version
         */
        inline ReverseRowConstIterator row_rend()const {
            return cells.crend();
        }

        /*!
         * \fn inline ReverseColterator col_rbegin(const ReverseRowIterator &it)
         * \brief reverse column begin iterator
         * \param it row iterator representing the current queried row
         */
        inline ReverseColIterator col_rbegin(const ReverseRowIterator &it) {
            return it->rbegin();
        }

        /*!
         * \fn inline ReverseColterator col_rend(const ReverseRowIterator &it)
         * \brief reverse column end iterator
         * \param it row iterator representing the current queried row
         */
        inline ReverseColIterator col_rend(const ReverseRowIterator &it) {
            return it->rend();
        }

        /*!
         * \fn inline ReverseColConstIterator col_rbegin(const ReverseRowConstIterator &it) const
         * \brief reverse column begin iterator, const version
         * \param it row iterator representing the current queried row
         */
        inline ReverseColConstIterator col_rbegin(const ReverseRowConstIterator &it)const {
            return it->crbegin();
        }

        /*!
         * \fn inline ReverseColConstIterator col_end(const ReverseRowConstIterator &it) const
         * \brief reverse column end iterator, const version
         * \param it row iterator representing the current queried row
         */
        inline ReverseColConstIterator col_rend(const ReverseRowConstIterator &it)const {
            return it->crend();
        }

        /*!
         * \fn inline BfsIterator bfs_begin()
         * \brief bfs iterator representing the beginning
         */
        inline BfsIterator bfs_begin() {
            return cells_in_bfs.begin();
        }

        /*!
         * \fn inline BfsIterator bfs_end()
         * \brief bfs iterator representing the ending
         */
        inline BfsIterator bfs_end() {
            return cells_in_bfs.end();
        }

        /*!
         * \fn inline BfsConstIterator bfs_begin()
         * \brief const bfs iterator representing the beginning
         */
        inline BfsConstIterator bfs_begin() const {
            return cells_in_bfs.begin();
        }

        /*!
         * \fn inline BfsConstIterator bfs_end()
         * \brief const bfs iterator representing the ending
         */
        inline BfsConstIterator bfs_end() const {
            return cells_in_bfs.end();
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

        PolarizationList &&get_input_polarizations() const;

        PolarizationList &&get_output_polarizations() const;
        /////////////////////////////////////////////////////////////////////////////////////////////

    private:
        std::vector<std::vector<std::shared_ptr<QCACell>>>  cells;//!< the internal data structure for circuit layout
        std::vector<std::shared_ptr<QCACell>> cells_in_bfs;

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
