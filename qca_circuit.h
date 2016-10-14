//
// Created by Administrator on 2016/10/13.
//

#ifndef QCASIM_QCACIRCUIT_H
#define QCASIM_QCACIRCUIT_H

#include <iostream>
#include <vector>
#include <memory>

#include "qca_cell.h"

namespace hfut {

    class QCACircuit {
    public:
        typedef std::vector<std::vector<std::shared_ptr<QCACell>>>::iterator RowIterator;
        typedef std::vector<std::vector<std::shared_ptr<QCACell>>>::const_iterator RowConstIterator;
        typedef std::vector<std::shared_ptr<QCACell>>::iterator ColIterator;
        typedef std::vector<std::shared_ptr<QCACell>>::const_iterator ColConstIterator;

        void populate_cells(const std::vector<std::vector<int>> &cell_structure_matrix);
        std::shared_ptr<QCACell> get_cell(int i, int j);

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

//        std::size_t row_size;
//        std::size_t col_size;



    private:
        std::vector<std::vector<std::shared_ptr<QCACell>>>  cells;

        friend class SimEngine;
        friend std::ostream & operator<<(std::ostream &out, const QCACircuit &circuit);
    };

    std::ostream & operator<<(std::ostream &out, const QCACircuit &circuit);

}

#endif //QCASIM_QCACIRCUIT_H
