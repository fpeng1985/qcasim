/*!
 * \file sim_manager.h
 * \author Peng Fei
 * \date 2016/10/14
 * \brief SimManager definition
 */

#ifndef QCASIM_SIMMANAGER_H
#define QCASIM_SIMMANAGER_H

#include "qca_circuit.h"
#include "sim_engine.h"

#include <vector>
#include <tuple>
#include <map>
#include <set>
#include <string>
#include <memory>
#include <iostream>
#include <cassert>

#ifndef NDBUG
#define private   public
#define protected public
#endif

namespace hfut {

    //! helper template class for generating all the combinations of interger sequence [0,1,...,n-1]
    template<typename T>
    class CombinationGenerator {
    public:
        /*!
         * \fn void generate_combination(const std::vector<T> &a, int m, std::vector<std::vector<T>> &combinations) {
         * \brief the main generation routine
         * \param a the vector representing the actual values
         * \param m the size of the combination to be selected
         * \param combinations the data structure storing the selected combination results
         */
        void generate_combination(const std::vector<T> &a, int m, std::vector<std::vector<T>> &combinations) {
            int n = a.size();

            assert(m>=1 && m<=n);
            assert(combinations.empty());

            std::vector<int> b;
            b.resize(a.size());

            combine(a, n, b, m, m, combinations);
        }

    private:
        /*!
         * \fn void combine(const std::vector<T> &a, int n, std::vector<int> &b, int m, const int M, std::vector<std::vector<T>> &combinations) {
         * \brief the actual generation algorithm implemented in a recursive manner
         * \param a the vector representing the actual values
         * \param n the number of the current vector length in a recursive call
         * \param b the vector stroring the results in a recursive call
         * \param m the nubmer to be selected in a recursive call
         * \param M the most top level number to be selected from the vector a
         * \param combinations the data structure storing the generated results
         */
        void combine(const std::vector<T> &a, int n, std::vector<int> &b, int m, const int M, std::vector<std::vector<T>> &combinations) {
            for(int i=n; i>=m; i--) {
                b[m-1] = i - 1;
                if (m > 1)
                    combine(a, i-1, b, m-1, M, combinations);
                else {
                    std::vector<T> tmp;
                    for(int j=M-1; j>=0; j--)
                        tmp.push_back(a[b[j]]);
                    combinations.push_back(tmp);
                }
            }
        }
    };

    /*!
     * \typedef std::tuple<int, int, int> QCATruthValue;
     * \brief truth value type definition
     */
    typedef std::tuple<int, int, int> QCATruthValue;

    /*!
     * \ typedef std::set<QCATruthValue> QCATruthValueSet;
     * \brief a set of QCATruthValue representing a group of input or output truth value
     */
    typedef std::set<QCATruthValue> QCATruthValueSet;

    class QCATruthTable  {
    public:
        //! QCATruthTable default constructor
        QCATruthTable() = default;

        /*!
         * \brief QCATruthTable constructor, used only in the original circuit structure's initialization
         * \param input_idx the indices of the input
         * \param output_idx the indices of the output
         */
        QCATruthTable(const std::vector<std::pair<int, int>> &input_idx, const std::vector<std::pair<int, int>> &output_idx);

        /*!
         * \fn inline void set_value(const QCATruthValueSet &input_values, const QCATruthValueSet &output_values)
         * \brief set the input output mapping
         * \param input_values the value of the input which is meant to be in the table already
         * \param output_values the value of the output to be updated
         */
        inline void set_value(const QCATruthValueSet &input_values, const QCATruthValueSet &output_values) {
            assert(table.count(input_values) == 1);

            assert(input_values.size()  == input_size);
            assert(output_values.size() == output_size);

            table[input_values] = output_values;
        }

        inline size_t size() const {
            return table.size();
        }

        ///////////////////////////////Iterators///////////////////////////////////////////
        typedef std::map<QCATruthValueSet, QCATruthValueSet>::iterator QCATruthTableIterator;
        typedef std::map<QCATruthValueSet, QCATruthValueSet>::const_iterator QCATruthTableConstIterator;

        inline QCATruthTableIterator begin() {
            return table.begin();
        }

        inline QCATruthTableIterator end() {
            return table.end();
        }

        inline QCATruthTableConstIterator begin() const {
            return table.begin();
        }

        inline QCATruthTableConstIterator end() const {
            return table.end();
        }
        ///////////////////////////////Iterators///////////////////////////////////////////

    private:

#ifndef NDBUG
        size_t input_size;//!< the size of the input cells
        size_t output_size;//!< the size of the output cells
#endif

        std::map< QCATruthValueSet, QCATruthValueSet> table;//!< map from input values to output values

        friend std::ostream &operator<<(std::ostream &os, const QCATruthTable &truth_table);//!< friend declareation for operator<<
        friend bool operator==(const QCATruthTable &lhs, const QCATruthTable &rhs);//!< friend declaration for operator==
        friend bool operator!=(const QCATruthTable &lhs, const QCATruthTable &rhs);//!< friend declaration for operator!=
    };

    std::ostream &operator<<(std::ostream &os, const QCATruthTable &truth_table);
    bool operator==(const QCATruthTable &lhs, const QCATruthTable &rhs);
    bool operator!=(const QCATruthTable &lhs, const QCATruthTable &rhs);

    struct SimResult {
        SimResult(const QCACircuit::CircuitStructure &structure, const QCATruthTable &table,
                  bool sflag=true, bool cflag=true) : success(sflag), correct(cflag) {
            circuit_structure = structure;
            truth_table = table;
        }

        inline void update_sim_result(const QCATruthValueSet &input_values, const QCATruthValueSet &output_values) {
            truth_table.set_value(input_values, output_values);
        }

        bool success;
        bool correct;

        QCACircuit::CircuitStructure circuit_structure;
        QCATruthTable truth_table;
    };

    class SimResultGroup {
    public:
        SimResultGroup(const QCACircuit::CircuitStructure &structure, const QCATruthTable &table,
                       const std::vector<std::vector<std::pair<int, int>>> &index_combinations);

        inline double correction_ratial() const {
            return _correction_ratial;
        }

        inline void set_correction_ratial(double ratial) {
            _correction_ratial = ratial;
        }

        inline size_t size() const {
            return results.size();
        }

        ///////////////////////////////Iterators///////////////////////////////////////////
        typedef std::vector<SimResult>::iterator SimResultIterator;
        typedef std::vector<SimResult>::const_iterator SimResultConstIterator;

        inline SimResultIterator begin() {
            return results.begin();
        }

        inline SimResultIterator end() {
            return results.end();
        }

        inline SimResultConstIterator begin() const{
            return results.begin();
        }

        inline SimResultConstIterator end() const{
            return results.end();
        }
        ///////////////////////////////Iterators///////////////////////////////////////////
    private:
        std::vector<SimResult> results;
        double _correction_ratial;
    };

    //! the simulation management class
    class SimManager {
    public:
        //! SimManager constructor, create memory for simulation engine and QCA circuit
        SimManager();

        /*!
         * \fn void load_benchmark(const std::string &path)
         * \brief setup the all circuit structures to be simulated from one benchmark file
         * \param path the file path of the benchmark
         */
        void load_benchmark(const std::string &path);

        /*!
         * \fn void test_benchmark()
         * \brief for each circuit structure run its simulation and record the result into truth table data structures
         */
        void test_benchmark();

    protected:
        static inline long double convert_logic_to_polarization(int logic_val) {
            return 2*logic_val - 1;
        }

        static inline int convert_polarization_to_logic(long double pola_val) {

            if (pola_val >= 0.5 && pola_val <= 1) {
                return 1;
            } else if (pola_val <= -0.5 && pola_val >= -1){
                return 0;
            } else {
                return -1;
            }
        }

        void test_circuit_structure(SimResult &result);

        void compute_truth_table(const QCACircuit::CircuitStructure &structure, QCATruthTable &table);

    private:
        std::shared_ptr<SimEngine>  engine;//!< pointer to the simulation engine, initialized in constructor
        std::shared_ptr<QCACircuit> circuit;//!< pointer to the QCA circuit, initialized in constructor

        QCACircuit::CircuitStructure benchmark_circuit_structure;//!< original circuit structure, initialized in loading process
        QCATruthTable benchmark_truth_table;//!< original circuit's truth table, initialized in loading process

        size_t input_cell_size;//!< the number of the input cells, initialized in loading process
        size_t normal_cell_size;//!< the number of the normal cells, initialized in loading process
        size_t output_cell_size;//!< the number of the output cells, initialized in loading process

        typedef std::pair<int, int> CellIndex;
        std::vector<CellIndex> input_idx;//!< the input cells' r-c index, initialized in loading process
        std::vector<CellIndex> normal_idx;//!< the normal cells' r-c index, initialized in loading process
        std::vector<CellIndex> output_idx;//!< the output cells' r-c index, initialized in loading process

        std::map<size_t, SimResultGroup> results;//!< the simulation result data structure, computed in test process
    };
}

#endif //QCASIM_SIMMANAGER_H
