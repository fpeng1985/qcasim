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
#include <string>
#include <memory>
#include <cassert>

namespace hfut {

    //! the simulation management class
    class SimManager {
    public:

        /*!
         * \typedef std::tuple<int, int, int> QCATruthValue;
         * \brief truth value type definition
         */
        typedef std::tuple<int, int, int> QCATruthValue;

        /*!
         * \typedef std::vector<QCATruthValue> QCATruthValueList;
         * \brief vector of QCATruthValue representing a group of input or output truth value
         */
        typedef std::vector<QCATruthValue> QCATruthValueList;

        /*!
         * \typedef std::map< QCATruthValueList, QCATruthValueList > QCATruthTable;
         * \brief truth table type definition
         */
        typedef std::map< QCATruthValueList, QCATruthValueList > QCATruthTable;

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
            assert( (0.5 < pola_val && pola_val <= 1) || (-1 <= pola_val && pola_val < -0.5) );

            if (pola_val > 0) {
                return 1;
            } else {
                return 0;
            }
        }

    private:
        std::shared_ptr<SimEngine>  engine;//!< pointer to the simulation engine
        std::shared_ptr<QCACircuit> circuit;//!< pointer to the QCA circuit

        std::vector<QCACircuit::CircuitStructure> structures;//!< all the circuit structures generated from one benchmark
        std::vector<QCATruthTable> tables;//! all the simulation results in logic form(instead of polarization form)

        unsigned int input_size;

    public:
        //! helper class for generating all the combinations of interger sequence [0,1,...,n-1]
        class CombinationGenerator {
        public:
         /*!
          * \fn void generate_combination(int n, std::vector<std::vector<int>> &combinations)
          * \brief the main generation routine
          * \param n the number of integers
          * \param combinations the data structure storing the generated results
          */
            void generate_combination(int n, std::vector<std::vector<int>> &combinations);

        private:
            /*!
             * \fn void combine(const std::vector<int> &a, int n, std::vector<int> &b, int m, const int M, std::vector<std::vector<int>> &combinations)
             * \brief the actual generation algorithm implemented in a recursive manner
             * \param a the vector representing the integers [0,1,...,n-1]
             * \param n the number of the current vector length in a recursive call
             * \param b the vector stroring the results in a recursive call
             * \param m the nubmer to be selected in a recursive call
             * \param M the most top level number to be selected from the vector a
             * \param combinations the data structure storing the generated results
             */
            void combine(const std::vector<int> &a, int n, std::vector<int> &b, int m, const int M, std::vector<std::vector<int>> &combinations);
        };
    };
}

#endif //QCASIM_SIMMANAGER_H
