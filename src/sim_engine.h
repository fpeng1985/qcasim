/*!
 * \file sim_engine.h
 * \author Peng Fei
 * \date 2016/10/13
 * \brief SimEngine definition
 */

#ifndef QCASIM_SIMENGINE_H
#define QCASIM_SIMENGINE_H

#include "qca_circuit.h"

#include <fstream>
#include <vector>
#include <map>
#include <memory>

namespace hfut {
        /*!
         * \typedef std::map<std::pair<int, int>, long double> Polarization
         * \brief data structure mapping from (row, col) index to polarization
         */
    typedef std::map<std::pair<int, int>, long double> Polarization;

    //! the simulated anealing simulation engine class
    class SimEngine {
    public:
        //! SimEngine constructor setting the sa parameters etc
        SimEngine();

        /*!
         * \fn void set_circuit(std::shared_ptr<QCACircuit> circuit)
         * \brief set the QCA circuit to be simulated
         * \param circuit a pointer to a QCACircuit
         */
        void set_circuit(std::shared_ptr<QCACircuit> circuit);

        /*!
         * \fn void run_simulation(const Polarization &input_p)
         * \brief the main routine to simuate the circuit
         * \param input_o a map of Polarization representing the circuit input
         */
        void run_simulation(const Polarization &input_p);

    private:
        /////////////////////////////////////////////////////////////////////////////////////////////
        //////////////////////////high level routine in sa algorithm/////////////////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////////

        /*!
         * \fn void set_input_polarization(const Polarization &pola) const
         * \brief set the circuit input in polarization format
         * \param pola a map of Polarization representing the circuit input
         */
        void set_input_polarization(const Polarization &pola) const;

        /*!
         * \fn void set_output_polarization_randomly() const
         * \brief set the circuit's normal and ouput cell's polarization randomly
         */
        void set_non_input_polarization_randomly() const;

        /*!
         * \fn void setup_runtime_states()
         * \brief set the simulated anealing algorithm runtime variables
         */
        void setup_runtime_states();

        /*!
         * \fn void neighbour()
         * \brief the main routine in sa algorithm, generating a neighbour state
         */
        void neighbour();

        /*!
         * \fn void accept()
         * \brief the main routine in sa algorithm, judge the acceptance of the neighbour states
         */
        void accept();
        /////////////////////////////////////////////////////////////////////////////////////////////

        /////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////low level routin in sa algorithm/////////////////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////////

        long double compute_polarization_energy(const Polarization &output_p) const;

        /*!
         * \fn long double compute_polarization_energy(const Polarization &old_pola) const
         * \brief compute the polarization energy(difference between the inputed polarization to it's synthesised counterpart)
         * \param old_pola the inputed Polarization
         */
        long double compute_polarization_difference(const Polarization &old_pola) const;

        /*!
         * \fn bool compare_energy(long double circuit_diff, long double neighbour_diff) const
         * \brief compare two circuit's energy value and deside wether to accept the neighbour state
         * \param circuit_energy the energy of the current Polarization
         * \param neighbour_energy the energy of the neighbour Polarization
         */
        bool compare_energy(long double circuit_energy, long double neighbour_energy) const;
        /////////////////////////////////////////////////////////////////////////////////////////////

        //simulated anealing algorithm parameters
        long double sa_temp;//!< the simulated anealing algorithm parameter representing current temperature
        long double cooling_rate; //!< the simulated anealing algorithm parameter representing cooling rate
        long double terminate_temp;//!< the simulated anealing algorithm parameter representing terminating temperature
        long double convergence_factor;//!< when the current energy is smaller than the convergence factor, the algorithm will be terminated

        std::shared_ptr<QCACircuit> circuit;//!< the pointer to the QCA circuit to be simulated

        //run time states
        Polarization input_p;//!< input cells's polarizations
        Polarization output_p;//!< normal and output cells' polarizations
        long double  output_energy;//!< the energy of the non input cells
        Polarization neighbour_p;//!< the neighbour polarization generated by neighbour method
        long double  neighbour_energy;//!< the energy of the neighbour polarization
        Polarization best_p;//!< the currently best polarization
        long double  best_diff;//!< the energy of the currently best polarization

#ifndef NDBUG
        std::ofstream fs;//!< debug infomation file stream
#endif
    };

}

#endif //QCASIM_SIMENGINE_H
