//
// Created by Administrator on 2016/10/13.
//

#ifndef QCASIM_SIMENGINE_H
#define QCASIM_SIMENGINE_H

#include "qca_circuit.h"

#include <fstream>
#include <vector>
#include <map>
#include <memory>

namespace hfut {

    typedef std::map<std::pair<int, int>, long double> Polarization;

    class SimEngine {
    public:
        SimEngine();

        void set_circuit(std::shared_ptr<QCACircuit> circuit);
        void run_simulation(const Polarization &input_p);

    private:
        //high level routine in sa algorithm
        void set_input_polarization(const Polarization &pola) const;
        void set_output_polarization_randomly() const;
        void setup_runtime_states();
        void neighbour();
        void accept();

        //low level routin in sa algorithm
        long double compute_polarization_energy(const Polarization &old_pola) const;
        bool compare_energy(long double circuit_diff, long double neighbour_diff) const;

        //simulated anealing algorithm parameters
        long double sa_temp;
        long double cooling_rate;
        long double terminate_temp;
        long double convergence_factor;

        std::shared_ptr<QCACircuit> circuit;

        //run time states
        Polarization output_p;
        long double  output_diff;
        Polarization neighbour_p;
        long double  neighbour_diff;
        Polarization best_p;
        long double  best_diff;

#ifndef NDBUG
        std::ofstream fs;
#endif
    };

}

#endif //QCASIM_SIMENGINE_H
