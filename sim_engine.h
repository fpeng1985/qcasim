//
// Created by Administrator on 2016/10/13.
//

#ifndef QCASIM_SIMENGINE_H
#define QCASIM_SIMENGINE_H

#include "qca_circuit.h"

#include <vector>
#include <map>
#include <memory>

namespace hfut {

    typedef std::map<std::pair<size_t, size_t>, double> Polarization;

    class SimEngine {

    public:
        SimEngine(std::shared_ptr<QCACircuit> circuit, const Polarization &input_p);

        void generate_neighbour_polarization();
        long double accept_circuit_polarization();

    private:
        static const long double QCA_TEMPERATURE;
        static const long double RI;
        static const long double BOLTZMANN;

        static const long double EK1;
        static const long double EK2;
        static const long double EK3;

        long double compute_new_polarization_from_old(const Polarization &old_pola, Polarization &new_pola);
        long double compute_acceptance_probability(long double circuit_diff, long double neighbour_diff);

        std::shared_ptr<QCACircuit> circuit;
        Polarization input_p;
        Polarization output_p;
        Polarization neighbour_p;
        Polarization best_p;

        long double sa_temp;
        long double terminate_temp;
        long double cooling_rate;
        long double convergence_factor;
    };

}

#endif //QCASIM_SIMENGINE_H
