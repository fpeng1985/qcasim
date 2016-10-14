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

    typedef std::map<std::pair<std::size_t, std::size_t>, long double> Polarization;

    class SimEngine {
    public:
        SimEngine(std::shared_ptr<QCACircuit> circuit, const Polarization &input_p, const Polarization &output_p=Polarization());

        void generate_neighbour_polarization();
        long double accept_circuit_polarization();
        void run_simulation();

        static const long double QCA_TEMPERATURE;
        static const long double RI;
        static const long double BOLTZMANN;

        static const long double EK1;
        static const long double EK2;
        static const long double EK3;

    private:
        long double compute_new_polarization_from_old(const Polarization &old_pola, Polarization &new_pola);
        long double compute_acceptance_probability(long double circuit_diff, long double neighbour_diff);

        std::shared_ptr<QCACircuit> circuit;

        Polarization output_p;
        long double  output_diff;

        Polarization neighbour_p;

        Polarization best_p;
        long double  best_diff;

        long double sa_temp;
        long double cooling_rate;
        long double terminate_temp;
        long double convergence_factor;
    };

}

#endif //QCASIM_SIMENGINE_H
