//
// Created by Administrator on 2016/10/13.
//

#ifndef QCASIM_SIMENGINE_H
#define QCASIM_SIMENGINE_H

#include "qca_circuit.h"

#include <vector>
#include <tuple>
#include <memory>

namespace hfut {

    typedef std::tuple<size_t, size_t, double> Polarization;

    class SimEngine {

    public:
        SimEngine(std::shared_ptr<QCACircuit> circuit, const std::vector<Polarization> &input_p);

        //void initialize_output_polarization();
        void generate_neighbour_polarization();

    private:
        static const double QCA_TEMPERATURE;

        std::shared_ptr<QCACircuit> circuit;
        std::vector<Polarization> input_p;
        std::vector<Polarization> output_p;
        std::vector<Polarization> neighbour_p;
        std::vector<Polarization> best_p;

        double sa_temp;
        double terminate_temp;
        double cooling_rate;
        double convergence_factor;
    };

}

#endif //QCASIM_SIMENGINE_H
