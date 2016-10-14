//
// Created by Administrator on 2016/10/13.
//

#ifndef QCASIM_SIMENGINE_H
#define QCASIM_SIMENGINE_H

#include "qca_circuit.h"

#include <vector>
#include <tuple>

namespace hfut {

    typedef std::tuple<size_t, size_t, double> Polarization;

    class SimEngine {

    public:
        SimEngine(QCACircuit *circuit, const std::vector<Polarization> &input_p);

        void initialize_output_polarization();
        void generate_neighbour_polarization();

    private:
        static const double TEMPERATURE;

        QCACircuit *circuit;
        std::vector<Polarization> input_p;
        std::vector<Polarization> output_p;
        std::vector<Polarization> neighbour_p;
        std::vector<Polarization> best_p;
    };

}

#endif //QCASIM_SIMENGINE_H
