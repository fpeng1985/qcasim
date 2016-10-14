//
// Created by Administrator on 2016/10/14.
//

#ifndef QCASIM_SIMMANAGER_H
#define QCASIM_SIMMANAGER_H

#include "qca_circuit.h"
#include "sim_engine.h"

#include <vector>
#include <string>
#include <memory>

namespace hfut {

    class SimManager {
    public:
        void load_benchmark(const std::string &path);

    private:
        std::shared_ptr<SimEngine> engine;
        std::shared_ptr<QCACircuit> circuit;

        std::vector<QCACircuit::CircuitStructure> structures;

    };

}

#endif //QCASIM_SIMMANAGER_H
