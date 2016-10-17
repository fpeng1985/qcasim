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
        SimManager();

        void load_benchmark(const std::string &path);

    private:
        std::shared_ptr<SimEngine>  engine;
        std::shared_ptr<QCACircuit> circuit;

        std::vector<QCACircuit::CircuitStructure> structures;

    public:
        class CombinationGenerator {
        public:
            void generate_combination(int n, std::vector<std::vector<int>> &combinations);

        private:
            void combine(const std::vector<int> &a, int n, std::vector<int> &b, int m, const int M, std::vector<std::vector<int>> &combinations);
        };

    };
}

#endif //QCASIM_SIMMANAGER_H
