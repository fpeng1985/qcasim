//
// Created by Administrator on 2016/10/13.
//

#include "sim_engine.h"

namespace hfut {

    SimEngine::TEMPERATURE = 1;

    SimEngine::SimEngine(QCACircuit *c, const std::vector<Polarization> &in_p) :circuit(c) {
        input_p.assign(in_p.begin(), in_p.end());

    }

    void SimEngine::initialize_output_polarization() {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(-1, 1);

        for (auto &line : circuit->cells) {
            for (auto &cell : line) {
                if (cell != nullptr && cell->cell_type != CellType::Input) {
                    cell->polarization = dis(gen);
                }
            }
        }
    }

    void SimEngine::generate_neighbour() {

    }

}
