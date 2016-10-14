//
// Created by Administrator on 2016/10/13.
//

#include "sim_engine.h"

namespace hfut {

    SimEngine::QCATEMPERATURE = 1;

    SimEngine::SimEngine(QCACircuit *circuit, const std::vector<Polarization> &input_p) {
        //set simulated annealing algorithm parameters
        sa_temp = 1000;
        terminate_temp = 0.01;
        cooling_rate = 0.96;
        convergence_factor = 1e-8;

        //initialize circuit structure and input polarizations
        this->circuit = circuit;
        this->input_p.assign(input_p.begin(), input_p.end());

        //randomly generate output polarizations
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(-1, 1);

        for (size_t i=0; i<circuit->cells.size(); ++i) {
            for (size_t j=0; j<circuit->cells[i].size(); ++j) {
                if (cell != nullptr && cell->cell_type != CellType::Input) {
                    cell->polarization = dis(gen);
                    output_p.push_back(std::make_tuple(i, j, cell->polarization));
                }
            }
        }
    }

    /*
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
     */

    void SimEngine::generate_neighbour_polarization() {

    }

}
