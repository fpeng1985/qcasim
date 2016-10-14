//
// Created by Administrator on 2016/10/13.
//

#include "sim_engine.h"

#include <cmath>

namespace hfut {

    using namespace std;

    SimEngine::QCATEMPERATURE = 1;
    SimEngine::RI = 3.8e-10;
    SimEngine::BOLTZMANN = 1.38e-23;

    SimEngine::EK1 = 2.3637e-22;
    SimEngine::EK2 = 7.68e-24;
    SimEngine::EK3 = -5.1638e-23;

    SimEngine::SimEngine(shared_ptr<QCACircuit> circuit, const Polarization &input_p) {
        //set simulated annealing algorithm parameters
        sa_temp = 1000;
        terminate_temp = 0.01;
        cooling_rate = 0.96;
        convergence_factor = 1e-8;

        //initialize circuit structure and input polarizations
        this->circuit = circuit;
        this->input_p = input_p;

        //randomly generate output polarizations
        random_device rd;
        mt19937 gen(rd());
        uniform_real_distribution<> dis(-1, 1);

        for (size_t i=0; i<circuit->cells.size(); ++i) {
            for (size_t j=0; j<circuit->cells[i].size(); ++j) {
                if (cell != nullptr && cell->cell_type != CellType::Input) {
                    cell->polarization = dis(gen);
                    output_p.insert(make_pair(make_pair(i, j), cell->polarization));
                }
            }
        }

        best_p = output_p;
    }

    void SimEngine::generate_neighbour_polarization() {
        //initialize neighbour to the last iteration output polarization
        neighbour_p = output_p;

        //randomly select an output cell and change its polarization
        random_device rd;
        mt19937 gen(rd());

        uniform_int_distribution<> uidis(0, output_p.size()-1);
        uniform_real_distribution<> urdis(-1, 1);

        get<2>( neighbour_p[uidis(gen)] ) = urdis(gen);
    }

    long double SimEngine::accept_circuit_polarization() {
        //compute new polarization from current circuit polarization
        Polarization pola_from_circuit;
        long double diff_for_circuit = compute_new_polarization_from_old(output_p, pola_from_circuit);

        //compute new polarization from current neighbour polarization
        Polarization pola_from_neighbour;
        long double diff_for_neighbour = compute_new_polarization_from_old(neighbour_p, pola_from_neighbour);

        long double acceptance_prob = compute_acceptance_probability(diff_for_circuit, diff_for_neighbour);

        //judge the acceptance of the neighbour
        random_device rd;
        mt19937 gen(rd());
        uniform_real_distribution<> dis(0, 1);

        bool accept = acceptance_prob > dis(gen);
        if (accept) {
            output_p = neighbour_p;

            for (auto it=output_p.begin(); it!=output_p.end(); ++it) {
                const size_t &xidx = it->first.first;
                const size_t &yidx = it->first.second;
                circuit->cells[xidx][yidx] = it->second;
            }

            if (diff_for_neighbour < diff_for_circuit) {
                best_p = neighbour_p;
            }
        }

        //cooling the sa algorithm temperature
        sa_temp *= cooling_rate;

        if (accept) {
            return diff_for_circuit;
        } else {
            return diff_for_neighbour;
        }
    }

    long double SimEngine::compute_new_polarization_from_old(const Polarization &old_pola, Polarization &new_pola) {
        long double diff = 0;

        for (Polarization::iterator it=old_pola.begin(); it!=old_pola.end(); ++it) {
            int xidx = it->first->first;
            int yidx = it->first.second;
            const long double &cur_pola_val = it->second;

            QCACell *curcell;
            long double sigma = 0;

            //type 1 neighbour cells
            curcell = circuit->get_cell(xidx+1, yidx);
            if (curcell != nullptr) {
                sigma += EK1 * curcell->polarization;
            }

            curcell = circuit->get_cell(xidx-1, yidx);
            if (curcell != nullptr) {
                sigma += EK1 * curcell->polarization;
            }

            curcell = circuit->get_cell(xidx, yidx+1);
            if (curcell != nullptr) {
                sigma += EK1 * curcell->polarization;
            }

            curcell = circuit->get_cell(xidx, yidx-1);
            if (curcell != nullptr) {
                sigma += EK1 * curcell->polarization;
            }

            //type 2 neighbour cells
            curcell = circuit->get_cell(xidx+2, yidx);
            if (curcell != nullptr) {
                sigma += EK2 * curcell->polarization;
            }

            curcell = circuit->get_cell(xidx-2, yidx);
            if (curcell != nullptr) {
                sigma += EK2 * curcell->polarization;
            }

            curcell = circuit->get_cell(xidx, yidx+2);
            if (curcell != nullptr) {
                sigma += EK2 * curcell->polarization;
            }

            curcell = circuit->get_cell(xidx, yidx-2);
            if (curcell != nullptr) {
                sigma += EK2 * curcell->polarization;
            }

            //type 3 neighbour cells
            curcell = circuit->get_cell(xidx+1, yidx+1);
            if (curcell != nullptr) {
                sigma += EK3 * curcell->polarization;
            }

            curcell = circuit->get_cell(xidx+1, yidx-1);
            if (curcell != nullptr) {
                sigma += EK3 * curcell->polarization;
            }

            curcell = circuit->get_cell(xidx-1, yidx+1);
            if (curcell != nullptr) {
                sigma += EK3 * curcell->polarization;
            }

            curcell = circuit->get_cell(xidx-1, yidx-1);
            if (curcell != nullptr) {
                sigma += EK3 * curcell->polarization;
            }

            //sqrt of RI and sigma
            long double tmp = sqrt(4*RI*RI + sigma*sigma);
            //new polarization value
            long double new_pola_val = sigma/tmp * tanh(tmp/2*BOLTZMANN*QCA_TEMPERATURE);

            //update pola_from_neighbour and target_for_neighbour
            new_pola.insert(make_pair(make_pair(xidx, yidx), new_pola_val));
            diff += pow(new_pola_val-cur_pola_val, 2);
        }

        return sqrt(diff);
    }

    long double compute_acceptance_probability(long double circuit_diff, long double neighbour_diff) {
        long double probability = 1;

        if ( neighbour_diff >= circuit_diff ) {
            probability = exp((circuit_diff-neighbour_diff) / sa_temp);
        }

        return probability;
    }
}
