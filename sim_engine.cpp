//
// Created by Administrator on 2016/10/13.
//

#include "sim_engine.h"

#include <iostream>
#include <fstream>
#include <algorithm>
#include <cstdlib>
#include <cassert>
#include <ctime>
#include <cmath>

namespace hfut {

    using namespace std;

    const long double SimEngine::QCA_TEMPERATURE = 1.0L;
    const long double SimEngine::RI = 3.8e-10L;
    const long double SimEngine::BOLTZMANN = 1.38e-23L;

    const long double SimEngine::EK1 = 2.3637e-22L;
    const long double SimEngine::EK2 = 7.68e-24L;
    const long double SimEngine::EK3 = -5.1638e-23L;

    SimEngine::SimEngine(shared_ptr<QCACircuit> circuit, const Polarization &input_p, const Polarization &output_p) {
        //set simulated annealing algorithm parameters
        sa_temp = 1000;
        terminate_temp = 0.01;
        cooling_rate = 0.9999;
        convergence_factor = 1e-8;

        //initialize circuit structure
        this->circuit = circuit;

        //initialize input polarizations
        for (auto &it : input_p) {
            auto &xidx = it.first.first;
            auto &yidx = it.first.second;
            circuit->cells[xidx][yidx]->polarization = it.second;
        }

        srand(time(0));

        //initialize output polarizations
        this->output_p = output_p;
        if (output_p.empty()) {
            //randomly generate output polarizations

            shared_ptr<QCACell> cell;
            for (size_t i=0; i<circuit->cells.size(); ++i) {
                for (size_t j=0; j<circuit->cells[i].size(); ++j) {
                    cell = circuit->cells[i][j];

                    if (cell != nullptr && cell->cell_type != CellType::Input) {
                        cell->polarization = (rand() * 1.0 / RAND_MAX - 0.5)*2;//[-1, 1]
                        this->output_p.insert(make_pair(make_pair(i, j), cell->polarization));
                    }
                }
            }
        } else {
            for (auto &it : output_p) {
                auto &xidx = it.first.first;
                auto &yidx = it.first.second;
                circuit->cells[xidx][yidx]->polarization = it.second;
            }
        }

        Polarization tmp_pola;
        output_diff = compute_new_polarization_from_old(this->output_p, tmp_pola);

        best_p = this->output_p;
        best_diff = output_diff;
    }

    void SimEngine::generate_neighbour_polarization() {
        //initialize neighbour to the last iteration output polarization
        neighbour_p = output_p;

        //randomly select an output cell and change its polarization
//        random_device rd;
//        mt19937 gen(rd());

//        uniform_int_distribution<> uidis(0, output_p.size()-1);
        auto it = neighbour_p.begin();
        advance(it, rand()%output_p.size());

//        uniform_real_distribution<> urdis(-1, 1);
        it->second = (rand() * 1.0 / RAND_MAX - 0.5)*2;//[-1, 1]
    }

    long double SimEngine::accept_circuit_polarization() {
        //compute new polarization from current neighbour polarization
        Polarization pola_from_neighbour;
        long double neighbour_diff = compute_new_polarization_from_old(neighbour_p, pola_from_neighbour);

        long double acceptance_prob = compute_acceptance_probability(output_diff, neighbour_diff);

        //judge the acceptance of the neighbour
//        random_device rd;
//        mt19937 gen(rd());
//        uniform_real_distribution<> dis(0, 1);

        bool accept = acceptance_prob > rand()*1.0/RAND_MAX;//[0,1]
        if (accept) {
            output_p = neighbour_p;

            for (auto it=output_p.begin(); it!=output_p.end(); ++it) {
                const size_t &xidx = it->first.first;
                const size_t &yidx = it->first.second;
                circuit->cells[xidx][yidx]->polarization = it->second;
            }

            output_diff = neighbour_diff;

            if (neighbour_diff < best_diff) {
                best_p = neighbour_p;
                best_diff = neighbour_diff;
            }
        }

        //cooling the sa algorithm temperature
        sa_temp *= cooling_rate;

        return output_diff;
    }

    void SimEngine::run_simulation() {

        int i = 0;
        long double diff = 0;

        ofstream fs(getenv("HOME") + string("/output.txt"));

        for (;;++i) {
            generate_neighbour_polarization();
            diff = accept_circuit_polarization();

            fs << i  << "th iteration, diff is " << diff <<  endl;
            fs << *circuit << endl;

            if (sa_temp < terminate_temp || diff < convergence_factor)
                break;
        }
    }

    long double SimEngine::compute_new_polarization_from_old(const Polarization &old_pola, Polarization &new_pola) {
        assert( new_pola.empty() );

        long double diff = 0;

        for (Polarization::const_iterator it=old_pola.begin(); it!=old_pola.end(); ++it) {
            int xidx = int(it->first.first);
            int yidx = int(it->first.second);
            const long double &cur_pola_val = it->second;

            shared_ptr<QCACell> curcell;
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

    long double SimEngine::compute_acceptance_probability(long double circuit_diff, long double neighbour_diff) {
        long double probability = 1;

        if ( neighbour_diff >= circuit_diff ) {
            probability = exp((circuit_diff-neighbour_diff) / sa_temp);
        }

        return probability;
    }

}
