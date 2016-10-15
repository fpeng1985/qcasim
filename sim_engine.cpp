//
// Created by Administrator on 2016/10/13.
//

#include "sim_engine.h"

#include <algorithm>

#include <cstdlib>
#include <ctime>
#include <cmath>
#include <cassert>

namespace hfut {

    using namespace std;

    const long double QCA_TEMPERATURE = 1.0L;
    const long double RI = 3.8e-10L;
    const long double BOLTZMANN = 1.38e-23L;

    const long double EK1 = 2.3637e-22L;
    const long double EK2 = 7.68e-24L;
    const long double EK3 = -5.1638e-23L;

    SimEngine::SimEngine() {
        srand(time(0));

        //set simulated annealing algorithm parameters
        sa_temp = 1000;
        cooling_rate = 0.9999;
        terminate_temp = 0.01;
        convergence_factor = 1e-8;

#ifndef NDBUG
        fs.open(getenv("HOME") + string("/output.txt"));
#endif
    }

    void SimEngine::set_circuit(std::shared_ptr<QCACircuit> circuit) {
        this->circuit = circuit;
    }

    void SimEngine::run_simulation(const Polarization &input_p) {
        assert(circuit != nullptr);

        set_input_polarization(input_p);
        set_output_polarization_randomly();
        setup_runtime_states();

        int i=0;
        while ( sa_temp > terminate_temp && output_diff> convergence_factor ) {
            make_anealing_iteration();

#ifndef NDBUG
            fs << i  << "th iteration, diff is " << output_diff <<  endl;
            fs << *circuit << endl;
#endif
        }
    }

    void SimEngine::set_input_polarization(const Polarization &pola) const {
        for (auto &it : pola) {
            auto &xidx = it.first.first;
            auto &yidx = it.first.second;

            assert(circuit->get_cell(xidx, yidx) != nullptr);
            assert(circuit->get_cell(xidx, yidx)->cell_type == CellType::Input);
            circuit->get_cell(xidx, yidx)->polarization = it.second;
        }
    }

    void SimEngine::set_output_polarization_randomly() const {
        shared_ptr<QCACell> cell;
        for (auto rit=circuit->row_begin(); rit!=circuit->row_end(); ++rit) {
            for (auto cit=circuit->col_begin(rit); cit!=circuit->col_end(rit); ++cit) {
                cell = *cit;
                if (cell != nullptr && cell->cell_type != CellType::Input) {
                    cell->polarization = (rand() * 1.0 / RAND_MAX - 0.5)*2;//[-1, 1]
                }
            }
        }
    }

    void SimEngine::setup_runtime_states() {
        //compute output_p
        shared_ptr<QCACell> cell;
        int i=0;
        for (auto rit=circuit->row_begin(); rit!=circuit->row_end(); ++rit, ++i) {
            int j=0;
            for (auto cit=circuit->col_begin(rit); cit!=circuit->col_end(rit); ++cit, ++j) {
                cell = *cit;
                if (cell != nullptr && cell->cell_type != CellType::Input) {
                    output_p.insert(make_pair(make_pair(i, j), cell->polarization));
                }
            }
        }

        for (auto it : output_p) {
            cout << it.first.first << " " << it.first.second << " " << it.second << endl;
        }

        //compute output_diff
        output_diff = compute_polarization_energy(output_p);

        //initialize best state to current state
        best_p    = output_p;
        best_diff = output_diff;
    }

    void SimEngine::make_anealing_iteration() {
        neighbour();
        accept();
    }

    void SimEngine::neighbour() {
        //initialize neighbour to the last iteration output polarization
        neighbour_p = output_p;

        //randomly select an output cell and change its polarization
        auto it = neighbour_p.begin();
        advance(it, rand()%output_p.size());

        it->second += (rand() * 1.0 / RAND_MAX - 0.5)*2  * sa_temp / 1000;//[-1, 1]
        it->second = fmod(it->second, 1);

        //compute diff value for neighbour
        neighbour_diff = compute_polarization_energy(neighbour_p);
    }

    void SimEngine::accept() {
        //judge the acceptance of the neighbour
        bool accept = compare_energy(output_diff, neighbour_diff);

        if (accept) {
            output_p = neighbour_p;
            output_diff = neighbour_diff;

            for (auto it=output_p.begin(); it!=output_p.end(); ++it) {
                auto &xidx = it->first.first;
                auto &yidx = it->first.second;

                assert(circuit->get_cell(xidx, yidx) != nullptr);
                circuit->get_cell(xidx, yidx)->polarization = it->second;
            }

            //check the quality of the accepted neighbour
            if (neighbour_diff < best_diff) {
                best_p    = neighbour_p;
                best_diff = neighbour_diff;
            }
        }

        //cooling the sa algorithm temperature
        sa_temp *= cooling_rate;
    }


    long double SimEngine::compute_polarization_energy(const Polarization &pola) const {

        long double diff = 0;

        for (Polarization::const_iterator it=pola.begin(); it!=pola.end(); ++it) {
            auto &xidx = it->first.first;
            auto &yidx = it->first.second;
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
            long double new_pola_val = sigma/tmp * tanh(tmp/(2*BOLTZMANN*QCA_TEMPERATURE));

            //update pola_from_neighbour and target_for_neighbour
            diff += pow(new_pola_val-cur_pola_val, 2);
        }

        return sqrt(diff);
    }

    bool SimEngine::compare_energy(long double circuit_diff, long double neighbour_diff) const {
        if (neighbour_diff < circuit_diff) {
            return true;
        } else {
            long double prob = exp((-1*output_diff) / sa_temp);
#ifndef NDBUG
            const_cast<ofstream&>(fs) << string("enter judge mode") << endl;
            const_cast<ofstream&>(fs) << string("accept probability is ") << prob << endl;
#endif
            if (prob > rand()*1.0/RAND_MAX) {//[0,1]
#ifndef NDBUG
                const_cast<ofstream&>(fs) << string("accept neighbour") << endl;
#endif
                return true;
            }

#ifndef NDBUG
            const_cast<ofstream&>(fs) << string("not accept neighbour") << endl;
#endif
        }
        return false;
    }

}
