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
//        convergence_factor = 1e-8;

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
        set_non_input_polarization_randomly();
        setup_runtime_states();

        int i=0;
        while ( sa_temp > terminate_temp ) {
            neighbour();
            accept();

#ifndef NDBUG
            fs << i  << "th iteration, diff is " << output_energy <<  endl;
            fs << *circuit << endl;
#endif
        }
    }

    void SimEngine::set_input_polarization(const Polarization &pola) const {
        for (auto &it : pola) {
            auto &ridx = it.first.first;
            auto &cidx = it.first.second;

            assert(circuit->get_cell(ridx, cidx) != nullptr);
            assert(circuit->get_cell(ridx, cidx)->cell_type == CellType::Input);
            circuit->get_cell(ridx, cidx)->polarization = it.second;
        }
    }

    void SimEngine::set_non_input_polarization_randomly() const {
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
        shared_ptr<QCACell> cell;
        int r=0;
        for (auto rit=circuit->row_begin(); rit!=circuit->row_end(); ++rit, ++r) {
            int c=0;
            for (auto cit=circuit->col_begin(rit); cit!=circuit->col_end(rit); ++cit, ++c) {
                cell = *cit;
                if (cell != nullptr) {
                    if(cell->cell_type == CellType::Input) {
                        input_p.insert(make_pair(make_pair(r, c), cell->polarization));
                    } else {
                        output_p.insert(make_pair(make_pair(r, c), cell->polarization));
                    }
                }
            }
        }

        //compute output_energy
        output_energy = compute_polarization_energy(output_p);

        //initialize best state to current state
        best_p    = output_p;
        best_diff = output_energy;
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
        neighbour_energy = compute_polarization_energy(neighbour_p);
    }

    void SimEngine::accept() {
        //judge the acceptance of the neighbour
        bool accept = compare_energy(output_energy, neighbour_energy);

        if (accept) {
            output_p = neighbour_p;
            output_energy = neighbour_energy;

            for (auto it=output_p.begin(); it!=output_p.end(); ++it) {
                auto &ridx = it->first.first;
                auto &cidx = it->first.second;

                assert(circuit->get_cell(ridx, cidx) != nullptr);
                circuit->get_cell(ridx, cidx)->polarization = it->second;
            }

            //check the quality of the accepted neighbour
            if (neighbour_energy < best_diff) {
                best_p    = neighbour_p;
                best_diff = neighbour_energy;
            }
        }

        //cooling the sa algorithm temperature
        sa_temp *= cooling_rate;
    }

    long double SimEngine::compute_polarization_energy(const Polarization &output_p) const {
        static constexpr long double pi = 3.14159265359;
        static constexpr long double eps0 = 8.854e-12;
        static constexpr long double epsr = 12.9;
        static constexpr long double e = 1.602e-19;
        static constexpr long double base_factor = 1/(4*pi*epsr*eps0) *10e30;

        static const long double individual_factor = base_factor *(sqrt(2)-4)/9.0 * pow((e/2),2);
        static const long double mutal_factor_1 = base_factor *
                (1.0/5.0 + 2.0/sqrt(202) + 2.0/sqrt(922) -4.0/sqrt(481) -2.0/11 - 2.0/29) * pow((e/2), 2);
        static const long double mutal_factor_2 = base_factor *
                (1.0/10 + 2.0/sqrt(979) + 2.0/sqrt(2482) - 4/sqrt(1618) - 2.0/31 - 2.0/49) * pow((e/2), 2);
        static const long double mutal_factor_3 = base_factor *
                (1.0/sqrt(50) + 2.0/sqrt(962) + 1.0/sqrt(242) + 1.0/sqrt(1682) - 4.0/sqrt(1241) - 4.0/sqrt(521)) * pow((e/2), 2);

        long double individual_energy = 0;
        long double mutal_energy = 0;

        //fill the input polarizations
        Polarization pola = output_p;
        for (auto &p : input_p) {
            pola.insert(p);
        }

        for (Polarization::const_iterator it=pola.begin(); it!=pola.end(); ++it) {
            auto &ridx = it->first.first;
            auto &cidx = it->first.second;
            const long double &cur_pola_val = it->second;

            shared_ptr<QCACell> curcell;

            //individual energy
            curcell = circuit->get_cell(ridx, cidx);
            assert(curcell != nullptr);
//            cout << ridx << " " << cidx << " " << cur_pola_val << " " << circuit->get_cell(ridx, cidx)->polarization << endl;
//            assert(cur_pola_val == circuit->get_cell(ridx, cidx)->polarization);

            individual_energy += individual_factor * pow(cur_pola_val, 2);

            //type 1 mutal energy
            curcell = circuit->get_cell(ridx+1, cidx);
            if (curcell != nullptr) {
                mutal_energy += mutal_factor_1 * cur_pola_val * curcell->polarization;
            }

            curcell = circuit->get_cell(ridx-1, cidx);
            if (curcell != nullptr) {
                mutal_energy += mutal_factor_1 * cur_pola_val * curcell->polarization;
            }

            curcell = circuit->get_cell(ridx, cidx+1);
            if (curcell != nullptr) {
                mutal_energy += mutal_factor_1 * cur_pola_val * curcell->polarization;
            }

            curcell = circuit->get_cell(ridx, cidx-1);
            if (curcell != nullptr) {
                mutal_energy += mutal_factor_1 * cur_pola_val * curcell->polarization;
            }

            //type 2 mutal energy
            curcell = circuit->get_cell(ridx+2, cidx);
            if (curcell != nullptr) {
                mutal_energy += mutal_factor_2 * cur_pola_val * curcell->polarization;
            }

            curcell = circuit->get_cell(ridx-2, cidx);
            if (curcell != nullptr) {
                mutal_energy += mutal_factor_2 * cur_pola_val * curcell->polarization;
            }

            curcell = circuit->get_cell(ridx, cidx+2);
            if (curcell != nullptr) {
                mutal_energy += mutal_factor_2 * cur_pola_val * curcell->polarization;
            }

            curcell = circuit->get_cell(ridx, cidx-2);
            if (curcell != nullptr) {
                mutal_energy += mutal_factor_2 * cur_pola_val * curcell->polarization;
            }

            //type 3 mutal energy
            curcell = circuit->get_cell(ridx+1, cidx+1);
            if (curcell != nullptr) {
                mutal_energy += mutal_factor_3 * cur_pola_val * curcell->polarization;
            }

            curcell = circuit->get_cell(ridx+1, cidx-1);
            if (curcell != nullptr) {
                mutal_energy += mutal_factor_3 * cur_pola_val * curcell->polarization;
            }

            curcell = circuit->get_cell(ridx-1, cidx+1);
            if (curcell != nullptr) {
                mutal_energy += mutal_factor_3 * cur_pola_val * curcell->polarization;
            }

            curcell = circuit->get_cell(ridx-1, cidx-1);
            if (curcell != nullptr) {
                mutal_energy += mutal_factor_3 * cur_pola_val * curcell->polarization;
            }
        }

        long double total_energy = individual_energy + mutal_energy/2.0;
        return total_energy;
    }

    long double SimEngine::compute_polarization_difference(const Polarization &pola) const {

        long double diff = 0;

        for (Polarization::const_iterator it=pola.begin(); it!=pola.end(); ++it) {
            auto &ridx = it->first.first;
            auto &cidx = it->first.second;
            const long double &cur_pola_val = it->second;

            shared_ptr<QCACell> curcell;
            long double sigma = 0;

            //type 1 neighbour cells
            curcell = circuit->get_cell(ridx+1, cidx);
            if (curcell != nullptr) {
                sigma += EK1 * curcell->polarization;
            }

            curcell = circuit->get_cell(ridx-1, cidx);
            if (curcell != nullptr) {
                sigma += EK1 * curcell->polarization;
            }

            curcell = circuit->get_cell(ridx, cidx+1);
            if (curcell != nullptr) {
                sigma += EK1 * curcell->polarization;
            }

            curcell = circuit->get_cell(ridx, cidx-1);
            if (curcell != nullptr) {
                sigma += EK1 * curcell->polarization;
            }

            //type 2 neighbour cells
            curcell = circuit->get_cell(ridx+2, cidx);
            if (curcell != nullptr) {
                sigma += EK2 * curcell->polarization;
            }

            curcell = circuit->get_cell(ridx-2, cidx);
            if (curcell != nullptr) {
                sigma += EK2 * curcell->polarization;
            }

            curcell = circuit->get_cell(ridx, cidx+2);
            if (curcell != nullptr) {
                sigma += EK2 * curcell->polarization;
            }

            curcell = circuit->get_cell(ridx, cidx-2);
            if (curcell != nullptr) {
                sigma += EK2 * curcell->polarization;
            }

            //type 3 neighbour cells
            curcell = circuit->get_cell(ridx+1, cidx+1);
            if (curcell != nullptr) {
                sigma += EK3 * curcell->polarization;
            }

            curcell = circuit->get_cell(ridx+1, cidx-1);
            if (curcell != nullptr) {
                sigma += EK3 * curcell->polarization;
            }

            curcell = circuit->get_cell(ridx-1, cidx+1);
            if (curcell != nullptr) {
                sigma += EK3 * curcell->polarization;
            }

            curcell = circuit->get_cell(ridx-1, cidx-1);
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

    bool SimEngine::compare_energy(long double circuit_energy, long double neighbour_energy) const {
        if (neighbour_energy < circuit_energy) {
            return true;
        } else {
            long double prob = exp((output_energy) / sa_temp);//please note output_energy is a negative value!!!
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
