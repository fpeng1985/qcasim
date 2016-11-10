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

    void SimEngine::set_circuit(std::shared_ptr<QCACircuit> circuit) {
        this->circuit = circuit;
    }

    void SimEngine::run_simulation(const PolarizationList &input_p) {
        assert(circuit != nullptr);

        set_input_polarization(input_p);
        set_non_input_polarization_zero();

        static const long double convergence_factor = 1e-10;

        shared_ptr<QCACell> cell = nullptr;
        PolarizationList old_pola;
        PolarizationList new_pola;
        size_t cell_cnt = 0;
        long double convergence_val = 100;//any positive value greater than convergence_factor is OK!

        for (auto rit = circuit->row_begin(); rit != circuit->row_end(); ++rit) {
            for (auto cit = circuit->col_begin(rit); cit != circuit->col_end(rit); ++cit) {
                cell = *cit;
                if (cell != nullptr) {
                    if (cell->cell_type == CellType::Input) continue;

                    new_pola.push_back(make_tuple(cell->r_index, cell->c_index, cell->polarization));
                    ++cell_cnt;
                }
            }
        }

#ifndef NDBUG
        int loop_cnt = 0;
        static int ii = 0;
        string filename = "/iterative_output_" + to_string(ii++) + ".txt";
        fs.open(getenv("HOME") + filename);
#endif

        do {
            //[1]save current circuit state into old_pola
            old_pola = new_pola;

            //[2]compute new polarization according to the formula, and save the results into new_pola
            new_pola.clear();
            for (auto rit = circuit->row_begin(); rit != circuit->row_end(); ++rit) {
                for (auto cit = circuit->col_begin(rit); cit != circuit->col_end(rit); ++cit) {
                    cell = *cit;
                    if (cell != nullptr) {
                        if (cell->cell_type == CellType::Input) continue;

                        auto &ridx = cell->r_index;
                        auto &cidx = cell->c_index;

                        auto pola_val = compute_polarization_from_neighbour_cells(ridx, cidx);

                        new_pola.push_back(make_tuple(ridx, cidx, pola_val));
                        cell->polarization = pola_val;//write back the newly computed polarization to the circuit
                    }
                }
            }
            assert(old_pola.size() == new_pola.size());

            //[3]update convergence_val
            long double diff_1 = 0;

            for (size_t i = 0; i < old_pola.size(); ++i) {
                auto &old_pola_val = get<2>(old_pola[i]);
                auto &new_pola_val = get<2>(new_pola[i]);

                diff_1 += pow(old_pola_val - new_pola_val, 2);
            }

            diff_1 = sqrt(diff_1);

            convergence_val = diff_1;

#ifndef NDBUG
            fs << loop_cnt++ << "th iteration" << endl;
            fs << "convergence_val : " << convergence_val << endl;
            fs << *circuit << endl;
#endif

        } while (convergence_val > convergence_factor);//end while loop
    }

    void SimEngine::set_input_polarization(const PolarizationList &input_p) const {
        int ridx, cidx;
        long double pola_val;

        for (auto &p : input_p) {
            tie(ridx, cidx, pola_val) = p;

            assert(circuit->get_cell(ridx, cidx) != nullptr);
            assert(circuit->get_cell(ridx, cidx)->cell_type == CellType::Input);
            circuit->get_cell(ridx, cidx)->polarization = pola_val;
        }
    }

    void SimEngine::set_non_input_polarization_zero() const {
        shared_ptr<QCACell> cell;
        for (auto rit=circuit->row_begin(); rit!=circuit->row_end(); ++rit) {
            for (auto cit=circuit->col_begin(rit); cit!=circuit->col_end(rit); ++cit) {
                cell = *cit;
                if (cell != nullptr && cell->cell_type != CellType::Input) {
                    cell->polarization = 0;
                }
            }
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

    long double SimEngine::compute_polarization_from_neighbour_cells(int ridx, int cidx) const {

        static const long double QCA_TEMPERATURE = 1.0;
        static const long double RI = 3.8e-23;
        static const long double BOLTZMANN = 1.3806505e-23;

        static const long double EK1 = 2.3637e-22;
        static const long double EK2 = 7.068e-24;
        static const long double EK3 = -5.1638e-23;

        assert(circuit->get_cell(ridx, cidx) != nullptr);

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

        return new_pola_val;
    }

/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////

    //Simulated Anealing Algorithm
    const long double SimulatedAnealingSimEngine::MAX_TEMP = 1000;
    const long double SimulatedAnealingSimEngine::MIN_TEMP = 0.01;
    const long double SimulatedAnealingSimEngine::cooling_rate = 0.999;

    SimulatedAnealingSimEngine::SimulatedAnealingSimEngine() : SimEngine() {
        srand(time(0));
    }

    void SimulatedAnealingSimEngine::run_simulation(const PolarizationList &input_p) {
        assert(circuit != nullptr);

        set_input_polarization(input_p);
        set_non_input_polarization_randomly();
        setup_runtime_states();

        int i=0;
        while ( sa_temp > MIN_TEMP ) {
            neighbour();
            accept();
            cooling();

#ifndef NDBUG
            fs << i++  << "th iteration, energy is " << output_energy <<  endl;
            fs << *circuit << endl;
#endif
        }
    }

    void SimulatedAnealingSimEngine::setup_runtime_states() {
#ifndef NDBUG
        static int ii = 0;
        string filename = "/output" + to_string(ii++) + ".txt";
        fs.open(getenv("HOME") + filename);
#endif

        shared_ptr<QCACell> cell;

        //using the polarization formula to set the initial state of the circuit
        for (size_t i=0; i<5; ++i) {
            for (auto rit=circuit->row_begin(); rit!=circuit->row_end(); ++rit) {
                for (auto cit=circuit->col_begin(rit); cit!=circuit->col_end(rit); ++cit) {
                    cell = *cit;
                    if (cell != nullptr && cell->cell_type != CellType::Input) {
                        cell->polarization = compute_polarization_from_neighbour_cells(cell->r_index, cell->c_index);
                    }
                }
            }
        }

        //setup input_p and output_p states
        int r=0;
        for (auto rit=circuit->row_begin(); rit!=circuit->row_end(); ++rit, ++r) {
            int c=0;
            for (auto cit=circuit->col_begin(rit); cit!=circuit->col_end(rit); ++cit, ++c) {
                cell = *cit;
                if (cell != nullptr) {
                    if(cell->cell_type == CellType::Input) {
                        input_p.push_back(make_tuple(r,c,cell->polarization));
                    } else {
                        output_p.push_back(make_tuple(r,c,cell->polarization));
                    }
                }
            }
        }

        //compute output_energy
        output_energy = compute_polarization_energy(output_p);

        //initialize best state to current state
        best_p    = output_p;
        best_energy = output_energy;

        //set simulated annealing algorithm parameters
        sa_temp = MAX_TEMP;
    }

    void SimulatedAnealingSimEngine::neighbour() {
        //initialize neighbour to the last iteration output polarization
        neighbour_p = output_p;

        //randomly select an output cell and change its polarization
        auto it = neighbour_p.begin();
        advance(it, rand()%output_p.size());

        long double rand_val = (rand() * 1.0 /RAND_MAX - 0.5) * 2;

        int ridx, cidx;
        long double pola_val;
        tie(ridx, cidx, pola_val) = *it;

        long double avg_val = compute_polarization_from_neighbour_cells(ridx, cidx);

        get<2>(*it) = rand_val * sa_temp / MAX_TEMP +  avg_val * (1000 - sa_temp) / MAX_TEMP;

        //compute energy value for neighbour
        neighbour_energy = compute_polarization_energy(neighbour_p);
    }

    void SimulatedAnealingSimEngine::accept() {
        //judge the acceptance of the neighbour
        bool accept = compare_energy(output_energy, neighbour_energy);

        if (accept) {
            output_p = neighbour_p;
            output_energy = neighbour_energy;

            for (auto it=output_p.begin(); it!=output_p.end(); ++it) {
                int ridx, cidx;
                long double pola_val;
                tie(ridx, cidx, pola_val) = *it;

                assert(circuit->get_cell(ridx, cidx) != nullptr);
                circuit->get_cell(ridx, cidx)->polarization = pola_val;
            }

            //check the quality of the accepted neighbour
            if (neighbour_energy < best_energy) {
                best_p    = neighbour_p;
                best_energy = neighbour_energy;
            }
        }
    }

    void SimulatedAnealingSimEngine::cooling() {
        //cooling the sa algorithm temperature
        sa_temp *= cooling_rate;
    }

    long double SimulatedAnealingSimEngine::compute_polarization_energy(const PolarizationList &output_p) const {
        static constexpr long double pi = 3.14159265359;
        static constexpr long double eps0 = 8.854e-12;
        static constexpr long double epsr = 12.9;
        static constexpr long double e = 1.602e-19;
        static constexpr long double base_factor = 1/(4*pi*epsr*eps0);

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
        PolarizationList pola = output_p;
        for (auto &p : input_p) {
            pola.push_back(p);
        }

        int ridx, cidx;
        long double cur_pola_val;
        for (PolarizationList::const_iterator it=pola.begin(); it!=pola.end(); ++it) {
            tie(ridx, cidx, cur_pola_val) = *it;

            shared_ptr<QCACell> curcell;

            //individual energy
            curcell = circuit->get_cell(ridx, cidx);
            assert(curcell != nullptr);

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

        long double total_energy = (individual_energy + mutal_energy/2.0);
        return total_energy;
    }

    bool SimulatedAnealingSimEngine::compare_energy(long double circuit_energy, long double neighbour_energy) const {
        static const long double FACTOR = 10e50;

        if (neighbour_energy < circuit_energy) {
#ifndef NDBUG
            const_cast<ofstream&>(fs) << string("directly accept neighbour") << endl;
#endif
            return true;
        } else {
            long double prob = exp( (circuit_energy*FACTOR - output_energy*FACTOR) / sa_temp );//please note output_energy is a negative value!!!
#ifndef NDBUG
            const_cast<ofstream&>(fs) << string("enter judge mode") << endl;
            const_cast<ofstream&>(fs) << string("energy diff : ") << (circuit_energy*FACTOR - output_energy*FACTOR) << endl;
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
