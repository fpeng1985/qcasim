//
// Created by Administrator on 2016/10/15.
//

#define CATCH_CONFIG_MAIN
#include "catch.hpp"

using namespace std;

#include "sim_manager.h"
using namespace hfut;

SCENARIO("majority gate 1", "[majority_gate_1]") {
    SimManager manager;

    SimManager::CombinationGenerator comgen;

    GIVEN("A circuit majority gate 1") {
        string mg1 = getenv("COMMON") + string("/benchmark/qcasim/majority_gate_1.txt");

        WHEN("we load the benchmark") {

            vector<vector<int>> combinations;
            comgen.generate_combination(10, combinations);

            for (auto &comb : combinations) {
                for (auto id : comb) {
                    cout << id << " ";
                }
                cout << endl;
            }

            manager.load_benchmark(mg1);
            THEN("we get the corresponding circuit structures") {
            }
        }
    }
}
