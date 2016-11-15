//
// Created by Administrator on 2016/10/15.
//

#define CATCH_CONFIG_MAIN
#include "catch.hpp"

using namespace std;

#include "sim_manager.h"
using namespace hfut;

SCENARIO("combination generation", "[combination generation]") {

    vector<int> a;
    for (int i=0; i<10; ++i) {
        a.push_back(i);
    }

    vector<vector<int>> combinations;

    CombinationGenerator<int> comgen;

    GIVEN("A size number m = 1") {
        int m = 1;

        WHEN("we generate the combinations") {
            comgen.generate_combination(a, m, combinations);

            THEN("we get the 1 item combinations") {
                for (auto &comb : combinations) {
                    for (auto id : comb) {
                        cout << id << " ";
                    }
                    cout << endl;
                }
            }
        }
    }

    GIVEN("A size number m = 3") {
        int m = 3;

        WHEN("we generate the combinations") {
            comgen.generate_combination(a, m, combinations);

            THEN("we get the 3 item combinations") {
                for (auto &comb : combinations) {
                    for (auto id : comb) {
                        cout << id << " ";
                    }
                    cout << endl;
                }
            }
        }
    }
}
