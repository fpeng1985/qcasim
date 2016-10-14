//
// Created by Administrator on 2016/10/13.
//

#ifndef QCASIM_SIMENGINE_H
#define QCASIM_SIMENGINE_H

#include "qca_circuit.h"


namespace hfut {

    class SimEngine {
    public:
        void initialize_circuit();
        void create_neighbour();

    private:
        QCACircuit *circuit;

    };

}

#endif //QCASIM_SIMENGINE_H
