include_directories(SYSTEM ${PROJECT_SOURCE_DIR})

set(QCA_SRC ../qca_cell.cpp ../qca_circuit.cpp)

add_executable(test_qca_circuit test_qca_circuit.cpp ${QCA_SRC})
add_executable(test_simulated_anealing_sim_engine test_simulated_anealing_sim_engine.cpp  ../sim_engine.cpp ${QCA_SRC})
add_executable(test_iterative_sim_engine          test_iterative_sim_engine.cpp           ../sim_engine.cpp ${QCA_SRC})
add_executable(test_sim_manager test_sim_manager.cpp ../sim_manager.cpp ../sim_engine.cpp ${QCA_SRC})