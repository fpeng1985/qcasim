include_directories(SYSTEM ${PROJECT_SOURCE_DIR})

add_executable(test_qca_circuit test_qca_circuit.cpp ../qca_circuit.cpp ../qca_cell.cpp)