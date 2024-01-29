##################################################
### for test
##################################################
enable_testing (test)
### for tutorial
add_test (NAME step2-2d COMMAND asfem "-i" "${CMAKE_CURRENT_SOURCE_DIR}/examples/tutorial/step2-2d.json")
add_test (NAME step2-3d COMMAND asfem "-i" "${CMAKE_CURRENT_SOURCE_DIR}/examples/tutorial/step2-2d.json")
add_test (NAME step3-2d COMMAND asfem "-i" "${CMAKE_CURRENT_SOURCE_DIR}/examples/tutorial/step3-2d.json")
add_test (NAME step3-3d COMMAND asfem "-i" "${CMAKE_CURRENT_SOURCE_DIR}/examples/tutorial/step3-3d.json")
### for test_input
add_test (NAME mesh1d COMMAND asfem "-i" "${CMAKE_CURRENT_SOURCE_DIR}/test_input/mesh/mesh1d.json --read-only")
add_test (NAME mesh2d COMMAND asfem "-i" "${CMAKE_CURRENT_SOURCE_DIR}/test_input/mesh/mesh2d.json --read-only")
add_test (NAME mesh3d COMMAND asfem "-i" "${CMAKE_CURRENT_SOURCE_DIR}/test_input/mesh/mesh3d.json --read-only")
add_test (NAME bcs COMMAND asfem "-i" "${CMAKE_CURRENT_SOURCE_DIR}/test_input/bcs/poisson-2d-mixed.json")
add_test (NAME importmesh2 COMMAND asfem "-i" "${CMAKE_CURRENT_SOURCE_DIR}/test_input/importmesh/poisson-msh2-2d.json")
add_test (NAME importmesh4 COMMAND asfem "-i" "${CMAKE_CURRENT_SOURCE_DIR}/test_input/importmesh/poisson-msh4-2d.json")
add_test (NAME postprocess COMMAND asfem "-i" "${CMAKE_CURRENT_SOURCE_DIR}/test_input/postprocess/poisson-2d-pps.json")
add_test (NAME postprocess-sinxy COMMAND asfem "-i" "${CMAKE_CURRENT_SOURCE_DIR}/test_input/postprocess/sinxy-integration.json")