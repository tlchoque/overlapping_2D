########### next target ###############

add_executable(MHM main.cpp TPZMatLaplacianLagrange.h TPZMatLaplacianLagrange.cpp meshgen.h meshgen.cpp)

target_link_libraries(MHM pz)

add_executable(MHMElast mainElast.cpp TPZMatLaplacianLagrange.h TPZMatLaplacianLagrange.cpp meshgen.h meshgen.cpp)

target_link_libraries(MHMElast pz)

configure_file(porous.txt porous.txt COPYONLY)
