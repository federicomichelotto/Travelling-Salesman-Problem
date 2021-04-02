# $ cmake --version
cmake_minimum_required(VERSION 2.8)

# Setting the project name
project(tsp)

# Create a sources variable with a link to all cpp files to compile
set(SOURCE_FILES
        src/main.c
        src/tsp.c
        src/utils.c)

# Add an executable with the above sources
add_executable(${PROJECT_NAME} ${SOURCE_FILES})

message(STATUS "You're using ${CMAKE_SYSTEM_NAME}")

if (WIN32)

    # put here some stuff related to WIN32

elseif (UNIX AND NOT APPLE)

    set(CPLEX_HOME /opt/ibm/ILOG/CPLEX_Studio201/cplex)     # mandatory, if you're using Linux set the right path to CPLEX
    target_link_libraries(${PROJECT_NAME} -lcplex -lm -lpthread -ldl)
    target_link_libraries(${PROJECT_NAME} -L${CPLEX_HOME}/lib/x86-64_linux/static_pic)
    target_include_directories(${PROJECT_NAME} PRIVATE ${CPLEX_HOME}/include/ilcplex)

else (APPLE)

    set(CPLEX_HOME /Applications/CPLEX_Studio201/cplex)     # mandatory, if you're using MacOS set the right path to CPLEX
    target_link_libraries(${PROJECT_NAME} -lcplex -lm -lpthread -ldl)
    target_link_libraries(${PROJECT_NAME} -L${CPLEX_HOME}/lib/x86-64_osx/static_pic)
    target_include_directories(${PROJECT_NAME} PRIVATE ${CPLEX_HOME}/include/ilcplex)

endif ()