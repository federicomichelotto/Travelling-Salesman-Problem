# $ cmake --version
cmake_minimum_required(VERSION 2.8...3.21)

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
message(STATUS "You're using ${CMAKE_VERSION}")

if (WIN32)

    # put here some stuff related to WIN32

elseif (UNIX AND NOT APPLE)

    # To link CPLEX libraries and directories stored inside Blade, you must uncomment the following line
    # the latest version installed in the cluster is 12.10 (but obviously you can choose another one)
#    set(CPLEX_HOME /nfsd/opt/CPLEX12.10/cplex)
#    set(CONCORDELIB /home/miliamikel/concorde)

    # To link CPLEX libraries and directories stored inside your personal Linux distro just set the right path
    set(CPLEX_HOME /opt/ibm/ILOG/CPLEX_Studio201/cplex)
    set(CONCORDELIB /opt/concorde)

    # To link and include useful CPLEX libraries and directories, play with the following lines
    target_link_libraries(${PROJECT_NAME} -lcplex -lm -lpthread -ldl)
    target_link_libraries(${PROJECT_NAME} -L${CPLEX_HOME}/lib/x86-64_linux/static_pic)
    target_link_libraries(${PROJECT_NAME} ${CONCORDELIB}/concorde.a)
    target_include_directories(${PROJECT_NAME} PRIVATE ${CPLEX_HOME}/include/ilcplex ${CONCORDELIB} )
    

else (APPLE)

    # To link CPLEX libraries and directories stored inside your personal MacOS distro just set the right path
    set(CPLEX_HOME /Applications/CPLEX_Studio201/cplex)

    # To link and include useful CPLEX libraries and directories, play with the following lines
    target_link_libraries(${PROJECT_NAME} -lcplex -lm -lpthread -ldl)
    target_link_libraries(${PROJECT_NAME} -L${CPLEX_HOME}/lib/x86-64_osx/static_pic)
    target_include_directories(${PROJECT_NAME} PRIVATE ${CPLEX_HOME}/include/ilcplex)

endif ()
