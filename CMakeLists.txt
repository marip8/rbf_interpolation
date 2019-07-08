cmake_minimum_required(VERSION 2.8.3)
project(rbf_interpolation)

add_compile_options(-std=c++11)

find_package(catkin REQUIRED COMPONENTS
  roscpp
)

find_package(Eigen3 REQUIRED)

catkin_package(
  INCLUDE_DIRS
    include
  CATKIN_DEPENDS
    roscpp
)

###########
## Build ##
###########

include_directories(
  include
  ${catkin_INCLUDE_DIRS}
  ${EIGEN3_INCLUDE_DIRS}
)

#############
## Testing ##
#############
if(CATKIN_ENABLE_TESTING)
  find_package(rostest REQUIRED)
  catkin_add_gtest(${PROJECT_NAME}_test
    test/utest.cpp
  )
  target_link_libraries(${PROJECT_NAME}_test
    ${catkin_LIBRARIES}
  )
endif()

#############
## Install ##
#############

install(DIRECTORY include/${PROJECT_NAME}/
  DESTINATION ${CATKIN_PACKAGE_INCLUDE_DESTINATION}
)