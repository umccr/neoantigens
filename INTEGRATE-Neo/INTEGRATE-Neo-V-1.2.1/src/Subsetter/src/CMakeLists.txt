cmake_minimum_required(VERSION 2.8)

set(LIB_SOURCES
    Util.cpp
    Util.h
    Bedpe.cpp
    Bedpe.h
    CutterByRule.cpp
    CutterByRule.h
    MyTypes.h
)

add_library(fusionBedpeSubsetter ${LIB_SOURCES})
add_executable(fusionBedpeSubsetter-bin main.cpp)
target_link_libraries(fusionBedpeSubsetter-bin fusionBedpeSubsetter)
set_target_properties(fusionBedpeSubsetter-bin PROPERTIES OUTPUT_NAME fusionBedpeSubsetter)
