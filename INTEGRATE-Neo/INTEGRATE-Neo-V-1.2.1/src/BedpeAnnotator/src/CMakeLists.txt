cmake_minimum_required(VERSION 2.8)

set(LIB_SOURCES
    Util.cpp
    Util.h
    Bedpe.cpp
    Bedpe.h
    Gene2.cpp
    Gene2.h
    Reference.cpp
    Reference.h
    Annotate.cpp
    Annotate.h
    Dinucleo.cpp
    Dinucleo.h
    AnnotateDi.cpp
    AnnotateDi.h
    CMakeLists.txt
    MyTypes.h
)

add_library(bedpeAnnotator ${LIB_SOURCES})
add_executable(bedpeAnnotator-bin main.cpp)
target_link_libraries(bedpeAnnotator-bin bedpeAnnotator)
set_target_properties(bedpeAnnotator-bin PROPERTIES OUTPUT_NAME fusionBedpeAnnotator)
