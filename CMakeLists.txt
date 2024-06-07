cmake_minimum_required(VERSION 3.15...3.26)

project(${SKBUILD_PROJECT_NAME} VERSION ${SKBUILD_PROJECT_VERSION} LANGUAGES CXX)
set(CMAKE_INSTALL_MESSAGE NEVER)

find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
find_package(pybind11 CONFIG REQUIRED)
find_package(gemmi 0.6 CONFIG REQUIRED)
find_package(LAPACK REQUIRED)

include(CTest)

macro(add_common_props TARGET_NAME)
  # using c++ 17 standard
	set_property(TARGET ${TARGET_NAME} PROPERTY CXX_STANDARD 17)

	# generic compile options
	target_compile_options(${TARGET_NAME} PRIVATE -fPIC -Wall -fvisibility=hidden -Werror)

	# needed because of some errors in gemmi
	if ("${IS_CLANG_COMPILER}")
		target_compile_options(${TARGET_NAME} PRIVATE -Wno-range-loop-analysis)
	endif()
endmacro()

# DEPENDENCIES
include(FetchContent)

# BLAZE
FetchContent_Declare(
	blaze_proj
  GIT_REPOSITORY https://bitbucket.org/blaze-lib/blaze.git
  GIT_TAG ee64f719aef27cb8db2f27c2dc5f7de14a24980a
)
FetchContent_MakeAvailable(blaze_proj)

# CATCH
FetchContent_Declare(
  catch_proj
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG        v2.9.2
)
FetchContent_MakeAvailable(catch_proj)

# JSON
FetchContent_Declare(
	json_proj
  GIT_REPOSITORY https://github.com/ArthurSonzogni/nlohmann_json_cmake_fetchcontent
  GIT_TAG v3.11.2
)
FetchContent_MakeAvailable(json_proj)

######################################
#  Add libraries, modules, and tests #
######################################

set(CMAKE_INSTALL_RPATH "$ORIGIN")

# LIBRARIES
file(GLOB SRC "prebio/cpp/src/*.cpp")
add_library(prebiocommon SHARED ${SRC})
target_include_directories(prebiocommon PUBLIC prebio/cpp/include ${blaze_proj_SOURCE_DIR})
target_link_libraries(prebiocommon PUBLIC gemmi::gemmi_cpp nlohmann_json::nlohmann_json ${BLAS_LIBRARIES})
add_common_props(prebiocommon)

# # TESTS FOR LIBRARIES
# set(CPP_TESTS test/cpp)
# set(CPP_TEST_NAME inpro_cpp_tests)
# add_executable(${CPP_TEST_NAME} ${CPP_TESTS}/testCif.cpp)
# target_link_libraries(${CPP_TEST_NAME} PRIVATE Catch2::Catch2WithMain cifinpro)
# target_include_directories(${CPP_TEST_NAME} PRIVATE inpro/cpp/include)
# target_compile_definitions(${CPP_TEST_NAME} PRIVATE DATADIR="${CMAKE_CURRENT_SOURCE_DIR}/test/data/")
# add_common_props(${CPP_TEST_NAME})

# add_test(NAME ${CPP_TEST_NAME} COMMAND ${CPP_TEST_NAME})

# add_custom_command(TARGET ${CPP_TEST_NAME} POST_BUILD COMMAND ctest)

# MODULES
pybind11_add_module(prebioimpl MODULE prebio/cpp/bindings/bindings.cpp)
target_link_libraries(prebioimpl PUBLIC pybind11::headers prebiocommon)
target_compile_definitions(prebioimpl PUBLIC MOD_NAME=prebioimpl)
add_common_props(prebioimpl)

# install prebio
install(TARGETS prebiocommon prebioimpl LIBRARY DESTINATION ${SKBUILD_PROJECT_NAME})