#!/bin/sh
current_directory="$(cd "$(dirname "$0")" && pwd)"
project_root_dir=$current_directory

echo "Hello Friend! Let's try and build this application together :D"

mkdir -p "${project_root_dir}/third_party"
cd "${project_root_dir}/third_party"
#compile MATPLOTPLUSPLUS
if [ ! -d "matplotplusplus" ]; then
    git clone https://github.com/alandefreitas/matplotplusplus
    cd "${project_root_dir}/third_party/matplotplusplus"
    rm -rf build install
    mkdir -p build && mkdir -p install && cd build
    cmake -DCMAKE_INSTALL_PREFIX="${project_root_dir}/third_party/matplotplusplus/install" -DCMAKE_BUILD_TYPE=Release -DMATPLOTPP_BUILD_TESTS:BOOL=OFF -DMATPLOTPP_BUILD_EXAMPLES:BOOL=OFF ..
    cmake --build . -j 4
    cmake --build . --config debug
    cmake --build . --config release
    cmake --install .
fi

# compile the rest of the application
cd "${project_root_dir}"

# remove the build directory that has the current code in it
echo "deleting the BUILD directory"
rm -rf "${project_root_dir}/build"

echo "make a new BUILD directory to start the compiling process"
mkdir -p "${project_root_dir}/build"
cd "${project_root_dir}/build"

echo "cmake engage!"
cmake ..

echo "convert this to an executable application -- let's go!!"
cmake --build . -j 4
cmake --build . --target install --config debug
cd "${project_root_dir}"
echo "declare success -- hooray!"

echo "running the executable with some default parameters"
echo "./build/debug/matplot++ 2>&1"
echo "  the 2>&1 redirects the stderr to a 1 so we don't see the gnuplot problems"
./build/main 2>&1
