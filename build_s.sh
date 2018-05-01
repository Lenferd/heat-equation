#!/bin/bash

while [[ $# -ge 1 ]]
do
    key=$1
    case $key in
        -h|--help)
            printf "Help
$./build.sh -m debug
or  $./build.sh , and this should be release conf
"
            exit
            ;;
        -m|--mode)
            BUILD_MODE=$2
            shift
            ;;
        -c|--clear)
            Clear
            ;;
        *)
            ;;
    esac

    shift
done

Build () {
    root_dir="$( pwd )"
    mkdir build
    cd build
    echo $root_dir
    if [ "$BUILD_MODE" == "debug" ]; then
        BUILD_MODE="Debug"
    else
        BUILD_MODE="Release"
    fi

    cmake -DCMAKE_BUILD_TYPE=$BUILD_MODE $root_dir
    make -j4
}

Clear() {
  if [ -d build ]
    then
    rm -rf ./build
  fi
}

Build
