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
    mkdir _build
    cd _build
    echo $root_dir
    if [ "$BUILD_MODE" == "debug" ]; then
        BUILD_MODE="Debug"
    else
        BUILD_MODE="Release"
    fi

    cmake -DCMAKE_BUILD_TYPE=$BUILD_MODE $root_dir
    make -j4

    cd ..
    SetLink
}

Clear() {
  if [ -d _build ]
    then
    rm -rf ./_build
  fi
}

SetLink() {
        ln -sf "$(pwd)/initial" "$(pwd)/_build/Sergey_$BUILD_MODE/"
        ln -sf $(pwd)/initial/Sergey/setting2.ini $(pwd)/_build/Sergey_$BUILD_MODE/setting2.ini
        ln -sf "$(pwd)/result/Sergey" "$(pwd)/_build/Sergey_$BUILD_MODE/"
}

Build
# ln -s script/Sergey
