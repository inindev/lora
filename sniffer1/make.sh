#!/bin/sh

# Copyright (C) 2024, John Clark <inindev@gmail.com>

set -e

PROJ="$(dirname "$(realpath "$0")")"

# install
#   arduino-cli board listall
#   arduino-cli core list --all --additional-urls https://raw.githubusercontent.com/HelTecAutomation/CubeCell-Arduino/master/package/package_CubeCell_index.json
#   arduino-cli core install CubeCell:CubeCell:CubeCell-Board-V2 --additional-urls https://raw.githubusercontent.com/HelTecAutomation/CubeCell-Arduino/master/package/package_CubeCell_index.json

echo "building $PROJ"
$HOME/bin/arduino-cli compile -b CubeCell:CubeCell:CubeCell-Board-V2 "$PROJ"

echo
PORT="$(ls /dev/cu.usbserial* 2>/dev/null || ls /dev/cu.usbmodem* 2>/dev/null || ls /dev/ttyUSB* 2>/dev/null)"
if lsof -t -S 2 -O "$PORT" >/dev/null; then
    echo "\033[31merror: cannot program, serial port $PORT is in use\033[m\n"
    exit 1
fi
$HOME/bin/arduino-cli upload -b CubeCell:CubeCell:CubeCell-Board-V2 -p "$PORT" "$PROJ"

