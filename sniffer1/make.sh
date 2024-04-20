#!/bin/sh

# Copyright (C) 2024, John Clark <inindev@gmail.com>

set -e

PROJ="$(dirname "$(realpath "$0")")"

# install
#   arduino-cli board listall
#   arduino-cli core list --all --additional-urls https://raw.githubusercontent.com/HelTecAutomation/CubeCell-Arduino/master/package/package_CubeCell_index.json
#   arduino-cli core install CubeCell:CubeCell:CubeCell-Board-V2 --additional-urls https://raw.githubusercontent.com/HelTecAutomation/CubeCell-Arduino/master/package/package_CubeCell_index.json

echo "building $PROJ"
arduino-cli compile -b CubeCell:CubeCell:CubeCell-Board-V2 "$PROJ"
arduino-cli upload -b CubeCell:CubeCell:CubeCell-Board-V2 -p /dev/cu.usbserial-0001 "$PROJ"
