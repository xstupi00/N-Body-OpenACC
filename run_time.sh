#!/bin/bash

# NBODY parameters for running the binary
DT=0.01f
STEPS=500

INPUT_DIR="inputs"
OUTPUT_DIR="outputs"
GEN_BIN="./commons/gen"
NBODY_BIN="./step3/nbody"
RESULT_FILE="step3/times.txt"

# Check whether exists OUTPUT FILE to store the results
if [ ! -d "$RESULT_FILE" ]; then
    touch ${RESULT_FILE}
fi
# Delete the OUT FILE to store the results
> ${RESULT_FILE}

# Check whether exists INPUT directory
if [ ! -d "$INPUT_DIR" ]; then
  mkdir "$INPUT_DIR"
fi
# Check whether exists OUTPUT directory
if [ ! -d "$OUTPUT_DIR" ]; then
  mkdir "$OUTPUT_DIR"
fi
# Check whether INPUT directory is empty
if [ "$(ls -A $INPUT_DIR)" ]; then
  rm -- "${INPUT_DIR}"/*
fi
# Check whether OUTPUT directory is empty
if [ "$(ls -A $OUPUT_DIR)""step0/times.txt" ]; then
  rm -- "${OUTPUT_DIR}"/*
fi

# Generate INPUT files with available generator
for i in {1..30}; do
  ${GEN_BIN} $((i*1024)) "${INPUT_DIR}/$((i*1024)).h5"
done

# Run the individual INPUT files and write the measured time to the file
for i in {1..30}; do
  IN_FILE="${INPUT_DIR}/$((i*1024)).h5"
  OUT_FILE="${OUTPUT_DIR}/$((i*1024)).h5"
  OUT=$(${NBODY_BIN} $((i*1024)) ${DT} ${STEPS} 128 ${IN_FILE} ${OUT_FILE})
  TIME="$(echo "$OUT" |  sed -n '4 p' | cut -d' ' -f2)"
  RESULT="$(printf "$((i)) * 1024 \t\t\t %g\n" ${TIME})"
  echo "$RESULT" >> "${RESULT_FILE}"
done