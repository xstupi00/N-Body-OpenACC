#!/bin/bash

KERNELS_NUM=4
INPUT_DIR="../inputs"


# Check whether exists INPUT directory
if [ ! -d "$INPUT_DIR" ]; then
  mkdir "$INPUT_DIR"
fi
# Check whether INPUT directory is empty
if [ "$(ls -A $INPUT_DIR)" ]; then
  rm -- "${INPUT_DIR}"/*
fi
# Generate inputs files
for i in {7..17};
do
  ./../commons/gen $((2**i)) ../inputs/$((2**i)).h5
done

times=()

echo "Running measuring run-times ..."

for i in {7..17};
do
  OUT=$(./nbody $((2**i)) 0.01f 500 20 ../inputs/$((2**i)).h5 step4Output.h5)
  TIME="$(echo "$OUT" |  sed -n '4 p' | cut -d' ' -f2)"
  times[${i}]=${TIME}
  echo "${TIME}"
done

echo "... done measuring run-times."

echo "Running for global memory throughput ..."

for i in {7..17};
do
  sum=0
  N=$((2**i))
  make profile_mem N=${N} > prof.out 2>&1
  stats="$(tail -n $((KERNELS_NUM*5)) prof.out | sed '1d;6d;11d;16d' | awk 'NF>1{print $NF}')"
  for stat in ${stats};
  do
    if [ "$(echo "${stat}" | tail -c5)" == "MB/s" ]; then
      sum=$(echo "${sum} + $(echo "${stat}" | head -c6)" | bc)
    elif [ "$(echo "${stat}" | tail -c5)" == "GB/s" ]; then
      gb_to_mb=$(echo "$(echo "${stat}" | head -c6) * 1000" | bc)
      sum=$(echo "${sum} +  ${gb_to_mb}" | bc)
    elif [ "$(echo "${stat}" | tail -c4)" == "B/s" ]; then
      sum=${sum}
    else
      echo "${"$(echo "${stat}" | tail -c5)"}"
      echo "Unknown units in the stats!"
      break
    fi
  done
  echo "${sum}"
done

echo "... done running for global memory throughput."

echo "Running for performance in MFLOPS ..."

for i in {7..17};
do
  sum=0
  N=$((2**i))
  make profile_flops N=${N} > prof.out 2>&1
  stats="$(tail -n $((KERNELS_NUM*3)) prof.out | sed '1d;4d;7d;10d' | awk 'NF>1{print $NF}')"
  for stat in ${stats};
  do
    sum=$(echo "${sum} + ${stat}" | bc)
  done
  echo "$(echo "scale=4; ((${sum} * 500.0) / ${times[${i}]}) / 1000000.0" | bc)"
done

echo "... done running for performance in MFLOPS."
