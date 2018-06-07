#!/bin/bash
# use predefined variables to access passed arguments
#echo arguments to the shell

DB_NAME=random-walk
DB_TRAIN=home/aocsa/datasets/$DB_NAME/$DB_NAME.TRAIN
DB_TEST=home/aocsa/datasets/$DB_NAME/$DB_NAME.TEST
DB_SIZE=10000
QUERIES=10
DIM=64
SUCCESS=0.99
nNN=10
DISTANCE_TYPE=1
LEVELS=10
THREADS=256
MAX_TOPK=1024
LAMBDA_STEP=1000
ARCH=sm_61
MAIN=cuda_lsh.cu
EXE=cudalsh-$DB_NAME.exe

echo "Arguments to TopKMotifs"
echo "DB_NAME = $DB_NAME"
echo "DB_TRAIN = $DB_TRAIN"
echo "DB_TEST = $DB_TEST"
echo "DB_SIZE = $DB_SIZE"
echo "EXE = $EXE"

nvcc $MAIN -O3 -w -o $EXE -arch=$ARCH -DN_THREADS=$THREADS -DMAX_TOPK=$MAX_TOPK -DLAMBDA_STEP=$LAMBDA_STEP

echo  "Execution"
nvcc -run ./$EXE -run-args $DB_SIZE,$QUERIES,$DIM,$SUCCESS,$nNN,$DB_TRAIN,$DB_TEST,$DISTANCE_TYPE,$LEVELS