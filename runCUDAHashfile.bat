:: CUDA Pair Motif

:: 1. Compilation
@ECHO OFF
set SOURCE=src
set BIN=bin
IF NOT EXIST %BIN% MKDIR %BIN%
@ECHO ON

@echo ===========================================================
set DB_NAME=rwalk
set DB_TRAIN=E:/datasets/%DB_NAME%/%DB_NAME%.TRAIN
set DB_TEST=E:/datasets/%DB_NAME%/%DB_NAME%.TEST
set DB_SIZE=100
set QUERIES=10
set DIM=1024
set SUCCESS=0.99
set nNN=10
set DISTANCE_TYPE=1
set LEVELS=10
set THREADS=256
set MAX_TOPK=1024
set LAMBDA_STEP=1000
set ARCH=sm_12
set MACHINE=64
set MAIN=cuda_lsh.cu
set EXE=%BIN%/cudalsh-%DB_NAME%.exe
@echo ===========================================================

nvcc %SOURCE%/%MAIN% -O3 -w -o %EXE% -arch=%ARCH% -m%MACHINE% -DN_THREADS=%THREADS% -DMAX_TOPK=%MAX_TOPK% -DLAMBDA_STEP=%LAMBDA_STEP%

:: 2. Execution
nvcc -run %EXE% -run-args %DB_SIZE%,%QUERIES%,%DIM%,%SUCCESS%,%nNN%,%DB_TRAIN%,%DB_TEST%,%DISTANCE_TYPE%,%LEVELS%