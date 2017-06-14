#!/bin/bash
#PBS -N pp2g_block
#PBS -e /gpfs_scratch/ncasler/logs/pp2g/pp2g-block.err
#PBS -o /gpfs_scratch/ncasler/logs/pp2g/pp2g-block.out
#PBS -S /bin/bash
#PBS -l walltime=1:00:00
#PBS -l nodes=5:ppn=20

module load boost shapelib gdal2-stack MPICH

BIN="/home/ncasler/app-dev/p_points2grid/bin/test/BlockLas"
INPUT="/projects/isgs/lidar/champaign/las"
#INPUT="/projects/isgs/lidar/cumberland/las"
SCRATCH="/gpfs_scratch/ncasler/data/tmp"
BUF_SIZE=50000000
OUTPUT="/projects/isgs/output/pp2g/blocks"
RES=5
PROCS=100

START=`date +%s`

mpirun -np $PROCS $BIN -i "${INPUT}" -s "${SCRATCH}" -o "${OUTPUT}" -b $BUF_SIZE -r $RES

END=`date +%s`

echo "Exec time: $((END - $START))"
