#################################################
#  simplifed version
#################################################
#!/bin/sh

#SBATCH --job-name=1sk_long
#SBATCH --time=48:00:00
#SBATCH --ntasks=32

mpirun ./OSHUN-2D-im.e

