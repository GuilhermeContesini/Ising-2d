# Evolutionary Strategies in Game Theory


Read and use Make!

Example of execution
./main.out 256 12 100000 16 123456

Example of Parallel use (require bash-parallel):
parallel --bar --jobs 50% ./main.out 1024 12 1000000 32 $1 ::: {123456..654321..50000}
