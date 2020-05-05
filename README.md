# Evolutionary Strategies in Game Theory

Example of execution
./main.out 256 12 100000 16 123456


parallel --bar --jobs 50% ./main.out 1024 12 1000000 32 $1 ::: {123456..654321..50000}