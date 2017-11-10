# DiphantineSolverMk3
Diophantine Equation Solver

This is functionally the same as the Mk2 version. It now uses the Boost Multi-precision library on top of the MPIR library. 
This allows the bigintegers to be handled more or less like normal integers. (The native MPIR functions are a bit like using assembler,
Each operation would need a separate function call, whereas with Boost normal expressions can be used)
