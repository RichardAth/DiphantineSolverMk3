# DiphantineSolverMk3
Diophantine Equation Solver

This is functionally the same as the Mk2 version. It now uses the Boost Multi-precision library on top of the MPIR library. 
This allows the bigintegers to be handled more or less like normal integers. (The native MPIR functions are a bit like using assembler,
Each operation would need a separate function call, whereas with Boost library normal expressions can be used, similar to Python)

A certain amount of tidying up has been done to make the code a bit less obscure e.g. reduce the usage of global variables to pass data around.

A useful feature of the Boost multi-precision library is that integers are automatically extended  to multi-precision when necessary. Conversion the other way requires a special function I wrote, which checks that the number will actually fit into 64 bits, otherwise the program aborts.
