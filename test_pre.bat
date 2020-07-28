:: Script to presolve problems in .mps format and write them to .pre files
:: Author: Ivan Viviani

@echo off

setlocal enabledelayedexpansion
for %%f in (test\*.mps) do ( 
	::echo "fullname: %%f"
    ::echo "name: %%~nf"
	cplex -c "read %%f" "write test_pre\%%~nf.pre"
)